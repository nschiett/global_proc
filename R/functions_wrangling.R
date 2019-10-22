
## combine parameters ##

load("data/fishtree_glob.RData")


##### sptl #####
get_sptl <- function(x){
  sptl <- read.csv("data/species_list_tl.csv")
  sptl$Species <- gsub("_", " ", sptl$species)
  tax <- rfishbase::load_taxa()
  sptl <- right_join(dplyr::select(tax, Species, Family), sptl)
  sptl[sptl$species == "Rhinesomus_triqueter", "Family"] <- "Ostraciidae"
  sptl[sptl$species == "Abudefduf_luridus", "Family"] <- "Pomacentridae"
  sptl[sptl$Family == "Scaridae", "Family"] <- "Labridae"
  return(sptl)
}

##### metabolic parameters #####
get_metpar <- function(kmax, sptl){
  metpar <- read.csv("data/metpar_fam_smr.csv")
  B0_mean <- mean(metpar$B0)
  B0_sd <- sd(metpar$B0)
  alpha_mean <- mean(metpar$alpha)
  alpha_sd <- sd(metpar$alpha)
  theta_mean <- mean(metpar$theta)
  
  metpar <- metpar %>% full_join(unique(select(sptl, Family)))
  metpar[is.na(metpar$alpha), "alpha"] <- alpha_mean
  metpar[is.na(metpar$alpha_sd), "alpha_sd"] <- alpha_sd
  metpar[is.na(metpar$B0), "B0"] <- B0_mean
  metpar[is.na(metpar$B0_sd), "B0_sd"] <- B0_sd
  metpar[is.na(metpar$theta), "theta"] <- theta_mean
  
  spt <- unique(dplyr::select(kmax, sst))
  spt <-  expand.grid(spt$sst, unique(sptl$Family))
  colnames(spt) <- c("sst", "Family")
  
  param <- left_join(spt, metpar)
  param$B0_adj <- 
    param$B0 * exp(0.59 / 8.62e-5 * (1 / (28 + 273.15) - 1 / (param$sst + 273.15)))
  metpar <- dplyr::select(param, Family, sst,
                   alpha_m = alpha, alpha_sd,
                   f0_m = B0_adj, f0_sd = B0_sd, theta_m = theta)
  
  return(metpar)
}


##### cnp #####
get_cnp <- function(){
  cnp <- read.csv("data/cnp_all_extrapolated.csv")
  return(cnp)
}

##### cnpdiet #####
get_cnpdiet <- function(){
  cnpdiet <- read.csv("data/cnpdiet_glob.csv")
  return(cnpdiet)
}

##### growth #####
get_kmax <- function(){
  kmax <- read.csv("data/kmax_predict.csv") %>% 
    select(species, sizemax = max_size, sst, k_m = kmax_m, k_sd = kmax_sd)
  return(kmax)
}

##### find lw #####
get_lw <- function(sptl){
lwl <- map(as.list(sptl$Species), possibly(fishflux::find_lw, otherwise =  NA_real_))
lwl2 <- lapply(lwl, function(x){
   if (is.na(x)){
     res <- data.frame(
       species = NA, lwa_m = NA, lwa_sd = NA, lwb_m = NA, lwb_sd = NA
     )
     return(res)
   } else{
     return(x)
   }
 }) 

 lwd <- plyr::ldply(lwl2) %>% drop_na()
 lwd$species <- gsub(" ", "_", lwd$species)
 traits <-left_join(select(sptl, species), lwd) %>% dplyr::select( species, lwa_m, lwb_m, lwa_sd, lwb_sd) %>% ungroup()
 treesp <- 
   data.frame(species = fishtree[[1]]$tip.label)
 traits <- dplyr::left_join(treesp, traits)
 
 test <- phylopars(traits, fishtree[[1]])
 recon <- test$anc_recon %>% as.data.frame()
 recon$species <- rownames(recon)
 recon <- recon[1:1110,]
 lw_glob <- recon
 return(lw_glob)
}

##### ar and troph #####

get_artr <- function(sptl){
 ar <- map(as.list(sptl$Species), possibly(fishflux::aspect_ratio, otherwise =  NA_real_))
 ar2 <- lapply(ar, function(x){x$aspect_ratio}) %>% unlist()
 tr <- map(as.list(sptl$Species), possibly(fishflux::trophic_level, otherwise =  NA_real_))
 tr2 <- lapply(tr, function(x){x$trophic_level}) %>% unlist()

 res <- sptl %>% mutate(h_m = tr2, r_m = ar2)

 # Family level average for missing
 ar_fam <- group_by(red, Family) %>% 
   summarize(ar_fam = mean(r_m, na.rm = TRUE), tr_fam = mean(h_m, na.rm = TRUE))
 res <- left_join(res, ar_fam)

 res[is.na(res$r_m), "r_m"] <- res[is.na(res$r_m), "ar_fam"]
 res[is.na(res$h_m), "h_m"] <- res[is.na(res$h_m), "tr_fam"]

 res <- dplyr::select(res, Species, Family, r_m, h_m)
 
 # add weight prop
 wp <- lapply(res$Family, fishflux::wprop) %>% plyr::ldply()
 res$mdwm <- wp$ww
 
 return(res)
}

##### combine all so far #####
combine_params <- function(sptl, kmax, lw, cnp, cnpdiet, metpar, artr){
  params <- dplyr::full_join(sptl, kmax) %>% 
    dplyr::left_join(cnp) %>% 
    dplyr::left_join(cnpdiet) %>% 
    dplyr::left_join(metpar) %>%
    dplyr::left_join(lw) %>% 
    dplyr::left_join(artr) %>%
    mutate(v_m = sst) %>% select(-sst)
  return(params)
}

get_tfish_unique <- function(){
  tfish_unique <- read.csv("data/tfish_size_temp.csv")
  tfish_unique$sst <- round(tfish_unique$sst)
  tfish_unique <- unique(tfish_unique) 
  return(tfish_unique)
}

##### run_fishflux
run_fishflux <- function(tfish_unique, params){

  tfish_unique$v_m <- round(tfish_unique$sst)
  tfish_unique <- unique(tfish_unique) %>% select(-sst)
  
  data <- left_join(select(tfish_unique, species, size_cm, v_m), params)
  
  cnpflux <- lapply(1:nrow(data), function(x){
    dt <- data[1,]
    mod <- fishflux::cnp_model_mcmc(TL = data[x, "size_cm"],
                                    param = list(
                                      Qc_m = data[x, "C"],
                                      Qn_m = data[x, "N"],
                                      Qp_m = data[x, "P"],
                                      Qc_sd = data[x, "C_sd"],
                                      Qn_sd = data[x, "N_sd"],
                                      Qp_sd = data[x, "P_sd"],
                                      Fc_m = data[x, "fc"],
                                      Fn_m = data[x, "fn"],
                                      Fp_m = data[x, "fp"],
                                      Fc_sd = data[x, "fc_sd"],
                                      Fn_sd = data[x, "fn_sd"],
                                      Fp_sd = data[x, "fp_sd"],
                                      Linf_m = data[x, "max_size"],
                                      k_m = data[x, "kmax"],
                                      k_sd = data[x, "kmax_sd"],
                                      lwa_m = data[x, "lwa_m"],
                                      lwa_sd = data[x, "lwa_sd"],
                                      lwb_m = data[x, "lwb_m"],
                                      lwb_sd = data[x, "lwb_sd"],
                                      asp_m = data[x, "ar"],
                                      troph_m = data[x, "troph"],
                                      temp_m = data[x, "sst"],
                                      B0_m = data[x, "B0_adj"],
                                      B0_sd = data[x, "B0_sd"],
                                      a_m = data[x, "alpha"],
                                      a_sd = data[x, "alpha_sd"],
                                      f_m = 2, 
                                      Tn_m = 0.0037,
                                      Tn_sd = 0.0046,
                                      Tp_m = 0.00037,
                                      Tp_sd = 0.00051,
                                      w_prop_m = data[x, "wprop"],
                                      AEc_m = 0.8, AEn_m = 0.8, AEp_m = 0.7
                                    )
    )
    extr <- fishflux::extract(mod, par = c("N_ex", "P_ex", "C_in", "N_in", "P_in","C_g", "C_eg", "N_eg", "P_eg", "C_r"))
    return(extr)
  }) %>% plyr::ldply()
  
  return(cnpflux)
}

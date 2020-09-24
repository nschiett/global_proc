
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
  cnpdiet[,3:8] <- cnpdiet[,3:8] * 100
  colnames(cnpdiet) <-
    c("family", "species", "Dc_m", "Dc_sd", 
      "Dn_m", "Dn_sd", "Dp_m", "Dp_sd", "diet_cat")
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
 ar2 <- lapply(ar, function(x){
   if(is.data.frame(x)){
     return(x$aspect_ratio)
   } else{return(NA)}
   }) %>% unlist()
 tr <- map(as.list(sptl$Species), possibly(fishflux::trophic_level, otherwise =  NA_real_))
 tr2 <- lapply(tr, function(x){
   if(is.data.frame(x)){
     return(x$trophic_level)
   } else{return(NA)}
 }) %>% unlist()
 res <- sptl %>% mutate(h_m = tr2, r_m = ar2)

 # Family level average for missing
 ar_fam <- group_by(res, Family) %>% 
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
  params <- dplyr::full_join(select(sptl, Family, Species, species), kmax) %>% 
    dplyr::left_join(cnp) %>% 
    dplyr::left_join(select(cnpdiet, -family)) %>% 
    dplyr::left_join(metpar) %>%
    dplyr::left_join(lw) %>% 
    dplyr::left_join(artr) %>%
    mutate(v_m = sst, linf_m = sizemax) %>% select(-sst, -sizemax) %>% 
    mutate( F0nz_m = 0.0037, F0nz_sd = 0.0046, F0pz_m = 0.00037, F0pz_sd = 0.00051,
            ac_m = 0.8, an_m = 0.8, ap_m = 0.7)
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
  
  cnpflux <- parallel::mclapply(1:nrow(tfish_unique), function(x){
    print(x)
    
    dt <- data[x,] 
    par <- dt %>% select(-species, - Family, - size_cm, -Species, -diet_cat) %>% as.list()
    mod <- fishflux::cnp_model_mcmc(TL = dt$size_cm,
                                    param = par, iter = 2000, seed = x)
    

    extr <- fishflux::extract(mod, par = c("F0c", "F0n", "F0p", "Gc", "Gn", "Gp", "Sc", "Sn", "Sp", 
                                           "Ic", "In", "Ip", "Wc", "Wn", "Wp", "Fc", "Fn", "Fp"))
    extr <- cbind(dt[,1:5], extr) 
    lim <- fishflux::limitation(mod, plot = FALSE)
    extr$limitation <-first(lim[lim$prop_lim == max(lim$prop_lim), "nutrient"])
    
     return(extr)
  }, mc.cores = 50) %>% plyr::ldply()
  
  
  return(cnpflux)
}

##### summarise ecosystem services per transect #####
summarise_pertransect <- function(tfish, cnpflux, params){

sst <- read.csv("data/avSst.csv")
tfish <- left_join(tfish, sst)
tfish$v_m <- round(tfish$mean)

tfi <- left_join(tfish, params)

tfi <- mutate(tfi, biomass = abun * (lwa_m * size_cm ^ lwb_m))

tf_sum <- tfi %>% group_by(sites, transect_id, area) %>%
  dplyr::summarise(biomass_tot = sum(biomass),
                   nspec = length(unique(species)),
                   biomass_c = sum(biomass[which(Dp_m > 0.33)]), 
                   biomass_h = sum(biomass[which(Dp_m < 0.33)]),
                   abu_tot = sum(abun)) %>% ungroup() %>%
  mutate(biomass_tot = biomass_tot/area,
         biomass_c = biomass_c/area,
         biomass_h = biomass_h/area,
         abu_tot = abu_tot/area)

nt <- table(tf_sum$sites) %>% as.data.frame()

lengths <- lapply(1:nrow(tfi), FUN = function(x){
  res <- data.frame(
    transect_id = rep(tfi$transect_id[x],tfi$abun[x]),
    species = rep(tfi$species[x],tfi$abun[x]),
    size_cm = rep(tfi$size_cm[x], tfi$abun[x]),
    size_max = rep(tfi$max_size[x], tfi$abun[x])
  )
  return(res)
})

lengths <- plyr::ldply(lengths) 
ll <- group_by(lengths, transect_id) %>% summarise(logtl_m = mean(log(size_cm)), size_median = median(size_cm), sizemax_median = median(size_max))
tf_sum <- left_join(tf_sum, ll)

##### nutrient fluxes #####

tfish <- left_join(tfish, cnpflux) %>% left_join(select(params, Species, v_m, diet_cat)) %>% mutate(
  Fc = Fc_median * abun,
  Fn = Fn_median * abun,
  Fp = Fp_median * abun,
  Ic = Ic_median * abun,
  In = In_median * abun,
  Ip = Ip_median * abun,
  Gc = Gc_median * abun,
  Gn = Gn_median * abun,
  Gp = Gp_median * abun,
  Wc = Wc_median * abun,
  Wn = Wn_median * abun,
  Wp = Wp_median * abun)

tfs <- group_by(tfish, studyName, region, locality, sites, area, transect_id, lon, lat) %>% 
  dplyr::summarise(
    Fc = sum(Fc)/unique(area),
    Fn = sum(Fn)/unique(area),
    Fp = sum(Fp)/unique(area),
    Ic = sum(Ic)/unique(area),
    In = sum(In)/unique(area),
    Ip = sum(Ip)/unique(area),
    Gc = sum(Gc)/unique(area),
    Gn = sum(Gn)/unique(area),
    Gp = sum(Gp)/unique(area),
    Wc = sum(Wc)/unique(area),
    Wn = sum(Wn)/unique(area),
    Wp = sum(Wp)/unique(area)
  ) %>% ungroup()


hum <- read.csv("data/humanPopulation.csv")
geo <- read.csv("data/geogrVar.csv")

tfs <- ungroup(tfs) %>% left_join(geo) 

tfs <- left_join(tfs, tf_sum)

tfs <- filter(tfs, nspec > 5)

regions <- data.frame(
  region = unique(tfs$region),
  bioregion = c("c_indopacific", "c_pacific", "w_atlantic", "e_atlantic", "w_atlantic", "e_pacific", "w_atlantic",
                "c_pacific", "c_pacific", "c_pacific", "c_pacific", "w_indian", "e_pacific", "w_indian"))

tfs <- left_join(tfs, regions)

# community variables

flux_tall <- left_join(tfish, select(params, Species, species, k_m, linf_m, h_m))

# immaturity
flux_tall$imm <- flux_tall$k_m * (flux_tall$linf_m - flux_tall$size_cm)
flux_tall$imm[flux_tall$imm<0] <- 0

# community variables per transect
comm <- 
  lapply(unique(flux_tall$transect_id), function(i){
    sub <-  filter(flux_tall, transect_id == i) %>% tidyr::uncount(abun)
    imm_wci <- quantile(sub$imm, 0.975) - quantile(sub$imm, 0.025)
    imm_q1 <- quantile(sub$imm, 0.025)
    imm_q3 <- quantile(sub$imm, 0.975)
    imm_m <- median(sub$imm)
    size_wci <- quantile(sub$size_cm, 0.975) - quantile(sub$size_cm, 0.025)
    size_q3 <- quantile(sub$size_cm, 0.975) 
    size_q1 <- quantile(sub$size_cm, 0.025) 
    size_m <- median(sub$size_cm)
    sizemax_wci <- quantile(sub$linf_m, 0.975) - quantile(sub$linf_m, 0.025)
    sizemax_q3 <- quantile(sub$linf_m, 0.975) 
    sizemax_q1 <- quantile(sub$linf_m, 0.025) 
    sizemax_m <- median(sub$linf_m)
    troph_wci <- quantile(sub$h_m, 0.975) - quantile(sub$h_m, 0.025)
    troph_q3 <- quantile(sub$h_m, 0.975) 
    troph_q1 <- quantile(sub$h_m, 0.025) 
    troph_m <- median(sub$h_m)
    imm_kur <- kurtosis(sub$imm)
    imm_ske <- skewness(sub$imm)
    size_kur <- kurtosis(sub$size_cm)
    size_ske <- skewness(sub$size_cm)
    sizemax_kur <- kurtosis(sub$linf_m)
    sizemax_ske <- skewness(sub$linf_m)
    troph_kur <- kurtosis(sub$h_m)
    troph_ske <- skewness(sub$h_m)
    result <- data.frame(transect_id = i, imm_m, imm_wci, size_m, size_wci, troph_m, troph_wci,
                         imm_kur, imm_ske, size_kur, size_ske, troph_kur, troph_ske,
                         imm_q1, troph_q3, size_q3, imm_q3, troph_q1, size_q1, sizemax_m, 
                         sizemax_q1, sizemax_q3, sizemax_wci, sizemax_kur, sizemax_ske)
    return(result)
  }) %>% plyr::ldply()

## combine all variables per transect
tfs <- left_join(tfs, comm)

return(tfs)
}

##### Proportions contributions #####

get_contributions <- function(tfish, cnpflux, params, summary_transect){
  
  tfish <- left_join(tfish, cnpflux) %>% 
    left_join(select(params, Species, lwa_m, lwb_m, diet_cat)) %>% 
    mutate(
      Fc = Fc_median * abun,
      Fn = Fn_median * abun,
      Fp = Fp_median * abun,
      Ic = Ic_median * abun,
      In = In_median * abun,
      Ip = Ip_median * abun,
      Gc = Gc_median * abun,
      Gn = Gn_median * abun,
      Gp = Gp_median * abun,
      Wc = Wc_median * abun,
      Wn = Wn_median * abun,
      Wp = Wp_median * abun,
      biomass = abun * lwa_m * size_cm^lwb_m)
  
  tfs <- group_by(tfish, studyName, region, locality, sites, area, transect_id, lon, lat) %>% 
    dplyr::summarise(
      Fn_s = sum(Fn),
      Fp_s = sum(Fp),
      Ic_s = sum(Ic),
      In_s = sum(In),
      Ip_s = sum(Ip),
      Gc_s = sum(Gc),
      Gn_s = sum(Gn),
      Gp_s = sum(Gp),
      Fc_s = sum(Fc),
      Wc_s = sum(Wc),
      Wn_s = sum(Wn),
      Wp_s = sum(Wp),
      biomass_s = sum(biomass),
      biomass_h_s = sum(biomass[diet_cat == 2 & Family != "Mugilidae"]),
      biomass_p_s = sum(biomass[diet_cat == 4])
    ) %>% right_join(tfish) %>% 
    group_by(studyName, region, locality, sites, area, transect_id, lon, lat, 
             Family, species, diet_cat) %>%
    dplyr::summarise(
      Fn_p = sum(Fn)/unique(Fn_s),
      Fp_p = sum(Fp)/unique(Fp_s),
      Ic_p = sum(Ic)/unique(Ic_s),
      In_p = sum(In)/unique(In_s),
      Ip_p = sum(Ip)/unique(Ip_s),
      Gc_p = sum(Gc)/unique(Gc_s),
      Gn_p = sum(Gn)/unique(Gn_s),
      Gp_p = sum(Gp)/unique(Gp_s),
      Fc_p = sum(Fc)/unique(Fc_s),
      Wc_p = sum(Wc)/unique(Wc_s),
      Wn_p = sum(Wn)/unique(Wn_s),
      Wp_p = sum(Wp)/unique(Wp_s),
      biomass_p = sum(biomass)/unique(biomass_s),
      biomass_herb_p = sum(biomass)/unique(biomass_h_s),
      biomass_pisc_p = sum(biomass)/unique(biomass_p_s)
    ) %>%
    # 
    mutate(
      biomass_herb_p = case_when(
        diet_cat == 2 & Family != "Mugilidae" ~ biomass_herb_p,
        TRUE ~ 0),
      biomass_pisc_p = case_when(
        diet_cat == 4 ~ biomass_pisc_p,
        TRUE ~ 0))
  
  tfss <- left_join(tfs, select(summary_transect, transect_id, nspec, bioregion)) %>% 
    filter(!is.na(nspec))
  
  return(tfss)
}


##### supplement herbivory carnivory #####

get_herb_pisc <- function(tfish, cnpflux, params, summary_transect){

  
  sst <- read.csv("data/avSst.csv")
  tfish <- left_join(tfish, sst)
  tfish$v_m <- round(tfish$mean)
  
tfish <- left_join(tfish, unique(dplyr::select(cnpflux, species, v_m, size_cm, (ends_with("_median"))))) %>% 
  left_join(select(params, Species, v_m, diet_cat, Dc_m)) %>% 
  mutate(Ic = Ic_median * abun) %>%
  mutate(I = Ic/area) 


tfish$herb <- "no"
tfish[tfish$diet_cat == 2, "herb"] <- "yes"
tfish[tfish$Family == "Mugilidae", "herb"] <- "no"
tfish$pisc <- "no"
tfish[tfish$diet_cat == 4, "pisc"] <- "yes"

contr <- dplyr::group_by(tfish, region, locality, sites, area, transect_id) %>% 
  dplyr::summarise(
    I_herb_s = sum(I[herb == "yes"], na.rm = TRUE),
    I_pisc_s = sum(I[pisc == "yes"], na.rm = TRUE)
  ) %>% dplyr::ungroup() %>% unique() %>% dplyr::right_join(tfish)  %>%
  dplyr::group_by(region, locality, sites, area, transect_id, species, I_herb_s, I_pisc_s) %>%
  dplyr::summarise(
    I_herb_p = sum(I[herb == "yes"], na.rm = TRUE)/unique(I_herb_s),
    I_pisc_p = sum(I[pisc == "yes"], na.rm = TRUE)/unique(I_pisc_s)
  ) %>% ungroup()

tfs <- dplyr::group_by(tfish, region, locality, sites, area, transect_id) %>% 
  dplyr::summarise(
    I_herb = sum(I[herb == "yes"], na.rm = TRUE),
    I_pisc = sum(I[pisc == "yes"], na.rm = TRUE)
  ) %>% ungroup()

contr <- left_join(contr, select(summary_transect, transect_id, nspec, bioregion)) %>% 
  filter(!is.na(nspec)) %>% 
  mutate(I_herb_p = replace_na(I_herb_p, 0)) %>%
  mutate(I_pisc_p = replace_na(I_pisc_p, 0))


return(list(contributions_herb_pisc = contr, summary_herb_pisc = tfs))
}

##### final transect summary #####
combine_summary <- function(summary_transect, herb_pisc){
  
  flux <- summary_transect
  hp <- herb_pisc$summary_herb_pisc
  sst <- read.csv("data/avSst.csv")
  flux <- left_join(flux, hp) %>% left_join(sst)
  
  return(flux)
  
}


##### fit contribution models #####
fit_contribution_Fp <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("Fp_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    return(fit_new)
  }, mc.cores = 7)
  
  return(Fp = list)
}

fit_contribution_Fn <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("Fn_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    return(fit_new)
  }, mc.cores = 7)
  
  return(Fn = list)
}

fit_contribution_Fc <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("Fc_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    return(fit_new)
  }, mc.cores = 7)
  
  return(Fc = list)
}

fit_contribution_Ic <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("Ic_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    return(fit_new)
  }, mc.cores = 7)
  
  return(Ic = list)
}

fit_contribution_In <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("In_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    return(fit_new)
  }, mc.cores = 7)
  
  return(In = list)
}

fit_contribution_Ip <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("Ip_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    return(fit_new)
  }, mc.cores = 7)
  
  return(Ip = list)
}




fit_contribution_Gn <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("Gn_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  subsp[1:10, "Gn_p"] <- 0
  fitzi <- brm("Gn_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "zero_inflated_beta", data = subsp)
  
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    if(min(subsp$Gn_p) == 0){
      fit_new <- update(fitzi, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    }else{
      fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    }
    return(fit_new)
  }, mc.cores = 6)
  
  return(Gn = list)
}


fit_contribution_Gc <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("Gc_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  subsp[1:10, "Gc_p"] <- 0
  fitzi <- brm("Gc_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "zero_inflated_beta", data = subsp)
  
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    if(min(subsp$Gc_p) == 0){
      fit_new <- update(fitzi, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    }else{
      fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    }
    return(fit_new)
  }, mc.cores = 6)
  
  return(Gc = list)
}

fit_contribution_Gp <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("Gp_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  subsp[1:10, "Gp_p"] <- 0
  fitzi <- brm("Gp_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "zero_inflated_beta", data = subsp)
  
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    if(min(subsp$Gp_p) == 0){
      fit_new <- update(fitzi, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    }else{
      fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    }
    return(fit_new)
  }, mc.cores = 6)
  
  return(Gp = list)
}


fit_contribution_Wp <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("Wp_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    return(fit_new)
  }, mc.cores = 7)
  return(Wp = list)
}

fit_contribution_Wn <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("Wn_p ~1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    return(fit_new)
  }, mc.cores = 7)
  return(Wn = list)
}

fit_contribution_Wc <- function(contributions, sp_loc){
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("Wc_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    return(fit_new)
  }, mc.cores = 7)
  return(Wc = list)
}

fit_contribution_I_herb <- function(herb_pisc, sp_loc){
  
  # get herbivores
  contributions <- herb_pisc$contributions_herb_pisc %>%
    filter(I_herb_p > 0) %>%
    mutate(I_herb_p = case_when(I_herb_p == 1 ~ 0.999999, !I_herb_p == 1 ~ I_herb_p))

  sp_loc <- filter(sp_loc, species %in% unique(contributions$species))
    
  ## compile model once 
  subsp <- filter(contributions, bioregion == "c_indopacific", species == "Acanthurus_leucocheilus")
  fit <- brm("I_herb_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    return(fit_new)
  }, mc.cores = 7)
  
  return(I_herb = list)
}

fit_contribution_I_pisc <- function(herb_pisc, sp_loc){
  
  # get carnivores
  contributions <- herb_pisc$contributions_herb_pisc %>%
    filter(I_pisc_p > 0) %>%
    mutate(I_pisc_p = case_when(I_pisc_p == 1 ~ 0.999999, !I_pisc_p == 1 ~ I_pisc_p))
  
  sp_loc <- filter(sp_loc, species %in% unique(contributions$species))
  
  ## compile model once 
  subsp <- filter(contributions, bioregion == "w_atlantic", species == "Cephalopholis_cruentata")
  fit <- brm("I_pisc_p ~ 1 + (1|locality) + (1|sites)", chains = 1, family = "beta", data = subsp)
  
  # run same model on all species
  list <- parallel::mclapply(1:nrow(sp_loc), function(x){
    print(x)
    subsp <- filter(contributions, bioregion == sp_loc$bioregion[x], species == sp_loc$species[x])
    fit_new <- update(fit, newdata = subsp, recompile = FALSE, chains = 2, control = list(adapt_delta = 0.99), save_ranef = FALSE)
    return(fit_new)
  }, mc.cores = 7)
  
  return(I_pisc = list)
}


####### extract contributions from models ######
expit <- function(x){exp(x)/(1 + exp(x))}

extract_contr <- function(list){
summary <- parallel::mclapply(list, function(x){
  post <- expit(brms::posterior_samples(x, "b_Intercept"))
  phi <- brms::posterior_samples(x, "phi")
  result <- data.frame(
    contr_m = median(post$b_Intercept),
    contr_sd = sd(post$b_Intercept),
    contr_q1 = quantile(post$b_Intercept, 0.025),
    contr_q3 = quantile(post$b_Intercept, 0.975),
    phi_m = median(phi$phi),
    phi_sd = sd(phi$phi)
  )
  
  return(result)
}, mc.cores = 50) %>% plyr::ldply()

n <- names(list[[1]]$data)[1]

colnames(summary) <- c(
  paste0(n, "_m"), 
  paste0(n, "_sd"), 
  paste0(n, "_q1"), 
  paste0(n, "_q3"),
  paste0(n, "phi_m"),
  paste0(n, "phi_sd"))

return(summary)
}

combine_contributions <- function(sp_loc, contributions_Fn, 
                                  contributions_Fp, contributions_Gc, 
                                  contributions_I_herb, contributions_I_pisc, herb_pisc){
  
  herb <-  herb_pisc$contributions_herb_pisc %>%
    filter(I_herb_p > 0) 
  
  sploc_h <- filter(sp_loc, species %in% unique(herb$species)) %>%
    cbind(contributions_I_herb)
  
  pisc <-  herb_pisc$contributions_herb_pisc %>%
    filter(I_pisc_p > 0) 
  
  sploc_p <- filter(sp_loc, species %in% unique(pisc$species)) %>%
    cbind(contributions_I_pisc)
  
  comb <- cbind(
    sp_loc, contributions_Fn, contributions_Fp,
    contributions_Gc
  ) %>% left_join(sploc_h) %>% left_join(sploc_p) %>% replace(., is.na(.), 0)
  
  return(comb)
  
}

add_occurence <- function(contributions, contributions_sp_loc, herb_pisc){
  tot <- contributions %>%
    group_by(bioregion) %>%
    summarise(
      tot = n()) 
  occ <- contributions %>%
    group_by(bioregion, species) %>%
    summarise(occurence = n()) %>%
    left_join(tot) %>% mutate(rel_occ = occurence/tot) %>% 
    dplyr::select(-tot)
  result <- left_join(contributions_sp_loc, occ)
  
  correct <-  select(result, 1:2) %>%
    inner_join(select(ungroup(contributions), bioregion, species, Fn_p, Fp_p, Gc_p)) %>%
    left_join(herb_pisc$contributions_herb_pisc) %>%
    group_by(bioregion, species) %>%
    summarise(
      # Fc_p_m = median(Fc_p),
      # Fc_p_sd = sd(Fc_p),
      # Fc_p_q1 = quantile(Fc_p, 0.025),
      # Fc_p_q3 = quantile(Fc_p, 0.975),
      Fn_p_m = median(Fn_p),
      Fn_p_sd = sd(Fp_p),
      Fn_p_q1 = quantile(Fn_p, 0.025),
      Fn_p_q3 = quantile(Fn_p, 0.975),
      Fp_p_m = median(Fp_p),
      Fp_p_sd = sd(Fp_p),
      Fp_p_q1 = quantile(Fp_p, 0.025),
      Fp_p_q3 = quantile(Fp_p, 0.975),  
      # Ic_p_m = median(Ic_p),
      # Ic_p_sd = sd(Ic_p),
      # Ic_p_q1 = quantile(Ic_p, 0.025),
      # Ic_p_q3 = quantile(Ic_p, 0.975),
      Gc_p_m = median(Gc_p),
      Gc_p_sd = sd(Gc_p),
      Gc_p_q1 = quantile(Gc_p, 0.025),
      Gc_p_q3 = quantile(Gc_p, 0.975),
      I_herb_p_m = median(I_herb_p),
      I_herb_p_sd = sd(I_herb_p),
      I_herb_p_q1 = quantile(I_herb_p, 0.025),
      I_herb_p_q3 = quantile(I_herb_p, 0.975),
      I_pisc_p_m = median(I_pisc_p),
      I_pisc_p_sd = sd(I_pisc_p),
      I_pisc_p_q1 = quantile(I_pisc_p, 0.025),
      I_pisc_p_q3 = quantile(I_pisc_p, 0.975)
    ) %>% left_join(select(result, species, bioregion, occurence, rel_occ)) 
  
  
  return(correct)
}

# res <- sp_loc %>% mutate(con_Fp = unlist(l3))

# 
# # 
# summary(fit_new)
# 
# me <- marginal_effects(fit_new)
# 
# np <- nuts_params(fit_new)
# sum(subset(np, Parameter == "divergent__")$Value)

# l2 <- parallel::mclapply(models_contributions_Fn, function(x){
#   np <- nuts_params(x)
#   return(sum(subset(np, Parameter == "divergent__")$Value))
# }, mc.cores = 50)
# hist(unlist(l2))




# test <- pivot_longer(low, cols = c(Fn_p_m, Fp_p_m, Gc_p_m, I_herb_p_m, I_pisc_p_m)) %>%
#   filter(!value == 0)
# 
# ggplot(test, aes(x = (rel_occ), y = (value), color = bioregion, fill = bioregion)) +
#   geom_point(size = 0.5, alpha = 0.5) +
#   geom_smooth(formula = y ~ x + I(x^2), method = "lm", alpha = 0.2, se = FALSE) +
#   facet_wrap(~name, scales = "free") +
#   theme_bw() +
#   scale_color_fish_d(option = "Scarus_quoyi") 
# 
# ggplot(low, aes(y = log(Fp_p_m), x = log(rel_occ*100))) +
#   geom_point() + geom_smooth(formula = y ~ 1 + x + I(x^2)+ I(x^3) , method = "lm")
# 
# fit <- brm(data = contributions_sp_loc_occ, log(Fp_p_m) ~ 1 + log(occurence) + I(log(occurence)^2))  
# fit2 <- brm(data = contributions_sp_loc_occ, (Fp_p_m) ~ 1 + log(occurence) + I(log(occurence)^2) +  I(log(occurence)^3))  
# 
# marginal_effects(fit) +
#   scale_x_continuous(trans = "log")
# p <- fitted(fit2) %>% as.data.frame()
# 
# ggplot(contributions_sp_loc_occ, aes(y = (Fp_p_m), x = log(occurence))) +
#   geom_point() + 
#   geom_line(aes(y = (p$Estimate)), color = "red")
# 
# bayes_R2(fit2)
# marginal_effects(fit2) 
#   

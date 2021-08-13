
## combine parameters ##

load("data/fishtree_glob.RData")


##### sptl #####
get_sptl <- function(x){
  sptl <- read.csv("data/species_list_tl.csv")
  sptl$Species <- gsub("_", " ", sptl$species)
  tax <- dplyr::collect(rfishbase::load_taxa())
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
  # temperature adjustment of metabolic constant 
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
  
  cnpdiet <- read.csv("data/cnp_per_trophic_guild_1000draws.csv") %>%
    rename(trophic_guild_predicted = pred_cat)
  diet <- read.csv("data/extrapolation_trophic_guilds.csv") 
  
  
  diets <- diet %>%
    select(family, species, intersect(ends_with("_m"), starts_with("p"))) %>%
    pivot_longer(cols = intersect(ends_with("_m"), starts_with("p")),
                 names_to = "trophic_guild_predicted", values_to = "probability") %>%
    mutate(trophic_guild_predicted = as.factor(trophic_guild_predicted)) %>%
    mutate(trophic_guild_predicted = fct_recode(trophic_guild_predicted, 
                                 "1" = "p1_m",
                                 "2" = "p2_m",
                                 "3" = "p3_m",
                                 "4" = "p4_m",
                                 "5" = "p5_m",
                                 "6" = "p6_m",
                                 "7" = "p7_m",
                                 "8" = "p8_m")) %>%
    right_join(mutate(cnpdiet, trophic_guild_predicted = as.factor(trophic_guild_predicted)))
  
  result <- diets %>%
    group_by(family, species, draw) %>%
    summarise(c = sum(probability * c)/100,
              n = sum(probability * n)/100,
              p = sum(probability * p)/100) %>%
    ungroup() %>%
    group_by(family, species) %>%
    summarise(c_m = mean(c), c_sd = sd(c),
              n_m = mean(n), n_sd = sd(n),
              p_m = mean(p), p_sd = sd(p)) 
  
  result <- ungroup(result) %>% 
    left_join(select(diet, species, trophic_guild_predicted)) %>%
    ungroup()
  
  result[,3:8] <- result[,3:8] * 100
  colnames(result) <-
    c("family", "species", "Dc_m", "Dc_sd", 
      "Dn_m", "Dn_sd", "Dp_m", "Dp_sd", "trophic_guild_predicted")
  
  return(result)
}

##### growth #####
get_kmax <- function(){
  kmax <- read.csv("data/kmax_predict.csv") %>% 
    select(species, sizemax = max_size, sst, k_m = kmax_m, k_sd = kmax_sd)
  return(kmax)
}

##### get lw #####
get_lw <- function(sptl){
  readr::read_csv("data/length_weight_fishbase.csv")
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
    mutate(F0nz_m = 0.0037, F0nz_sd = 0.0046, F0pz_m = 0.00037, F0pz_sd = 0.00051,
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
    par <- dt %>% select(-species, - Family, - size_cm, -Species, -trophic_guild_predicted) %>% as.list()
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

tfi <- inner_join(tfish, params)

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

tfish <- left_join(tfish, cnpflux) %>% left_join(select(params, Species, v_m, trophic_guild_predicted)) %>% mutate(
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
    left_join(select(params, Species, lwa_m, lwb_m, trophic_guild_predicted)) %>% 
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
      biomass_h_s = sum(biomass[trophic_guild_predicted == 2 & Family != "Mugilidae"]),
      biomass_p_s = sum(biomass[trophic_guild_predicted == 4])
    ) %>% right_join(tfish) %>% 
    group_by(studyName, region, locality, sites, area, transect_id, lon, lat, 
             Family, species, trophic_guild_predicted) %>%
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
        trophic_guild_predicted == 2 & Family != "Mugilidae" ~ biomass_herb_p,
        TRUE ~ 0),
      biomass_pisc_p = case_when(
        trophic_guild_predicted == 4 ~ biomass_pisc_p,
        TRUE ~ 0))
  
  tfss <- left_join(tfs, select(summary_transect, transect_id, nspec, bioregion)) %>% 
    filter(!is.na(nspec))
  
  return(tfss)
}


##### herbivory carnivory #####

get_herb_pisc <- function(tfish, cnpflux, params, summary_transect){

  
  sst <- read.csv("data/avSst.csv")
  tfish <- left_join(tfish, sst)
  tfish$v_m <- round(tfish$mean)

diet <- read.csv("data/extrapolation_trophic_guilds.csv") 

diets <- diet %>%
  select(family, species, p2_m, p4_m, trophic_guild_predicted = trophic_guild_predicted) %>% 
  mutate(herb = trophic_guild_predicted == 2 & p2_m > 0.5,
         pisc = trophic_guild_predicted == 4 & p4_m > 0.5)

tfish <- left_join(tfish, unique(dplyr::select(cnpflux, species, v_m, size_cm, (ends_with("_median"))))) %>% 
  left_join(diets) %>%
  mutate(Ic = Ic_median * abun) %>%
  mutate(I = Ic/area) 

#tfish[tfish$Family == "Mugilidae", "herb"] <- FALSE

tfs <- dplyr::group_by(tfish, region, locality, sites, area, transect_id) %>% 
  dplyr::summarise(
    I_herb = sum(I[herb == TRUE], na.rm = TRUE),
    I_pisc = sum(I[pisc == TRUE], na.rm = TRUE),
  ) %>% ungroup() %>%
  mutate(I_herb = case_when(I_herb == 0 ~ NA_real_,
                            TRUE ~ I_herb),
         I_pisc = case_when(I_pisc == 0 ~ NA_real_,
                            TRUE ~ I_pisc))


contr <- dplyr::group_by(tfish, region, locality, sites, area, transect_id) %>% 
  dplyr::summarise(
    I_herb_s = sum(I[herb == TRUE], na.rm = TRUE),
    I_pisc_s = sum(I[pisc == TRUE], na.rm = TRUE)
  ) %>% dplyr::ungroup() %>% unique() %>% dplyr::right_join(tfish)  %>%
  dplyr::group_by(region, locality, sites, area, transect_id, species, I_herb_s, I_pisc_s) %>%
  dplyr::summarise(
    I_herb_p = sum(I[herb == TRUE], na.rm = TRUE)/unique(I_herb_s),
    I_pisc_p = sum(I[pisc == TRUE], na.rm = TRUE)/unique(I_pisc_s)
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

##### coordinates per location and site #####
get_coords_siteloc <- function(data) {
  data %>%
    group_by(bioregion, locality, sites) %>%
    summarize_at(vars(lat, lon, mean), median) %>%
    group_by(locality) %>%
    mutate(lon_loc = median(lon), lat_loc = median(lat)) 
    
}


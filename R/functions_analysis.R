standard <- function(x){
  m <- mean(x, na.rm = TRUE)
  sd <- sd(x, na.rm = TRUE)
  st <- sapply(x, FUN = function(i){
    st <- (i - m)/sd
    return(st)
  })
  return(st)
}

normalize <- function(x){
  100 * (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

##### models with biomass, sst, locality, and site effects #######
# models with only biomass and sst
run_bmmodels <- function(summary_transect_complete){
  
  flux <- summary_transect_complete
  
  # run models
  fit_Fn <- brm(log(Fn) ~ (log(biomass_tot)) + mean, 
                data = flux, cores = 4)
  
  fit_Fp <- brm(log(Fp) ~ log(biomass_tot) + mean, 
                data = flux, cores = 4)
  
  fit_Gc <- brm(log(Gc) ~ log(biomass_tot) + mean, 
                data = flux[flux$Gc > 0, ], cores = 4)
  
  fit_I_herb <- brm(log(I_herb) ~ log(biomass_tot) + mean, 
                    data = flux[flux$I_herb > 0,], cores = 4)  
  
  fit_I_pisc <- brm(log(I_pisc) ~ log(biomass_tot) + mean, 
                    data = flux[flux$I_pisc > 0,], cores = 4)
  
  return(list(fit_Fn = fit_Fn, fit_Fp = fit_Fp, fit_Gc = fit_Gc, 
              fit_I_herb = fit_I_herb, fit_I_pisc = fit_I_pisc))
}


run_procmodels <- function(summary_transect_complete){
  
  flux <- summary_transect_complete
  
  # run models
  fit_Fn <- brm(standard(log(Fn)) ~ (log(biomass_tot)) + mean + (1|locality) + (1|sites), 
                data = flux, cores = 4)
  
  fit_Fp <- brm(standard(log(Fp)) ~ log(biomass_tot) + mean + (1|locality) + (1|sites), 
              data = flux, cores = 4)
  
  fit_Gc <- brm(standard(log(Gc)) ~ log(biomass_tot) + mean + (1|locality) + (1|sites), 
                data = flux[flux$Gc > 0, ], cores = 4)
  
  fit_I_herb <- brm(standard(log(I_herb)) ~ log(biomass_tot) + mean + (1|locality) + (1|sites), 
                    data = flux[flux$I_herb > 0,], cores = 4)  
  
  fit_I_pisc <- brm(standard(log(I_pisc)) ~ log(biomass_tot) + mean + (1|locality) + (1|sites), 
                    data = flux[flux$I_pisc > 0,], cores = 4)
  
  return(list(fit_Fn = fit_Fn, fit_Fp = fit_Fp, fit_Gc = fit_Gc, 
              fit_I_herb = fit_I_herb, fit_I_pisc = fit_I_pisc))
}

get_residuals <- function(summary_transect_complete, procmodels){
  
  flux <- summary_transect_complete
  
  res1 <- residuals(procmodels[[1]], re_formula = "standard(log(Fn)) ~ log(biomass_tot) + mean")
  res2 <- residuals(procmodels[[2]], re_formula = "standard(log(Fp)) ~ log(biomass_tot) + mean")
  res3 <- residuals(procmodels[[3]], re_formula = "standard(log(Gc)) ~ log(biomass_tot) + mean")
  res4 <- residuals(procmodels[[4]], re_formula = "standard(log(I_herb)) ~ log(biomass_tot) + mean")
  res5 <- residuals(procmodels[[5]], re_formula = "standard(log(I_pisc)) ~ log(biomass_tot) + mean")

  res <- data.frame(
    transect_id = flux$transect_id,
    Fn_r = res1[,1],
    Fp_r = res2[,1],
    Gc_r = res3[,1]
  )
  
  res_h <- data.frame(
    transect_id = flux[flux$I_herb>0, ]$transect_id,
    I_herb_r = res4[,1]
  )
  
  res_p <- data.frame(
    transect_id = flux[flux$I_pisc>0, ]$transect_id,
    I_pisc_r = res5[,1]
  )
  
  
  res <- res %>% left_join(res_h) %>% left_join(res_p) 
  
  return(res)
}

get_location_effect <- function(model, summary_transect){
  result <- ranef(model)$locality[,1,] %>% as.data.frame() %>%
    rownames_to_column("locality")

  colnames(result) <- c("locality", "r_loc_Fn", "r_loc_Fp", "r_loc_Gc", 
                        "r_loc_I_herb", "r_loc_I_pisc")
  
  # add coordinates
  coord <- select(summary_transect, bioregion, locality, lat, lon) %>%
    unique() %>%
    group_by(bioregion, locality) %>%
    summarize_all(mean)
  
  result <- coord %>% left_join(result) %>% ungroup()
 
  return(result)
    
}

get_location_effect2 <- function(model, summary_transect){
  result <- ranef(model)$locality %>% as.data.frame() %>%
    rownames_to_column("locality") %>%
    select(1:2)
  
  colnames(result) <- c("locality", "r_loc_multi")
  
  # add coordinates
  coord <- select(summary_transect, bioregion, locality, lat, lon) %>%
    unique() %>%
    group_by(bioregion, locality) %>%
    summarize_all(mean)
  
  result <- coord %>% left_join(result) %>% ungroup()
  
  return(result)
  
}

####### community models ######
run_commodels <- function(summary_transect_complete){
  
  flux <- summary_transect_complete
  
  flux$logbiomass <- log(flux$biomass_tot)
  
  fit_Fn_st <- brm(standard(log(Fn)) ~ standard(mean) + standard(logbiomass) + standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
                     standard(size_q3) + standard(troph_q3) + standard(imm_q1) + standard(imm_q3) + standard(troph_q1) + standard(size_q1) ,
                   data = flux, chain = 3, cores = 3)
  
  fit_Fp_st <- brm(standard(log(Fp)) ~ standard(mean) + standard(logbiomass) + standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
                     standard(size_q3) + standard(troph_q3) + standard(imm_q1) + standard(imm_q3) + standard(troph_q1) + standard(size_q1) ,
                   data = flux, chain = 3, cores = 3)
  
  fit_Gc_st <- brm(standard(log(Gc)) ~ standard(mean) + standard(logbiomass) + standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
                     standard(size_q3) + standard(troph_q3) + standard(imm_q1) + standard(imm_q3) + standard(troph_q1) + standard(size_q1) ,
                   data = flux, chain = 3, cores = 3)
  
  fit_Iherb_st <- brm(standard(log(I_herb)) ~standard(mean) + standard(logbiomass) + standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
                        standard(size_q3) + standard(troph_q3) + standard(imm_q1) + standard(imm_q3) + standard(troph_q1) + standard(size_q1) ,
                      data = flux[flux$I_herb>0,], chain = 3, cores = 3)
  
  fit_Ipisc_st <- brm(standard(log(I_pisc)) ~ standard(mean) + standard(logbiomass) + standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
                        standard(size_q3) + standard(troph_q3) + standard(imm_q1) + standard(imm_q3) + standard(troph_q1) + standard(size_q1) ,
                      data = flux[flux$I_pisc>0,], chain = 3, cores = 3)
  
  result <- list(fit_Fn_st, fit_Fp_st, fit_Gc_st, fit_Iherb_st, fit_Ipisc_st)
  
  lapply(result, function(x){
    brms::bayes_R2(x)
  })
  
  return(result)
}

####### community models absolute values ######
run_commodels_abs <- function(summary_transect_complete){
  
  flux <- summary_transect_complete
  
  flux$logbiomass <- log(flux$biomass_tot)
  
  fit_Fn <- brm(log(Fn) ~ mean + logbiomass + nspec + size_m + troph_m + imm_m +
                     size_q3 + troph_q3 + imm_q1 + imm_q3 + troph_q1 + size_q1 ,
                   data = flux, chain = 3, cores = 1)
  
  fit_Fp <- brm(log(Fp) ~ mean + logbiomass + nspec + size_m + troph_m + imm_m +
                  size_q3 + troph_q3 + imm_q1 + imm_q3 + troph_q1 + size_q1 ,
                   data = flux, chain = 3, cores = 1)
  
  fit_Gc <- brm(log(Gc)  ~ mean + logbiomass + nspec + size_m + troph_m + imm_m +
                  size_q3 + troph_q3 + imm_q1 + imm_q3 + troph_q1 + size_q1 ,
                   data = flux, chain = 3, cores = 1)
  
  fit_Iherb <- brm(log(I_herb) ~ mean + logbiomass + nspec + size_m + troph_m + imm_m +
                     size_q3 + troph_q3 + imm_q1 + imm_q3 + troph_q1 + size_q1 ,
                      data = flux[flux$I_herb>0,], chain = 3, cores = 1)
  
  fit_Ipisc <- brm(log(I_pisc) ~ mean + logbiomass + nspec + size_m + troph_m + imm_m +
                     size_q3 + troph_q3 + imm_q1 + imm_q3 + troph_q1 + size_q1 ,
                      data = flux[flux$I_pisc>0,], chain = 3, cores = 1)
  
  result <- list(fit_Fn, fit_Fp, fit_Gc, fit_Iherb, fit_Ipisc)
  
  lapply(result, function(x){
    brms::bayes_R2(x)
  })
  
  return(result)
}

####### contributions ######

####### vulnerability ######

get_vuln <- function(sptl){
  # fshing
  vuln_fi <- data.frame(
    species = unique(sptl$Species),
    vuln_fi = rfishbase::species(gsub("_", " ", sptl$Species), field = "Vulnerability")[[1]]/100
  )
  # climate
  covu <- read.csv("data/Link_to_corals.csv") %>% dplyr::select(1,2,3,4,7, 8, 11, 12)
  colnames(covu) <- c("Gaspar_code" , "Genus_and_species" ,  "FB_Spec_Code" ,"Family" ,"Size","Size_type","diet_spec", "habitat_spec")
  covu[covu$diet_spec == "ns", "diet_spec"] <- "NS"
  covu <-  covu %>% mutate(diet_spec, diet_spec = as.character(diet_spec))
  
  graham <- read.csv("data/graham_tolerance.csv") %>% mutate(Genus_and_species = as.character(Species)) %>%
    dplyr::left_join(covu) %>% drop_na(diet_spec)
  graham[graham$diet_spec == "ns", "diet_spec"] <- "NS"
  
  # build model with diet spec, habitat_spec and size and Family as random effect
  fit_cv <- brm(data = graham,
                Climate_change_vulnerability ~ diet_spec + habitat_spec + log(Size) +  (1|Family),
                control = list(adapt_delta = 0.95))
  
  # new data
  nd <- covu
  nd[!nd$Family %in% unique(graham$Family), "Family"] <- NA #family na if not in dataframe
  nd$Family <- as.character(nd$Family)
  
  bayes_R2(fit_cv)
  # Estimate  Est.Error      Q2.5    Q97.5
  # R2 0.7866488 0.01811666 0.7439637 0.815099
  
  predict <- fitted(fit_cv, newdata = nd)
  covu$vuln_climate <- predict[,1]
  
  covu$species <- covu$Genus_and_species
  covu <- left_join(covu, vuln_fi) %>% filter(species %in% sptl$Species) %>%
    mutate(species = gsub(" ", "_", species))
  
}


##### contributions ######

get_dd <- function(contributions, herb_pisc){
  
  # combine data
  con <- left_join(contributions, herb_pisc$contributions_herb_pisc) %>%
    ungroup()
  
  ### function to calculate area of ranked function contribution !! 
  # data has to be ordered
  # y = cumsum!!
  rank_area <- function(y) {
    n <- length(y)
    a <- y[2:n]
    b <- y[1:(n-1)]
    area <- sum(a + b)/2
    return(area)
  }
  
  # dominance function
  ks_area <- function(data, var){
    
    if (var %in% c("I_herb_p", "I_pisc_p")){
      data <- data %>%
        filter(.data[[var]] > 0)
    }
    
    if (nrow(data) == 0){
      return(NA)
    }else if (nrow(data) == 1){
      return(1)
    }else{
      drank <- data %>% 
        mutate(rank = dense_rank(desc(.data[[var]]))) %>%
        dplyr::arrange(rank) %>%
        mutate(cumsum = cumsum(.data[[var]]))
      
      A <- rank_area(drank$cumsum)
      
      r <- nrow(data) # species richness
      A_max <- 1 * (r - 1) # maximum if first species performs 100% (rectangle)
      A_min <- (r^2-1)/(2*r) # minimum if all sp contribute equally (trapezium)
      
      dd <- (A - A_min)/(A_max - A_min)
      
      return(dd)
    }
  }
  
  
  key <- parallel::mclapply(unique(con$transect_id), function(x){
    
    print(x/9118)
    t <- dplyr::filter(con, transect_id == x)
    
    result <- data.frame(
      transect_id = x,
      dd_Fn = ks_area(t, "Fn_p"),
      dd_Fp = ks_area(t, "Fp_p"),
      dd_Gc = ks_area(t, "Gc_p"),
      dd_I_herb = ks_area(t, "I_herb_p"),
      dd_I_pisc = ks_area(t, "I_pisc_p")
    )
    return(result)
  }, mc.cores = 40) %>%plyr::ldply()
  
  return(key)
}

# median contributions families 
get_cf <- function(contributions, herb_pisc){
  
  left_join(contributions, herb_pisc$contributions_herb_pisc) %>%
    group_by(bioregion, locality, sites, transect_id, Family) %>%
    summarize_at(.vars = vars(ends_with("_p")), sum, na.rm = TRUE) %>%
    group_by(bioregion, locality, sites, Family) %>%
    summarize_at(.vars = vars(ends_with("_p")), median, na.rm = TRUE) %>%
    group_by(bioregion, locality, Family) %>%
    summarize_at(.vars = vars(ends_with("_p")), median, na.rm = TRUE) %>%
    group_by(Family) %>%
    summarize_at(.vars = vars(ends_with("_p")), median, na.rm = TRUE) %>%
    ungroup() %>%
    group_by() %>%
    mutate_at(.vars = vars(ends_with("_p")), function(x){x/sum(x)}) %>%
    ungroup() 
  
}

get_importance <- function(contributions, herb_pisc){
  
  # if con > 1/N !!
  
  con <- left_join(contributions, herb_pisc$contributions_herb_pisc)
  
  parallel::mclapply(unique(con$transect_id), function(x){
    sub <- filter(con, transect_id == x)
    nherb <- sum(sub$I_herb_p>0)
    npisc <- sum(sub$I_pisc_p>0)
    sub <- sub %>%
      mutate(Gc_i = Gc_p > 1/nspec,
             Fn_i = Fn_p > 1/nspec,
             Fp_i = Fp_p > 1/nspec,
             I_pisc_i = I_pisc_p > 1/npisc,
             I_herb_i = I_herb_p > 1/nherb) %>%
      mutate(I_pisc_i = case_when(I_pisc_p == 0 ~ NA, TRUE ~ I_pisc_i),
             I_herb_i = case_when(I_herb_p == 0 ~ NA, TRUE ~ I_herb_i)) %>%
      dplyr::select(bioregion, transect_id, species, Gc_i, Fn_i, Fp_i, I_pisc_i, I_herb_i)
    return(sub)
  }, mc.cores = 30) %>% plyr::ldply()
  
}
  
get_fd <- function(sp_importance){
  # relative frequency when important 
  sp_importance %>%
    dplyr::select(bioregion, locality, sites, species, Gc_i, Fn_i, Fp_i, I_pisc_i, I_herb_i) %>%
    group_by(species) %>%
    mutate(occ = length(Gc_i)) %>%
    group_by(species, occ) %>%
    dplyr::summarise_if(is.logical, function(x){sum(x, na.rm = TRUE)/length(x[!is.na(x)])})
}
 
get_spi_vuln <- function(sp_importance, vulnerability){
  
  sp_importance %>%
    dplyr::select(species, Gc_i, Fn_i, Fp_i, I_pisc_i, I_herb_i) %>%
    unique() %>%
    group_by(species) %>%
    dplyr::summarise_if(is.logical, function(x){sum(x)>0}) %>%
    left_join(vulnerability)%>% 
    drop_na(vuln_climate, vuln_fi) %>%
    ungroup() %>%
    mutate(vulncat = case_when(
      (vuln_climate > median(vuln_climate) & vuln_fi > median(vuln_fi)) ~ "vuln_high",
      (vuln_climate > median(vuln_climate) & vuln_fi <= median(vuln_fi)) ~ "vuln_cli",
      (vuln_climate <= median(vuln_climate) & vuln_fi > median(vuln_fi)) ~ "vuln_fi",
      (vuln_climate <= median(vuln_climate) & vuln_fi <= median(vuln_fi)) ~ "vuln_low")) %>%
    dplyr::select(species, Family, vulncat, Gc = Gc_i, Fn = Fn_i, Fp = Fp_i, I_pisc = I_pisc_i, I_herb = I_herb_i) %>%
    pivot_longer(cols = c(Gc, Fn, Fp, I_pisc, I_herb), names_to = "name", values_to = "ep") %>%
    group_by(name) %>%
    mutate(n_spec_tot = sum(!is.na(ep))) %>%
    dplyr::group_by(vulncat, name, n_spec_tot) %>%
    dplyr::summarize(ep_n = sum(ep, na.rm = TRUE), 
                     n = (sum(ep, na.rm = TRUE) + sum(ep == FALSE, na.rm = TRUE))) %>%
    mutate(ep_prop = ep_n/n_spec_tot, # proportion of important species per process and category
           n_prop = n/n_spec_tot) %>% # proportion of species in vuln category
    mutate(perc = 100 * ep_prop/n_prop)
  
}

############# multivariate ################
stand <- function(x){
  max <- max(x, na.rm = TRUE)
  st <- sapply(x, FUN = function(i){
    st <- (i )/max
    return(st)
  })
}

geomean <- function(x, ...) {
  prod(x, ...)^(1/length(x))
}



##### functions for imputation #####

impute <- function(summary_transect_complete){
  fitimp1 <- brm(log(I_herb) ~ standard(mean) + standard(log(biomass_tot)) +
                   standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
                   standard(size_q3) + standard(troph_q3) + standard(imm_q1) +
                   standard(imm_q3) + standard(troph_q1) + standard(size_q1) +
                   (1|s|sites) + (1|p|locality),
                 data = summary_transect_complete[summary_transect_complete$I_herb>0,],
                 cores = 4,
                 backend = "cmdstanr",
                 threads = threading(10))
  
  imp1 <- fitted(fitimp1,
                 newdata = summary_transect_complete[summary_transect_complete$I_herb==0,],
                 allow_new_levels = TRUE)
  
  
  #### imputation ####
  fitimp2 <- brm(log(I_pisc) ~ standard(mean) + standard(log(biomass_tot)) +
                   standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
                   standard(size_q3) + standard(troph_q3) + standard(imm_q1) +
                   standard(imm_q3) + standard(troph_q1) + standard(size_q1) +
                   (1|s|sites) + (1|p|locality),
                 data = summary_transect_complete[summary_transect_complete$I_pisc>0,],
                 cores = 4,
                 backend = "cmdstanr",
                 threads = threading(10))
  
  imp2 <- fitted(fitimp2,
                 newdata = summary_transect_complete[summary_transect_complete$I_pisc==0,],
                 allow_new_levels = TRUE)
  
  # replace
  flux <- summary_transect_complete
  flux[flux$I_herb == 0, "I_herb"] <- exp(imp1[,1])
  flux[flux$I_pisc == 0, "I_pisc"] <- exp(imp2[,1])
  
  # get multifunction
  
  result <- flux %>%
    mutate(Fn_st = stand(Fn),
           Fp_st = stand(Fp),
           Gc_st = stand(Gc),
           I_herb_st = stand(I_herb),
           I_pisc_st = stand(I_pisc)) %>%
    rowwise() %>%
    mutate(multi = geomean(Fn_st, Fp_st, Gc_st, I_herb_st, I_pisc_st)) %>%
    ungroup() %>%
    mutate(multi = stand(multi))
  
  return(list(summary_transect_imp = result, 
              models = list(fitimp1, fitimp2)))
}

subset_data <- function(data) {
  data %>%
    dplyr::filter(I_herb > 0, I_pisc > 0) %>%
    mutate(Fn_st = stand(Fn),
           Fp_st = stand(Fp),
           Gc_st = stand(Gc),
           I_herb_st = stand(I_herb),
           I_pisc_st = stand(I_pisc)) %>%
    rowwise() %>%
    mutate(multi = geomean(Fn_st, Fp_st, Gc_st, I_herb_st, I_pisc_st)) %>%
    ungroup() %>%
    mutate(multi = stand(multi))
  
}


validate_model <- function(fit, save_plots = TRUE, name = "Title", pars = NULL){
  
  p1 <- 
    brms::pp_check(fit)  +
    labs(title = paste0(name))
  
  p2 <- 
    mcmc_plot(fit, type = "trace", pars = pars) +
    labs(title = paste0(name))
  
  p3 <- 
    mcmc_plot(fit, type = "rhat_hist", pars = pars) +
    labs(title = paste0(name))
  
  if (save_plots == TRUE) {
    ggsave(plot = p1, 
           filename = paste0("output/plots/__modeval_", name, "_ppcheck.pdf"), 
           width = 8, height = 6)
    ggsave(plot = p2, 
           filename = paste0("output/plots/__modeval_", name, "_trace.pdf"), 
           width = 8, height = 6)
    ggsave(plot = p3, 
           filename = paste0("output/plots/__modeval_", name, "_Rhat.pdf"), 
           width = 8, height = 6)
  }
  
  R2 <- bayes_R2(fit)
  return(list(pp_plot = p1, par_plot = p2, R2 = R2))
}


##### Models #####
##### 1) only site and location 
fit_mvfun_siteloc <- function(data){
  brm(mf(mvbind(log(Fn), log(Fp), log(Gc), log(I_herb), log(I_pisc))|mi() ~ 
        (1|s|sites) + (1|p|locality)) + 
        set_rescor(rescor = TRUE),
      data = data, cores = 4,
      backend = "cmdstanr",
      threads = threading(12))
  
}
fit_mf_siteloc <- function(data){
  brm(log(multi) ~ (1|s|sites) + (1|p|locality),
      data = data, cores = 4,
      backend = "cmdstanr",
      threads = threading(10))
}



##### 2) add biomass and sst
fit_mvfun_bm <- function(data){
  prior <- 
    set_prior(prior = "normal(1,1)", class = "b", coef = "logbiomass_tot", resp = "logFn") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "mean", resp = "logFn") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "logbiomass_tot", resp = "logFp") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "mean", resp = "logFp") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "logbiomass_tot", resp = "logGc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "mean", resp = "logGc") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "logbiomass_tot", resp = "logIherb") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "mean", resp = "logIherb") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "logbiomass_tot", resp = "logIpisc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "mean", resp = "logIpisc") 
  
  brm(bf(mvbind(log(Fn), log(Fp), log(Gc), log(I_herb), log(I_pisc))|mi() ~ 
        log(biomass_tot) + mean + (1|s|sites) + (1|p|locality)) + 
        set_rescor(rescor = TRUE),
      data = data, cores = 4,
      backend = "cmdstanr",
      threads = threading(10),
      prior = prior)
}

fit_mf_bm <- function(data){
  
  prior <- 
    set_prior(prior = "normal(1,1)", class = "b", coef = "logbiomass_tot", ) +
    set_prior(prior = "normal(0,1)", class = "b", coef = "mean")
  
  brm(log(multi) ~ log(biomass_tot) + mean + (1|s|sites) + (1|p|locality),
      data = data, cores = 4,
      backend = "cmdstanr",
      threads = threading(10),
      prior = prior)
}


fit_mvfun_com <- function(data){
  data <- data %>%
    dplyr::mutate(Fn_st = standard(log(Fn)),
                  Fp_st = standard(log(Fp)),
                  Gc_st = standard(log(Gc)),
                  I_herb_st = standard(log(I_herb)),
                  I_pisc_st = standard(log(I_pisc)))
  
  prior <- 
    set_prior(prior = "normal(1,1)", class = "b", coef = "standardlogbiomass_tot", resp = "Fnst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardmean", resp = "Fnst") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "standardlogbiomass_tot", resp = "Fpst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardmean", resp = "Fpst") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "standardlogbiomass_tot", resp = "Gcst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardmean", resp = "Gcst") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "standardlogbiomass_tot", resp = "Iherbst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardmean", resp = "Iherbst") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "standardlogbiomass_tot", resp = "Ipiscst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardmean", resp = "Ipiscst") +
    
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_m", resp = "Fnst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_q1", resp = "Fnst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_q3", resp = "Fnst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardnspec", resp = "Fnst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_m", resp = "Fnst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_q1", resp = "Fnst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_q3", resp = "Fnst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_m", resp = "Fnst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_q1", resp = "Fnst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_q3", resp = "Fnst") +
    
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_m", resp = "Fpst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_q1", resp = "Fpst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_q3", resp = "Fpst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardnspec", resp = "Fpst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_m", resp = "Fpst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_q1", resp = "Fpst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_q3", resp = "Fpst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_m", resp = "Fpst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_q1", resp = "Fpst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_q3", resp = "Fpst") +
    
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_m", resp = "Gcst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_q1", resp = "Gcst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_q3", resp = "Gcst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardnspec", resp = "Gcst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_m", resp = "Gcst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_q1", resp = "Gcst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_q3", resp = "Gcst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_m", resp = "Gcst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_q1", resp = "Gcst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_q3", resp = "Gcst") +
    
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_m", resp = "Iherbst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_q1", resp = "Iherbst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_q3", resp = "Iherbst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardnspec", resp = "Iherbst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_m", resp = "Iherbst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_q1", resp = "Iherbst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_q3", resp = "Iherbst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_m", resp = "Iherbst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_q1", resp = "Iherbst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_q3", resp = "Iherbst") +
    
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_m", resp = "Ipiscst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_q1", resp = "Ipiscst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardimm_q3", resp = "Ipiscst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardnspec", resp = "Ipiscst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_m", resp = "Ipiscst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_q1", resp = "Ipiscst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardsize_q3", resp = "Ipiscst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_m", resp = "Ipiscst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_q1", resp = "Ipiscst") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "standardtroph_q3", resp = "Ipiscst") 
    
  
  brm(bf(mvbind(Fn_st, Fp_st, Gc_st, I_herb_st, I_pisc_st)|mi() ~ 
        standard(mean) + standard(log(biomass_tot)) +
        standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
        standard(size_q3) + standard(troph_q3) + standard(imm_q1) +
        standard(imm_q3) + standard(troph_q1) + standard(size_q1) + 
        (1|s|sites) + (1|p|locality)) + 
        set_rescor(rescor = TRUE),
      data = data, cores = 4,
      backend = "cmdstanr",
      threads = threading(10),
      prior = prior)
}


fit_mvfun_com2 <- function(data){
  
  prior <- 
    set_prior(prior = "normal(1,1)", class = "b", coef = "logbiomass_tot", resp = "logFn") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "mean", resp = "logFn") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "logbiomass_tot", resp = "logFp") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "mean", resp = "logFp") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "logbiomass_tot", resp = "logGc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "mean", resp = "logGc") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "logbiomass_tot", resp = "logIherb") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "mean", resp = "logIherb") +
    set_prior(prior = "normal(1,1)", class = "b", coef = "logbiomass_tot", resp = "logIpisc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "mean", resp = "logIpisc") +
    
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_m", resp = "logFn") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_q1", resp = "logFn") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_q3", resp = "logFn") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "nspec", resp = "logFn") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_m", resp = "logFn") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_q1", resp = "logFn") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_q3", resp = "logFn") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_m", resp = "logFn") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_q1", resp = "logFn") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_q3", resp = "logFn") +
    
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_m", resp = "logFp") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_q1", resp = "logFp") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_q3", resp = "logFp") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "nspec", resp = "logFp") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_m", resp = "logFp") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_q1", resp = "logFp") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_q3", resp = "logFp") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_m", resp = "logFp") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_q1", resp = "logFp") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_q3", resp = "logFp") +
    
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_m", resp = "logGc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_q1", resp = "logGc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_q3", resp = "logGc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "nspec", resp = "logGc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_m", resp = "logGc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_q1", resp = "logGc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_q3", resp = "logGc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_m", resp = "logGc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_q1", resp = "logGc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_q3", resp = "logGc") +
    
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_m", resp = "logIherb") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_q1", resp = "logIherb") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_q3", resp = "logIherb") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "nspec", resp = "logIherb") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_m", resp = "logIherb") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_q1", resp = "logIherb") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_q3", resp = "logIherb") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_m", resp = "logIherb") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_q1", resp = "logIherb") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_q3", resp = "logIherb") +
    
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_m", resp = "logIpisc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_q1", resp = "logIpisc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "imm_q3", resp = "logIpisc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "nspec", resp = "logIpisc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_m", resp = "logIpisc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_q1", resp = "logIpisc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "size_q3", resp = "logIpisc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_m", resp = "logIpisc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_q1", resp = "logIpisc") +
    set_prior(prior = "normal(0,1)", class = "b", coef = "troph_q3", resp = "logIpisc") 
  

  brm(bf(mvbind(log(Fn), log(Fp), log(Gc), log(I_herb), log(I_pisc))|mi()  ~ 
        (mean) + (log(biomass_tot)) +
        (nspec) + (size_m) + (troph_m) + (imm_m) +
        (size_q3) + (troph_q3) + (imm_q1) +
        (imm_q3) + (troph_q1) + (size_q1) + 
        (1|s|sites) + (1|p|locality)) + 
        set_rescor(rescor = TRUE),
      data = data, cores = 4,
      prior = prior,
      backend = "cmdstanr",
      threads = threading(10))
}



fit_mf_com <- function(data){
  data <- data %>%
    dplyr::mutate(Fn_st = standard(log(Fn)),
                  Fp_st = standard(log(Fp)),
                  Gc_st = standard(log(Gc)),
                  I_herb_st = standard(log(I_herb)),
                  I_pisc_st = standard(log(I_pisc)),
                  multi_st = standard(log(multi)))
  brm(multi_st ~ 
        standard(mean) + standard(log(biomass_tot)) +
        standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
        standard(size_q3) + standard(troph_q3) + standard(imm_q1) +
        standard(imm_q3) + standard(troph_q1) + standard(size_q1) + 
        (1|s|sites) + (1|p|locality),
      data = data, cores = 4,
      backend = "cmdstanr",
      threads = threading(10))
}


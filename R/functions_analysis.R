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

##### calculate multifunctionality #####
add_multi <- function(summary_transect, herb_pisc){
  
  flux <- summary_transect
  hp <- herb_pisc$summary_herb_pisc
  sst <- read.csv("data/avSst.csv")
  flux <- left_join(flux, hp) %>% left_join(sst)
  
  # transform into centiles then divide by 100
  scores <- select(flux, transect_id, Fn, Fp, Gc, I_herb, I_pisc)  %>% 
    filter(I_herb>0, I_pisc>0) %>%
    mutate_at(2:6, function(x){normalize(log(x))}) %>%
    drop_na()
  # correlation for weighing
  m <-  cor(scores[2:6])
  m <- 1-m
  we <-  rowSums((m))
  we <- we/sum(we)
  

  multi <- scores %>%
    rowwise() %>% 
    mutate(multi = mean(c(Fn, Fp, Gc, I_herb, I_pisc)), 
           var = var(c( Fn, Fp, Gc, I_herb, I_pisc))) %>%
    group_by() %>%
    mutate(multi = multi*(1-(var/3000))/100) # 3000 = maximum theoretical variance
  
  flux <- flux %>% left_join(select(multi, transect_id, multi))
  
  return(flux)
  
}

flux <- left_join(flux, res)
flux$good <- ( flux$biomass_tot< 50 & flux$multi < quantile(flux$multi, 0.05, na.rm = TRUE))

flux <- flux %>%
  mutate(cat = case_when(
    biomass_tot < 50 & 
      multi < quantile(flux$multi, 0.1, na.rm = TRUE) ~ "bb",
    biomass_tot < 50 & 
      multi > quantile(flux$multi, 0.9, na.rm = TRUE) ~ "bg",
    biomass_tot > 150 & 
      multi < quantile(flux$multi, 0.1, na.rm = TRUE) ~ "gb",
    biomass_tot > 150 & 
      multi > quantile(flux$multi, 0.9, na.rm = TRUE) ~ "gg",
    TRUE~"other"
  ))

t <-flux %>%
  group_by(sites, locality) %>% summarise_if(is.numeric, mean)

t <- t %>%
  mutate(cat = case_when(
    biomass_tot < 50 & 
      multi < quantile(flux$multi, 0.1, na.rm = TRUE) ~ "bb",
    biomass_tot < 50 & 
      multi > quantile(flux$multi, 0.9, na.rm = TRUE) ~ "bg",
    biomass_tot > 150 & 
      multi < quantile(flux$multi, 0.1, na.rm = TRUE) ~ "gb",
    biomass_tot > 200 & 
      multi > quantile(flux$multi, 0.9, na.rm = TRUE) ~ "gg",
    TRUE~"other"
  ))
test <- table(flux$cat, flux$sites) %>% as.data.frame() %>%
  group_by(Var2) %>%
  mutate(tot = sum(Freq)) %>%
  mutate(Freq/tot) %>%
  filter(!(Var1 %in% c("bb", "bg") & Freq>1))

ggplot(flux) +
  geom_boxplot(aes(x = cat, y = imm_m))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = imm_q1))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = imm_q3))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = size_m))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = log(sizemax_m)))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = log(sizemax_q1)))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = log(sizemax_q3)))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = size_q1)) +
  scale_y_continuous(trans = "log")
ggplot(flux) +
  geom_violin(aes(x = cat, y = size_q1))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = size_q3))
ggplot(flux) +
  geom_violin(aes(x = cat, y = size_m))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = troph_m))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = troph_q3))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = troph_q1))

ggplot(flux) +
  geom_boxplot(aes(x = cat, y = multi))

ggplot(flux) +
  geom_boxplot(aes(x = cat, y = log(nspec)))
ggplot(flux) +
  geom_boxplot(aes(x = cat, y = mean))

ggplot(flux, aes(x = nspec, y = multi)) +
  geom_point() +
  geom_smooth() 

fitm <- brm(standard(multi) ~ standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
              standard(size_q3) + standard(troph_q3) + standard(imm_q1) + standard(imm_q3) + standard(troph_q1) + standard(size_q1) ,
            data = flux[!is.na(flux$multi),], chain = 3, cores = 3)

summary(fitm)

t <-flux %>%
  group_by(sites) %>% summarise_if(is.numeric, mean)

##### models with biomass, sst, locality, and site effects #######

run_procmodels <- function(summary_transect_complete){
  
  flux <- summary_transect_complete
  
  # run models
  fit_Fn <- brm(standard(log(Fn)) ~ (log(biomass_tot)) + (mean), data = flux, cores = 4)
  
  fit_Fp <- brm(standard(log(Fp)) ~ log(biomass_tot) + mean + (1|locality) + (1|sites), 
              data = flux, cores = 4)
  
  fit_Gc <- brm(standard(log(Gc)) ~ log(biomass_tot) + mean + (1|locality) + (1|sites), 
                data = flux[flux$Gc > 0, ], cores = 4)
  
  fit_I_herb <- brm(standard(log(I_herb)) ~ log(biomass_tot) + mean + (1|locality) + (1|sites), 
                    data = flux[flux$I_herb > 0,], cores = 4)  
  
  fit_I_pisc <- brm(standard(log(I_pisc)) ~ log(biomass_tot) + mean + (1|locality) + (1|sites), 
                    data = flux[flux$I_pisc > 0,], cores = 4)
  
  fit_multi <- brm(standard(multi) ~ log(biomass_tot) + mean + (1|locality) + (1|sites), 
                    data = flux[flux$multi > 0,], cores = 4)
  
  return(list(fit_Fn = fit_Fn, fit_Fp = fit_Fp, fit_Gc = fit_Gc, 
              fit_I_herb = fit_I_herb, fit_I_pisc = fit_I_pisc, fit_multi = fit_multi))
}

get_residuals <- function(summary_transect_complete, procmodels){
  
  flux <- summary_transect_complete
  
  res1 <- residuals(procmodels[[1]], re_formula = "standard(log(Fn)) ~ log(biomass_tot) + mean")
  res2 <- residuals(procmodels[[2]], re_formula = "standard(log(Fp)) ~ log(biomass_tot) + mean")
  res3 <- residuals(procmodels[[3]], re_formula = "standard(log(Gc)) ~ log(biomass_tot) + mean")
  res4 <- residuals(procmodels[[4]], re_formula = "standard(log(I_herb)) ~ log(biomass_tot) + mean")
  res5 <- residuals(procmodels[[5]], re_formula = "standard(log(I_pisc)) ~ log(biomass_tot) + mean")
  res6 <- residuals(procmodels[[6]], re_formula = "standard(multi) ~ log(biomass_tot) + mean")
  
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
  
  res_m <- data.frame(
    transect_id = flux[!is.na(flux$multi), ]$transect_id,
    multi_r = res6[,1]
  )
  
  res <- res %>% left_join(res_h) %>% left_join(res_p) %>% left_join(res_m)
  
  #add multifunctionality
  scores <- select(res, transect_id, Fn = Fn_r, Fp = Fp_r, Gc = Gc_r, I_herb = I_herb_r, I_pisc = I_pisc_r)  %>% 
    drop_na() %>%
    mutate_at(2:6, function(x){normalize(x)}) 
  
  multi <- scores %>%
    rowwise() %>% 
    mutate(multi = mean(c(Fn, Fp, Gc, I_herb, I_pisc)), 
           var = var(c( Fn, Fp, Gc, I_herb, I_pisc))) %>%
    group_by() %>%
    mutate(multi2 = multi*(1-(var/3000))/100) # 3000 = maximum theoretical variance
  
  res <- left_join(res, select(multi, transect_id, multi, multi2))
  
  return(res)
}

get_location_effect <- function(procmodels, summary_transect){
  result <- lapply(procmodels, function(x){
      x %>%
      spread_draws(r_locality[locality, Intercept]) %>%
      median_qi() %>%
      select(locality, effect = r_locality)
  }) %>% purrr::reduce(dplyr::left_join, by = "locality")
  
  colnames(result) <- c("locality", "r_loc_Fn", "r_loc_Fp", "r_loc_Gc", 
                        "r_loc_I_herb", "r_loc_I_pisc", "r_loc_multi")
  
  # add coordinates
  coord <- select(summary_transect, bioregion, locality, lat, lon) %>%
    unique() %>%
    group_by(bioregion, locality) %>%
    summarize_all(mean)
  
  result <- coord %>% left_join(result) %>% ungroup()
  
  #add multifunctionality
  scores <- select(result, locality, Fn = r_loc_Fn, Fp = r_loc_Fp, Gc = r_loc_Gc, I_herb = r_loc_I_herb, I_pisc = r_loc_I_pisc)  %>% 
    drop_na() %>%
    mutate_at(2:6, function(x){normalize(x)}) 
  
  multi <- scores %>%
    rowwise() %>% 
    mutate(multi = mean(c(Fn, Fp, Gc, I_herb, I_pisc)), 
           var = var(c( Fn, Fp, Gc, I_herb, I_pisc))) %>%
    group_by() %>%
    mutate(multi = multi*(1-(var/3000))/100) # 3000 = maximum theoretical variance
  
  result <- left_join(result, select(multi, locality, multi))
  
  return(result)
    
}


# summary(fitN)
# 
# res_Fn <- 
#   fitN %>%
#   spread_draws(r_sites[sites, Intercept]) %>%
#   median_qi() %>%
#   select(sites, Fn_r = r_sites)
# 
# res_Fp <- 
#   fitP %>%
#   spread_draws(r_sites[sites, Intercept]) %>%
#   median_qi() %>%
#   select(sites, Fp_r = r_sites)
# 
# res <- left_join(res_Fn, res_Fp)
# 
# ggplot(res) +
# geom_point(aes(x = Fn_r, y = Fp_r))
# 
# summary(cor(res$Fn_r, res$Fp_r))




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


plot_con_vuln_all <- function(contributions_sp_loc_occ, vulnerability){
  
  con <- left_join(contributions_sp_loc_occ, vulnerability) %>% drop_na(vuln_fi, vuln_climate)
  
  plot_con_vuln <- function(var, ylab){
    
    data <- con
    data <- data %>% filter(!.data[[var]] == 0)
    probs <- c(0.05, 0.25, 0.75, 0.95)
    quantiles <- quantile(data[[var]], prob = probs, na.rm = TRUE)
    data$quant <- factor(findInterval(data[[var]], quantiles))
    
    data$density1 <- fields::interp.surface(
      MASS::kde2d(data[["vuln_fi"]], data[[var]]), as.data.frame(data[,c("vuln_fi", var)]))
    data$density2 <- fields::interp.surface(
      MASS::kde2d(data[["vuln_climate"]], data[[var]]), as.data.frame(data[,c("vuln_climate", var)]))
    
    data <- filter(data, .data[[var]] > 0.00001)
    
    p1 <-
      ggplot(data) +
      geom_point(aes(x = 100 * vuln_fi, y = 100 * .data[[var]], size = rel_occ, color = quant, alpha = 1/density1)) +
      #geom_density2d(aes(x = 100 * vuln_fi, y = 100 * .data[[var]]), color = "grey") +   
      geom_hline(aes(yintercept = 100 * (median(.data[[var]]))), linetype = 2, color = "grey", size = 1) +
      geom_vline(aes(xintercept =  100 *  median(con$vuln_fi, na.rm = TRUE) ), linetype = 2,  color = "grey", size = 1) +
      scale_x_continuous(trans = "log", breaks = c(15, 25, 50, 100)) +
      scale_y_continuous(trans = "log", breaks = c(0.0009, 0.04, 2, 100)) +
      scale_alpha(range = c(.25, .75) ) +
      xlab("Vulnerability to fishing") +
      ylab(ylab) +
      labs(color = "Quantile", size = "Relative occurrence") +
      scale_color_fish_d("Hypsypops_rubicundus", begin = 0 ,
                         labels=c("[0,0.05[", "[0.05,0.25[", "[0.25,0.75[",  "[0.75, 0.95[", "[0.95,1[")) +
      theme_bw()  + theme(legend.position = "none")
    
    p1
    p2 <-
      ggplot(data) +
      geom_point(aes(x = 100 * vuln_climate, y = 100 * .data[[var]], size = rel_occ, color = quant, alpha = 1/density2)) +
      #geom_smooth(aes(x = 100 * vuln_climate, y = 100 * (Fp_p)), method = "lm", se = FALSE, color = "black") +
      geom_hline(aes(yintercept = 100 * (median(.data[[var]]))), linetype = 2, color = "grey", size = 1) +
      geom_vline(aes(xintercept =  100 *  median(con$vuln_climate, na.rm = TRUE) ), linetype = 2,  color = "grey", size = 1) +
      scale_x_continuous(trans = "log", breaks = c(15, 25, 50, 100)) +
      scale_y_continuous(trans = "log", breaks = c(0.0009, 0.04, 2, 100)) +
      scale_alpha(range = c(.25, .75) ) +
      xlab("Vulnerability to climate change") +
      ylab(ylab) +
      labs(color = "Quantile", size = "Relative occurrence") +
      scale_color_fish_d("Hypsypops_rubicundus", begin = 0 ,
                         labels=c("[0,0.05[", "[0.05,0.25[", "[0.25,0.75[",  "[0.75, 0.95[", "[0.95,1[")) +
      theme_bw()  + theme(legend.position = "none")
    
    plot <- cowplot::plot_grid(p1, p2, nrow = 1)
    return(plot)
    
  }
  
  pa1 <- plot_con_vuln("Fn_p_m", "Contribution Fn")
  pa2 <- plot_con_vuln("Fp_p_m", "Contribution Fp")
  pa3 <- plot_con_vuln("Gc_p_m", "Contribution Gc")
  pa4 <- plot_con_vuln("I_herb_p_m", "Contribution I_herb")
  pa5 <- plot_con_vuln("I_pisc_p_m", "Contribution I_pisc")
  
  layout <- "
AA
BB
CC
DD
EE
"
  
  pa1 + pa2 + pa3 + pa4 + pa5 + plot_layout(nrow = 5)
  ggsave("output/plots/con_vuln.pdf", width = 8, height = 16)
  
}



##### dominance ######

get_dominance <- function(contributions, herb_pisc){
  
  # combine data
  con <- left_join(contributions, herb_pisc$contributions_herb_pisc)
  
  # dominance function
  dominance <- function(data, var){
    test <- data %>% mutate(rank = rank(- .data[[var]], ties.method = c("first"))) %>%
      filter(.data[[var]] > 0) %>%
      dplyr::arrange(rank) %>%
      mutate(cumsum = cumsum(.data[[var]]))
    
    if (nrow(test) < 4){
      return(list(slope = NA, dom = NA))
    }else{
    
    dom <- (length(test$cumsum[test$cumsum > 0.5])/nrow(test))
    
    fit <- lm(log(test[[var]])~log(test[["rank"]]))
    slope <-  fit$coefficients[2]
    return(list(slope = slope, dom = dom))
    }
  }
  
  dom <- parallel::mclapply(unique(con$transect_id), function(x){
    
    print(x)
    t <- dplyr::filter(con, transect_id == x)
    
    d_Fn <- dominance(t, "Fn_p")
    d_Fp <- dominance(t, "Fp_p")
    d_Gc <- dominance(t, "Gc_p")
    d_I_herb <- dominance(t, "I_herb_p")
    d_I_pisc <- dominance(t, "I_pisc_p")
    
    result <- data.frame(
      transect_id = x,
      Fn_d = d_Fn[[2]],
      Fn_dslope = d_Fn[[1]],
      Fp_d = d_Fp[[2]],
      Fp_dslope = d_Fp[[1]],
      Gc_d = d_Gc[[2]],
      Gc_dslope = d_Gc[[1]],
      I_herb_d = d_I_herb[[2]],
      I_herb_dslope = d_I_herb[[1]],
      I_pisc_d = d_I_pisc[[2]],
      I_pisc_dslope = d_I_pisc[[1]]
    )
    return(result)
  }, mc.cores = 40) %>%plyr::ldply()
  
  return(dom)
}

### calculate area of ranked function contribution !! data has to be ordered
# y = cumsum!!
rank_area <- function(y) {
  n <- length(y)
  a <- y[2:n]
  b <- y[1:(n-1)]
  area <- sum(a + b)/2
  return(area)
}


get_keystoneness <- function(contributions, herb_pisc){
  
  # combine data
  con <- left_join(contributions, herb_pisc$contributions_herb_pisc)
  
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
      test <- data %>% mutate(rank = rank(- .data[[var]], ties.method = c("first"))) %>%
        dplyr::arrange(rank) %>%
        mutate(cumsum = cumsum(.data[[var]]))
      
      A <- rank_area(test$cumsum)
      
      r <- nrow(data) # species richness
      A_max <- 1 * (r - 1) # maximum if first species performs 100% (rectangle)
      A_min <- (r^2-1)/(2*r) # minimum if all sp contribute equally (trapezium)
      
      ks <- (A - A_min)/(A_max - A_min)
      
      return(ks)
    }
  }
  
  
  key <- parallel::mclapply(unique(con$transect_id), function(x){
    
    print(x)
    t <- dplyr::filter(con, transect_id == x)
    
    result <- data.frame(
      transect_id = x,
      ks_Fn = ks_area(t, "Fn_p"),
      ks_Fp = ks_area(t, "Fp_p"),
      ks_Gc = ks_area(t, "Gc_p"),
      ks_I_herb = ks_area(t, "I_herb_p"),
      ks_I_pisc = ks_area(t, "I_pisc_p")
    )
    return(result)
  }, mc.cores = 40) %>%plyr::ldply()
  
  return(key)
}







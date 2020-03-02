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







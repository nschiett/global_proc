
make_fig1 <- function(location_effect){
  
  library("rnaturalearth")
  library("rnaturalearthdata")
  library(rgeos)
  library(sf)
  world <- ne_countries(scale = "medium", returnclass = "sf")
  class(world)
  
  
  
  theme_worldmap <- function(){
    theme_bw()+
      theme(legend.position = "none", 
            panel.grid = element_blank(),
            panel.border = element_rect(fill = NA, color = "grey", size = 1),
            axis.title = element_blank(),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 6),
            plot.title = element_text(size = 6, hjust = 0, face = "bold"), 
            plot.margin = unit(c(0.002,0.002,0.002,0.002), units = , "cm"),
            panel.spacing = unit(0, "cm")
      ) 
  }
  
  col <- fish(n = 5, option = "Callanthias_australis")
  
  get_scale <- function(f){
    quant <- quantile(f, c(0,0.25,0.75,1), na.rm = TRUE)
    quant[1] <- quant[1] - 0.5   # adjust lower bound
    factor(cut(f, breaks = quant), labels = c( "low", "medium", "high"))
  }
  
  get_scores <- function(proc){
    probs <- seq(0,1,0.01)
    quantiles <- quantile(proc, prob = probs, na.rm = TRUE)
    quant <- as.numeric(factor(findInterval(proc, quantiles)))/100
    return(quant)
  }
  

  getcol <- colorRampPalette(c("white", col[1]))
  pal2 <- getcol(n = 3)
  getcol <- colorRampPalette(c("black","white"))
  pal1 <- getcol(n = 3)
  pal <- c(pal1[1], pal2[2:3])
  pal <- c("black", "grey50", col[1])
  
  
  
  a <-
    ggplot(location_effect) + 
    geom_sf(data = world, color = NA, fill = "lightgrey") +
    geom_point(aes(x = lon, y = lat, color = get_scale(r_loc_Fn), 
                   size = get_scores(r_loc_Fn)),  alpha = 0.8,
               position = position_jitter(0,0), shape = 16) +
    geom_point(aes(x = lon, y = lat),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(location_effect, r_loc_Fn > quantile(location_effect$r_loc_Fn, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    geom_text(aes(x = -175, y = 30, label = "N excretion (g N/m²day)"), size = 3, hjust = 0) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.000,0.000,0.000,0.000), units = , "cm"))
  a
  
  
  getcol <- colorRampPalette(c("white", col[2]))
  pal2 <- getcol(n = 3)
  getcol <- colorRampPalette(c("black","white"))
  pal1 <- getcol(n = 3)
  pal <- c("black", "grey50", col[2])
  
  
  b <-
    ggplot(location_effect) + 
    geom_sf(data = world, color = NA, fill = "lightgrey") +
    geom_point(aes(x = lon, y = lat, color = get_scale(r_loc_Fp), 
                   size = get_scores(r_loc_Fp)),  alpha = 0.8,
               position = position_jitter(0,0), shape = 16) +
    geom_point(aes(x = lon, y = lat),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(location_effect, r_loc_Fp > quantile(location_effect$r_loc_Fp, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    geom_text(aes(x = -175, y = 30, label = "P excretion (g P/m²day)"), size = 3, hjust = 0) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.0001,0.0001,0.0001,0.0001), units = , "cm"))
  b
  
  
  
  getcol <- colorRampPalette(c("white", col[3]))
  pal2 <- getcol(n = 3)
  getcol <- colorRampPalette(c("black","white"))
  pal1 <- getcol(n = 3)
  pal <- c("black", "grey50", col[3])
  
  c <-
    ggplot(location_effect) + 
    geom_sf(data = world, color = NA, fill = "lightgrey") +
    geom_point(aes(x = lon, y = lat, color = get_scale(r_loc_Gc), 
                   size = get_scores(r_loc_Gc)),  alpha = 0.8,
               position = position_jitter(0,0), shape = 16) +
    geom_point(aes(x = lon, y = lat),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(location_effect, r_loc_Gc > quantile(location_effect$r_loc_Gc, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    geom_text(aes(x = -175, y = 30, label = "Production (g C/m²day)"), size = 3, hjust = 0) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.00,0.00,0.00,0.00), units = , "cm"))
  c
  
  
  
  
  
  
  getcol <- colorRampPalette(c("lightgrey", col[4]))
  pal2 <- getcol(n = 3)
  getcol <- colorRampPalette(c("black","lightgrey"))
  pal1 <- getcol(n = 3)
  pal <- c(pal1[1:2], pal2[2:3])
  pal <- c("black", "grey50", col[4])
  
  
  d <-
    ggplot(location_effect) + 
    geom_sf(data = world, color = NA, fill = "lightgrey") +
    geom_point(aes(x = lon, y = lat, color = get_scale(r_loc_I_herb), 
                   size = get_scores(r_loc_I_herb)),  alpha = 0.8,
               position = position_jitter(0,0), shape = 16) +
    geom_point(aes(x = lon, y = lat),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(location_effect, r_loc_I_herb > quantile(location_effect$r_loc_I_herb, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    geom_text(aes(x = -175, y = 30, label = "Herbivory (g C/m²day)"), size = 3, hjust = 0) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.00,0.00,0.00,0.00), units = , "cm"))
  d
  
  
  
  getcol <- colorRampPalette(c("lightgrey", col[5]))
  pal2 <- getcol(n = 3)
  getcol <- colorRampPalette(c("black","lightgrey"))
  pal1 <- getcol(n = 3)
  pal <- c(pal1[1:2], pal2[2:3])
  
  pal <- c("black", "grey50", col[5])
  
  e <-
    ggplot(location_effect) + 
    geom_sf(data = world, color = NA, fill = "lightgrey") +
    geom_point(aes(x = lon, y = lat, color = get_scale(r_loc_I_pisc), 
                   size = get_scores(r_loc_I_pisc)),  alpha = 0.8,
               position = position_jitter(0,0), shape = 16) +
    geom_point(aes(x = lon, y = lat),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(location_effect, r_loc_I_pisc > quantile(location_effect$r_loc_I_pisc, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "",
                       drop = TRUE, na.translate = F) +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    geom_text(aes(x = -175, y = 30, label = "Piscivory (g C/m²day)"), size = 3, hjust = 0) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
    theme_worldmap()+ 
    theme(legend.position = "none", 
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.0,0.00,0.00,0.00), units = , "cm"))
  e
  
  
  multimap <-
    a + b + c + d + e +
    plot_layout(ncol = 1)
  
  #multimap
  ggsave("output/plots/fig1_multimap.png", multimap, width = 8, height = 8)
}


make_fig2 <- function(commodels){
  
  standard <- function(x){
    m <- mean(x, na.rm = TRUE)
    sd <- sd(x, na.rm = TRUE)
    st <- sapply(x, FUN = function(i){
      st <- (i - m)/sd
      return(st)
    })
    return(st)
  }
  
  
  
  ## function to extract slope
  
  extract_b <- function(x, name){
    draws <- extract_draws(x)$dpars$mu$fe$b
    means <- apply(draws, 2, mean)
    lq <- apply(draws, 2, quantile, 0.025)
    uq <- apply(draws, 2, quantile, 0.975)
    result <- 
      data.frame(
        model = name,
        variable = names(means),
        mean = means,
        lq = lq,
        uq = uq
      )
  }
  
  # extract all slopes for each model
  result <- bind_rows(
    extract_b(commodels[[1]], "Fn"),
    extract_b(commodels[[2]], "Fp"),
    extract_b(commodels[[3]], "Gc"),
    extract_b(commodels[[4]], "I_herb"),
    extract_b(commodels[[5]], "I_pisc")) 
  
  # Indicate which ones overlap with 0
  result <- mutate(result, diff = (uq > 0 & lq < 0)) 
  
  # Specify importance of variables per slope values
  importance <- group_by(result, variable) %>%
    summarise(tot = max(abs(mean))) 
  
  # filter out variables that we don't need 
  result <- result %>%
    dplyr::filter(!variable == "b_Intercept") %>%
    left_join(importance) %>% 
    mutate(variable = fct_reorder(variable, tot)) %>%
    filter(diff == FALSE) %>%
    filter(!variable == "b_standardlogbiomass") %>%
    filter(!variable == "b_standardmean") 
    
  
  # plot
  slopes <-
    ggplot(result) +
    geom_vline(xintercept = seq(1, 10 ,1) + .5, color = "lightgrey", 
               size = 0.5, linetype = 1, alpha = 0.7) +
    geom_hline(yintercept = 0, size = 1, color = "black") +
    geom_linerange(aes(x = variable, ymin = 0, ymax = mean, color = model),
                   position = position_dodge(.7), size = 0.5, linetype = 2) + 
    geom_linerange(aes(x = variable, ymin = lq, ymax = uq, color = model),
                   position = position_dodge(.7), size = 1, linetype = 1) + 
    geom_point(aes(x = variable, y = mean, color = model), position = position_dodge(.7), size = 3) +
    coord_flip() +
    theme_minimal() +
    scale_color_fish_d(option = "Callanthias_australis", 
                       labels = c("N excretion", "P excretion", "Production", "Herbivory", "Piscivory"),
                       name = "")  + 
    labs(y = "Effect size", x = "") +
    scale_x_discrete(labels = c(  "Immaturity (2.5%)", "Size (median)",  "Immaturity (97.5%)", 
                                  "Size (2.5%)", "Richness",  "Immaturity (median)",
                                  "Size (97.5%)", "Trophic level (2.5%)",  "Trophic level (median)",
                                  "Trophic level (97.5%)")) +
   
    theme(legend.position = "bottom", axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12), legend.text = element_text(size = 12), 
          axis.line.x = element_line(), axis.ticks.x = element_line(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank()
    )   
  slopes
  ggsave("output/plots/fig2_com_slopes.png", slopes, width = 8, height = 9)
  
}

make_fig3 <- function(contributions, herb_pisc, degree_dominance, freq_dominance){
  
  # a
  con_fam <- left_join(contributions, herb_pisc$contributions_herb_pisc) %>%
    filter(!is.na(Family)) %>%
    mutate(I_herb_p = case_when(I_herb_p == 0 ~ NA_real_,
                                TRUE ~ I_herb_p),
           I_pisc_p = case_when(I_pisc_p == 0 ~ NA_real_,
                                TRUE ~ I_pisc_p)) %>%
    mutate(Fn_p = Fn_p - biomass_p,
           Fp_p = Fp_p - biomass_p,
           Gc_p = Gc_p - biomass_p,
           I_herb_p = I_herb_p - biomass_herb_p,
           I_pisc_p = I_pisc_p - biomass_pisc_p) %>%
    group_by(bioregion, locality, sites, transect_id, Family) %>%
    summarize_at(.vars = vars(ends_with("_p")), sum, na.rm = TRUE)  %>%
    group_by(bioregion, locality, sites, Family) %>%
    summarize_at(.vars = vars(ends_with("_p")), median, na.rm = TRUE)  %>%
    group_by(bioregion, locality, Family) %>%
    summarize_at(.vars = vars(ends_with("_p")), mean, na.rm = TRUE)  %>%
    ungroup() %>%
    pivot_longer(c(Fn_p, Fp_p, Gc_p, I_herb_p, I_pisc_p)) %>%
    mutate(value = case_when(value == 0 ~ NA_real_,
                             TRUE ~ value)) %>%
    group_by(name, Family) %>%
    summarize(biomass_p = median(biomass_p),
              value_m = median(value, na.rm = TRUE), 
              value_lq = quantile(value, 0.25, na.rm = TRUE),
              value_uq = quantile(value, 0.75, na.rm = TRUE)) %>%
    arrange(-biomass_p) %>%
    filter(biomass_p > 0.025) 
  
  con_fam <- con_fam %>%
    arrange(biomass_p) %>%
    filter(biomass_p > 0.025) %>%
    mutate(Family = factor(Family, unique(Family)))

  
  a <- 
    ggplot(con_fam) +
    geom_point(aes(xmin = value_lq, xmax = value_uq, x = value_m, 
                   y = Family, color = name),
               position = ggstance::position_dodgev(height = -0.7)) +
    geom_hline(yintercept = seq(1, 12 ,1) + .5, color = "lightgrey", 
               size = 0.5, linetype = 1, alpha = 0.7) +
    geom_vline(xintercept = 0) +
    geom_errorbarh(aes(xmin = value_lq, xmax = value_uq, x = value_m, 
                       y = Family, color = name),
                   position = ggstance::position_dodgev(height = -0.7), height = 0) +

    scale_color_fish_d(option = "Callanthias_australis",
                       labels = c("N excretion", "P excretion", "Production", "Herbivory", "Piscivory"),
                       name = "") +
    
    add_fishape(family = "Acanthuridae", option = NA, alpha = 0.7,
                xmin = 0.075, xmax = 0.11, ymin = 11.75, ymax = 12.25) +
    add_fishape(family = "Labridae", option = NA, alpha = 0.7,
                xmax = -0.07, xmin = -0.11, ymin = 10.75, ymax = 11.25) +
    add_fishape(family = "Pomacentridae", option = NA,  alpha = 0.7,
                xmin = 0.075, xmax = 0.11, ymin = 9.75, ymax = 10.25) +
    add_fishape(family = "Haemulidae", option = NA,  alpha = 0.7,
                xmax = -0.075, xmin = -0.11, ymin = 8.75, ymax = 9.25) +
    add_fishape(family = "Lutjanidae", option = NA,  alpha = 0.7,
                xmin = 0.075, xmax = 0.11, ymin = 7.75, ymax = 8.25) +
    add_fishape(family = "Serranidae", option = "Cephalopholis_argus",  alpha = 0.7,
                xmax = -0.075, xmin = -0.11, ymin = 6.75, ymax = 7.25) +
    add_fishape(family = "Balistidae", option = NA,  alpha = 0.7,
                xmin = 0.075, xmax = 0.11, ymin = 5.75, ymax = 6.25) +
    add_fishape(family = "Mullidae", option = NA,  alpha = 0.7,
                xmax = -0.075, xmin = -0.11, ymin = 4.75, ymax = 5.25) +
    add_fishape(family = "Kyphosidae", option = NA,  alpha = 0.7,
                xmin = 0.075, xmax = 0.11, ymin = 3.75, ymax = 4.25) +
    add_fishape(family = "Holocentridae", option = NA,  alpha = 0.7,
                xmax = -0.075, xmin = -0.11, ymin = 2.75, ymax = 3.25) +
    add_fishape(family = "Chaetodontidae", option = NA,  alpha = 0.7,
                xmin = 0.075, xmax = 0.11, ymin = 1.75, ymax = 2.25) +
    add_fishape(family = "Lethrinidae", option = NA,  alpha = 0.7,
                xmax = -0.075, xmin = -0.11, ymin = 0.75, ymax = 1.25) +
    
    xlim(c(-0.11, 0.11)) +

    labs(x = "contribution function - contribution biomass", y = "") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(), panel.border = element_blank(), 
          legend.position = "right",
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line(), axis.ticks.x = element_line(),
          axis.ticks.y = element_blank(), axis.text = element_text(color = "black"))
  
  # b
  
  dd_long <- degree_dominance %>%
    pivot_longer(c(dd_Fn, dd_Fp, dd_Gc, dd_I_herb, dd_I_pisc))
  
  b <- 
    ggplot(dd_long, alpha = 0.8) +
    geom_eyeh(aes(x = value, y = reorder(name, desc(name)), fill = name),
              alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
    scale_fill_fish_d(option = "Callanthias_australis",
                      labels = c("N excretion", "P excretion", "Production", "Herbivory", "Piscivory"),
                      name = "Function") +
    scale_x_continuous(limits = c(0,1)) +
    labs(x = "Degree of dominance (communities)") +
    theme_bw() +
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
          legend.position = "none", 
          panel.border = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line()
    )
  
  # c
  spi_sp_long <- freq_dominance %>%
    pivot_longer(c(Fn_i, Fp_i, Gc_i, I_herb_i, I_pisc_i))
  
  c <-
    ggplot(spi_sp_long, alpha = 0.8) +
    geom_eyeh(aes(x = value, y = reorder(name, desc(name)), fill = name),
              alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
    scale_fill_fish_d(option = "Callanthias_australis",
                      labels = c("N excretion", "P excretion", "Production", "Herbivory", "Piscivory"),
                      name = "Function") +
    scale_x_continuous(limits = c(0,1)) +
    labs(x = "Frequency of dominance (species)") +
    theme_bw() +
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
          legend.position = "none", 
          panel.border = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line()
    )
  
  # combine plots
  layout = 
    "AB
     AB
     AC
     AC"
  
  plot <- 
  a + b + c  +
    plot_layout(design = layout, guides = "collect") +
    plot_annotation(tag_levels = 'a', tag_suffix = ")") #& 
    #theme(plot.tag = element_text(size = 12)) 
  
  ggsave("output/plots/fig3_contributions.png", plot,  width = 9, height = 6)
  
}

make_fig4 <- function(spi_vuln){
  
  # barplots for each function with labels at angles
  
  bp1 <-
    ggplot(filter(spi_vuln, vulncat == "vuln_high") %>% drop_na()) +
    geom_bar(aes(x = name, y = n_prop, fill = name), stat = "identity", alpha = 0.5) +
    geom_bar(aes(x = name, y = ep_prop, fill = name), stat = "identity") +
    geom_text(aes(x = name, y = n_prop + 0.1, label = paste0(round(perc), "%")), angle = 45) +
    scale_fill_fish_d(option = "Callanthias_australis",
                      labels = c("N excretion", "P excretion",
                                 "Production", "Herbivory", "Piscivory"),
                      name = "Function") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_void() +
    theme(legend.position = "none",
          text = element_text(angle = 45), plot.margin = unit(c(0,0,0,0), "cm")
    )
  ggsave("output/plots/fig4_glob_vuln_high.png", bp1, width = 4, height = 4)
  
  bp2 <-
    ggplot(filter(spi_vuln, vulncat == "vuln_cli") %>% drop_na()) +
    geom_bar(aes(x = name, y = n_prop, fill = name), stat = "identity", alpha = 0.5) +
    geom_bar(aes(x = name, y = ep_prop, fill = name), stat = "identity") +
    geom_text(aes(x = name, y = n_prop + 0.1, label = paste0(round(perc), "%")), angle = 135) +
    scale_fill_fish_d(option = "Callanthias_australis",
                      labels = c("N excretion", "P excretion",
                                 "Production", "Herbivory", "Piscivory"),
                      name = "Function") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_void() +
    theme(legend.position = "none", plot.margin = unit(c(0,0,0,0), "cm")
    )
  
  ggsave("output/plots/fig4_glob_vuln_cc.png", bp2, width = 4, height = 4)
  
  bp3 <-
    ggplot(filter(spi_vuln, vulncat == "vuln_low")) +
    geom_bar(aes(x = name, y = n_prop, fill = name), stat = "identity", alpha = 0.5) +
    geom_bar(aes(x = name, y = ep_prop, fill = name), stat = "identity") +
    geom_text(aes(x = name, y = n_prop + 0.1, label = paste0(round(perc), "%")), angle = -135) +
    scale_fill_fish_d(option = "Callanthias_australis",
                      labels = c("N excretion", "P excretion",
                                 "Production", "Herbivory", "Piscivory"),
                      name = "Function") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_void() +
    theme(legend.position = "none",
          text = element_text(angle = 45), plot.margin = unit(c(0,0,0,0), "cm")
    )
  
  bp3
  ggsave("output/plots/fig4_glob_vuln_low.png", bp3, width = 4, height = 4)
  
  
  bp4 <-
    ggplot(filter(spi_vuln,  vulncat == "vuln_fi")) +
    geom_bar(aes(x = name, y = n_prop, fill = name), stat = "identity", alpha = 0.5) +
    geom_bar(aes(x = name, y = ep_prop, fill = name), stat = "identity") +
    geom_text(aes(x = name, y = n_prop + 0.1, label = paste0(round(perc), "%")), angle = -45) +
    scale_fill_fish_d(option = "Callanthias_australis",
                      labels = c("N excretion", "P excretion",
                                 "Production", "Herbivory", "Piscivory"),
                      name = "Function") +
    theme_void() +
    theme(legend.position = "none",
          text = element_text(angle = 45), plot.margin = unit(c(0,0,0,0), "cm")
    )
  
  bp4
  ggsave("output/plots/fig4_glob_vuln_fi.png", bp4, width = 4, height = 4)
  
  ## center
  
  middle <- ggplot() +
    geom_polygon(aes(x = c(0, 0, 0, 1), y = c(0, 1,1,0)), 
                 color = "black", fill = "grey60") +
    geom_polygon(aes(x = c(0, 0, 0, -1), y = c(0, 1,1,0)),  
                 color = "black", fill = "grey80") +
    geom_polygon(aes(x = c(0, 0, 0, 1), y = c(0, -1,-1,0)), 
                 color = "black", fill = "grey80") +
    geom_polygon(aes(x = c(0, 0, 0, -1), y = c(0, -1,-1,0)), 
                 color = "black", fill = "grey90") +
    scale_x_continuous(limits = c(-2,2), expand = c(0,0)) +
    scale_y_continuous(limits = c(-2,2), expand = c(0,0)) +
    geom_text(aes(x = -0.33, y = -0.33, label = "Low \n vulnerability \n to both"), angle = -45, size = 2.5) +
    geom_text(aes(x = 0.33, y = 0.33, label = "High \n vulnerability \n to both"), angle = -45, size = 2.5) +
    geom_text(aes(x = -0.33, y = 0.33, label = "High \n vulnerability \n to fishing"), angle = 45, size = 2.5) +
    geom_text(aes(x = 0.33, y = -0.33, label = "High \n vulnerability \n to climate change"), angle = 45, size = 2.5) +
    coord_fixed() +
    theme_void()

  ggsave("output/plots/fig4_middle.png", middle, width = 4, height = 4)
  
}

make_annex_fig1 <- function(summary_transect_complete, bmmodels){

  nd <- summary_transect_complete %>%
    select(mean) %>%
    unique() %>%
    mutate(biomass_tot = 100) 
  
  ndp <- nd %>%
    mutate(Fn_ref = exp(fitted(bmmodels[[1]], nd)[,1])) %>%
    mutate(Fp_ref = exp(fitted(bmmodels[[2]], nd)[,1])) %>%
    mutate(Gc_ref = exp(fitted(bmmodels[[3]], nd)[,1])) %>%
    mutate(I_herb_ref = exp(fitted(bmmodels[[4]], nd)[,1])) %>%
    mutate(I_pisc_ref = exp(fitted(bmmodels[[5]], nd)[,1]))  %>%
    select(-biomass_tot) %>%
    right_join(summary_transect_complete) %>%
    mutate(multif = as.numeric(Fn > Fn_ref & 
                                 Fp > Fp_ref & 
                                 Gc > Gc_ref & 
                                 I_herb > I_herb_ref & 
                                 I_pisc > I_pisc_ref)) %>%
    mutate(logbiomass = log(biomass_tot))

 fit_mf1 <- brm(multif ~ mean + logbiomass,
                 data = ndp, chain = 1, cores = 1, family = "bernoulli")

 
  me <- marginal_effects(fit_mf1, method = "fitted")

  plot <- 
  ggplot(me$logbiomass) +
    geom_ribbon(aes(x = logbiomass, ymin = lower__, ymax = upper__), fill = "grey60", alpha = 0.5) +
    geom_line(aes(x = logbiomass, y = estimate__), size = 1, color = "grey40") +
    geom_vline(xintercept = log(100), linetype = 3, alpha = 0.7) +
    geom_vline(xintercept = log(450), linetype = 2) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    #geom_vline(xintercept = log(1000), linetype = 3, alpha = 0.7)+
    labs(x = "log(biomass) (g/m2)", y = "fitted probability of MF") +
    theme_minimal() +
    theme(axis.line = element_line(), axis.ticks = element_line())
  
  ggsave("output/plots/annex_fig1_mf_logbiomass.png", plot, width = 8, height = 6)
  
  # sum(ndp$biomass_tot > 100)/nrow(ndp)
  # sum(ndp$biomass_tot > 450 & ndp$multif>0)/nrow(ndp)
  # sum(ndp$biomass_tot > 1000)/nrow(ndp)
  # sum(ndp$multif)/nrow(ndp)
  # test <- ndp %>%
  #   filter(multif>0)
  # 
  # f <- fitted(fit_mf1)
  # sum(f[,1]>0.5)/nrow(f)
  
  return(plot)
  
}

make_annex_fig2 <- function(summary_transect_complete, residuals){
 
   bins <- summary_transect_complete %>%
    mutate(bin = cut(biomass_tot, breaks = seq(0, 36200, by = 50))) %>%
    group_by(bin) %>%
    mutate(biomass_bin = mean(biomass_tot)) %>%
    group_by(bin, biomass_bin) %>%
    summarise_if(is.numeric, function(x){max(x, na.rm = TRUE)/min(x,na.rm = TRUE)}) %>%
    pivot_longer(names_to = "var", values_to = "value", cols = c(Fn, Fp, Gc, I_herb, I_pisc)) %>%
    filter(!value == 1)
   
   plot_biomass <-
     ggplot(bins, aes(x = var, y = value)) +
     geom_violin(size = 0.5, alpha = 0.7, fill = "darkgrey", color = "black", draw_quantiles = c(0.25, 0.5, 0.75)) +
     scale_y_continuous(trans = "log10", breaks = c(2,10,100,1000, 10000)) +
     theme_bw() +
     labs(x = "", y = "fold variation per biomass class") +
     scale_x_discrete(labels = c("N excretion", "P excretion",
                                 "Production", "Herbivory", "Piscivory"),
                      position = "top") +
     theme(axis.text.x = element_text(angle = 45, hjust = 0),
           panel.grid = element_blank(),
           panel.border = element_blank(),
           text = element_text(size = 10),
           axis.line.y = element_line(),
           axis.ticks.x = element_blank())
   
   fluxs <- residuals %>%
     drop_na() %>%
     left_join(summary_transect_complete)
   
   fcor <- cor(select(fluxs, Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r)) %>% 
     as.data.frame() 
   
   fcor <- fcor %>%
     mutate(x = rownames(fcor)) %>%
     pivot_longer(names_to = "y", values_to = "cor", -x) 
   
   cplot <-
     ggplot(fcor) +
     geom_tile(aes(x = x, y = y, fill = cor), color = "white", size = 2) +
     geom_text(aes(x = x, y = y, label = round(cor, 2))) +
     scale_fill_fish("Hypsypops_rubicundus", direction = -1, 
                     limits = c(-1, 1), name = "correlation \nof residuals") +
     scale_x_discrete(labels = c("N excretion", "P excretion", "Production",
                                 "Herbivory", "Piscivory"), 
                      position = "top") +
     scale_y_discrete(labels = c("N excretion", "P excretion", "Production",
                                 "Herbivory", "Piscivory")) +
     coord_equal() +
     theme_bw() +
     theme(panel.grid = element_blank(), panel.border = element_blank(), 
           axis.title = element_blank(), 
           axis.text.x = element_text(angle = 45, hjust = 0),
           axis.text.y = element_text(angle = 45, vjust = 0),
           text = element_text(size = 10),
           axis.ticks = element_blank())
   
   plot <-
     plot_biomass + cplot  +  
     plot_annotation(tag_levels = 'a', tag_suffix = ")")
   
   ggsave("output/plots/annex_fig2.png", plot, width = 10, height = 6)
   
   return(plot)
}

make_annex_fig3 <- function(summary_transect_complete, residuals){
  
  
  hum <- read_csv("data/humanPopulation.csv") 
  
  fluxs <- residuals %>%
    drop_na() %>%
    left_join(summary_transect_complete) %>%
    left_join(select(hum, - lat,- lon, - locality, - region))
    
  tall <- pivot_longer(fluxs, c(Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r)) %>%
    group_by(name) %>%
    mutate(value = standard(value))
  
  col <- fish(n = 5, option = "Callanthias_australis")
  
  plot <- 
  ggplot(tall, aes(x = log(grav_Nmarket), y = value), alpha = 0.2) +
    geom_smooth(method = "lm", aes(color = name, fill = name), alpha = 0.2) +
    labs(x = "log(gravity to markets)", y = "Standarized residuals") +
    scale_color_manual(values = c(col), name = "",
                       labels = c("N excretion", "P excretion",
                                  "Production", 
                                  "Herbivory", "Piscivory")) +
    scale_fill_manual(values = c(col),  name = "",
                      labels = c("N excretion", "P excretion",
                                 "Production", 
                                 "Herbivory", "Piscivory")) +
    theme_classic() +
    theme(legend.position = "bottom")
  
  ggsave("output/plots/annex_fig3_gravity_markets.png", plot, width = 8, height = 6)
  
  
}







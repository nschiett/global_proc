
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
  tot <- group_by(result, model) %>%
    summarise(tot = sum(abs(mean))) 
  
  # filter out variables that we don't need 
  result <- result %>%
    dplyr::filter(!variable == "b_Intercept") %>%
    left_join(tot) %>% 
    filter(diff == FALSE) %>%
    filter(!variable == "b_standardlogbiomass") %>%
    filter(!variable == "b_standardmean") %>%
    mutate(variable = fct_relevel(variable, 
                                  "b_standardimm_q1",
                                  "b_standardimm_m",
                                  "b_standardimm_q3",
                                  "b_standardsize_q1",
                                  "b_standardsize_m",
                                  "b_standardsize_q3",
                                  "b_standardtroph_q1",
                                  "b_standardtroph_m",
                                  "b_standardtroph_q3",
                                  "b_standardnspec"
                                  ))
    
  
  # plot
  slopes <-
    ggplot(result) +
    geom_vline(xintercept = seq(1, 10 ,1) + .5, color = "lightgrey",
               size = 0.5, linetype = 1, alpha = 0.7) +
    geom_hline(yintercept = 0, size = 1, color = "black") +
    geom_linerange(aes(x = variable, ymin = lq, ymax = uq, color = model),
                   position = position_dodge(.7), size = 1, linetype = 1) + 
    geom_point(aes(x = variable, y = mean, color = model), position = position_dodge(.7), size = 3) +
    coord_flip() +
    theme_minimal() +
    scale_color_fish_d(option = "Callanthias_australis", 
                       labels = c("N excretion", "P excretion", "Production", "Herbivory", "Piscivory"),
                       name = "")  + 
    labs(y = "Effect size", x = "") +
    scale_x_discrete(
      labels = rev(c( "Richness", "Trophic level (97.5%)", "Trophic level (median)", "Trophic level (2.5%)",
                  "Size (97.5%)", "Size (median)", "Size (2.5%)",
                  "Immaturity (97.5%)", "Immaturity (median)", "Immaturity (2.5%)"
                   ))
      ) +
   
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
    filter(!Family == "Fistulariidae") %>%
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
                xmin = -0.075, xmax = -0.11, ymin = 11.75, ymax = 12.25) +
    add_fishape(family = "Labridae", option = NA, alpha = 0.7,
                xmax = -0.075, xmin = -0.11, ymin = 10.75, ymax = 11.25) +
    add_fishape(family = "Pomacentridae", option = NA,  alpha = 0.7,
                xmin = -0.075, xmax = -0.11, ymin = 9.75, ymax = 10.25) +
    add_fishape(family = "Haemulidae", option = NA,  alpha = 0.7,
                xmax = -0.075, xmin = -0.11, ymin = 8.75, ymax = 9.25) +
    add_fishape(family = "Lutjanidae", option = NA,  alpha = 0.7,
                xmin = -0.075, xmax = -0.11, ymin = 7.75, ymax = 8.25) +
    add_fishape(family = "Serranidae", option = "Cephalopholis_argus",  alpha = 0.7,
                xmax = -0.075, xmin = -0.11, ymin = 6.75, ymax = 7.25) +
    add_fishape(family = "Balistidae", option = NA,  alpha = 0.7,
                xmin = -0.075, xmax = -0.11, ymin = 5.75, ymax = 6.25) +
    add_fishape(family = "Mullidae", option = NA,  alpha = 0.7,
                xmax = -0.075, xmin = -0.11, ymin = 4.75, ymax = 5.25) +
    add_fishape(family = "Kyphosidae", option = NA,  alpha = 0.7,
                xmin = -0.075, xmax = -0.11, ymin = 3.75, ymax = 4.25) +
    add_fishape(family = "Holocentridae", option = NA,  alpha = 0.7,
                xmax = -0.075, xmin = -0.11, ymin = 2.75, ymax = 3.25) +
    add_fishape(family = "Chaetodontidae", option = NA,  alpha = 0.7,
                xmin = -0.075, xmax = -0.11, ymin = 1.75, ymax = 2.25) +
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

make_fig4 <- function(contributions, vulnerability, herb_pisc, residuals){

  
  con <- left_join(contributions, vulnerability) %>%
    left_join(herb_pisc$contributions_herb_pisc) %>%
    group_by(transect_id) %>%
    dplyr::summarize(
      vuln_bm_fi = sum(biomass_p * vuln_fi, na.rm = TRUE),
      vuln_Fn_fi = sum(Fn_p * vuln_fi, na.rm = TRUE),
      vuln_Fp_fi = sum(Fp_p * vuln_fi, na.rm = TRUE),
      vuln_Gc_fi = sum(Gc_p * vuln_fi, na.rm = TRUE),
      vuln_Iherb_fi = sum(I_herb_p[I_herb_p > 0] * vuln_fi[I_herb_p > 0], na.rm = TRUE),
      vuln_Ipisc_fi = sum(I_pisc_p[I_pisc_p > 0] * vuln_fi[I_pisc_p > 0], na.rm = TRUE),
      vuln_bm_cl = sum(biomass_p * vuln_climate, na.rm = TRUE),
      vuln_Fn_cl = sum(Fn_p * vuln_climate, na.rm = TRUE),
      vuln_Fp_cl = sum(Fp_p * vuln_climate, na.rm = TRUE),
      vuln_Gc_cl = sum(Gc_p * vuln_climate, na.rm = TRUE),
      vuln_Iherb_cl = sum(I_herb_p[I_herb_p > 0] * vuln_climate[I_herb_p > 0], na.rm = TRUE),
      vuln_Ipisc_cl = sum(I_pisc_p[I_pisc_p > 0] * vuln_climate[I_pisc_p > 0], na.rm = TRUE)) %>%
    filter(vuln_Ipisc_cl > 0, vuln_Iherb_cl > 0) %>%
    mutate(dv_Fn_fi = vuln_bm_fi - vuln_Fn_fi,
           dv_Fp_fi = vuln_bm_fi - vuln_Fp_fi,
           dv_Gc_fi = vuln_bm_fi - vuln_Gc_fi,
           dv_Iherb_fi = vuln_bm_fi - vuln_Iherb_fi,
           dv_Ipisc_fi = vuln_bm_fi - vuln_Ipisc_fi,
           dv_Fn_cl = vuln_bm_cl - vuln_Fn_cl,
           dv_Fp_cl = vuln_bm_cl - vuln_Fp_cl,
           dv_Gc_cl = vuln_bm_cl - vuln_Gc_cl,
           dv_Iherb_cl = vuln_bm_cl - vuln_Iherb_cl,
           dv_Ipisc_cl = vuln_bm_cl - vuln_Ipisc_cl)
  
  con_long <- con %>%
    pivot_longer(cols = 14:23) %>%
    separate(name, into = c("v", "fun", "impact"), sep = "_", remove = TRUE) %>%
    pivot_wider(names_from = impact, values_from = value) %>%
    mutate(high_fi = fi < 0, 
           high_cl = cl < 0,
           high_both = cl < 0 & fi < 0)
  
  sum <- con_long %>%
    group_by(fun) %>%
    summarize(fi_high = sum(high_fi)/n(),
              cl_high = sum(high_cl)/n(),
              both_high = sum(high_both)/n()) %>%
    pivot_longer(c(fi_high, cl_high, both_high)) %>%
    mutate(name = fct_relevel(name, (unique(name))))
  
  # ggplot(sum) + 
  #   geom_bar(aes(name, value, fill = fun), stat = "identity",
  #            position = position_dodge(), alpha = 0.9) +
  #   labs(y = "Proportion communities with increased vulnerability", x = "") +
  #   scale_x_discrete(labels = c("Fishing", "Climate change", "Both")) +
  #   scale_fill_fish_d(option = "Callanthias_australis",
  #                     labels = c("N excretion", 
  #                                "P excretion", 
  #                                "Production",
  #                                "Herbivory",
  #                                "Piscivory")) +
  #   theme_classic() +
  #   theme(panel.grid.major.x = element_blank(),
  #         panel.grid.major.y = element_line(color = "lightgrey"),
  #         panel.grid.minor = element_blank(), legend.position = "none", 
  #         text = element_text(size = 12))  
  
  sum <- con_long %>%
    group_by(fun) %>%
    summarize(fi_high = sum(high_fi)/n(),
              cl_high = sum(high_cl)/n(),
              both_high = sum(high_both)/n()) %>%
    pivot_longer(c(fi_high, cl_high, both_high)) %>%
    mutate(name = fct_relevel(name, c("fi_high", "both_high", "cl_high")))
  
  plot <- 
  ggplot(sum) + 
    geom_hline(yintercept = 0.5, color = "lightgrey", size = 0.4) +
    geom_hline(yintercept = 0.75, color = "lightgrey", size = 0.4) +
    geom_hline(yintercept = 0.25, color = "lightgrey", size = 0.4) +
    annotate("text", x = -Inf, y = c(0.3, 0.55, 0.8), 
             label = c("0.25", "0.50", "0.75") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
    #geom_text(aes(x = "fi_high", y = 0.6, label = "0.5")) +
    geom_bar(aes(name, value, fill = fun), stat = "identity",
             position = position_dodge(), alpha = 0.9, width = 0.7) +
    labs(y = "Proportion communities with increased vulnerability", x = "") +
    annotate("text", x = c("fi_high", "both_high", "cl_high"), y = c(0.9, 0.9, 0.9), 
             label = c("Fishing", "Both", "Climate change"), angle = c(-60, 0, 60), size = 5) +
    scale_fill_fish_d(option = "Callanthias_australis",
                      labels = c("N excretion", 
                                 "P excretion", 
                                 "Production",
                                 "Herbivory",
                                 "Piscivory"),
                      name = "Proportion communities with increased functional vulnerability") +
    theme_void() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), legend.position = c(0.5, 0.95), 
          text = element_text(size = 12),
          legend.margin = margin(0,0,0,0), panel.spacing = unit(c(0,0,0,0), units = "cm"), 
          legend.direction = "horizontal",
          )  +
    scale_y_continuous(limits = c(-0.3, 0.95), breaks = c(0,0.5, 1)) +
    guides(fill = guide_legend(title.position="top", legend.title.align = 0.5, title.vjust = 1,)) +
    coord_polar()
  
  
  # 
  # res_long <- pivot_longer(residuals, 2:6, values_to = "residual") %>%
  #   mutate(fun = case_when(name == "Fn_r" ~ "Fn",
  #                          name == "Fp_r" ~ "Fp",
  #                          name == "Gc_r" ~ "Gc",
  #                          name == "I_herb_r" ~ "Iherb",
  #                          name == "I_pisc_r" ~ "Ipisc")) %>%
  #   select(-name)
  # 
  # vu <- data.frame(
  #   fun = c("Fn", "Fp", "Gc", "Iherb", "Ipisc"),
  #   m_fi = c(
  #     median(con$vuln_Fn_fi),
  #     median(con$vuln_Fp_fi),
  #     median(con$vuln_Gc_fi),
  #     median(con[con$vuln_Iherb_fi>0,]$vuln_Iherb_fi),
  #     median(con[con$vuln_Ipisc_fi>0,]$vuln_Ipisc_fi)
  #   ),
  #   m_cl = c(
  #     median(con$vuln_Fn_cl),
  #     median(con$vuln_Fp_cl),
  #     median(con$vuln_Gc_cl),
  #     median(con[con$vuln_Iherb_cl>0,]$vuln_Iherb_cl),
  #     median(con[con$vuln_Ipisc_cl>0,]$vuln_Ipisc_cl)
  #   )
  # ) %>%
  #   mutate(m_fi = mean(m_fi),
  #          m_cl = mean(m_cl))
  # 
  # 
  # q_cl <- quantile(c(
  #   (con$vuln_Fn_cl),
  #   (con$vuln_Fp_cl),
  #   (con$vuln_Gc_cl),
  #   (con[con$vuln_Iherb_cl>0,]$vuln_Iherb_cl),
  #   (con[con$vuln_Ipisc_cl>0,]$vuln_Ipisc_cl)
  # ), c(0.25, 0.4, 0.6, 0.75))
  # q_fi <- quantile(c(
  #   (con$vuln_Fn_fi),
  #   (con$vuln_Fp_fi),
  #   (con$vuln_Gc_fi),
  #   (con[con$vuln_Iherb_fi>0,]$vuln_Iherb_fi),
  #   (con[con$vuln_Ipisc_fi>0,]$vuln_Ipisc_fi)
  # ), c(0.25, 0.4, 0.6, 0.75))
  # 
  # 
  # 
  # prop_fi <- con_long %>%
  #   left_join(res_long) %>%
  #   left_join(vu) %>%
  #   mutate(cat_fi = case_when(fi > q_fi[4] ~ "High",
  #                             fi > q_fi[3] ~ "Medium",
  #                             fi > q_fi[2] ~ "Medium",
  #                             fi > q_fi[1] ~ "Medium",
  #                             TRUE ~ "Low")) %>%
  #   mutate(high = residual > 0) %>%
  #   group_by(cat_fi, fun) %>%
  #   dplyr::summarise(n_high_fi = sum(high), n_fi = length(unique(transect_id))) %>%
  #   ungroup() %>% dplyr::group_by(fun) %>%
  #   mutate(ntot_fi = sum(n_fi)) %>%
  #   mutate(prop_high_fi = n_high_fi / ntot_fi, prop_fi = n_fi/ntot_fi) %>%
  #   mutate(cat = factor(cat_fi, levels = c("Very low", "Low", "Medium", "High", "Very high")))
  # 
  # prop_cl <- con_long %>%
  #   left_join(res_long) %>%
  #   left_join(vu) %>%
  #   mutate(cat_cl = case_when(cl > q_cl[4] ~ "High",
  #                             cl > q_cl[3] ~ "Medium",
  #                             cl > q_cl[2] ~ "Medium",
  #                             cl > q_cl[1] ~ "Medium",
  #                             TRUE ~ "Low")) %>%
  #   mutate(high = residual > 0) %>%
  #   group_by(cat_cl, fun) %>%
  #   dplyr::summarise(n_high_cl = sum(high), n_cl = length(unique(transect_id))) %>%
  #   ungroup() %>% dplyr::group_by(fun) %>%
  #   mutate(ntot_cl = sum(n_cl)) %>%
  #   mutate(prop_high_cl = n_high_cl / ntot_cl, prop_cl = n_cl/ntot_cl) %>%
  #   mutate(cat = factor(cat_cl, levels = c("Very low", "Low", "Medium", "High", "Very high")))
  # 
  # p1 <- ggplot(prop_fi[!prop_fi$cat_fi == "Medium",]) +
  #   geom_bar(aes(x = cat, y = prop_fi, fill = fun), 
  #            alpha = 0.9, stat = "identity", position = "dodge", width = 0.9) +
  #   scale_fill_fish_d(option = "Callanthias_australis",
  #                     labels = c("N excretion", 
  #                                "P excretion", 
  #                                "Production",
  #                                "Herbivory",
  #                                "Piscivory")) +
  #   labs(y = "Proportion of communities", x = "Vulnerability to fishing", fill = "") +
  #   theme_bw() +
  #   theme(panel.grid.major.x = element_blank(), 
  #         panel.grid.minor = element_blank(), legend.position = "top",
  #         text = element_text(size = 12))  
  # 
  # p2 <- ggplot(prop_cl[!prop_cl$cat_cl == "Medium",]) +
  #   geom_bar(aes(x = cat, y = prop_cl, fill = fun), 
  #            alpha = 0.9, stat = "identity", position = "dodge", width = 0.9) +
  #   scale_fill_fish_d(option = "Callanthias_australis",
  #                     labels = c("N excretion", 
  #                                "P excretion", 
  #                                "Production",
  #                                "Herbivory",
  #                                "Piscivory")) +
  #   labs(y = "Proportion of communities", x = "Vulnerability to climate change", fill = "") +
  #   theme_bw() +
  #   theme(panel.grid.major.x = element_blank(), 
  #         panel.grid.minor = element_blank(), legend.position = "none",
  #         text = element_text(size = 12))  
  # 
  # double <- 
  #   con_long %>%
  #   left_join(res_long) %>%
  #   left_join(vu) %>%
  #   mutate(cat_double = case_when(fi > q_fi[4] & cl > q_cl[4] ~ "High",
  #                                 fi < q_fi[1] & cl < q_cl[1] ~ "Low",
  #                                 TRUE ~ "Medium")) %>%
  #   group_by(cat_double, fun) %>%
  #   dplyr::summarise( n = length(unique(transect_id))) %>%
  #   ungroup() %>% dplyr::group_by(fun) %>%
  #   mutate(ntot = sum(n)) %>%
  #   mutate(prop = n/ntot) %>%
  #   mutate(cat_double = factor(cat_double, levels = c("Low", "Medium", "High")))
  # 
  # 
  # p3 <- ggplot(double[!double$cat_double == "Medium",]) +
  #   geom_bar(aes(x = cat_double, y = prop, fill = fun), 
  #            alpha = 0.9, stat = "identity", position = "dodge", width = 0.9) +
  #   scale_fill_fish_d(option = "Callanthias_australis",
  #                     labels = c("N excretion", 
  #                                "P excretion", 
  #                                "Production",
  #                                "Herbivory",
  #                                "Piscivory")) +
  #   labs(y = "Proportion of communities", x = "Vulnerability to both", fill = "") +
  #   theme_bw() +
  #   theme(panel.grid.major.x = element_blank(), 
  #         panel.grid.minor = element_blank(), 
  #         legend.position = "none",
  #         text = element_text(size = 12))  
  # 
  # p3
  # 
  # p1 + p2 + p3 + plot_layout(ncol = 1)
  
  ggsave("output/plots/figure_4_vuln_com.png", plot,  width = 6, height = 6)
  
}

# make_annex_fig1 <- function(summary_transect_complete, bmmodels){
# 
#   nd <- summary_transect_complete %>%
#     select(mean) %>%
#     unique() %>%
#     mutate(biomass_tot = 100) 
#   
#   ndp <- nd %>%
#     mutate(Fn_ref = exp(fitted(bmmodels[[1]], nd)[,1])) %>%
#     mutate(Fp_ref = exp(fitted(bmmodels[[2]], nd)[,1])) %>%
#     mutate(Gc_ref = exp(fitted(bmmodels[[3]], nd)[,1])) %>%
#     mutate(I_herb_ref = exp(fitted(bmmodels[[4]], nd)[,1])) %>%
#     mutate(I_pisc_ref = exp(fitted(bmmodels[[5]], nd)[,1]))  %>%
#     select(-biomass_tot) %>%
#     right_join(summary_transect_complete) %>%
#     mutate(multif = as.numeric(Fn > Fn_ref & 
#                                  Fp > Fp_ref & 
#                                  Gc > Gc_ref & 
#                                  I_herb > I_herb_ref & 
#                                  I_pisc > I_pisc_ref)) %>%
#     mutate(logbiomass = log(biomass_tot))
# 
#  fit_mf1 <- brm(multif ~ mean + logbiomass,
#                  data = ndp, chain = 1, cores = 1, family = "bernoulli")
# 
#  
#   me <- marginal_effects(fit_mf1, method = "fitted")
# 
#   plot <- 
#   ggplot(me$logbiomass) +
#     geom_ribbon(aes(x = logbiomass, ymin = lower__, ymax = upper__), fill = "grey60", alpha = 0.5) +
#     geom_line(aes(x = logbiomass, y = estimate__), size = 1, color = "grey40") +
#     geom_vline(xintercept = log(100), linetype = 3, alpha = 0.7) +
#     geom_vline(xintercept = log(450), linetype = 2) +
#     geom_hline(yintercept = 0.5, linetype = 2) +
#     #geom_vline(xintercept = log(1000), linetype = 3, alpha = 0.7)+
#     labs(x = "log(biomass) (g/m2)", y = "fitted probability of MF") +
#     theme_minimal() +
#     theme(axis.line = element_line(), axis.ticks = element_line())
#   
#   ggsave("output/plots/annex_fig1_mf_logbiomass.png", plot, width = 8, height = 6)
#   
#   # sum(ndp$biomass_tot > 100)/nrow(ndp)
#   # sum(ndp$biomass_tot > 450 & ndp$multif>0)/nrow(ndp)
#   # sum(ndp$biomass_tot > 1000)/nrow(ndp)
#   # sum(ndp$multif)/nrow(ndp)
#   # test <- ndp %>%
#   #   filter(multif>0)
#   # 
#   # f <- fitted(fit_mf1)
#   # sum(f[,1]>0.5)/nrow(f)
#   
#   return(plot)
#   
# }

########### Figures supplemental ###########

make_annex_fig1 <- function(summary_transect_complete, residuals, bmmodels){
  
  fold <- function(x){max(x, na.rm = TRUE)/min(x,na.rm = TRUE)}
  ma <- function(x){max(x , na.rm = TRUE)}
  mi <- function(x){min(x , na.rm = TRUE)}
  
   bins <- summary_transect_complete %>%
    select(biomass_tot, Fn, Fp, Gc, I_herb, I_pisc) %>%
    mutate(bin = cut(biomass_tot, breaks = seq(0, 36200, by = 100))) %>%
    group_by(bin) %>%
    mutate(biomass_bin = mean(biomass_tot)) %>%
    group_by(bin, biomass_bin) %>%
    summarise_if(is.numeric, funs(fold, ma, mi))%>%
    filter(biomass_bin < 4000)
   
   # predicted relationship

   pdata <- data.frame(
     biomass_tot = seq(50, 4000, by = 50)
   ) %>%
     mutate(mean = 26)
   
   pred1 <- predict(bmmodels[[1]], newdata = pdata) %>%
     as.data.frame() %>%
     cbind(pdata)
   pred2 <- predict(bmmodels[[2]], newdata = pdata) %>%
     as.data.frame() %>%
     cbind(pdata)
   pred3 <- predict(bmmodels[[3]], newdata = pdata) %>%
     as.data.frame() %>%
     cbind(pdata)
   pred4 <- predict(bmmodels[[4]], newdata = pdata) %>%
     as.data.frame() %>%
     cbind(pdata)
   pred5 <- predict(bmmodels[[5]], newdata = pdata) %>%
     as.data.frame() %>%
     cbind(pdata)
   
   p1 <- ggplot(bins) +
     geom_ribbon(aes(x = biomass_tot, ymin = exp(`Q2.5`), ymax = exp(`Q97.5`)), 
                 data = pred1, alpha = 0.3, fill = "#003366") +
     geom_line(aes(x = biomass_tot, y = exp(Estimate)), 
                 data = pred1, alpha = 0.8, color = "#003366") +
     # geom_point(aes(x = biomass_tot, y = Fn), 
     #            size = 0.5, alpha = 0.1,
     #            data = filter(summary_transect_complete, biomass_tot < 4000)) +
     geom_linerange(aes(x = biomass_bin, ymin = Fn_mi, ymax = Fn_ma),
                    size = 1.5 , alpha = 0.7) +
     labs(x = "bins of biomass (g/m2)", y = "N excretion (gN/m2/day)") +
     theme_bw()
   
   
   p2 <- ggplot(bins) +
     geom_ribbon(aes(x = biomass_tot, ymin = exp(`Q2.5`), ymax = exp(`Q97.5`)), 
                 data = pred2, alpha = 0.3, fill = "#003366") +
     geom_line(aes(x = biomass_tot, y = exp(Estimate)), 
               data = pred2, alpha = 0.8, color = "#003366") +
     geom_linerange(aes(x = biomass_bin, ymin = Fp_mi, ymax = Fp_ma),
                    size = 1.5 , alpha = 0.7) +
     labs(x = "bins of biomass (g/m2)", y = "P excretion (gP/m2/day)") +
     theme_bw()
   p3 <- ggplot(bins) +
     geom_ribbon(aes(x = biomass_tot, ymin = exp(`Q2.5`), ymax = exp(`Q97.5`)), 
                 data = pred3, alpha = 0.3, fill = "#003366") +
     geom_line(aes(x = biomass_tot, y = exp(Estimate)), 
               data = pred3, alpha = 0.8, color = "#003366") +
     geom_linerange(aes(x = biomass_bin, ymin = Gc_mi, ymax = Gc_ma),
                    size = 1.5 , alpha = 0.7) +
     labs(x = "bins of biomass (g/m2)", y = "Production (gC/m2/day)") +
     theme_bw()
   p4 <- ggplot(bins) +
     geom_ribbon(aes(x = biomass_tot, ymin = exp(`Q2.5`), ymax = exp(`Q97.5`)), 
                 data = pred4, alpha = 0.3, fill = "#003366") +
     geom_line(aes(x = biomass_tot, y = exp(Estimate)), 
               data = pred4, alpha = 0.8, color = "#003366") +
     geom_linerange(aes(x = biomass_bin, ymin = I_herb_mi, ymax = I_herb_ma),
                    size = 1.5 , alpha = 0.7) +
     labs(x = "bins of biomass (g/m2)", y = "Herbivory (gC/m2/day)") +
     theme_bw()
   p5 <- ggplot(bins) +
     geom_ribbon(aes(x = biomass_tot, ymin = exp(`Q2.5`), ymax = exp(`Q97.5`)), 
                 data = pred5, alpha = 0.3, fill = "#003366") +
     geom_line(aes(x = biomass_tot, y = exp(Estimate)), 
               data = pred5, alpha = 0.8, color = "#003366") +
     geom_linerange(aes(x = biomass_bin, ymin = I_pisc_mi, ymax = I_pisc_ma),
                    size = 1.5 , alpha = 0.7) +
     labs(x = "bins of biomass (g/m2)", y = "Piscivory (gC/m2/day)") +
     theme_bw()
   
   
   bins_long <- bins %>%
    pivot_longer(names_to = "var", values_to = "value", cols = c(Fn_fold, Fp_fold, Gc_fold, I_herb_fold, I_pisc_fold)) %>%
    filter(!value == 1)
   
   
   plot_biomass <-
     ggplot(bins_long, aes(x = var, y = value)) +
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
   
   layout <-
   "AAAA
   BBBB
   CCCC
   DDDD
   EEEE
   FFGG
   FFGG"
   
   plot <- p1 + p2 + p3 + p4 + p5 +
     plot_biomass + cplot  +  
     plot_annotation(tag_levels = 'a', tag_suffix = ")") +
     plot_layout(design = layout)
   
   ggsave("output/plots/annex_fig1.png", plot, width = 10, height = 16)
   
   return(plot)
}

make_annex_fig2 <- function(summary_transect_complete, residuals){
  
  
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
  
  ggsave("output/plots/annex_fig2_gravity_markets.png", plot, width = 8, height = 6)
  
  
}



make_rank_plots <- function(contributions, herb_pisc){

  con <- left_join(contributions, herb_pisc$contributions_herb_pisc) %>%
    ungroup() %>%
    dplyr::group_by(bioregion) %>%
    dplyr::mutate(n_transect = length(unique(as.character(transect_id)))) %>%
    group_by(bioregion, Family, species) %>%
    summarize(rel_occ = length(Fn_p)/unique(n_transect),
              Fn_p_m = median(Fn_p),
              Fp_p_m = median(Fp_p),
              Gc_p_m = median(Gc_p),
              I_herb_p_m = median(I_herb_p),
              I_pisc_p_m = median(I_pisc_p))
  
  main_fam <- dplyr::select(con, species, Family, bioregion) %>%
    group_by(Family) %>%
    summarize(n = n()) %>%
    mutate(
      rank = row_number(desc(n) )
    ) %>%
    filter(rank < 12) %>%
    filter(!Family == "Pomacanthidae")
  
  
  rank <- con %>%
    filter(Family %in% main_fam$Family) %>%
    filter(rel_occ > 0.01) %>%
    group_by(bioregion) %>%
    mutate(
      Fn_r = dense_rank(desc(Fn_p_m)),
      Fp_r = dense_rank(desc(Fp_p_m)),
      Gc_r = dense_rank(desc(Gc_p_m)),
      I_herb_r = dense_rank(desc(I_herb_p_m)),
      I_pisc_r = dense_rank(desc(I_pisc_p_m))
    ) %>%
    ungroup() %>%
    filter(
      Fn_r < 11 |
        Fp_r < 11 |
        Gc_r < 11 |
        I_herb_r < 11|
        I_pisc_r < 11
    ) %>% pivot_longer(
      c(Fn_p_m, Fp_p_m, Gc_p_m, I_herb_p_m, I_pisc_p_m),
      names_to = "name",
      values_to = "value"
    ) %>%
    group_by(name, bioregion) %>%
    mutate(
      rank = row_number(desc(value) )
    ) %>%
    filter(rank < 11) %>%
    mutate(Family = as_factor(as.character(Family))) %>%
    filter(value> 0) %>%
    group_by(bioregion, name) %>%
    dplyr::mutate(value2 = value/max(value))%>%
    ungroup() %>%
    mutate(Family = as_factor(as.character(Family)))
  
  rank_fam <- con %>%
    group_by(bioregion, Family) %>%
    summarise_if(is.numeric, sum
    ) %>%
    ungroup() %>%
    pivot_longer(
      c(Fn_p_m, Fp_p_m, Gc_p_m, I_herb_p_m, I_pisc_p_m),
      names_to = "name",
      values_to = "value"
    ) %>%
    group_by(name, bioregion) %>%
    mutate(
      rank = dense_rank(desc(value))
    ) %>%
    filter(rank < 11) %>%
    mutate(Family = as_factor(as.character(Family))) %>%
    filter(value> 0) %>%
    group_by(bioregion, name) %>%
    dplyr::mutate(value2 = value/max(value))  %>%
    ungroup() %>%
    mutate(Family = as_factor(as.character(Family)))
  
 
  ggplot(rank) +
    geom_bar(aes(y = value2, x = as.character(1/rank), fill = Family), stat = "identity", alpha = 0.8) +
    geom_text(aes(y = 0.5, x = as.character(1/rank), label = species), size = 2) +
    facet_grid(bioregion~name, scales="free", space="free",
               labeller = labeller(
                 name = c(Fn_p_m = "N excretion",
                          Fp_p_m = "P excretion",
                          Gc_p_m = "Production",
                          I_herb_p_m = "Herbovory",
                          I_pisc_p_m = "Piscivory"),
                 bioregion = c(c_indopacific = "CIP",
                               c_pacific = "CP",
                               e_atlantic = "EA",
                               e_pacific = "EP",
                               w_atlantic = "WA",
                               w_indian = "WI")
               )) + coord_flip() +
    scale_fill_fish_d(option = "Pseudocheilinus_tetrataenia") +
    labs(y = "Rescaled contribution") +
    theme_bw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), panel.grid = element_blank())
  ggsave("output/plots/con_species_rank.png",  width = 8, height = 8)
  
  
  ggplot(rank_fam) +
    geom_bar(aes(y = value2, x = as.character(1/rank), fill = Family), alpha = 0.8, stat = "identity") +
    geom_text(aes(y = 0.5, x = as.character(1/rank), label = Family), size = 2) +
    facet_grid(bioregion~name, scales="free", space="free",
               labeller = labeller(
                 name = c(Fn_p_m = "N excretion",
                          Fp_p_m = "P excretion",
                          Gc_p_m = "Production",
                          I_herb_p_m = "Herbovory",
                          I_pisc_p_m = "Piscivory"),
                 bioregion = c(c_indopacific = "CIP",
                               c_pacific = "CP",
                               e_atlantic = "EA",
                               e_pacific = "EP",
                               w_atlantic = "WA",
                               w_indian = "WI")
               )) + coord_flip() +
    scale_fill_fish_d(option = "Pseudocheilinus_tetrataenia") +
    labs(y = "Rescaled contribution") +
    theme_bw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          panel.grid = element_blank())
  
  ggsave("output/plots/con_family_rank.png",  width = 8, height = 8)
}


make_pp_plots <- function(bmmodels, procmodels, commodels){
  
  a1 <- pp_check(bmmodels[[1]], nsamples = 50)
  #ggsave("output/plots/pp_plot_bmmodel_Fn.png", width = 8, height = 6)
  a2 <- pp_check(bmmodels[[2]], nsamples = 50)
  #ggsave("output/plots/pp_plot_bmmodel_Fn.png", width = 8, height = 6)
  a3 <- pp_check(bmmodels[[3]], nsamples = 50)
  #ggsave("output/plots/pp_plot_bmmodel_Gc.png", width = 8, height = 6)
  a4 <- pp_check(bmmodels[[4]], nsamples = 50)
  #ggsave("output/plots/pp_plot_bmmodel_herb.png", width = 8, height = 6)
  a5 <- pp_check(bmmodels[[5]], nsamples = 50)
  #ggsave("output/plots/pp_plot_bmmodel_pisc.png", width = 8, height = 6)
  

  
  # pp_check(procmodels[[1]], nsamples = 50)
  # ggsave("output/plots/pp_plot_procmodel_Fn.png", width = 8, height = 6)
  # pp_check(procmodels[[2]], nsamples = 50)
  # ggsave("output/plots/pp_plot_procmodel_Fp.png", width = 8, height = 6)
  # pp_check(procmodels[[3]], nsamples = 50)
  # ggsave("output/plots/pp_plot_procmodel_Gc.png", width = 8, height = 6)
  # pp_check(procmodels[[4]], nsamples = 50)
  # ggsave("output/plots/pp_plot_procmodel_herb.png", width = 8, height = 6)
  # pp_check(procmodels[[5]], nsamples = 50)
  # ggsave("output/plots/pp_plot_procmodel_pisc.png", width = 8, height = 6)
  
  b1 <- pp_check(commodels[[1]], nsamples = 50)
  #ggsave("output/plots/pp_plot_commodel_Fn.png", width = 8, height = 6)
  b2 <- pp_check(commodels[[2]], nsamples = 50)
  #ggsave("output/plots/pp_plot_commodel_Fp.png", width = 8, height = 6)
  b3 <- pp_check(commodels[[3]], nsamples = 50)
  #ggsave("output/plots/pp_plot_commodel_Gc.png", width = 8, height = 6)
  b4 <- pp_check(commodels[[4]], nsamples = 50)
  #ggsave("output/plots/pp_plot_commodel_herb.png", width = 8, height = 6)
  b5 <- pp_check(commodels[[5]], nsamples = 50)
  #ggsave("output/plots/pp_plot_commodel_pisc.png", width = 8, height = 6)

  pp_plots <- a1 + a2 + a3 + a4 + a5 + b1 + b2 + b3 + b4 + b5 + 
    plot_layout(ncol = 2) +
    plot_annotation(tag_levels = 'a', tag_suffix = ")")
  
  ggsave("output/plots/pp_plots.png", pp_plots, width = 6, height = 10)
  
}

####### extra analysis diets #######
alt_diet <- function(tfish, cnpflux){
  
  diet <- read.csv("data/extrapolation_trophic_guilds.csv") 
  
  diet_experts <- read.csv("data/original_experts_classification.csv")%>%
    select(species = Genus_and_species_, Gaspar, Siquiera)
  
  diets <- diet %>%
    select(family, species, p2_m, p4_m, trophic_guild_predicted) %>% 
    left_join(diet_experts) %>%
    mutate(herb = trophic_guild_predicted == 2,
           herb_gas = Gaspar == "HD",
           herb_siq = Siquiera == "HD",
           herb_t1 = trophic_guild_predicted == 2 & p2_m > 0.5,
           pisc = trophic_guild_predicted == 4,
           pisc_gas = Gaspar == "FC",
           pisc_siq = Siquiera == "GC",
           pisc_t1 = trophic_guild_predicted == 4 & p4_m > 0.5
    )
  
  
  tfish <- left_join(tfish, unique(dplyr::select(cnpflux, species, v_m, size_cm, (ends_with("_median"))))) %>% 
    left_join(diets) %>%
    mutate(Ic = Ic_median * abun) %>%
    mutate(I = Ic/area) 
  
  
  tfs <- dplyr::group_by(tfish, region, locality, sites, area, transect_id) %>% 
    dplyr::summarise(
      I_herb = sum(I[herb == TRUE], na.rm = TRUE),
      I_pisc = sum(I[pisc == TRUE], na.rm = TRUE),
      I_herb_t1 = sum(I[herb_t1 == TRUE], na.rm = TRUE),
      I_pisc_t1 = sum(I[pisc_t1 == TRUE], na.rm = TRUE),
      I_herb_gas= sum(I[herb_gas == TRUE], na.rm = TRUE),
      I_pisc_gas = sum(I[pisc_gas == TRUE], na.rm = TRUE),
      I_herb_siq = sum(I[herb_siq == TRUE], na.rm = TRUE),
      I_pisc_siq = sum(I[pisc_siq == TRUE], na.rm = TRUE)
    ) %>% ungroup() %>%
    mutate(d_herb1 =  (I_herb_t1 - I_herb),
           d_pisc1 =  (I_pisc_t1 - I_pisc),
           d_herb_gas = (I_herb_gas - I_herb),
           d_pisc_gas = (I_pisc_gas - I_pisc),
           d_herb_siq = (I_herb_siq - I_herb),
           d_pisc_siq = (I_pisc_siq - I_pisc))
  
  b1 <- ggplot(tfs) +
    geom_histogram(aes(x = (d_herb_gas), y = stat(count) / sum(count)), bins = 30) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_continuous(breaks = c(-5, - 2.5, 0, 2.5, 5), limits = c(-5, 5)) +
    labs(y = "Percent of total", x = "Difference herbivory rate", title = "C) Comparison Gaspar") +
    theme_bw()
  b2 <- ggplot(tfs) +
    geom_histogram(aes(x = (d_pisc_gas), y = stat(count) / sum(count)), bins = 30) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_continuous(breaks = c(-1, - 0.5, 0, 0.5, 1), limits = c(-1, 1)) +
    labs(y = "Percent of total", x = "Difference piscivory rate", title = "D) Comparison Gaspar") +
    theme_bw()
  c1 <- ggplot(tfs) +
    geom_histogram(aes(x = (d_herb_siq), y = stat(count) / sum(count)), bins = 30) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_continuous(breaks = c(-5, - 2.5, 0, 2.5, 5), limits = c(-5, 5)) +
    labs(y = "Percent of total", x = "Difference herbivory rate", title = "E) Comparison Siqueira") +
    theme_bw()
  c2 <- ggplot(tfs) +
    geom_histogram(aes(x = (d_pisc_siq), y = stat(count) / sum(count)), bins = 30) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    scale_x_continuous(breaks = c(-1, - 0.5, 0, 0.5, 1), limits = c(-1, 1)) +
    labs(y = "Percent of total", x = "Difference piscivory rate", title = "F) Comparison Siqueira") +
    theme_bw()
  
  p <-  b1 + b2 + c1 + c2 +plot_layout(ncol = 2)
  
  ggsave("output/plots/Supplemental_herbivory_piscivory_alternative_classification.png", 
         plot = p, 
         height = 6, width = 8)
  
}





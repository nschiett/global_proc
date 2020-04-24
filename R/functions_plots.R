
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
            plot.title = element_text(size = 8, hjust = 0, face = "bold"), 
            plot.margin = unit(c(0.01,0.01,0.01,0.01), units = , "cm"),
            panel.spacing = unit(0, "cm")
      ) 
  }
  
  col <- fish(n = 5, option = "Callanthias_australis")
  
  get_scale <- function(f){
    quant1 <- quantile(f[f < 0], c(0, 0.5, 1), na.rm = TRUE)
    quant2 <- quantile(f[f > 0], c(0, 0.5, 1), na.rm = TRUE)
    quant <- c(quant1[1:2], 0, quant2[2:3])
    quant <- quantile(f, c(0,0.2,0.8,1), na.rm = TRUE)
    quant[1] <- quant[1] - 0.5   # adjust lower bound
    factor(cut(f, breaks = quant), labels = c( "low", "medium", "high"))
  }
  
  get_scale(location_effect$r_loc_Fn)
  
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
               shape = 0, size = 4,
               data = filter(location_effect, r_loc_multi > quantile(location_effect$r_loc_multi, 0.95, na.rm = TRUE))) +
    geom_point(aes(x = lon, y = lat),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(location_effect, r_loc_Fn > quantile(location_effect$r_loc_Fn, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    labs(title = "N excretion") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "right", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.001,0.001,0.001,0.001), units = , "cm"))
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
               shape = 0, size = 4,
               data = filter(location_effect, r_loc_multi > quantile(location_effect$r_loc_multi, 0.95, na.rm = TRUE))) +
    geom_point(aes(x = lon, y = lat),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(location_effect, r_loc_Fp > quantile(location_effect$r_loc_Fp, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    labs(title = "P excretion") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "right", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.001,0.001,0.001,0.001), units = , "cm"))
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
               shape = 0, size = 4,
               data = filter(location_effect, r_loc_multi > quantile(location_effect$r_loc_multi, 0.95, na.rm = TRUE))) +
    geom_point(aes(x = lon, y = lat),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(location_effect, r_loc_Gc > quantile(location_effect$r_loc_Gc, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    labs(title = "Production") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "right", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.001,0.001,0.001,0.001), units = , "cm"))
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
               shape = 0, size = 4,
               data = filter(location_effect, r_loc_multi > quantile(location_effect$r_loc_multi, 0.95, na.rm = TRUE))) +
    geom_point(aes(x = lon, y = lat),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(location_effect, r_loc_I_herb > quantile(location_effect$r_loc_I_herb, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "") +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    labs(title = "Herbivory") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
    theme_worldmap() + 
    theme(legend.position = "right", 
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
               shape = 0, size = 4,
               data = filter(location_effect, r_loc_multi > quantile(location_effect$r_loc_multi, 0.95, na.rm = TRUE))) +
    geom_point(aes(x = lon, y = lat),  alpha = 0.7,
               shape = 1, size = 4,
               data = filter(location_effect, r_loc_I_pisc > quantile(location_effect$r_loc_I_pisc, 0.95, na.rm = TRUE))) +
    scale_color_manual(values = pal, name = "",
                       drop = TRUE, na.translate = F) +
    coord_sf(ylim = c(-35, 35), expand = FALSE) +
    labs(title = "Piscivory") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
    scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
    theme_worldmap()+ 
    theme(legend.position = "right", 
          plot.title = element_text(vjust = -8, hjust = 0.005),
          plot.margin = unit(c(0.0,0.00,0.00,0.00), units = , "cm"))
  e
  
  
  multimap <-
    a + b + c + d + e +
    plot_layout(ncol = 1)
  
  multimap
  ggsave("output/plots/fig1_multimap.png", multimap, width = 8, height = 8)
}




tt <- flux %>%
  group_by(locality, sites) %>%
  summarise(n = length(transect_id)) %>%
  group_by(locality) %>%
  summarise(n = length(sites))

test <- flux %>%
  group_by(bioregion, locality, sites) %>%
  summarise_if(is.numeric, median, na.rm = TRUE) %>%
  ungroup() %>%
  group_by(bioregion, locality) %>%
  summarise_if(is.numeric, median, na.rm = TRUE)

test <- filter(flux,
               Fn_r > quantile(Fn_r, 0.5) &
                 Fp_r > quantile(Fp_r, 0.5) &
                 Gc_r > quantile(Gc_r, 0.5) &
                 I_herb_r > quantile(I_herb_r, 0.5, na.rm = TRUE) &
                 I_pisc_r > quantile(I_pisc_r, 0.5, na.rm = TRUE) 
)


test2 <- flux %>%
  mutate(
    f1 = case_when(Fn_r < quantile(Fn_r, 0.33) ~ 0,
                   Fn_r > quantile(Fn_r, 0.66) ~ 2, 
                   TRUE ~ 1),
    f2 = case_when(Fp_r < quantile(Fp_r, 0.33) ~ 0,
                   Fp_r > quantile(Fp_r, 0.66) ~ 2, 
                   TRUE ~ 1),
    f3 = case_when(Gc_r < quantile(Gc_r, 0.33) ~ 0,
                   Gc_r > quantile(Gc_r, 0.66) ~ 2, 
                   TRUE ~ 1),
    f4 = case_when(I_herb_r < quantile(I_herb_r, 0.33, na.rm = TRUE) ~ 0,
                   I_herb_r > quantile(I_herb_r, 0.66, na.rm = TRUE) ~ 2, 
                   is.na(I_herb_r) ~ NA_real_,
                   TRUE ~ 1),
    f5 = case_when(I_pisc_r < quantile(I_pisc_r, 0.33, na.rm = TRUE) ~ 0,
                   I_pisc_r > quantile(I_pisc_r, 0.66, na.rm = TRUE) ~ 2, 
                   is.na(I_herb_r) ~ NA_real_,
                   TRUE ~ 1)) %>%
  rowwise() %>%
  mutate(mf = sum(f1, f2, f3, f4, f5))  %>%
  filter(mf %in% c(0,1,2,8,9,10))

  test2$mf2 <- normalize(log(test2$biomass_tot)) * test2$mf/100
  
  loc <- test2%>%
    group_by(locality, sites) %>%
    summarize_if(is.numeric, median) %>%
    group_by(locality) %>%
    summarize_if(is.numeric, median)

ggplot(test2, aes(x = as.factor(mf), y = imm_m)) +
  geom_boxplot()
ggplot(test2, aes(x = as.factor(mf), y = imm_q1)) +
  geom_boxplot()
ggplot(test2, aes(x = as.factor(mf), y = imm_q3)) +
  geom_boxplot()
ggplot(test2, aes(x = as.factor(mf), y = size_m)) +
  geom_boxplot()
ggplot(test2, aes(x = as.factor(mf), y = size_q3)) +
  geom_boxplot()
ggplot(test2, aes(x = as.factor(mf), y = size_q1)) +
  geom_boxplot()
ggplot(test2, aes(x = as.factor(mf), y = troph_m)) +
  geom_boxplot()
ggplot(test2, aes(x = as.factor(mf), y = troph_q3)) +
  geom_boxplot()
ggplot(test2, aes(x = as.factor(mf), y = troph_q1)) +
  geom_boxplot()
ggplot(test2, aes(x = as.factor(mf), y = log(nspec))) +
  geom_boxplot()
ggplot(test2, aes(x = as.factor(mf), y = log(biomass_tot))) +
  geom_boxplot()

ggplot(test2, aes(y = (mf2), x = (nspec))) +
  geom_point() + geom_smooth()
ggplot(test2, aes(y = (mf2), x = log(biomass_tot))) +
  geom_point()

cor(select(flux, imm_m, imm_q1, imm_q3, troph_m, troph_q1, 
                          troph_q3, size_m, size_q1, size_q3, nspec, biomass_tot))






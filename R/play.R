loadd(contributions_sp_loc_occ)
loadd(vulnerability)

con <- left_join(contributions_sp_loc_occ, vulnerability) %>% 
  drop_na(vuln_fi, vuln_climate) 
  


tax <- rfishbase::load_taxa()

## functions

quant <- function(x){
  probs <- c(0.05, 0.25, 0.75, 0.95, 0.99)
  quantiles <- quantile(x, prob = probs, na.rm = TRUE)
  quant <- factor(findInterval(x, quantiles))
  return(quant)
}

get_scores <- function(proc){
  probs <- seq(0,1,0.01)
  quantiles <- quantile(proc, prob = probs, na.rm = TRUE)
  quant <- as.numeric(factor(findInterval(proc, quantiles)))/100
  return(quant)
}

normalize <- function(x){
  100 * (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

################ c indopacific ###########

sub <- filter(contributions_sp_loc_occ, occurence > 1, bioregion == "c_indopacific")

tree1 <- fishtree::fishtree_phylogeny(species = sub$species)

tree1t <- as_tibble(tree1)
tax$label <- gsub(" ", "_", tax$Species)
tree1t <- left_join(tree1t, tax)
#tree1t[!is.na(tree1t$Family) & tree1t$Family == "Scaridae", "Family"] <- "Labridae"

tree1d <- as.treedata(tree1t)


data <- filter(sub, species%in% tree1$tip.label) %>% 
  mutate(I_herb_p_m = case_when(I_herb_p_m == 0 ~ NA_real_, 
                                      TRUE ~ I_herb_p_m),
         I_pisc_p_m = case_when(I_pisc_p_m == 0 ~ NA_real_, 
                                TRUE ~ I_pisc_p_m)) %>%
  mutate(Fp_q = quant(Fp_p_m/rel_occ),
         Fn_q = quant(Fn_p_m/rel_occ),
         Gc_q = quant(Gc_p_m/rel_occ),
         I_herb_q = quant(I_herb_p_m/rel_occ),
         I_pisc_q = quant(I_pisc_p_m/rel_occ),
         Fp_n = normalize(Fp_p_m/rel_occ),
         Fn_n = normalize(Fn_p_m/rel_occ),
         Gc_n = normalize(Gc_p_m/rel_occ),
         I_herb_n = normalize(I_herb_p_m/rel_occ),
         I_pisc_n = normalize(I_pisc_p_m/rel_occ)
         ) %>%
  as.data.frame()

rownames(data) <- data$species


p <- ggtree(tree1d, layout = "circular")

p

p1 <- 
  gheatmap(p, data[,c("Fn_q", "Fp_q","Gc_q", "I_herb_q", "I_pisc_q"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish_d("Hypsypops_rubicundus", begin = 0, 
                     labels=c("[0,0.05[", "[0.05,0.25[", "[0.25,0.75[",  "[0.75, 0.95[", "[0.95,0.99[", "[0.99,1[")) + 
  geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)

ggsave("output/plots/phyloplot_con_c_indopacific.pdf", p1, height = 12, width = 12)

p2 <- 
  gheatmap(p, data[,c("Fn_n", "Fp_n","Gc_n", "I_herb_n", "I_pisc_n"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish("Hypsypops_rubicundus", begin = 0, trans = "log1p", breaks = c(0,5,25,50,100)) + 
  geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)

ggsave("output/plots/phyloplot_con_c_indopacific_2.pdf", p2, height = 12, width = 12)

################ c pacific ###########
sub <- filter(contributions_sp_loc_occ, occurence > 1,bioregion == "c_pacific")

tree1 <- fishtree::fishtree_phylogeny(species = sub$species)

tree1t <- as_tibble(tree1)
tax$label <- gsub(" ", "_", tax$Species)
tree1t <- left_join(tree1t, tax)
#tree1t[!is.na(tree1t$Family) & tree1t$Family == "Scaridae", "Family"] <- "Labridae"

tree1d <- as.treedata(tree1t)


data <- filter(sub, species%in% tree1$tip.label) %>% 
  mutate(I_herb_p_m = case_when(I_herb_p_m == 0 ~ NA_real_, 
                                TRUE ~ I_herb_p_m),
         I_pisc_p_m = case_when(I_pisc_p_m == 0 ~ NA_real_, 
                                TRUE ~ I_pisc_p_m)) %>%
  mutate(Fp_q = quant(Fp_p_m/rel_occ),
         Fn_q = quant(Fn_p_m/rel_occ),
         Gc_q = quant(Gc_p_m/rel_occ),
         I_herb_q = quant(I_herb_p_m/rel_occ),
         I_pisc_q = quant(I_pisc_p_m/rel_occ),
         Fp_n = normalize(Fp_p_m/rel_occ),
         Fn_n = normalize(Fn_p_m/rel_occ),
         Gc_n = normalize(Gc_p_m/rel_occ),
         I_herb_n = normalize(I_herb_p_m/rel_occ),
         I_pisc_n = normalize(I_pisc_p_m/rel_occ)
  ) %>%
  as.data.frame()

rownames(data) <- data$species


p <- ggtree(tree1d, layout = "circular")

p

p1 <- 
  gheatmap(p, data[,c("Fn_q", "Fp_q","Gc_q", "I_herb_q", "I_pisc_q"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish_d("Hypsypops_rubicundus", begin = 0, 
                    labels=c("[0,0.05[", "[0.05,0.25[", "[0.25,0.75[",  "[0.75, 0.95[", "[0.95,0.99[", "[0.99,1[")) 

p1 <- p1 + geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)

ggsave("output/plots/phyloplot_con_c_pacific.pdf", p1, height = 12, width = 12)

p2 <- 
  gheatmap(p, data[,c("Fn_n", "Fp_n","Gc_n", "I_herb_n", "I_pisc_n"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish("Hypsypops_rubicundus", begin = 0, trans = "log1p", breaks = c(0,5,25,50,100)) + 
  geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)

ggsave("output/plots/phyloplot_con_c_pacific_2.pdf", p2, height = 12, width = 12)

################ e atlantic ###########
sub <- filter(contributions_sp_loc_occ, occurence > 1, bioregion == "e_atlantic")

tree1 <- fishtree::fishtree_phylogeny(species = sub$species)

tree1t <- as_tibble(tree1)
tax$label <- gsub(" ", "_", tax$Species)
tree1t <- left_join(tree1t, tax)
#tree1t[!is.na(tree1t$Family) & tree1t$Family == "Scaridae", "Family"] <- "Labridae"

tree1d <- as.treedata(tree1t)


data <- filter(sub, species%in% tree1$tip.label) %>% 
  mutate(I_herb_p_m = case_when(I_herb_p_m == 0 ~ NA_real_, 
                                TRUE ~ I_herb_p_m),
         I_pisc_p_m = case_when(I_pisc_p_m == 0 ~ NA_real_, 
                                TRUE ~ I_pisc_p_m)) %>%
  mutate(Fp_q = quant(Fp_p_m/rel_occ),
         Fn_q = quant(Fn_p_m/rel_occ),
         Gc_q = quant(Gc_p_m/rel_occ),
         I_herb_q = quant(I_herb_p_m/rel_occ),
         I_pisc_q = quant(I_pisc_p_m/rel_occ),
         Fp_n = normalize(Fp_p_m/rel_occ),
         Fn_n = normalize(Fn_p_m/rel_occ),
         Gc_n = normalize(Gc_p_m/rel_occ),
         I_herb_n = normalize(I_herb_p_m/rel_occ),
         I_pisc_n = normalize(I_pisc_p_m/rel_occ)
  ) %>%
  as.data.frame()

rownames(data) <- data$species


p <- ggtree(tree1d, layout = "circular")

p

p1 <- 
  gheatmap(p, data[,c("Fn_q", "Fp_q","Gc_q", "I_herb_q", "I_pisc_q"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish_d("Hypsypops_rubicundus", begin = 0, 
                    labels=c("[0,0.05[", "[0.05,0.25[", "[0.25,0.75[",  "[0.75, 0.95[", "[0.95,0.99[", "[0.99,1[")) 

p1 <- p1 + geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)


ggsave("output/plots/phyloplot_con_e_atlantic.pdf", p1, height = 12, width = 12)

p2 <- 
  gheatmap(p, data[,c("Fn_n", "Fp_n","Gc_n", "I_herb_n", "I_pisc_n"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish("Hypsypops_rubicundus", begin = 0, trans = "log1p", breaks = c(0,5,25,50,100)) + 
  geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)

ggsave("output/plots/phyloplot_con_e_atlantic_2.pdf", p2, height = 12, width = 12)

################ e pacific ###########
sub <- filter(contributions_sp_loc_occ, occurence > 1,  bioregion == "e_pacific")

tree1 <- fishtree::fishtree_phylogeny(species = sub$species)

tree1t <- as_tibble(tree1)
tax$label <- gsub(" ", "_", tax$Species)
tree1t <- left_join(tree1t, tax)
#tree1t[!is.na(tree1t$Family) & tree1t$Family == "Scaridae", "Family"] <- "Labridae"

tree1d <- as.treedata(tree1t)

data <- filter(sub, species%in% tree1$tip.label) %>% 
  mutate(I_herb_p_m = case_when(I_herb_p_m == 0 ~ NA_real_, 
                                TRUE ~ I_herb_p_m),
         I_pisc_p_m = case_when(I_pisc_p_m == 0 ~ NA_real_, 
                                TRUE ~ I_pisc_p_m)) %>%
  mutate(Fp_q = quant(Fp_p_m/rel_occ),
         Fn_q = quant(Fn_p_m/rel_occ),
         Gc_q = quant(Gc_p_m/rel_occ),
         I_herb_q = quant(I_herb_p_m/rel_occ),
         I_pisc_q = quant(I_pisc_p_m/rel_occ),
         Fp_n = normalize(Fp_p_m/rel_occ),
         Fn_n = normalize(Fn_p_m/rel_occ),
         Gc_n = normalize(Gc_p_m/rel_occ),
         I_herb_n = normalize(I_herb_p_m/rel_occ),
         I_pisc_n = normalize(I_pisc_p_m/rel_occ)
  ) %>%
  as.data.frame()

rownames(data) <- data$species


p <- ggtree(tree1d, layout = "circular")

p

p1 <- 
  gheatmap(p, data[,c("Fn_q", "Fp_q","Gc_q", "I_herb_q", "I_pisc_q"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish_d("Hypsypops_rubicundus", begin = 0, 
                    labels=c("[0,0.05[", "[0.05,0.25[", "[0.25,0.75[",  "[0.75, 0.95[", "[0.95,0.99[", "[0.99,1[")) 

p1 <- p1 + geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)


ggsave("output/plots/phyloplot_con_e_pacific.pdf", p1, height = 12, width = 12)

p2 <- 
  gheatmap(p, data[,c("Fn_n", "Fp_n","Gc_n", "I_herb_n", "I_pisc_n"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish("Hypsypops_rubicundus", begin = 0, trans = "log1p", breaks = c(0,5,25,50,100)) + 
  geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)

ggsave("output/plots/phyloplot_con_e_pacific_2.pdf", p2, height = 12, width = 12)

################ w atlantic ###########
sub <- filter(contributions_sp_loc_occ, occurence > 1,  bioregion == "w_atlantic")

tree1 <- fishtree::fishtree_phylogeny(species = sub$species)

tree1t <- as_tibble(tree1)
tax$label <- gsub(" ", "_", tax$Species)
tree1t <- left_join(tree1t, tax)
#tree1t[!is.na(tree1t$Family) & tree1t$Family == "Scaridae", "Family"] <- "Labridae"

tree1d <- as.treedata(tree1t)


data <- filter(sub, species%in% tree1$tip.label) %>% 
  mutate(I_herb_p_m = case_when(I_herb_p_m == 0 ~ NA_real_, 
                                TRUE ~ I_herb_p_m),
         I_pisc_p_m = case_when(I_pisc_p_m == 0 ~ NA_real_, 
                                TRUE ~ I_pisc_p_m)) %>%
  mutate(Fp_q = quant(Fp_p_m/rel_occ),
         Fn_q = quant(Fn_p_m/rel_occ),
         Gc_q = quant(Gc_p_m/rel_occ),
         I_herb_q = quant(I_herb_p_m/rel_occ),
         I_pisc_q = quant(I_pisc_p_m/rel_occ),
         Fp_n = normalize(Fp_p_m/rel_occ),
         Fn_n = normalize(Fn_p_m/rel_occ),
         Gc_n = normalize(Gc_p_m/rel_occ),
         I_herb_n = normalize(I_herb_p_m/rel_occ),
         I_pisc_n = normalize(I_pisc_p_m/rel_occ)
  ) %>%
  as.data.frame()
rownames(data) <- data$species


p <- ggtree(tree1d, layout = "circular")

p

p1 <- 
  gheatmap(p, data[,c("Fn_q", "Fp_q","Gc_q", "I_herb_q", "I_pisc_q"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish_d("Hypsypops_rubicundus", begin = 0, 
                    labels=c("[0,0.05[", "[0.05,0.25[", "[0.25,0.75[",  "[0.75, 0.95[", "[0.95,0.99[", "[0.99,1[")) 

p1 <- p1 + geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)


ggsave("output/plots/phyloplot_con_w_atlantic.pdf", p1, height = 12, width = 12)

p2 <- 
  gheatmap(p, data[,c("Fn_n", "Fp_n","Gc_n", "I_herb_n", "I_pisc_n"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish("Hypsypops_rubicundus", begin = 0, trans = "log1p", breaks = c(0,5,25,50,100)) + 
  geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)

ggsave("output/plots/phyloplot_con_w_atlantic_2.pdf", p2, height = 12, width = 12)

################ w indian ###########
sub <- filter(contributions_sp_loc_occ, occurence > 1,  bioregion == "w_indian")

tree1 <- fishtree::fishtree_phylogeny(species = sub$species)

tree1t <- as_tibble(tree1)
tax$label <- gsub(" ", "_", tax$Species)
tree1t <- left_join(tree1t, tax)
#tree1t[!is.na(tree1t$Family) & tree1t$Family == "Scaridae", "Family"] <- "Labridae"

tree1d <- as.treedata(tree1t)

data <- filter(sub, species%in% tree1$tip.label) %>% 
  mutate(I_herb_p_m = case_when(I_herb_p_m == 0 ~ NA_real_, 
                                TRUE ~ I_herb_p_m),
         I_pisc_p_m = case_when(I_pisc_p_m == 0 ~ NA_real_, 
                                TRUE ~ I_pisc_p_m)) %>%
  mutate(Fp_q = quant(Fp_p_m/rel_occ),
         Fn_q = quant(Fn_p_m/rel_occ),
         Gc_q = quant(Gc_p_m/rel_occ),
         I_herb_q = quant(I_herb_p_m/rel_occ),
         I_pisc_q = quant(I_pisc_p_m/rel_occ),
         Fp_n = normalize(Fp_p_m/rel_occ),
         Fn_n = normalize(Fn_p_m/rel_occ),
         Gc_n = normalize(Gc_p_m/rel_occ),
         I_herb_n = normalize(I_herb_p_m/rel_occ),
         I_pisc_n = normalize(I_pisc_p_m/rel_occ)
  ) %>%
  as.data.frame()

rownames(data) <- data$species


p <- ggtree(tree1d, layout = "circular")

p

p1 <- 
  gheatmap(p, data[,c("Fn_q", "Fp_q","Gc_q", "I_herb_q", "I_pisc_q"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish_d("Hypsypops_rubicundus", begin = 0, 
                    labels=c("[0,0.05[", "[0.05,0.25[", "[0.25,0.75[",  "[0.75, 0.95[", "[0.95,0.99[", "[0.99,1[")) 

p1 <- p1 + geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)

ggsave("output/plots/phyloplot_con_w_indian.pdf", p1, height = 12, width = 12)

p2 <- 
  gheatmap(p, data[,c("Fn_n", "Fp_n","Gc_n", "I_herb_n", "I_pisc_n"), drop = FALSE], 
           offset = 0.2, width=1, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish("Hypsypops_rubicundus", begin = 0, trans = "log1p", breaks = c(0,5,25,50,100)) + 
  geom_tiplab(aes(angle = angle, label = tree1d@data$Family), size = 2, offset = 150)

ggsave("output/plots/phyloplot_con_c_indopacific_2.pdf", p2, height = 12, width = 12)

ggplot(sub) +
  geom_point(aes(x = Gc_p_m, y = Fp_p_m, color = Family)) +
  geom_abline(slope = 1)+ scale_x_continuous(limits = c(0.01, 0.5))+ 
  scale_y_continuous(limits = c(0.01, 0.5))

############## vulnerability ##########

data <- as.data.frame(vulnerability) %>% drop_na(vuln_climate, vuln_fi)
tree1 <- fishtree::fishtree_phylogeny(species = data$species)

tree1t <- as_tibble(tree1)
tax$label <- gsub(" ", "_", tax$Species)
tree1t <- left_join(tree1t, tax)
#tree1t[!is.na(tree1t$Family) & tree1t$Family == "Scaridae", "Family"] <- "Labridae"

tree1d <- as.treedata(tree1t)

data <- data %>% filter(species %in% tree1$tip.label) 
rownames(data) <- data$species


p <- ggtree(tree1d, layout = "circular")

p

p1 <- 
  gheatmap(p, data[,c("vuln_fi", "vuln_climate"), drop = FALSE], 
           offset = 0.2, width=0.8, 
           colnames_position = "top", font.size=2, color = "white") + 
  scale_fill_fish("Hypsypops_rubicundus", begin = 0)+ 
  geom_tiplab(aes(angle = angle, label = tree1d@data$Species), size = 1, offset = 150)
p1

ggsave("output/plots/phyloplot_vuln.pdf", p1, height = 12, width = 12)


####### rank barplots #####
con[con$Family == "Scaridae", "Family"] <- "Labridae"

main_fam <- dplyr::select(con, species, Family, bioregion) %>%
  group_by(Family) %>%
  summarize(n())

rank <- con %>%
  filter(rel_occ > 0.001) %>%
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
  dplyr::mutate(value2 = value/max(value))

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
  dplyr::mutate(value2 = value/max(value))

freq  <- table(rank_fam$Family) %>% as.data.frame()
library(forcats)

plots <- lapply(unique(rank$bioregion), function(i){
  sub <- subset(rank, bioregion == i)%>%
    arrange(desc(rank)) %>%                          
    mutate(spname = factor(paste(species, name, sep = "__"), 
                           levels = (paste(species, name, sep = "__")))) %>%
    subset(value>0)
  print(i)
  plot <-
  ggplot(sub) +
    geom_bar(aes(y = value, x = spname, fill = Family), stat = "identity") +
    geom_text(aes(y = value*0.75, x = spname, label = species), size = 2) +
    facet_wrap(~name, scales = "free", ncol = 5) +coord_flip() +
    scale_fill_fish_d(option = "Pseudocheilinus_tetrataenia") +
    labs(y = i) +
    theme_bw()+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +theme(legend.position = "none")
  return(plot)
})

leg <- cowplot::get_legend(ggplot(rank) +
                              geom_bar(aes(y = value, x = species, fill = Family), stat = "identity") +
                              scale_fill_fish_d(option = "Pseudocheilinus_tetrataenia") +
                             theme(legend.position = "bottom"))
                              
comb <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] +  leg +
  plot_layout(nrow = 7, heights = c(1,1,1,1,1,1,0.2))
ggsave("output/plots/con_species_rank.pdf", comb, width = 8, height = 12)

ggplot(rank) +
  geom_bar(aes(y = value2, x = as.character(1/rank), fill = Family), stat = "identity") +
  geom_text(aes(y = 0.5, x = as.character(1/rank), label = species), size = 2) +
  facet_grid(bioregion~name, scales="free", space="free") +coord_flip() +
  scale_fill_fish_d(option = "Pseudocheilinus_tetrataenia") +
  labs(y = "Rescaled contribution") +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("output/plots/con_species_rank.pdf",  width = 8, height = 8)


ggplot(rank_fam) +
  geom_bar(aes(y = value2, x = as.character(1/rank), fill = Family), stat = "identity") +
  geom_text(aes(y = 0.5, x = as.character(1/rank), label = Family), size = 2) +
  facet_grid(bioregion~name, scales="free", space="free") +coord_flip() +
  scale_fill_fish_d(option = "Pseudocheilinus_tetrataenia") +
  labs(y = "Rescaled contribution") +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave("output/plots/con_family_rank.pdf",  width = 8, height = 8)


######## entropy of contribution #####
entropy <- function(prob) {
  
  k              <- length(prob)                       #number of states
  prob[prob>1/k] <- prob[prob>1/k]/(1-k) - 1/(1-k)   #state entropies
  tent           <- sum(prob)               #node entropy
  
  return(tent)
}

loadd(contributions)
loadd(herb_pisc)

summary <- contributions %>% left_join(herb_pisc$contributions_herb_pisc) %>%
  group_by(bioregion, transect_id, nspec) %>%
  summarise_at(vars(Fn_p, Fp_p, Gc_p, I_pisc_p, I_herb_p), entropy) %>%
  pivot_longer(cols = ends_with("_p")) %>%
  filter(nspec >10) %>%
  filter(value>0)

ggplot(summary) +
  geom_point(aes(x = log(nspec), y = (value), color = name), size = 0.5) +
  geom_smooth(aes(x = log(nspec), y = (value), color = name), size = 0.5, method = "lm") +
  facet_wrap(~bioregion, scales = "free")

#fit <- brm(value ~ nspec:name + (1|bioregion), data = summary, family = "beta")  


####### fig 1 ######


## load data 
flux <- readd(summary_transect)
hp <- readd(herb_pisc)$summary_herb_pisc
sst <- read.csv("data/avSst.csv")

flux <- left_join(flux, hp) %>% left_join(sst)

flux <- flux %>% mutate(
  Fc = Fc/biomass_tot,
  Fn = Fn/biomass_tot,
  Fp = Fp/biomass_tot,
  Gc = Gc/biomass_tot,
  Wc = Wc/biomass_tot,
  Wn = Wn/biomass_tot,
  Wp = Wp/biomass_tot,
  I_herb = I_herb/biomass_tot,
  I_pisc = I_pisc/biomass_tot)

#multi
get_scores <- function(proc){
  probs <- seq(0,1,0.01)
  quantiles <- quantile(proc, prob = probs)
  quant <- as.numeric(factor(findInterval(proc, quantiles)))/100
  return(quant)
}
test <- select(flux, biomass_tot, Fc, Fn, Fp, Gc, Wc, Wn, Wp, I_herb, I_pisc)  %>% 
  dplyr::select(-biomass_tot)  %>%
  mutate_all(function(x){x/flux$biomass_tot})%>%
  mutate_all(function(x){get_scores(x)})

test[test<0] <-0
test[test>1] <-1

tt <- get_scores(test$Fc)
hist(tt)

library(ggplot2)
ggplot(aes(x = Gc, y = Fp), data = test) + 
  geom_point() +
  xlim(c(0, 1)) +ylim(c(0, 1)) +
  geom_vline(xintercept = quantile(test$Gc, 0.75)) +
  geom_hline(yintercept = quantile(test$Fp, 0.75)) +
  geom_vline(xintercept = quantile(test$Gc, 0.5)) +
  geom_hline(yintercept = quantile(test$Fp, 0.5))


m <-  cor(test[,c(2,3,4,7,8,9)])
m <- 1-m
we <-  rowSums((m))
we <- we/sum(we)

test <- test %>%
  rowwise() %>%
  do({
    result <- (as_data_frame(.))
    result$var = var(c(result$Fn, result$Fp, result$Gc, result$Wp, result$I_herb, result$I_pisc))
    result$multi1 <- mean(c(result$Fn, result$Fp, result$Gc, result$Wp, result$I_herb, result$I_pisc))
    result$multi1w <- weighted.mean(c(result$Fn, result$Fp, result$Gc, result$Wp, result$I_herb, result$I_pisc), w = we)
    result
  }) 

flux <- flux %>% mutate(
  Fc_r = residuals(lm(log(Fc) ~ log(biomass_tot))),
  Fn_r = residuals(lm(log(Fn) ~  log(biomass_tot))),
  Fp_r = residuals(lm(log(Fp) ~  log(biomass_tot))),
  Gc_r = residuals(lm(log(Gc) ~  log(biomass_tot))),
  Wc_r = residuals(lm(log(Wc) ~  log(biomass_tot))),
  Wn_r = residuals(lm(log(Wn) ~  log(biomass_tot))),
  Wp_r = residuals(lm(log(Wp) ~  log(biomass_tot))),
  I_herb_r = residuals(lm(log1p(I_herb) ~  log(biomass_tot))),
  I_pisc_r = residuals(lm(log1p(I_pisc) ~  log(biomass_tot))))



cor(select(flux, troph_m ,
           imm_m, nspec, troph_q3, size_q3, troph_wci, size_wci, imm_wci, imm_q3, imm_q1) )

fsel <- select(flux, biomass_tot, Fc, Fn, Fp, Ic, In, Ip, Gc, Wc, Wn, Wp, I_herb, I_pisc, nspec) %>%
  mutate_all(function(x){x/flux$biomass_tot}) %>% select(-biomass_tot)

fluxs_sel <-  mutate(flux, logbiomass = log(biomass_tot), logabu = log(abu_tot)) %>%
  select(nspec, mean, troph_m, troph_q3, size_q3, size_m,
         imm_m, imm_q1) #%>%
 # mutate_all(log1p)
library(corrplot)
corrplot(cor(fluxs_sel), method="number")
hist(fluxs_sel$size_m)


pca <- prcomp(fluxs_sel, center = TRUE,scale. = TRUE)
summary(pca)


dpca <- pca$x %>% as.data.frame() %>% cbind(flux)

rpca <- pca$rotation %>% as.data.frame()

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circ <- circleFun(c(0,0),2,npoints = 500)

arrows_12 <- 
  ggplot() +
  geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), data = rpca, arrow = arrow(length=unit(0.30,"cm"))) +
  geom_text(data = rpca, aes(x =  PC1, 
                             y = PC2, 
                             label = rownames(rpca)), size = 4, vjust = 0, colour="black",
            angle = (180/pi) * atan(rpca$PC2/rpca$PC1), hjust = (1 - 2 * sign(rpca$PC1)) / 2) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  labs(x = "PC1 (23.41%)", y = "PC2 (20.34%)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))
arrows_12


arrows_34 <- 
  ggplot() +
  geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_segment(aes(x = 0, y = 0, xend = PC3, yend = PC4), data = rpca, arrow = arrow(length=unit(0.30,"cm"))) +
  geom_text(data = rpca, aes(x =  PC3, 
                             y = PC4, 
                             label = rownames(rpca)), size = 4, vjust = 0, colour="black",
            angle = (180/pi) * atan(rpca$PC4/rpca$PC3), hjust = (1 - 2 * sign(rpca$PC3)) / 2) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  labs(x = "PC3 (15.55%)", y = "PC4 (13.8%)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))
arrows_34

library(gghighlight)
library(RColorBrewer)
library(ggnewscale)

library(fishualize)
fishualize("Hypsypops_rubicundus", n = 5, begin = 0, end = 0.8)


# Add bivariate density for each point
dpca$density <- fields::interp.surface(
  MASS::kde2d(dpca$PC1, dpca$PC2), dpca[,c("PC1", "PC2")])

dpca$density2 <- fields::interp.surface(
  MASS::kde2d(dpca$PC1, dpca$PC2), dpca[,c("PC3", "PC4")])

dpca$multi <- test$multi1

pca_plot <- function(proc, title){
  
  probs <- c(0.05, 0.25, 0.75, 0.95, 0.99)
  quantiles <- quantile(proc, prob = probs)
  dpca$quant <- factor(findInterval(proc, quantiles))
  
  p1 <-
    ggplot(dpca) +
    geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_point(aes(x = PC1, y = PC2, color = quant, alpha = 1/density), shape = 16, size = 2) +
    geom_point(data = dpca[dpca$quant == "9",], aes(x = PC1, y = PC2, color = quant, alpha = 1/density), shape = 16, size = 2) +
    scale_color_fish_d("Hypsypops_rubicundus", begin = 0, labels=c("[0,0.05[", "[0.05,0.25[", "[0.25,0.75[",  "[0.75, 0.95[", "[0.95,1[")) +
    labs(x = "PC1", y = "PC2", title = title) +
    scale_alpha(range = c(.25, .6)) +
    theme_minimal() +
    guides(alpha = "none") +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(fill= "transparent"), legend.position = "none")
  
  p2 <-
    ggplot(dpca) +
    geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_point( aes(x = PC3, y = PC4, color = quant , alpha = 1/density), shape = 16, size = 2) +
    geom_point(data = dpca[dpca$quant == "9",], aes(x = PC3, y = PC4, color = quant, alpha = 1/density2), shape = 16, size = 2) +
    scale_alpha(range = c(.25, .6) ) +
    scale_color_fish_d("Hypsypops_rubicundus", begin = 0, labels=c("[0,0.05[", "[0.05,0.10[", "[0.10,0.25[", "[0.25,0.50[", "[0.50,0.75[", "[0.75, 0.90[", "[0.90,0.95[", "[0.95,1[")) +
    labs(x = "PC3", y = "PC4", title = title) +
    theme_minimal() +
    guides(alpha = "none") +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(fill= "transparent"), legend.position = "none")
  
  return(list(p1, p2))
  
}



dpca$multi <- test$multi1w
p1 <- pca_plot(dpca$Fc_r, "Fc")
p2 <- pca_plot(dpca$Fn_r, "Fn")
p3 <- pca_plot(dpca$Fp_r, "Fp")
p4 <- pca_plot(dpca$Wc_r, "Wc")
p5 <- pca_plot(dpca$Wn_r, "Wn")
p6 <- pca_plot(dpca$Wp_r, "Wp")
p7 <- pca_plot(dpca$Gc_r, "Gc")
p8 <- pca_plot(dpca$I_herb_r, "I_herb")
p9 <- pca_plot(dpca$I_pisc_r, "I_pisc")
p10 <- pca_plot(dpca$multi_r, "multi")

p7[[1]]
p2[[1]]
p3[[1]]

p10[[1]]
p10[[2]]

p2[[2]]
p3[[2]]
p7[[2]]
arrows_12
arrows_34

library(patchwork)
layout <- "
AABB
AABB
CDIJ
EFKL
GHMN
"
arrows_12 + arrows_34 + p2[[1]]+ p3[[1]]+ p7[[1]]+ p8[[1]]+ p9[[1]]+ p10[[1]] +
  p2[[2]]+ p3[[2]]+ p7[[2]]+ p8[[2]]+ p9[[2]]+ p10[[2]] +
  plot_layout(design = layout) 
ggsave("output/plots/pca_plot.pdf" ,height = 15, width = 12)

####### ppt #######


pca_plot <- function(proc, title){
  
  probs <- c(0.05, 0.25, 0.75, 0.95, 0.99)
  quantiles <- quantile(proc, prob = probs)
  dpca$quant <- factor(findInterval(proc, quantiles))
  
  
  
  p1 <-
    ggplot(dpca) +
    geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_point(aes(x = PC1, y = PC2, color = quant, alpha = 1/density), shape = 16, size = 2) +
    geom_point(data = dpca[dpca$quant == "9",], aes(x = PC1, y = PC2, color = quant, alpha = 1/density), shape = 16, size = 2) +
    scale_color_fish_d("Hypsypops_rubicundus", begin = 0, labels=c("[0,0.05[", "[0.05,0.10[", "[0.10,0.25[", "[0.25,0.50[", "[0.50,0.75[", "[0.75, 0.90[", "[0.90,0.95[", "[0.95,1[")) +
    labs(x = "PC1", y = "PC2", title = title) +
    scale_alpha(range = c(.25, .75)) +
    theme_bw() +
    guides(alpha = "none") +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(fill= "transparent"), legend.position = "none") +
    theme(panel.background = element_rect(color = "#041e3d", fill =  "#041e3d"), 
          axis.text = element_text(color = "white", size = 13), text = element_text(color = "white", size = 12), 
          plot.background = element_rect(color = "#041e3d", fill =  "#041e3d"), axis.title = element_text(size = 18),
          legend.position = "none"  )
  
  p2 <-
    ggplot(dpca) +
    geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_point( aes(x = PC3, y = PC4, color = quant , alpha = 1/density), shape = 16, size = 2) +
    geom_point(data = dpca[dpca$quant == "9",], aes(x = PC3, y = PC4, color = quant, alpha = 1/density2), shape = 16, size = 2) +
    scale_alpha(range = c(.25, .75) ) +
    scale_color_fish_d("Hypsypops_rubicundus", begin = 0, labels=c("[0,0.05[", "[0.05,0.10[", "[0.10,0.25[", "[0.25,0.50[", "[0.50,0.75[", "[0.75, 0.90[", "[0.90,0.95[", "[0.95,1[")) +
    labs(x = "PC3", y = "PC4", title = title) +
    theme_bw() +
    guides(alpha = "none") +
    theme(panel.grid = element_blank(), 
          panel.border = element_rect(fill= "transparent"), legend.position = "none") +
    theme(panel.background = element_rect(color = "#041e3d", fill =  "#041e3d"), 
          axis.text = element_text(color = "white", size = 13), text = element_text(color = "white", size = 12), 
          plot.background = element_rect(color = "#041e3d", fill =  "#041e3d"), axis.title = element_text(size = 18),
          legend.position = "none"  )
  
  return(list(p1, p2))
  
}


arrows_12 <- 
  ggplot() +
  geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_segment(color = "white", aes(x = 0, y = 0, xend = PC1, yend = PC2), data = rpca, arrow = arrow(length=unit(0.30,"cm"))) +
  geom_text(data = rpca, aes(x =  PC1, 
                             y = PC2, 
                             label = rownames(rpca)), size = 4, vjust = 0, colour="white",
            angle = (180/pi) * atan(rpca$PC2/rpca$PC1), hjust = (1 - 2 * sign(rpca$PC1)) / 2) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  labs(x = "PC1 (23.41%)", y = "PC2 (20.34%)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(panel.background = element_rect(color = "#041e3d", fill =  "#041e3d"), 
        axis.text = element_text(color = "white", size = 13), text = element_text(color = "white", size = 12), 
        plot.background = element_rect(color = "#041e3d", fill =  "#041e3d"), axis.title = element_text(size = 18),
        legend.position = "none"  )
arrows_12


arrows_34 <- 
  ggplot() +
  geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_segment(color = "white", aes(x = 0, y = 0, xend = PC3, yend = PC4), data = rpca, arrow = arrow(length=unit(0.30,"cm"))) +
  geom_text(data = rpca, aes(x =  PC3, 
                             y = PC4, 
                             label = rownames(rpca)), size = 4, vjust = 0, colour="white",
            angle = (180/pi) * atan(rpca$PC4/rpca$PC3), hjust = (1 - 2 * sign(rpca$PC3)) / 2) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  labs(x = "PC3 (15.55%)", y = "PC4 (13.8%)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(panel.background = element_rect(color = "#041e3d", fill =  "#041e3d"), 
        axis.text = element_text(color = "white", size = 13), text = element_text(color = "white", size = 12), 
        plot.background = element_rect(color = "#041e3d", fill =  "#041e3d"), axis.title = element_text(size = 18),
        legend.position = "none"  )
arrows_34


dpca$multi <- test$multi1w
p1 <- pca_plot(dpca$Fc_r, "Fc")
p2 <- pca_plot(dpca$Fn_r, "Fn")
p3 <- pca_plot(dpca$Fp_r, "Fp")
p4 <- pca_plot(dpca$Wc_r, "Wc")
p5 <- pca_plot(dpca$Wn_r, "Wn")
p6 <- pca_plot(dpca$Wp_r, "Wp")
p7 <- pca_plot(dpca$Gc_r, "Gc")
p8 <- pca_plot(dpca$I_herb_r, "I_herb")
p9 <- pca_plot(dpca$I_pisc_r, "I_pisc")
p10 <- pca_plot(dpca$multi, "multi")
p10


library(patchwork)
layout <- "
AABB
AABB
CDIJ
EFKL
GHMN
"
arrows_12 + arrows_34 + p2[[1]]+ p3[[1]]+ p7[[1]]+ p8[[1]]+ p9[[1]]+ p10[[1]] +
  p2[[2]]+ p3[[2]]+ p7[[2]]+ p8[[2]]+ p9[[2]]+ p10[[2]] +
  plot_layout(design = layout) 
ggsave("output/plots/ppt_pca_plot.pdf" ,height = 15, width = 12)

layout <- "
ABC
DEF
"
p2[[1]]+ p3[[1]]+ p7[[1]]+ p8[[1]]+ p9[[1]]+ p10[[1]] +
  plot_layout(design = layout) 
ggsave("output/plots/ppt_pca_plot_clouds.pdf")

layout <- "
AB
"
arrows_12 + arrows_34 +  plot_layout(design = layout) 
ggsave("output/plots/ppt_pca_plot_arrows.pdf" )


########## human impact #########

hum <- read_csv("data/humanPopulation.csv")

hum <- select(flux, sites, bioregion) %>% unique() %>%
  left_join(hum) %>% right_join(select(flux, sites, transect_id))

summary(hum)

hum <- lapply(unique(flux$bioregion), function(x){
  sub <- filter(hum, bioregion == x) %>%
    dplyr::mutate(hum_centile = ntile(tt_Nmarket, n = 100)) %>%
    dplyr::select(bioregion, sites, hum_centile, pop_Npop, transect_id, tt_Npop, grav_tot, tt_Nmarket)
  return(sub)
}) %>% plyr::ldply() 

summary(hum)

dpca2 <- left_join(dpca, hum) #%>% filter(!is.na(tt_Nmarket))



pca_plot((dpca2$tt_Nmarket), "test")
pca_plot((1/dpca2$pop_Npop), "test")
pca_plot((dpca2$biomass_tot), "test")


summary(dpca2$tt_Nmarket)

ggplot(dpca2) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_point(aes(x = PC1, y = PC2, color = hum_centile, alpha = 1/density), shape = 16, size = 3) +
  geom_point(data = dpca2[dpca2$hum_centile < 3,], aes(shape = bioregion, x = PC1, y = PC2, color = hum_centile, alpha = 1/density), size = 3) +
  geom_point(data = dpca2[dpca2$hum_centile > 97,], aes(shape = bioregion, x = PC1, y = PC2, color = hum_centile, alpha = 1/density), size = 3) +
  scale_color_fish("Hypsypops_rubicundus", begin = 0) + 
  # facet_wrap(~bioregion) +
  labs(x = "PC1", y = "PC2", title = "human") +
  scale_alpha(range = c(.25, .75)) +
  theme_bw() +
  guides(alpha = "none") +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) 

ggplot(dpca2) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  # geom_point(aes(x = PC1, y = PC2, color = hum_centile, alpha = 1/density), shape = 16, size = 3) +
  geom_point(data = dpca2[dpca2$hum_centile == 1,], aes(shape = bioregion, x = PC3, y = PC4, color = hum_centile, alpha = 1/density), size = 3) +
  geom_point(data = dpca2[dpca2$hum_centile == 100,], aes(shape = bioregion, x = PC3, y = PC4, color = hum_centile, alpha = 1/density), size = 3) +
  scale_color_fish("Hypsypops_rubicundus", begin = 0) + 
  facet_wrap(~bioregion) +
  labs(x = "PC1", y = "PC2", title = "human") +
  scale_alpha(range = c(.25, .75)) +
  theme_bw() +
  guides(alpha = "none") +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) 



###### 
test <- contributions_sp_loc_occ# %>% filter(rel_occ > 0.001)
ggplot(test) +
  geom_histogram(aes(x = Fn_pphi_m * Fn_p_m)) +
  scale_x_continuous(trans = "log10") +
  facet_wrap(~bioregion, scales = "free_y")

ggplot(test) +
  geom_point(aes(x = Gc_p_m, y = Gc_pphi_m , size = rel_occ, color = Gc_pphi_m * Gc_p_m), alpha = 0.4) +
  scale_x_continuous(trans = "log10", limits = c(0.001, 1)) +
  scale_color_fish(trans = "log10") +
  scale_y_continuous(trans = "log10", limits = c(0.01, 1000)) +
  geom_hline(aes(yintercept = 50)) +
  geom_hline(aes(yintercept = 10)) +
  geom_hline(aes(yintercept = 5)) +
  geom_hline(aes(yintercept = 1)) +
  geom_vline(aes(xintercept = 0.1)) +
  geom_vline(aes(xintercept = 0.05)) +
  theme_bw()


phi = 5
mu = 0.02
a = phi * mu
b = phi * (1 - mu)
ggplot() + geom_density(aes(x = (rbeta(1000, a, b)))) 

##########dominance##########
dom <- readd(dominance)


hist(dom$I_pisc_d)

################keystoneness
key <- inner_join(con, key)

#fit_k_Fp <- brm(ks_Fp ~ 1 + (1|bioregion) + (1|sites), data = key)
hist(key$ks_I_herb)
hist(key$ks_I_pisc)
hist(key$ks_Gc)
hist(key$ks_Fn)
hist(key$ks_Fp)
plot(key$ks_Fp, key$ks_Fn)

keys <- key %>% left_join(flux) %>% group_by(bioregion, sites, locality) %>%
  summarise_if(is.numeric, median, na.rm = TRUE) %>%
  group_by(bioregion, locality) %>%
  summarise_if(is.numeric, median, na.rm = TRUE) 

ggplot(keys)+
  geom_point(aes(x = ks_I_herb, y = ks_Fp))
cor(keys$ks_Fn, keys$ks_Fp)
cor(keys$ks_Gc, keys$ks_Fp)
cor(keys$ks_Gc, keys$ks_Fn)




library("rnaturalearth")
library("rnaturalearthdata")
library(rgeos)
library(sf)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)


ggplot(keys) + mapworld + coord_quickmap() +
  geom_point(aes(x = lon, y = lat, color = ks_I_herb), alpha = 0.5) +
  scale_color_fish(option = "Hypsypops_rubicundus", limits = c(0,1)) 

theme_worldmap <- function(){
  theme_bw()+
    theme(legend.position = "none", 
          panel.grid = element_blank(),
          panel.border = element_rect(fill = NA, color = "grey", size = 1),
          axis.title = element_blank(),
          axis.text = element_text(size = 6),
          plot.title = element_text(size = 8, hjust = 0, face = "bold"), 
          plot.margin = unit(c(0.1,0.1,0.1,0.1), units = , "cm")
          ) 
}


breaks <- c(0, 0.6, 0.7, 0.8, 0.9, 1)

keys$cut <- (cut(keys$ks_Gc, breaks = breaks))
summary(keys$cut)

a <-
  ggplot(keys) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = cut, size = (Gc)),  alpha = 0.7,
             position = position_jitter(0,0), shape = 16) +
  scale_color_fish_d(option = "Hypsypops_rubicundus",  
                     name = "Keystoneness", drop = FALSE)  +
  coord_sf(ylim = c(-35, 35), expand = FALSE) +
  labs(title = "Production") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
  scale_size_continuous(range = c(0.5, 3)) +
  theme_worldmap() 
a

keys$cut <- cut(keys$ks_Fn, breaks = breaks)

b <-
  ggplot(keys) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = cut, size = (Fn)), alpha = 0.7,
             position = position_jitter(0,0), shape = 16) +
  scale_color_fish_d(option = "Hypsypops_rubicundus", 
                     name = "Keystoneness", drop = FALSE)  +
  coord_sf(ylim = c(-35, 35), expand = FALSE) +
  labs(title = "N cycling") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
  scale_size_continuous(range = c(0.5, 3)) +
  theme_worldmap() 
b
keys$cut <- cut(keys$ks_Fp, breaks = breaks)

c <-
  ggplot(keys) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = cut, size = (Fp)), alpha = 0.7,
             position = position_jitter(0,0), shape = 16) +
  scale_color_fish_d(option = "Hypsypops_rubicundus", 
                     name = "Keystoneness", drop = FALSE)  +
  coord_sf(ylim = c(-35, 35), expand = FALSE) +
  labs(title = "P cycling") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
  scale_size_continuous(range = c(0.5, 3)) +
  theme_worldmap() 
c
keys$cut <- cut(keys$ks_I_pisc, breaks =  breaks)

d <-
  ggplot(keys) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = cut, size = (I_pisc)), alpha = 0.7,
             position = position_jitter(0,0), shape = 16) +
  scale_color_fish_d(option = "Hypsypops_rubicundus", 
                     name = "Keystoneness", drop = FALSE)  +
  coord_sf(ylim = c(-35, 35), expand = FALSE) +
  labs(title = "Piscivory") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
  scale_size_continuous(range = c(0.5, 3)) +
  theme_worldmap() 
d
keys$cut <- cut(keys$ks_I_herb, breaks = breaks)

e <-
  ggplot(keys) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = cut, size = (I_herb)), alpha = 0.7,
             position = position_jitter(0,0), shape = 16) +
  scale_color_fish_d(option = "Hypsypops_rubicundus",  
                     name = "Degree of dominance", na.translate = F, drop = FALSE)  +
  coord_sf(ylim = c(-35, 35), expand = FALSE) +
  labs(title = "Herbivory") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  theme_worldmap() + theme(legend.position = "bottom", 
                           legend.title = element_text(size = 8),
                           legend.text = element_text(size = 8))

e 

multimap <-
a + b + c + d + e + plot_layout(ncol = 1)
#ggsave("output/plots/multimap.tiff", multimap, width = 8, height = 12)


key_long <- pivot_longer(key, cols = c(ks_Gc, ks_Fn, ks_Fp, ks_I_pisc, ks_I_herb))
ggplot(key_long) +
  geom_density(aes(x = value, fill = name, color = name), alpha = 0.5) +
  facet_grid(name~bioregion, scales = "free_y")

####### fig3 part 2######

con <- left_join(contributions, herb_pisc$contributions_herb_pisc)

con[con$Family == "Scaridae", "Family"] <- "Labridae"

main_fam <- dplyr::select(con, species, Family, bioregion) %>%
  group_by(Family) %>%
  summarize(n())

library(forcats)

tax2 <- select(tax, Species, Family) %>%
  mutate(species = gsub(" ", "_", Species))
tax2$Family[tax2$Family == "Scaridae"] <- "Labridae"

con_fam <- left_join(contributions, herb_pisc$contributions_herb_pisc) %>%
  left_join(tax2) %>%
  filter(!is.na(Family)) %>%
  group_by(bioregion, locality, sites, transect_id, Family) %>%
  summarize_at(.vars = vars(ends_with("_p")), sum, na.rm = TRUE) %>%
  group_by(bioregion, locality, sites, Family) %>%
  summarize_at(.vars = vars(ends_with("_p")), median, na.rm = TRUE) %>%
  group_by(bioregion, locality, Family) %>%
  summarize_at(.vars = vars(ends_with("_p")), median, na.rm = TRUE) %>%
  group_by(bioregion, Family) %>%
  summarize_at(.vars = vars(ends_with("_p")), median, na.rm = TRUE) %>%
  ungroup() 

famrank <- con_fam %>%
  group_by(Family, bioregion) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  group_by(bioregion) %>%
  mutate(
    Fn_r = dense_rank(desc(Fn_p)),
    Fp_r = dense_rank(desc(Fp_p)),
    Gc_r = dense_rank(desc(Gc_p)),
    I_herb_r = dense_rank(desc(I_herb_p)),
    I_pisc_r = dense_rank(desc(I_pisc_p))
  ) %>% rowwise() %>%
  mutate(rank = min(Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r)) %>%
  group_by(Family) %>%
  summarize(rank = quantile(rank, 0.25)) %>%
  filter(rank < 3) %>%
  filter(!Family == "Bothidae")
  
unique(famrank$Family)
  
fishualize(option = "Callanthias_australis") 
pal <- fish(option = "Callanthias_australis",n = 5) 

theme_heat <- function(){
  theme_bw() +
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(), 
          panel.grid = element_blank(),
          panel.border = element_blank(), 
          panel.background = element_rect(fill = "white"),
          legend.position = "none", legend.title = element_blank(),
          legend.text = element_text(size = 6),
          #axis.text.x = element_text(angle = 60, hjust = 1),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), units = , "cm"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6)
    ) 
}

a2 <-
ggplot(filter(con_fam, Family %in% unique(famrank$Family))) +
  geom_tile(aes(y = Family, x = bioregion, alpha = Gc_p), fill = pal[3]) +
  scale_y_discrete(position = "right") +
  scale_alpha_continuous(range = c(0,1)) +
  theme_heat()+
  theme(axis.text.x = element_text(angle = 90, size = 6))
a2
b2 <-
  ggplot(filter(con_fam, Family %in% unique(famrank$Family))) +
  geom_tile(aes(y = Family, x = bioregion, alpha = Fn_p), fill = pal[1]) +
  scale_alpha_continuous(range = c(0,1)) +
  scale_y_discrete(position = "right") +
  theme_heat()+
  theme(axis.text.x = element_text(angle = 90, size = 6))
b2
c2 <-
  ggplot(filter(con_fam, Family %in% unique(famrank$Family))) +
  geom_tile(aes(y = Family, x = bioregion, alpha = Fp_p), fill = pal[2]) +
  scale_alpha_continuous(range = c(0,1)) +
  scale_y_discrete(position = "right") +
  theme_heat()+
  theme(axis.text.x = element_text(angle = 90, size = 6))
c2
d2 <-
  ggplot(filter(con_fam, Family %in% unique(famrank$Family))) +
  geom_tile(aes(y = Family, x = bioregion, alpha = I_pisc_p), fill = pal[5]) +
  scale_y_discrete(position = "right") +
  scale_alpha_continuous(range = c(0,1)) +
  theme_heat()+
  theme(axis.text.x = element_text(angle = 90, size = 6))
d2
e2 <-
  ggplot(filter(con_fam, Family %in% unique(famrank$Family))) +
  geom_tile(aes(y = Family, x = bioregion, alpha = I_herb_p), fill = pal[4]) +
  scale_y_discrete(position = "right") +
  scale_alpha_continuous(range = c(0,1)) +
  theme_heat() +
  theme(axis.text.x = element_text(angle = 90, size = 6)) 
e2

p1 <- a +  b + c + d + e + plot_layout(ncol = 1)
p2 <- a2 + b2 + c2 + d2 + e2 + plot_layout(ncol = 1)
p2
p <- a + a2 + b + b2 + c + c2 + d  + d2 + e + e2 + plot_layout(ncol = 2, widths = c(7, 1))
p

ggsave("output/plots/multimap.tiff", p, height = 11, width = 8)


####### figure 4 #######

vuln <- as.data.frame(vulnerability) %>% drop_na(vuln_climate, vuln_fi) %>%
  mutate(Family = case_when(Family == "Scaridae" ~ "Labridae", 
                            TRUE ~ as.character(Family))) %>%
  dplyr::left_join(contributions_sp_loc_occ) %>%
  group_by(species, Family) %>%
  dplyr::summarise_if(is.numeric, max) %>%
  ungroup() %>%
  mutate(Gc = case_when(Gc_p_m > median(Gc_p_m, na.rm = TRUE) ~ "high", TRUE ~ "low"),
         Fn = case_when(Fn_p_m > median(Fn_p_m, na.rm = TRUE) ~ "high", TRUE ~ "low"),
         Fp = case_when(Fp_p_m > median(Fp_p_m, na.rm = TRUE) ~ "high", TRUE ~ "low"),
         I_pisc = case_when(I_pisc_p_m > median(I_pisc_p_m[I_pisc_p_m>0], na.rm = TRUE) ~ "high", 
                              I_pisc_p_m <= median(I_pisc_p_m[I_pisc_p_m>0], na.rm = TRUE) & I_pisc_p_m > 0 ~ "low",
                              TRUE ~ NA_character_),
         I_herb = case_when(I_herb_p_m > median(I_herb_p_m[I_herb_p_m>0], na.rm = TRUE) ~ "high", 
                              I_herb_p_m <= median(I_herb_p_m[I_herb_p_m>0], na.rm = TRUE) & I_herb_p_m > 0 ~ "low",
                              TRUE ~ NA_character_),
         vulncat = case_when(
           (vuln_climate > median(vuln_climate) & vuln_fi > median(vuln_fi)) ~ "vuln_high",
           (vuln_climate > median(vuln_climate) & vuln_fi <= median(vuln_fi)) ~ "vuln_cli",
           (vuln_climate <= median(vuln_climate) & vuln_fi > median(vuln_fi)) ~ "vuln_fi",
           (vuln_climate <= median(vuln_climate) & vuln_fi <= median(vuln_fi)) ~ "vuln_low")) %>%
  dplyr::select(species, Family, vulncat, Gc, Fn, Fp, I_pisc, I_herb) %>%
  pivot_longer(cols = c(Gc, Fn, Fp, I_pisc, I_herb), names_to = "name", values_to = "ep") %>%
  group_by(name) %>%
  mutate(n_spec_tot = sum(!is.na(ep))) %>%
  dplyr::group_by(vulncat, name, n_spec_tot) %>%
  dplyr::summarize(ep_n = sum(ep == "high", na.rm = TRUE), 
                   n = (sum(ep == "high", na.rm = TRUE) + sum(ep == "low", na.rm = TRUE))) %>%
  mutate(ep_prop = ep_n/n_spec_tot, n_prop = n/n_spec_tot) %>%
  mutate(perc = 100 * ep_prop/n_prop)

ggplot(vuln) +
  geom_bar(aes(x = name, y = n_prop, fill = name), stat = "identity", alpha = 0.5) +
  geom_bar(aes(x = name, y = ep_prop, fill = name), stat = "identity") +
  geom_text(aes(x = name, y = n_prop + 0.05, label = paste0(round(perc), "%"))) +
  scale_fill_fish_d(option = "Callanthias_australis",
                    labels = c("N excretion", "P excretion",
                               "Production", "Herbivory", "Piscivory"),
                    name = "Function") +
  facet_wrap(~vulncat) +
  theme_void()

ggsave("output/plots/vuln_bars.png", width = 6, height = 8)

# case 2
# idea: if con > 1/N !!
vuln <- con %>%
  group_by(sites, species) %>%
  dplyr::summarise_if(is.numeric, median, na.rm = TRUE) %>%
  group_by(species) %>%
  dplyr::summarise_if(is.numeric, max) %>% left_join(vulnerability)%>% 
  drop_na(vuln_climate, vuln_fi) %>%
  ungroup() %>%
  mutate(Gc = case_when(Gc_p > 0.05 ~ "high", TRUE ~ "low"),
                Fn = case_when(Fn_p > 0.05 ~ "high", TRUE ~ "low"),
                Fp = case_when(Fp_p > 0.05 ~ "high", TRUE ~ "low"),
                I_pisc = case_when(I_pisc_p > 0.1 ~ "high",
                                   I_pisc_p <= 0.1 & I_pisc_p > 0 ~ "low",
                                   TRUE ~ NA_character_),
                I_herb = case_when(I_herb_p > 0.1 ~ "high",
                                   I_herb_p <= 0.1 & I_herb_p > 0 ~ "low",
                                   TRUE ~ NA_character_),
         vulncat = case_when(
           (vuln_climate > median(vuln_climate) & vuln_fi > median(vuln_fi)) ~ "vuln_high",
           (vuln_climate > median(vuln_climate) & vuln_fi <= median(vuln_fi)) ~ "vuln_cli",
           (vuln_climate <= median(vuln_climate) & vuln_fi > median(vuln_fi)) ~ "vuln_fi",
           (vuln_climate <= median(vuln_climate) & vuln_fi <= median(vuln_fi)) ~ "vuln_low")) %>%
  dplyr::select(species, Family, vulncat, Gc, Fn, Fp, I_pisc, I_herb) %>%
  pivot_longer(cols = c(Gc, Fn, Fp, I_pisc, I_herb), names_to = "name", values_to = "ep") %>%
  group_by(name) %>%
  mutate(n_spec_tot = sum(!is.na(ep))) %>%
  dplyr::group_by(vulncat, name, n_spec_tot) %>%
  dplyr::summarize(ep_n = sum(ep == "high", na.rm = TRUE), 
                   n = (sum(ep == "high", na.rm = TRUE) + sum(ep == "low", na.rm = TRUE))) %>%
  mutate(ep_prop = ep_n/n_spec_tot, n_prop = n/n_spec_tot) %>%
  mutate(perc = 100 * ep_prop/n_prop)


ggplot(vuln) +
  geom_bar(aes(x = name, y = n_prop, fill = name), stat = "identity", alpha = 0.5) +
  geom_bar(aes(x = name, y = ep_prop, fill = name), stat = "identity") +
  geom_text(aes(x = name, y = n_prop + 0.05, label = paste0(round(perc), "%"))) +
  scale_fill_fish_d(option = "Callanthias_australis",
                    labels = c("N excretion", "P excretion",
                               "Production", "Herbivory", "Piscivory"),
                    name = "Function") +
  facet_wrap(~vulncat) +
  theme_void()

ggsave("output/plots/vuln_bars_option2.png", width = 6, height = 8)





vdata <- vulnerability %>% drop_na(vuln_climate, vuln_fi) 

fit_vuln <- brm(vuln_fi ~ vuln_climate + I(vuln_climate^2), data = vdata,
                family = "lognormal")

fitted_vuln <- fitted(fit_vuln, 
                      newdata = data.frame(vuln_climate = seq(0.01,1, 0.01))) %>%
  as_tibble() %>%
  mutate(vuln_climate = seq(0.01,1, 0.01))

ggplot(vdata) +
  # geom_rect(xmin = median(vdata$vuln_climate), ymin = median(vdata$vuln_fi),
  #           xmax = 1, ymax = 1, 
  #           fill = "red", alpha = 0.1) +
  # geom_rect(xmin = median(vdata$vuln_climate), ymin = -0.47,
  #           xmax = 1, ymax = median(vdata$vuln_fi), 
  #           fill = "orange", alpha = 0.1) +
  # geom_rect(xmin = -0.35, ymin = median(vdata$vuln_fi),
  #           xmax = median(vdata$vuln_climate), ymax = 1, 
  #           fill = "orange", alpha = 0.1) +
  # geom_rect(xmin = -0.35, ymin = -0.47,
  #           xmax = median(vdata$vuln_climate), ymax = median(vdata$vuln_fi), 
  #           fill = "lightgreen", alpha = 0.1) +
  
  geom_point(aes(x = vuln_climate, y = vuln_fi), 
             color = "black", size = 0.5, alpha = 0.1) +
  geom_hline(aes(yintercept = median(vuln_fi)), linetype = 1) +
  geom_vline(aes(xintercept = median(vuln_climate)), linetype = 1) +
  #geom_ribbon(aes(x = vuln_climate, ymin = Q2.5, ymax = Q97.5), 
  #            alpha = 0.2, size = 0.5,
  #            data = fitted_vuln) +
  geom_line(aes(x = vuln_climate, y = Estimate), 
              color = "black", alpha = 0.5, size = 0.5,
            linetype = 2,
              data = fitted_vuln) +
  scale_x_continuous(limits = c(-0.8424677,1.5)) +
  scale_y_continuous(limits = c(-0.9682,1.5)) +
  coord_equal() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white"))

ggsave("output/plots/vuln_cross.png", width = 6, height = 6)


median(vdata$vuln_fi) - (1.5 - median(vdata$vuln_fi))
median(vdata$vuln_climate) - (1.5 - median(vdata$vuln_climate))


# option 3
con <- left_join(contributions, herb_pisc$contributions_herb_pisc)

spi <- parallel::mclapply(unique(con$transect_id), function(x){
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
    dplyr::select(transect_id, species, Gc_i, Fn_i, Fp_i, I_pisc_i, I_herb_i)
  return(sub)
}, mc.cores = 50) %>% plyr::ldply()

sum(sub$Fp_i)

spi_unique <- spi %>%
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
  mutate(ep_prop = ep_n/n_spec_tot, n_prop = n/n_spec_tot) %>%
  mutate(perc = 100 * ep_prop/n_prop)


ggplot(spi_unique) +
  geom_bar(aes(x = name, y = n_prop, fill = name), stat = "identity", alpha = 0.5) +
  geom_bar(aes(x = name, y = ep_prop, fill = name), stat = "identity") +
  geom_text(aes(x = name, y = n_prop + 0.05, label = paste0(round(perc), "%"))) +
  scale_fill_fish_d(option = "Callanthias_australis",
                    labels = c("N excretion", "P excretion",
                               "Production", "Herbivory", "Piscivory"),
                    name = "Function") +
  facet_wrap(~vulncat) +
  theme_void()

ggsave("output/plots/vuln_bars_option3.png", width = 6, height = 8)













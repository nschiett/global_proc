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


ggplot(flux) +
  geom_point(aes(x = biomass_tot, y = Gc)) +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")
ggplot(flux) +
  geom_point(aes(x = mean, y = Gc)) +
  #scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")

# Get residuals
flux <- flux %>% 
  mutate(I_herb = case_when(I_herb == 0 ~ NA_real_, TRUE ~ I_herb),
         I_pisc = case_when(I_pisc == 0 ~ NA_real_, TRUE ~ I_pisc)) %>%
  mutate(
  Fn_r = residuals(lm(log(Fn) ~  log(biomass_tot) + mean)),
  Fp_r = residuals(lm(log(Fp) ~  log(biomass_tot) + mean)),
  Gc_r = residuals(lm(log(Gc) ~  log(biomass_tot) + mean)),
  I_herb_r = residuals(lm(log(I_herb) ~  log(biomass_tot) + mean, na.action=na.exclude)),
  I_pisc_r = residuals(lm(log(I_pisc) ~  log(biomass_tot) + mean, na.action=na.exclude)))

fit <- lm(log(Fn) ~  log(biomass_tot) + mean, data = flux)

#multi
get_scores <- function(proc){
  probs <- seq(0,1,0.01)
  quantiles <- quantile(proc, prob = probs, na.rm = TRUE)
  quant <- as.numeric(factor(findInterval(proc, quantiles)))/100
  return(quant)
}

# transform into centiles then divide by 100
scores <- select(flux, transect_id, Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r)  %>% 
  mutate_at(2:6, function(x){get_scores(x)}) %>%
  drop_na()
# correlation for weighing
m <-  cor(scores[2:6])
m <- 1-m
we <-  rowSums((m))
we <- we/sum(we)

multi <- scores %>%
  rowwise() %>% 
  mutate(multi = weighted.mean(c(Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r), w = we ),
         var = var(c(Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r))) %>%
  mutate(multi = 10*multi*(0.3-var)/3)


  
flux <- left_join(flux, select(multi, transect_id, multi))

## remote areas
hist(flux$biomass_tot)
wild <- read_csv("data/humanPopulation.csv") %>%
  left_join(unique(select(flux ,
                          transect_id, sites, bioregion, biomass_tot, 
                          Fn, Fp, Gc, I_herb, I_pisc))) %>%
  inner_join(select(flux, sites, transect_id, biomass_tot)) %>%
  group_by(bioregion) %>%
  filter(pop_100 == 0, 12 < tt_Nmarket) %>%
  select(transect_id, sites, bioregion, biomass_tot, Fn, Fp, Gc, I_herb, I_pisc) %>%
  filter(biomass_tot>100)
summary(wild)


fluxs_sel <-  mutate(flux, logbiomass = log(biomass_tot), logabu = log(abu_tot)) %>%
  select(nspec, troph_m, troph_q3, troph_q1, size_q3,
         imm_m)

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

library(ggrepel)

arrows_12 <- 
  ggplot() +
  geom_path(data = circ,aes(x,y), lty = 2, color = "grey", alpha = 0.7) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), data = rpca, arrow = arrow(length=unit(0.30,"cm"))) +
   geom_text_repel(data = rpca, aes(x =  PC1, 
                                    y = PC2, 
                                    label = c("richness", "trophic level (50%)", "trophic level (97.5%)", "trophic level (2.5%)",
                                              "size (97.5%)", "immaturity (50%)")), 
  size = 4, vjust = 0, colour="black",
                   angle = (180/pi) * atan(rpca$PC2/rpca$PC1), 
                   hjust = (1 - 2 * sign(rpca$PC1)) / 2, segment.alpha = 0) +
  coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
  labs(x = "PC1 (26.7%)", y = "PC2 (23.9%)") +
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
  labs(x = "PC3 (17%)", y = "PC4 (13%)") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent"))
arrows_34


# Add bivariate density for each point
dpca$density <- fields::interp.surface(
  MASS::kde2d(dpca$PC1, dpca$PC2), dpca[,c("PC1", "PC2")])

dpca$density2 <- fields::interp.surface(
  MASS::kde2d(dpca$PC1, dpca$PC2), dpca[,c("PC3", "PC4")])


## pca with wilderness
dpca_wild <- inner_join(wild, dpca)

pca_wild <-
  ggplot(dpca) +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  geom_vline(xintercept = 0, lty = 2, color = "grey") +
  geom_point(aes(x = PC1, y = PC2), color = "grey" ,
             data = dpca, alpha = 0.5, shape = 16)  +
  geom_point(aes(x = PC1, y = PC2, color = bioregion), 
             data = dpca_wild,  shape = 8, size = 0.8, alpha = 0.6)  +
  stat_density_2d(aes(x = PC1, y = PC2, fill = bioregion, color = bioregion),
                  alpha = 0.3 , geom = "polygon", data = dpca_wild, bins = 2) +
  labs(x = "PC1", y = "PC2") +
  
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  scale_color_fish_d(option = "Scarus_quoyi",
                     labels = c("CIP", "CP", "EP", "WA")) +
  scale_fill_fish_d(option = "Scarus_quoyi",
                    labels = c("CIP", "CP", "EP", "WA")) +
  coord_equal() +
  theme(
    legend.position = c(0.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )
pca_wild

## pca with hotspots
dpca_res<- inner_join(select(res2, transect_id, hot, cold, all_q), dpca)

pca <-
  ggplot(dpca_res) +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  geom_vline(xintercept = 0, lty = 2, color = "grey") +
  geom_point(aes(x = PC1, y = PC2), color = "grey" ,
             data = dpca, alpha = 0.5, shape = 16)  +
  geom_point(aes(x = PC1, y = PC2, color = "red"), 
             data = dpca_res[dpca_res$hot == TRUE ,],  shape = 8, size = 3, alpha = 0.6)  +
  geom_point(aes(x = PC1, y = PC2, color = "blue"), 
             data = dpca_res[dpca_res$cold == TRUE ,],  shape = 3, size = 3, alpha = 0.6)  +
  # stat_density_2d(aes(x = PC1, y = PC2, fill = bioregion, color = bioregion),
  #                 alpha = 0.3 , geom = "polygon", data = dpca_wild, bins = 2) +
  labs(x = "PC1", y = "PC2") +
  
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  scale_color_fish_d(option = "Scarus_quoyi",
                     labels = c("CIP", "CP", "EP", "WA")) +
  scale_fill_fish_d(option = "Scarus_quoyi",
                    labels = c("CIP", "CP", "EP", "WA")) +
  coord_equal() +
  theme(
    legend.position = c(0.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) +
  facet_wrap(~ bioregion)
pca
pca <-
  ggplot(dpca_res) +
  geom_hline(yintercept = 0, lty = 2, color = "grey") +
  geom_vline(xintercept = 0, lty = 2, color = "grey") +
  geom_point(aes(x = PC3, y = PC4), color = "grey" ,
             data = dpca, alpha = 0.5, shape = 16)  +
  geom_point(aes(x = PC3, y = PC4, color = "red"), 
             data = dpca_res[dpca_res$hot == TRUE ,],  shape = 8, size = 3, alpha = 0.6)  +
  geom_point(aes(x = PC3, y = PC4, color = "blue"), 
             data = dpca_res[dpca_res$cold == TRUE ,],  shape = 3, size = 3, alpha = 0.6)  +
  # stat_density_2d(aes(x = PC1, y = PC2, fill = bioregion, color = bioregion),
  #                 alpha = 0.3 , geom = "polygon", data = dpca_wild, bins = 2) +
  labs(x = "PC1", y = "PC2") +
  
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  scale_color_fish_d(option = "Scarus_quoyi",
                     labels = c("CIP", "CP", "EP", "WA")) +
  scale_fill_fish_d(option = "Scarus_quoyi",
                    labels = c("CIP", "CP", "EP", "WA")) +
  coord_equal() +
  theme(
    legend.position = c(0.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) +
  facet_wrap(~ bioregion)
pca

ggplot() +
  geom_histogram(aes(x = multi, y = stat(count) / sum(count)), 
                 bins = 20, data = filter(dpca, !transect_id %in% dpca_wild$transect_id), fill = "black", alpha = 0.6) +
  geom_histogram(aes(x = multi, y = stat(count) / sum(count)), 
                 bins = 20, data = dpca_wild, fill = "red", alpha = 0.6) +
  #facet_wrap(~bioregion) +
  theme_bw()

library(fishualize)
fishualize("Hypsypops_rubicundus", n = 5, begin = 0, end = 0.8)



pca_plot <- function(proc, title){
  
  probs <- c(0.05, 0.25, 0.75, 0.95)
  quantiles <- quantile(proc, prob = probs, na.rm = TRUE)
  dpca$quant <- factor(findInterval(proc, quantiles))
  
  p1 <-
    ggplot(dpca) +
    geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_point(aes(x = PC1, y = PC2, color = quant, alpha = 1/density), shape = 16, size = 2) +
    geom_point(data = dpca[dpca$quant == "9",], aes(x = PC1, y = PC2, color = quant, alpha = 1/density), shape = 16, size = 2) +
    scale_color_fish_d("Hypsypops_rubicundus", begin = 0, name = "",
                       labels=c("[0,0.05[", "[0.05,0.25[", "[0.25,0.75[",  "[0.75, 0.95[", "[0.95,1[")) +
    labs(x = "PC1", y = "PC2", title = title) +
    scale_alpha(range = c(.25, .6)) +
    theme_void() +
    guides(alpha = "none") +
    coord_equal() +
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(), legend.position = "none", 
          plot.title = element_text(hjust = 0.5))
  
  p2 <-
    ggplot(dpca) +
    geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
    geom_point( aes(x = PC3, y = PC4, color = quant , alpha = 1/density), shape = 16, size = 2) +
    geom_point(data = dpca[dpca$quant == "9",], aes(x = PC3, y = PC4, color = quant, alpha = 1/density2), shape = 16, size = 2) +
    scale_alpha(range = c(.25, .6) ) +
    scale_color_fish_d("Hypsypops_rubicundus", begin = 0, name = "",
                       labels=c("[0,0.05[", "[0.05,0.10[", "[0.10,0.25[", "[0.25,0.50[", "[0.50,0.75[", "[0.75, 0.90[", "[0.90,0.95[", "[0.95,1[")) +
    labs(x = "PC3", y = "PC4", title = title) +
    coord_equal() +
    theme_void() +
    guides(alpha = "none") +
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(), 
          legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  
  return(list(p1, p2))
  
}



p2 <- pca_plot(dpca$Fn_r, "N excretion")
p3 <- pca_plot(dpca$Fp_r, "P excretion")
p1 <- pca_plot(dpca$Gc_r, "Production")
p5 <- pca_plot(dpca$I_herb_r, "Herbivory")
p4 <- pca_plot(dpca$I_pisc_r, "Piscivory")
p6 <- pca_plot(dpca$multi, "Multifunctionality")

p6[[1]]

p2[[1]] <- p2[[1]] + theme(legend.position = "bottom") 

library(patchwork)
layout <- "
AABB
AABB
CDIJ
EFKL
GHMN
"
arrows_12 + arrows_34 + p1[[1]]+ p2[[1]]+ p3[[1]]+ p4[[1]]+ p5[[1]]+ p6[[1]] +
  p1[[2]]+ p2[[2]]+ p3[[2]]+ p4[[2]]+ p5[[2]]+ p6[[2]] +
  plot_layout(design = layout) 
ggsave("output/plots/pca_plot.pdf" ,height = 15, width = 12)


library(patchwork)
layout <- "
AAA
AAA
AAA
BCD
EFG
"
arrows_12 + p1[[1]]+ p2[[1]]+ p3[[1]]+ p4[[1]]+ p5[[1]] +p6[[1]]+
  plot_layout(design = layout) +   plot_annotation(tag_levels = 'a', tag_suffix = ")") 
ggsave("output/plots/fig1_pca_plot.png" ,height = 14, width = 8)


layout <- "
AAABBB
AAABBB
AAABBB
CDEFGH
"
arrows_12 + pca_wild +p1[[1]]+ p2[[1]]+ p3[[1]]+ p4[[1]]+ p5[[1]] +p6[[1]]+
  plot_layout(design = layout) +   plot_annotation(tag_levels = 'a', tag_suffix = ")") 

ggsave("output/plots/fig1_pca_plot_2.png" ,height = 10, width = 12)



ggplot(flux, aes(x = imm_m, y = multi)) +
  geom_point() +
  geom_smooth()
ggplot(flux, aes(x = size_q3, y = multi)) +
  geom_point() +
  geom_smooth()
ggplot(flux, aes(x = nspec, y = multi)) +
  geom_point() +
  geom_smooth()
ggplot(flux, aes(x = troph_m, y = multi)) +
  geom_point() +
  geom_smooth()



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
  mutate(multi = weighted.mean(c( Fn, Fp, Gc, I_herb, I_pisc), w = we ),
         var = var(c( Fn, Fp, Gc, I_herb, I_pisc))) %>%
  mutate(multi = multi*(1-(var/3000))/100)

hist(multi$multi)
fluxs <- filter(flux, transect_id %in% scores$transect_id) %>%
  mutate(multi = multi$multi) 


wild <- read_csv("data/humanPopulation.csv") %>%
  left_join(unique(select(fluxs ,
                          transect_id, sites, bioregion, biomass_tot, 
                          Fn, Fp, Gc, I_herb, I_pisc, 
                          Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r, multi, mean))) %>%
  inner_join(select(flux, sites, transect_id, biomass_tot)) %>%
  group_by(bioregion) %>%
  mutate(cat = case_when(pop_100 == 0 & 12 < tt_Nmarket & biomass_tot>100 ~ "wild",
                         pop_50 > 0 & 6 > tt_Nmarket ~ "not wild",
                         TRUE ~ "not wild")) 
summary(wild)


summary(lm(multi~standard(log(biomass_tot)), data = wild))
wild$multi_r <- residuals(lm(multi~normalize(log(biomass_tot)), data = wild))
fluxs$multi_r <- residuals(lm(multi~normalize(log(biomass_tot)), data = fluxs))


wilds <- wild %>%
  #filter(!cat == "int") %>%
  group_by(sites, bioregion, cat) %>%
  summarize_if(is.numeric, median, na.rm = TRUE)
%>%
  group_by(locality, bioregion, cat) %>%
  summarize_if(is.numeric, median, na.rm = TRUE)

ggplot(wild) +
  geom_boxplot(aes(x = bioregion, y = multi_r, color = cat))

test <- flux %>%
  mutate(bin = cut(biomass_tot, breaks = seq(0, 36200, by = 50))) %>%
  group_by(bin) %>%
  mutate(biomass_bin = mean(biomass_tot)) %>%
  group_by(bin, biomass_bin) %>%
  summarise_if(is.numeric, function(x){max(x, na.rm = TRUE)/min(x,na.rm = TRUE)}) %>%
  pivot_longer(names_to = "var", values_to = "value", cols = c(Fn, Fp, Gc, I_herb, I_pisc)) %>%
  filter(!value == 1)


ggplot(bins, aes(x = biomass_bin, y = value, fill = var, color = var)) +
  geom_point(size = 1, alpha = 0.8) +
  scale_y_continuous(trans = "log10", breaks = c(2,10,100,1000, 10000)) +
  scale_x_continuous(trans = "log10") +
  theme_bw() 

## multipanel
b1 <-
  ggplot(flux, aes(x = Fn_r, y = Fp_r)) +
    geom_point(alpha = 0.1, size = 0.5) +
    #geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
    theme_bw() +
  labs(y = "Residuals P excretion") +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_text(size = 8), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
b1
c1 <-
  ggplot(flux, aes(x = Fn_r, y = Gc_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  #geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw()+
  labs(y = "Residuals production") +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_text(size = 8), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
c1
d1 <-
  ggplot(flux, aes(x = Fn_r, y = I_herb_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  #geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw()+
  labs(y = "Residuals herbivory") +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_text(size = 8), 
        axis.title.x = element_blank(), 
        panel.border = element_blank(), 
        axis.text = element_text(size = 6),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
d1
e1 <-
  ggplot(flux, aes(x = Fn_r, y = I_pisc_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  #geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw()+
  labs(y = "Residuals piscivory") +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_text(size = 8), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
e1
c2 <-
  ggplot(flux, aes(x = Fp_r, y = Gc_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  #geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
c2
d2 <-
  ggplot(flux, aes(x = Fp_r, y = I_herb_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  #geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw()+
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
d2

e2 <-
  ggplot(flux, aes(x = Fp_r, y = I_pisc_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  #geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw()+
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
e2

d3 <-
  ggplot(flux, aes(x = Gc_r, y = I_herb_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  #geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw()+
  scale_x_continuous(limits = c(-3, 2.2)) +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 

d3
e3 <-
  ggplot(flux, aes(x = Gc_r, y = I_pisc_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  #geom_smooth( size = 1, linetype = 1, color = "orange") +
  theme_bw()+
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
e3
e4 <-
  ggplot(flux, aes(x = I_herb_r, y = I_pisc_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  #geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
e4
f1 <-
  ggplot(fluxs, aes(x = Fn_r, y = multi_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw() +
  labs(x = "Residuals N excretion", y = "Residuals multifunctionality") +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_text(size = 8), 
        axis.title.x = element_text(size = 8), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
f1
f2 <-
  ggplot(fluxs, aes(x = Fp_r, y = multi_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw()+
  labs(x = "Residuals P excretion") +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 8), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
f2
f3 <-
  ggplot(fluxs, aes(x = Gc_r, y = multi_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw()+
  labs(x = "Residuals production") +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 8), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
f3
f4 <-
  ggplot(fluxs, aes(x = I_herb_r, y = multi_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw()+
  labs(x = "Residuals herbivory") +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 8), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
f4

f5 <-
  ggplot(fluxs, aes(x = I_pisc_r, y = multi_r)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "loess", size = 1, linetype = 1, color = "orange") +
  theme_bw()+
  labs(x = "Residuals piscivory") +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 8), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
f5

a1 <-
  ggplot(aes(x = Fn_r), data = flux) +
    geom_histogram() +
    theme_bw() +
    theme(axis.line = element_line(color = "black"),
        axis.title.y = element_text(size = 8), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
a1
b2 <-
  ggplot(aes(x = Fp_r), data = flux) +
  geom_histogram()+
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
b2
c3 <-
  ggplot(aes(x = Gc_r), data = flux) +
  geom_histogram() +
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
c3
d4 <-
  ggplot(aes(x = I_herb_r), data = flux) +
  geom_histogram() +
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
d4
e5 <-
  ggplot(aes(x = I_pisc_r), data = flux) +
  geom_histogram() +
  theme_bw() +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
e5
f6 <-
  ggplot(aes(x = multi_r), data = fluxs) +
  geom_histogram() +
  theme_bw() +
  labs(x = "Residuals multifunctionality") +
  theme(axis.line = element_line(color = "black"),
        axis.title.y = element_blank(), 
        axis.title.x = element_text(size = 8), 
        axis.text = element_text(size = 6),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.spacing = unit(0, "cm")) 
f6

a1 + plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + 
  b1 + b2+ plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() +
  c1 + c2 + c3 + plot_spacer() + plot_spacer() + plot_spacer()  +
  d1 + d2 + d3 + d4 + plot_spacer() + plot_spacer()  +
  e1 + e2 + e3 + e4 + e5 + plot_spacer()  +
  f1 + f2 + f3 + f4 + f5 + f6 + plot_layout(ncol = 6)
ggsave("output/plots/residuals_multipanel.png", width = 10, height = 10)

##### plot 1 simple #####
plot_biomass <-
  ggplot(test, aes(x = var, y = value)) +
  geom_violin(size = 0.5, alpha = 0.7, fill = "darkgrey", color = "black", draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_y_continuous(trans = "log10", breaks = c(2,10,100,1000, 10000)) +
  theme_bw() +
  labs(x = "", y = "fold variation per biomass class") +
  # scale_fill_fish_d(option = "Callanthias_australis",
  #                   labels = c("N excretion", "P excretion",
  #                              "Production", "Herbivory", "Piscivory"),
  #                   name = "Function") +
  # scale_color_fish_d(option = "Callanthias_australis",
  #                   labels = c("N excretion", "P excretion",
  #                              "Production", "Herbivory", "Piscivory"),
  #                   name = "Function") +
  scale_x_discrete(labels = c("N excretion", "P excretion",
                              "Production", "Herbivory", "Piscivory"),
                   position = "top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        text = element_text(size = 10),
        axis.line.y = element_line(),
        axis.ticks.x = element_blank())
plot_biomass
#ggsave("output/plots/fold_variation.png")


fcor <- cor(select(fluxs, Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r)) %>% as.data.frame() 
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
cplot

plot_biomass + cplot  +  
  plot_annotation(tag_levels = 'a', tag_suffix = ")") 
ggsave("output/plots/fig1_simple.png", width = 10, height = 6)

fluxs_tall <- fluxs %>%
  pivot_longer(names_to = "fun", values_to = "value", 
               cols = c(Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r))

pmulti <- 
ggplot(fluxs_tall, aes(x = value, y = multi_r)) +
  geom_point(size = 0.1, alpha = 0.2) +
  geom_smooth(size = 1,
    color = fish(option = "Hypsypops_rubicundus", n = 3)[2]) +
  facet_grid(~fun, scales = "free", 
             labeller = as_labeller(c(
               "Fn_r" = "N excretion", 
               "Fp_r" = "P excretion", 
               "Gc_r" = "Production",
               "I_herb_r"  =  "Herbivory", 
               "I_pisc_r"  =  "Piscivory"))) +
  labs(x = "Residuals", y = "Residuals multifunctionality") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(), 
        text = element_text(size = 10),
        axis.line = element_line(),
        strip.background = element_rect(colour=NA, fill=NA))

layout <- "
AAAAABBBBB
AAAAABBBBB
CCCCCCCCCC
"
plot_biomass + cplot +
  pmulti +
  plot_layout(design = layout) +  
  plot_annotation(tag_levels = 'a', tag_suffix = ")") 
ggsave("output/plots/fig1_simple_v2.png", width = 10, height = 8)


#fit <- brm(multi_r ~ cat + bioregion + (1|sites) + (1|locality), data = wild, cores = 4)
summary(fit)
marginal_effects(fit, probs = c(0.2, 0.8))

test <- brms::posterior_samples(fit)

ggplot() +
  geom_histogram(aes(x = test[,1]), fill = "blue", alpha = 0.5) +
  geom_histogram(aes(x = test[,2] + test[,1]), fill = "red", alpha = 0.5) 
  


ggplot(wild, aes( y = (multi), color = cat, x = log(biomass_tot))) +
  geom_point( alpha = 0.5) +
  #geom_smooth() +
  facet_wrap(~bioregion)

ggplot(wilds, aes( y = multi_r, color = cat, x = cat)) +
  geom_boxplot( alpha = 0.1)  +
  geom_jitter(alpha = 0.01)# +
  facet_wrap(~bioregion)

ggplot(wilds, aes( y = Fn_r, x = cat, color = cat)) +
  geom_boxplot( alpha = 0.1)
ggplot(wilds, aes( y = Fp_r, x = cat, color = cat)) +
  geom_boxplot( alpha = 0.1)
ggplot(wilds, aes( y = Gc_r, x = cat, color = cat)) +
  geom_boxplot( alpha = 0.1)
ggplot(wilds, aes( y = I_herb_r, x = cat, color = cat)) +
  geom_boxplot( alpha = 0.1)
ggplot(wilds, aes( y = I_pisc_r, x = cat, color = cat)) +
  geom_boxplot( alpha = 0.1)
ggplot(wilds, aes( y = I_pisc_r, x = cat, color = cat)) +
  geom_boxplot( alpha = 0.1)

ggplot(wild, aes(x = log(tt_Nmarket), y = log(biomass_tot))) +
  geom_point() +
  geom_smooth(method = "lm")

fit_m <- brm(multi ~ log(grav_Nmarket) + log(biomass_tot) + mean + (1|sites), data = wild, cores = 4)
summary(fit_m)
marginal_effects(fit_m)
test <- brms::posterior_samples(fit_m)
hist(test[,2])
summary(lm(normalize(multi) ~ (log(tt_Nmarket)) + (log(biomass_tot) + mean), data = wild))
summary(lm(normalize(log(Fp)) ~ log(grav_Nmarket) + log(biomass_tot), data = wild))
summary(lm(normalize(log(Fn)) ~ log(grav_Nmarket) + log(biomass_tot), data = wild))
summary(lm(normalize(log(Gc)) ~ log(grav_Nmarket) + log(biomass_tot), data = wild))
summary(lm(normalize(log(I_herb)) ~ log(grav_Nmarket) + log(biomass_tot), data = wild))
summary(lm(normalize(log(I_pisc)) ~ log(grav_Nmarket) + log(biomass_tot), data = wild))


###### suppl. gravity markets ######
tall <- pivot_longer(fluxs, c(Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r)) %>%
  group_by(name) %>%
  left_join(hum) %>%
  mutate(value = standard(value))
a <-
ggplot(tall, aes(x = log(grav_Nmarket), y = value), alpha = 0.2) +
  geom_smooth(method = "lm", aes(color = name, fill = name), alpha = 0.2) +
  scale_color_manual(values = c(col), name = "",
                     labels = c("N excretion", "P excretion",
                               "Production", 
                               "Herbivory", "Piscivory")) +
  scale_fill_manual(values = c(col),  name = "",
                    labels = c("N excretion", "P excretion",
                               "Production", 
                               "Herbivory", "Piscivory")) +
  theme_classic()
a
b <-
  ggplot(wild, aes(x = log(grav_Nmarket), y = standard(Fp_r))) +
  #geom_point(size = 0.2,alpha = 0.1) +
  #geom_point(data = wilds, shape = 8, fill = "black", size = 0.5) +
  geom_smooth(method = "lm", data = wild, color = col[2], fill = col[2]) +
  theme_bw()
b
c <-
  ggplot(wild, aes(x = log(grav_Nmarket), y = standard(Gc_r))) +
  #geom_point(size = 0.2,alpha = 0.1) +
  #geom_point(data = wilds, shape = 8, fill = "black", size = 0.5) +
  geom_smooth(method = "lm", data = wild, color = col[3], fill = col[3]) +
  theme_bw()
c
d <-
  ggplot(wild, aes(x = log(grav_Nmarket), y = standard(I_herb_r))) +
  #geom_point(size = 0.2,alpha = 0.1) +
  #geom_point(data = wilds, shape = 8, fill = "black", size = 0.5) +
  geom_smooth(method = "lm", data = wild, color = col[4], fill = col[4]) +
  theme_bw()
d
e <-
  ggplot(wild, aes(x = log(grav_Nmarket), y = standard(I_pisc_r))) +
  #geom_point(size = 0.2,alpha = 0.1) +
  #geom_point(data = wilds, shape = 8, fill = "black", size = 0.5) +
  geom_smooth(method = "lm", data = wild, color = col[5], fill = col[5]) +
  theme_bw()
e

f <-
  ggplot(wild, aes(x = log(grav_Nmarket), y = standard(multi_r))) +
  #geom_point(size = 0.2,alpha = 0.1) +
  #geom_point(data = wilds, shape = 8, fill = "black", size = 0.5) +
  geom_smooth(method = "lm", data = wild, color = "black", fill = "black") +
  theme_bw()
f
layout <- 
  "ABC
   DEF"


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
ggsave("output/plots/gravity_markets.png", width = 8, height = 6)





# fit_multi <-brm(multi ~ mean + logbiomass + nspec + size_m + troph_m + imm_m +
#                   size_q3 + troph_q3 + imm_q1 + size_q1 + troph_q1 + imm_q3,
#                 data = flux, chain = 3, cores = 3)

###### supplement beta params #######
test <- contributions_sp_loc_occ# %>% filter(rel_occ > 0.001)
ggplot(test) +
  geom_histogram(aes(x = Fn_pphi_m * Fn_p_m)) +
  scale_x_continuous(trans = "log10") +
  facet_wrap(~bioregion, scales = "free_y")

a <- 
  ggplot(test) +
  geom_point(aes(x = Fn_p_m, y = Fn_pphi_m), alpha = 0.4, size = 1) +
  geom_smooth(aes(x = Fn_p_m, y = Fn_pphi_m), alpha = 0.4, method = "lm") +
  scale_x_continuous(trans = "log10", limits = c(0.001, 0.25)) +
  scale_y_continuous(trans = "log10", limits = c(1, 1000)) +
  labs(x = "mu", y = "phi", title = "N excretion") +
  theme_classic()
a
b <- 
  ggplot(test) +
  geom_point(aes(x = Fn_p_m, y = Fn_pphi_m), alpha = 0.4, size = 1) +
  geom_smooth(aes(x = Fn_p_m, y = Fn_pphi_m), alpha = 0.4, method = "lm") +
  scale_x_continuous(trans = "log10", limits = c(0.001, 0.25)) +
  scale_y_continuous(trans = "log10", limits = c(1, 1000)) +
  labs(x = "mu", y = "phi", title = "P excretion") +
  theme_classic()
b
c <- 
  ggplot(test) +
  geom_point(aes(x = Gc_p_m, y = Gc_pphi_m), alpha = 0.4, size = 1) +
  geom_smooth(aes(x = Gc_p_m, y = Gc_pphi_m), alpha = 0.4, method = "lm") +
  scale_x_continuous(trans = "log10", limits = c(0.001, 0.25)) +
  scale_y_continuous(trans = "log10", limits = c(1, 1000)) +
  labs(x = "mu", y = "phi", title = "Production") +
  theme_classic()
c

d <- 
  ggplot(test) +
  geom_point(aes(x = I_herb_p_m, y = I_herb_pphi_m), alpha = 0.4, size = 1) +
  geom_smooth(aes(x = I_herb_p_m, y = I_herb_pphi_m), alpha = 0.4, method = "lm") +
  scale_x_continuous(trans = "log10", limits = c(0.001, 1)) +
  scale_y_continuous(trans = "log10", limits = c(0.1, 1000)) +
  labs(x = "mu", y = "phi", title = "Herbivory") +
  theme_classic()
d
e <- 
  ggplot(test) +
  geom_point(aes(x = I_pisc_p_m, y = I_pisc_pphi_m), alpha = 0.4, size = 1) +
  geom_smooth(aes(x = I_pisc_p_m, y = I_pisc_pphi_m), alpha = 0.4, method = "lm") +
  scale_x_continuous(trans  = "log10", limits = c(0.05, 1)) +
  scale_y_continuous(trans = "log10", limits = c(0.1, 100)) +
  labs(x = "mu", y = "phi", title = "Piscivory") +
  theme_classic()
e

a + b + c + d + e + plot_layout(ncol = 1)

ggsave("output/plots/beta_params.png", width = 8, height = 12)


zs <- filter(contributions, species == "Zebrasoma_scopas")

ggplot(zs) +
  geom_smooth(aes(x = log(nspec), y = Fn_p), method = "lm")

#fit <- brm(Fn_p ~ log(nspec) , data = zs,chains = 1, family = "beta")

summary(fit)


phi = 1000
mu = 0.005
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

hist(key$ks_Fp)




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
# (0,0.6] (0.6,0.7] (0.7,0.8] (0.8,0.9]   (0.9,1] 
# 0        13        52        31         2 
# >0.7 --> 86.7%

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
summary(keys$cut)
# (0,0.6] (0.6,0.7] (0.7,0.8] (0.8,0.9]   (0.9,1] 
# 4        15        59        19         1 
# >0.7 --> 80.6%

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
summary(keys$cut)
# (0,0.6] (0.6,0.7] (0.7,0.8] (0.8,0.9]   (0.9,1] 
# 2        14        61        20         1 
#>0.7 --> 83.7%
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
summary(keys$cut)

# (0,0.6] (0.6,0.7] (0.7,0.8] (0.8,0.9]   (0.9,1]      NA's 
#        14        12        11        19        40         2 
#>0.7 72.9  >0.9 41.7

d <-
  ggplot(keys) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = cut, size = (I_pisc)), alpha = 0.7,
             position = position_jitter(0,0), shape = 16) +
  scale_color_fish_d(option = "Hypsypops_rubicundus",  
                     name = "Degree of dominance", na.translate = F, drop = FALSE)  +
  coord_sf(ylim = c(-35, 35), expand = FALSE) +
  labs(title = "Piscivory") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.8),
                              title.position = "top")) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  theme_worldmap() + 
  theme(legend.position = "bottom", 
                           legend.title = element_text(size = 6),
                           legend.text = element_text(size = 6),
        legend.title.align = 0.5)
d
keys$cut <- cut(keys$ks_I_herb, breaks = breaks)
summary(keys$cut)
# 20        42        20         9         7 
#>0.7 36.7,  <0.6 20.4
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
  theme_worldmap() 

e 

multimap <-
  a + b + c + d + e + plot_layout(ncol = 1)
#ggsave("output/plots/multimap.tiff", multimap, width = 8, height = 12)



key_long <- pivot_longer(key, cols = c(ks_Gc, ks_Fn, ks_Fp, ks_I_pisc, ks_I_herb))
ggplot(key_long) +
  geom_density(aes(x = value, fill = name, color = name), alpha = 0.5) +
  facet_grid(name~bioregion, scales = "free_y")

ggplot(key_long) +
  geom_point(aes(x = nspec, y = value), size = 0.1, alpha = 0.5)+
  geom_smooth(aes(x = nspec, y = value, color = name)) +
  facet_wrap(~name)

####### fig3 part 1 ######

## dom

ggplot(key, alpha = 0.8) +
  #geom_hline(aes(yintercept = 1), color = "grey80") +
  geom_violin(aes(y = ks_Fn, x = bioregion),
              fill = col[1], color = col[1],
              alpha = 0.5, size = 1, draw_quantiles = c(0.25,0.5,0.75)) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
  theme_classic() +
  theme(plot.margin = unit(c(0,0.2,0,0.2), "cm"),
        axis.ticks.y = element_line(size = 1),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(size = 1), 
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8)) 

key_long <- key %>%
  pivot_longer(c(ks_Fn, ks_Fp, ks_Gc, ks_I_herb, ks_I_pisc))

dd_com <- 
ggplot(key_long, alpha = 0.8) +
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
        axis.ticks.x = element_blank()
          )
# 
# a <-
# ggplot(key, alpha = 0.8) +
#   geom_hline(aes(yintercept = 1), color = "grey80") +
#   geom_eyeh(aes(x = ks_Fn, y = 1),
#             fill = col[1],
#             alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
#   scale_x_continuous(limits = c(0,1)) +
#   labs(title = "N excretion") +
#   theme_void() +
#   theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
#         plot.title = element_text(vjust = -1, hjust = 0.005,
#                                   size = 10))
# a
# b <-
#   ggplot(key, alpha = 0.8) +
#   geom_hline(aes(yintercept = 1), color = "grey80") +
#   geom_eyeh(aes(x = ks_Fp, y = 1),
#             fill = col[2],
#             alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
#   scale_x_continuous(limits = c(0,1)) +
#   labs(title = "P excretion") +
#   theme_void() +
#   theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
#         plot.title = element_text(vjust = -1, hjust = 0.005,
#                                   size = 10))
# b
# c <-
# ggplot(key, alpha = 0.8) +
#   geom_hline(aes(yintercept = 1), color = "grey80") +
#   geom_eyeh(aes(x = ks_Gc, y = 1),
#             fill = col[3],
#             alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
#   scale_x_continuous(limits = c(0,1)) +
#   labs(title = "Production") +
#   theme_void() +
#   theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
#         plot.title = element_text(vjust = -1, hjust = 0.005,
#                                   size = 10))
# c
# d <-
#   ggplot(key, alpha = 0.8) +
#   geom_hline(aes(yintercept = 1), color = "grey80") +
#   geom_eyeh(aes(x = ks_I_herb, y = 1),
#             fill = col[4],
#             alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
#   scale_x_continuous(limits = c(0,1)) +
#   labs(title = "Herbivory") +
#   theme_void() +
#   theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
#         plot.title = element_text(vjust = -1, hjust = 0.005,
#                                   size = 10))
# d
# e <-
#   ggplot(key, alpha = 0.8) +
#   geom_hline(aes(yintercept = 1), color = "grey80") +
#   geom_eyeh(aes(x = ks_I_pisc, y = 1),
#             fill = col[5],
#             alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
#   scale_x_continuous(limits = c(0,1)) +
#   labs(x = "Degree of dominance (communities)") +
#   labs(title = "Predation") +
#   theme_classic() +
#   theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.line.y = element_blank(),
#         axis.line.x = element_line(),
#         axis.ticks.x = element_line(size = 1), 
#         axis.text.x = element_text(size = 8),
#         axis.title.x = element_text(size = 8),
#         plot.title = element_text(vjust = -1, hjust = 0.005,
#                                   size = 10)) 
# e
# 
# 
# a + b + c + d + e + plot_layout(ncol = 1)
ggplot(key, alpha = 0.8) +
  #geom_hline(aes(yintercept = 1), color = "grey80") +
  geom_violin(aes(y = ks_Fn, x = bioregion),
              fill = col[1], color = col[1],
              alpha = 0.5, size = 1, draw_quantiles = c(0.25,0.5,0.75)) +
  scale_y_continuous(limits = c(0,1), breaks = c(0,1)) +
  theme_classic() +
  theme(plot.margin = unit(c(0,0.2,0,0.2), "cm"),
        axis.ticks.y = element_line(size = 1),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(size = 1), 
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8)) 

spi_sp_long <- spi_sp %>%
  pivot_longer(c(Fn_i, Fp_i, Gc_i, I_herb_i, I_pisc_i))

fd_sp <-
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
        axis.ticks.x = element_blank()
  )
# 
# a1 <- 
#   ggplot(spi_sp, alpha = 0.8) +
#   geom_hline(aes(yintercept = 1), color = "grey80") +
#   geom_eyeh(aes(x = Fn_i, y = 1), 
#             fill = col[1], 
#             alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
#   scale_x_continuous(limits = c(0,1)) +
#   theme_void() +
#   theme(plot.margin = unit(c(0,0.2,0,0.2), "cm"))
# a1
# b1 <- 
#   ggplot(spi_sp, alpha = 0.8) +
#   geom_hline(aes(yintercept = 1), color = "grey80") +
#   geom_eyeh(aes(x = Fp_i, y = 1), 
#             fill = col[2],
#             alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
#   scale_x_continuous(limits = c(0,1)) +
#   theme_void() +
#   theme(plot.margin = unit(c(0,0.2,0,0.2), "cm"))
# b1
# c1 <- 
#   ggplot(spi_sp, alpha = 0.8) +
#   geom_hline(aes(yintercept = 1), color = "grey80") +
#   geom_eyeh(aes(x = Gc_i, y = 1), 
#             fill = col[3], 
#             alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
#   scale_x_continuous(limits = c(0,1)) +
#   theme_void() +
#   theme(plot.margin = unit(c(0,0.2,0,0.2), "cm"))
# c1
# d1 <- 
#   ggplot(spi_sp, alpha = 0.8) +
#   geom_hline(aes(yintercept = 1), color = "grey80") +
#   geom_eyeh(aes(x = I_herb_i, y = 1), 
#             fill = col[4], 
#             alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
#   scale_x_continuous(limits = c(0,1)) +
#   theme_void() +
#   theme(plot.margin = unit(c(0,0.2,0,0.2), "cm"))
# d1
# e1 <- 
#   ggplot(spi_sp, alpha = 0.8) +
#   geom_hline(aes(yintercept = 1), color = "grey80") +
#   geom_eyeh(aes(x = I_pisc_i, y = 1), 
#             fill = col[5], 
#             alpha = 0.5, size = 3, .width = c(0.25,0.5,0.75)) +
#   scale_x_continuous(limits = c(0,1)) +
#   labs(x = "Frequency of importance (species)") +
#   theme_classic() +
#   theme(plot.margin = unit(c(0,0.2,0,0.2), "cm"),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.line.y = element_blank(),
#         axis.line.x = element_line(),
#         axis.ticks.x = element_line(size = 1), 
#         axis.text.x = element_text(size = 8),
#         axis.title.x = element_text(size = 8)) 
# e1
# 
# 
# 

####### fig3 part 2######

con <- left_join(contributions, herb_pisc$contributions_herb_pisc)

con[con$Family == "Scaridae", "Family"] <- "Labridae"

main_fam <- dplyr::select(con, species, Family, bioregion) %>%
  group_by(Family) %>%
  summarize(n())

library(forcats)

tax <- rfishbase::load_taxa()
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
  ungroup() %>%
  group_by(bioregion) %>%
  mutate_at(.vars = vars(ends_with("_p")), function(x){x/sum(x)}) %>%
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
  geom_tile(aes(y = Family, x = bioregion, alpha = Gc_p), fill = "black") +
  scale_y_discrete(position = "right") +
  scale_x_discrete(labels = c("CIP", "CP", "EA", "EP", "WA", "WI")) +
  scale_alpha_continuous(limits = c(0,1), range = c(0,1)) +
  theme_heat()+
  theme(axis.text.x = element_text(angle = 90, size = 6))
a2
b2 <-
  ggplot(filter(con_fam, Family %in% unique(famrank$Family))) +
  geom_tile(aes(y = Family, x = bioregion, alpha = Fn_p), fill = "black") +
  scale_alpha_continuous(limits = c(0,1), range = c(0,1)) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(labels = c("CIP", "CP", "EA", "EP", "WA", "WI")) +
  theme_heat()+
  theme(axis.text.x = element_text(angle = 90, size = 6))
b2
c2 <-
  ggplot(filter(con_fam, Family %in% unique(famrank$Family))) +
  geom_tile(aes(y = Family, x = bioregion, alpha = Fp_p), fill = "black") +
  scale_alpha_continuous(limits = c(0,1), range = c(0,1)) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(labels = c("CIP", "CP", "EA", "EP", "WA", "WI")) +
  theme_heat()+
  theme(axis.text.x = element_text(angle = 90, size = 6))
c2
d2 <-
  ggplot(filter(con_fam, Family %in% unique(famrank$Family))) +
  geom_tile(aes(y = Family, x = bioregion, alpha = I_pisc_p)) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(labels = c("CIP", "CP", "EA", "EP", "WA", "WI")) +
  scale_alpha_continuous(limits = c(0,1), breaks = c(0.25,0.5,0.75), 
                         range = c(0,1), name = "Proportion contribution") +
  guides(alpha = guide_legend(override.aes = list(size = 0.5),
                              title.position = "top")) +
  theme_heat() +
  theme(axis.text.x = element_text(angle = 90, size = 6),
        legend.position = "bottom", legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.title.align = 0.5) 
d2
e2 <-
  ggplot(filter(con_fam, Family %in% unique(famrank$Family))) +
  geom_tile(aes(y = Family, x = bioregion, alpha = I_herb_p)) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(labels = c("CIP", "CP", "EA", "EP", "WA", "WI")) +
  scale_alpha_continuous(limits = c(0,1), breaks = c(0.25,0.5,0.75), 
                         range = c(0,1)) +
  theme_heat() +
  theme(axis.text.x = element_text(angle = 90, size = 6)) 
e2

# p1 <- a +  b + c + d + e + plot_layout(ncol = 1)
# p2 <- a2 + b2 + c2 + d2 + e2 + plot_layout(ncol = 1)
# p2
p <- b + b2 + c + c2 + a + a2 + e + e2 + d  + d2 + 
  plot_layout(ncol = 2, widths = c(7, 1)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 8, hjust = 1, vjust = -3))
 #p



ggsave("output/plots/multimap.png", p, height = 10, width = 9)

####### family dens ######

#prop biomass
loadd(tfish)
loadd(params)

pb <- tfish %>%
  left_join(select(params, species, lwa_m, lwb_m, diet_cat)) %>%
  mutate(biomass = abun * lwa_m * size_cm ^ lwb_m) %>%
  mutate(biomass = biomass/area) %>%
  group_by(transect_id, species, diet_cat) %>%
  summarize(biomass = sum(biomass)) %>%
  group_by(transect_id) %>%
  mutate(biomass_tot = sum(biomass),
         biomass_herb = case_when(diet_cat == 2 ~ sum(biomass[diet_cat == 2]),
                                  TRUE ~ NA_real_),
         biomass_pisc = case_when(diet_cat == 4 ~ sum(biomass[diet_cat == 4]),
                                  TRUE ~ NA_real_)) %>%
  mutate(biomass_p = biomass/biomass_tot,
         biomass_pisc_p = biomass/biomass_pisc,
         biomass_herb_p = biomass/biomass_herb) 

ggplot(con_fam[con_fam$Family %in% famrank$Family,]) +
  geom_density(aes(x = Gc_p - biomass_p, color = Family)) +
  xlim(c(-0.1, 0.1))

con_fam <- left_join(contributions, herb_pisc$contributions_herb_pisc) %>%
  left_join(tax2) %>%
  left_join(pb) %>%
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

famrank <- con_fam %>%
  group_by(Family) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  ungroup() %>%
  mutate(
    Fn_r = as.numeric(dense_rank(desc(Fn_p))),
    Fp_r = as.numeric(dense_rank(desc(Fp_p))),
    Gc_r = as.numeric(dense_rank(desc(Gc_p))),
    I_herb_r = as.numeric(dense_rank(desc(I_herb_p))),
    I_pisc_r = as.numeric(dense_rank(desc(I_pisc_p)))
  ) %>% rowwise() %>%
  mutate(rank = mean(c(Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r))) %>%
  filter(rank < 9.5) 

con_fam <- left_join(contributions, herb_pisc$contributions_herb_pisc) %>%
  left_join(tax2) %>%
  left_join(pb) %>%
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
         value_uq = quantile(value, 0.75, na.rm = TRUE))

famo <- unique(con_fam[con_fam$Family %in% famrank$Family,c("Family", "biomass_p")])
famo <- famo[order(famo$biomass_p),]

factor()
con_fam <- con_fam[con_fam$Family %in% famrank$Family,]
con_fam$Family <- factor(con_fam$Family, famo$Family)
levels(con_fam$Family)

cfam <- 
ggplot(con_fam[con_fam$Family %in% famrank$Family,]) +
  geom_vline(xintercept = 0) +
  geom_errorbarh(aes(xmin = value_lq, xmax = value_uq, x = value_m, 
                     y = Family, color = reorder(name, rev(name))),
                 position = ggstance::position_dodgev(height = -0.7), height = 0) +
  geom_point(aes(xmin = value_lq, xmax = value_uq, x = value_m, y = Family, color = name),
             position = ggstance::position_dodgev(height = -0.7)) +
  scale_color_fish_d(option = "Callanthias_australis",
                     labels = c("N excretion", "P excretion", "Production", "Herbivory", "Piscivory"),
                     name = "Function") +
  labs(x = "contribution function - contribution biomass", y = "") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), panel.border = element_blank(), 
        legend.position = "right")
  ggsave("output/plots/fig3_family.png")
  
 
layout = 
  "AB
   AB
   AC
   AC"
  
cfam + dd_com + fd_sp  +
  plot_layout(design = layout, guides = "collect") +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 12)) 
ggsave("output/plots/figure3_take2.png", width = 9, height = 6)

ggplot(con_fam[con_fam$Family %in% famrank$Family,]) +
  ggridges::geom_density_ridges(aes(x = log(value), y = Family, fill = name, color = name), 
                                 alpha = 0.2, outlier.size = 0.01) +
  scale_fill_fish_d(option = "Callanthias_australis") +
  scale_color_fish_d(option = "Callanthias_australis") +
  theme_bw()
ggplot(con_fam[con_fam$Family %in% famrank$Family,]) +
  ggridges::geom_density_ridges(aes(x = (value), y = Family, fill = name, color = name), 
                                alpha = 0.4, outlier.size = 0.01) +
  scale_fill_fish_d(option = "Callanthias_australis") +
  scale_color_fish_d(option = "Callanthias_australis") +
  theme_bw()

con_fam <- left_join(contributions, herb_pisc$contributions_herb_pisc) %>%
  left_join(tax2) %>%
  filter(!is.na(Family)) %>%
  group_by(bioregion, locality, sites, transect_id, Family) %>%
  summarize_at(.vars = vars(ends_with("_p")), sum, na.rm = TRUE)  %>%
  group_by(bioregion, locality, sites, Family) %>%
  summarize_at(.vars = vars(ends_with("_p")), mean, na.rm = TRUE)  %>%
  #group_by(bioregion, locality, Family) %>%
  #summarize_at(.vars = vars(ends_with("_p")), mean, na.rm = TRUE)  %>%
  ungroup() %>%
  pivot_longer(c(Fn_p, Fp_p, Gc_p, I_herb_p, I_pisc_p)) %>%
  mutate(value = case_when(value == 0 ~ NA_real_,
                           TRUE ~ value)) %>%
  group_by(Family)

ggplot(con_fam[con_fam$Family %in% famrank$Family,]) +
  geom_density(aes(x = log(Fp_p), color = Family), alpha = 0.7)


# option 2
p <-
a + a1 + b2 + b + b1 + c2 +
  c + c1 + a2 + d + d1 + e2 + e + e1 + d2 +
  plot_layout(ncol = 3, widths = c(3,3,1)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 8, hjust = 1, vjust = -3))
ggsave("output/plots/figure3_option2.png", p, height = 8, width = 7)


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
# if con > 1/N !!

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
    dplyr::select(bioregion, transect_id, species, Gc_i, Fn_i, Fp_i, I_pisc_i, I_herb_i)
  return(sub)
}, mc.cores = 30) %>% plyr::ldply()


# proportion when important 
spi_sp <- spi %>%
  dplyr::select(bioregion, locality, sites, species, Gc_i, Fn_i, Fp_i, I_pisc_i, I_herb_i) %>%
  group_by(species) %>%
  mutate(occ = length(Gc_i)) %>%
  group_by(species, occ) %>%
  dplyr::summarise_if(is.logical, function(x){sum(x, na.rm = TRUE)/length(x[!is.na(x)])}) 

ggplot(spi_sp) +
  geom_jitter(aes(x = log(occ), Fn_i), alpha = 0.1, width = 0.01, height = 0.01) +
  geom_smooth(aes(x = log(occ), Fn_i), alpha = 0.1, width = 0.01, height = 0.01) 
hist(spi_sp$Gc_i)

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


# global

bp1 <-
  ggplot(filter(spi_unique, vulncat == "vuln_high") %>% drop_na()) +
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
  ggplot(filter(spi_unique, vulncat == "vuln_cli") %>% drop_na()) +
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
  ggplot(filter(spi_unique, vulncat == "vuln_low")) +
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
  ggplot(filter(spi_unique,  vulncat == "vuln_fi")) +
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





# per bioregion
spi_unique <- spi %>%
  dplyr::select(bioregion, species, Gc_i, Fn_i, Fp_i, I_pisc_i, I_herb_i) %>%
  unique() %>%
  group_by(species, bioregion) %>%
  dplyr::summarise_if(is.logical, function(x){sum(x)>0}) %>%
  left_join(vulnerability)%>% 
  drop_na(vuln_climate, vuln_fi) %>%
  ungroup() %>%
  mutate(vulncat = case_when(
    (vuln_climate > median(vuln_climate) & vuln_fi > median(vuln_fi)) ~ "vuln_high",
    (vuln_climate > median(vuln_climate) & vuln_fi <= median(vuln_fi)) ~ "vuln_cli",
    (vuln_climate <= median(vuln_climate) & vuln_fi > median(vuln_fi)) ~ "vuln_fi",
    (vuln_climate <= median(vuln_climate) & vuln_fi <= median(vuln_fi)) ~ "vuln_low")) %>%
  dplyr::select(species, Family, bioregion, vulncat, Gc = Gc_i, Fn = Fn_i, Fp = Fp_i, I_pisc = I_pisc_i, I_herb = I_herb_i) %>%
  pivot_longer(cols = c(Gc, Fn, Fp, I_pisc, I_herb), names_to = "name", values_to = "ep") %>%
  group_by(name, bioregion) %>%
  mutate(n_spec_tot = sum(!is.na(ep))) %>%
  dplyr::group_by(vulncat, name, bioregion, n_spec_tot) %>%
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
  facet_wrap(bioregion~vulncat, ncol = 4) +
  theme_void()


### find way to plot angles

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
middle  

ggsave("output/plots/fig4_middle.png", middle, width = 4, height = 4)

#wa

bp1 <-
  ggplot(filter(spi_unique, bioregion == "w_atlantic", vulncat == "vuln_high") %>% drop_na()) +
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
ggsave("output/plots/fig4_wa_vuln_high.png", bp1, width = 4, height = 4)

bp2 <-
  ggplot(filter(spi_unique, bioregion == "w_atlantic", vulncat == "vuln_cli") %>% drop_na()) +
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

ggsave("output/plots/fig4_wa_vuln_cc.png", bp2, width = 4, height = 4)

bp3 <-
  ggplot(filter(spi_unique, bioregion == "w_atlantic", vulncat == "vuln_low")) +
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
ggsave("output/plots/fig4_wa_vuln_low.png", bp3, width = 4, height = 4)


bp4 <-
  ggplot(filter(spi_unique, bioregion == "w_atlantic", vulncat == "vuln_fi")) +
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
ggsave("output/plots/fig4_wa_vuln_fi.png", bp4, width = 4, height = 4)

# CP


bp1 <-
  ggplot(filter(spi_unique, bioregion == "c_pacific", vulncat == "vuln_high") %>% drop_na()) +
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
ggsave("output/plots/fig4_cp_vuln_high.png", bp1, width = 4, height = 4)

bp2 <-
  ggplot(filter(spi_unique, bioregion == "c_pacific", vulncat == "vuln_cli") %>% drop_na()) +
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

ggsave("output/plots/fig4_cp_vuln_cc.png", bp2, width = 4, height = 4)

bp3 <-
  ggplot(filter(spi_unique, bioregion == "c_pacific", vulncat == "vuln_low")) +
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
ggsave("output/plots/fig4_cp_vuln_low.png", bp3, width = 4, height = 4)


bp4 <-
  ggplot(filter(spi_unique, bioregion == "c_pacific", vulncat == "vuln_fi")) +
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
ggsave("output/plots/fig4_cp_vuln_fi.png", bp4, width = 4, height = 4)


# dev.off()
# pushViewport(viewport(layout = grid.layout(2,2), width = unit(8, "cm"), height = unit(8, "cm")))
# grid.draw(middle)
# pushViewport(viewport(x = 0.65, y = 0.7, angle = -45, width = 0.2, height = 0.2))
# grid.draw(ggplotGrob(bp1))
# popViewport()
# pushViewport(viewport(x = 0.35, y = 0.7, angle = 45, width = 0.2, height = 0.2))
# grid.draw(ggplotGrob(bp4))
# popViewport()
# pushViewport(viewport(x = 0.67, y = 0.32, angle = -135, width = 0.2, height = 0.2))
# grid.draw(ggplotGrob(bp2))
# popViewport()
# pushViewport(viewport(x = 0.34, y = 0.32, angle = 135, width = 0.2, height = 0.2))
# grid.draw(ggplotGrob(bp3))
# 
# library("ggplotify")
# 
# vp <- pushViewport(viewport(layout = grid.layout(1,1), width = unit(8, "cm"), height = unit(8, "cm")))
# print(grid.draw(middle), vp = vp)
# vp1 <- 
#   pushViewport(viewport(x = 0.63, y = 0.63, angle = -45, width = 0.2, height = 0.2))
# bpv1 <- as.grob(~grid.draw(ggplotGrob(bp1)))
# popViewport()
# vp4 <-pushViewport(viewport(x = 0.27, y = 0.63, angle = 45, width = 0.2, height = 0.2))
# print(grid.draw(ggplotGrob(bp4)), vp = vp4)
# popViewport()
# vp2 <-pushViewport(viewport(x = 0.63, y = 0.27, angle = -135, width = 0.2, height = 0.2))
# print(grid.draw(ggplotGrob(bp2)), vp = vp2)
# popViewport()
# vp3 <- pushViewport(viewport(x = 0.27, y = 0.27, angle = 135, width = 0.2, height = 0.2))
# print(grid.draw(ggplotGrob(bp3)), vp = vp3)
# dev.off()

###### figure 2 #######


standard <- function(x){
  m <- mean(x, na.rm = TRUE)
  sd <- sd(x, na.rm = TRUE)
  st <- sapply(x, FUN = function(i){
    st <- (i - m)/sd
    return(st)
  })
  return(st)
}



## load data 
hp <- readd(herb_pisc)[[2]]
flux <- readd(summary_transect)
sst <- read.csv("data/avSst.csv")

flux <- left_join(flux, hp) %>% left_join(sst)


flux$logbiomass <- log(flux$biomass_tot)


fit_Fp_st <- brm(standard(log(Fp)) ~ standard(mean) + standard(logbiomass) + standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
                   standard(size_q3) + standard(troph_q3) + standard(imm_q1) + standard(imm_q3) + standard(troph_q1) + standard(size_q1) ,
                 data = flux, chain = 3, cores = 3)

fit_Fn_st <- brm(standard(log(Fn)) ~ standard(mean) + standard(logbiomass) + standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
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


## not standardized
fit_Fp <- brm((log(Fp)) ~ mean + logbiomass + nspec + size_m + troph_m + imm_m +
                size_q3 + troph_q3 + imm_q1 + size_q1 + troph_q1 + imm_q3,
              data = flux, chain = 3, cores = 3)

fit_Fn <- brm((log(Fn)) ~ mean + logbiomass + nspec + size_m + troph_m + imm_m +
                size_q3 + troph_q3 + imm_q1 + size_q1 + troph_q1 + imm_q3,
              data = flux, chain = 3, cores = 3)

fit_Gc <- brm((log(Gc)) ~ mean + logbiomass + nspec + size_m + troph_m + imm_m +
                size_q3 + troph_q3 + imm_q1,
              data = flux, chain = 3, cores = 3)

fit_Iherb <- brm((log(I_herb)) ~ mean + logbiomass + nspec + size_m + troph_m + imm_m +
                   size_q3 + troph_q3 + imm_q1 + size_q1 + troph_q1 + imm_q3,
                 data = flux[flux$I_herb > 0,], chain = 3, cores = 3)

fit_Ipisc <- brm((log(I_pisc)) ~ mean + logbiomass + nspec + size_m + troph_m + imm_m +
                   size_q3 + troph_q3 + imm_q1 + size_q1 + troph_q1 + imm_q3,
                 data = flux[flux$I_pisc > 0,], chain = 3, cores = 3)


extract_b <- function(x, name) {
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


result <- bind_rows(
  extract_b(fit_Fp_st, "Fp"),
  extract_b(fit_Fn_st, "Fn"),
  extract_b(fit_Gc_st, "Gc"),
  extract_b(fit_Iherb_st, "I_herb"),
  extract_b(fit_Ipisc_st, "I_pisc")) 

head(result)

result <- mutate(result, diff = (uq > 0 & lq < 0)) 

importance <- group_by(result, variable) %>%
  summarise(tot = max(abs(mean))) 

result <- result %>%
  dplyr::filter(!variable == "b_Intercept") %>%
  left_join(importance) %>% mutate(variable = fct_reorder(variable, tot)) %>%
  filter(diff == FALSE) %>%
  filter(!variable == "b_standardlogbiomass")

slopes <-
  ggplot(result) +
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
  labs(y = "slope", x = "") +
  scale_x_discrete(labels = c(  "immaturity (2.5%)", "size (median)", "sst", "immaturity (97.5%)", 
                                "size (2.5%)", "richness",  "immaturity (median)",
                                "size (97.5%)", "trophic level (2.5%)",  "trophic level (median)",
                                "trophic level (97.5%)")) +
  theme(legend.position = "bottom", axis.text = element_text(size = 10),
        axis.title = element_text(size = 10), legend.text = element_text(size = 10), 
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()
  )   

ggsave("output/plots/slope_regs.png", slopes)

fitlist <- list(fit_Fn = fit_Fn, fit_Fp = fit_Fp, fit_Gc = fit_Gc, 
                fit_Iherb = fit_Iherb, fit_Ipisc = fit_Ipisc)
row_rep <- function(df, n) {
  df[rep(1:nrow(df), times = n),]
}

pplot <- function(fitlist, var, xlab){
  
  newdata <- select(flux,  mean, logbiomass, nspec, size_m, troph_m, imm_m,
                    size_q3, troph_q3, imm_q1, size_q1, troph_q1, imm_q3) %>%
    summarise_all(median) %>%
    row_rep(n = 100) %>%
    mutate(!!var := seq(min(flux[[var]]), max(flux[[var]]), length.out = 100)) 
  
  mm <- lapply(fitlist, fitted, newdata = newdata)
  
  l <-lapply(1:length(fitlist), function(m){
    pred <- fitted(fitlist[[m]], newdata = newdata)
    pr <- exp(pred[,1])
    #minmax
    min <- min(exp(mm[[m]][,1]))
    max <- max(exp(mm[[m]][,1]))
    #score
    res <- data.frame(newdata[[var]], normalize(pr))
    colnames(res) <- c(var, names(fitlist)[m])
    return(res)
  }) %>% bind_cols() %>%
    select(-.data[[paste0(var, "1")]], -.data[[paste0(var, "2")]],
           -.data[[paste0(var, "3")]], -.data[[paste0(var, "4")]])
  #get multi
  #weigths
  cm <-  cor(l[,c(2:6)])
  cm <- 1 - cm
  we <-  rowSums((cm))
  we <- we/sum(we)
  out <- l %>%
    rowwise() %>%
    mutate(multi = weighted.mean(c(fit_Fn, fit_Fp, fit_Gc,  fit_Iherb, fit_Ipisc), w = we),
           var = var(c(fit_Fn, fit_Fp, fit_Gc,  fit_Iherb, fit_Ipisc))) %>%
    as.data.frame() %>%
    mutate(multi_norm = normalize(multi*(0.3-var))) %>%
    pivot_longer(2:6)
  
  ggplot(out) +
    geom_line(aes(x = .data[[var]], y = value, color = name), size = 1, linetype = 2) +
    geom_line(aes(x = .data[[var]], y = multi_norm), size = 1) +
    scale_color_fish_d(option = "Callanthias_australis") +
    labs(x = xlab, y = "function (%)") +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(), axis.line = element_line(),
          legend.position = "none")
  
}

p1 <- pplot(fitlist, "imm_m", "immaturity (median)")
p2 <- pplot(fitlist, "troph_m", "trophic level (median)")
p3 <- pplot(fitlist, "size_q3", "size (97.5%)")

layout <- "
AB
AB
AC
AC
AD
AD
EE
"
slopes + theme(legend.position = "bottom") + p2 + p1 + p3  + guide_area() +
  plot_layout(design = layout, guides = "collect")+   
  plot_annotation(tag_levels = 'a', tag_suffix = ")")

ggsave("output/plots/figure2.png", width = 8, height = 10)


############# worldmaps 2.0 #############

sumflux <- flux %>% 
  left_join(select(fluxs, locality,  multi)) %>%
  group_by(bioregion, sites, locality) %>%
  summarise_if(is.numeric, median, na.rm = TRUE) %>%
  group_by(bioregion, locality) %>%
  summarise_if(is.numeric, median, na.rm = TRUE) %>%
  ungroup()

# scores <- select(sumflux, locality, Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r)  %>% 
#   mutate_if(is.numeric, function(x) normalize(x)/100) %>%
#   drop_na()
# # correlation for weighing
# m <-  cor(scores[2:6])
# m <- 1-m
# we <-  rowSums((m))
# we <- we/sum(we)

multi <- scores %>%
  rowwise() %>% 
  mutate(multi = weighted.mean(c(Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r), w = we )) %>%
  group_by() %>%
  mutate(multi = normalize(multi)/100) %>%
  rowwise() %>%
  mutate(var = var(c(Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r))) %>%
  ungroup() %>%
  mutate(multi = multi*(0.3-(var)))


ggplot(sumflux)+
  geom_point(aes(x = I_pisc_r, y = multi))

cor(select(sumflux, Fn_r, Fp_r, Gc_r, I_herb_r, I_pisc_r, multi_r))




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

get_scale(sumflux$Fn_r)

getcol <- colorRampPalette(c("white", col[1]))
pal2 <- getcol(n = 3)
getcol <- colorRampPalette(c("black","white"))
pal1 <- getcol(n = 3)
pal <- c(pal1[1], pal2[2:3])
pal <- c("black", "grey50", col[1])



a <-
  ggplot(sumflux) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = get_scale(Fn_r), 
                 size = get_scores(Fn_r)),  alpha = 0.8,
             position = position_jitter(0,0), shape = 16) +
  geom_point(aes(x = lon, y = lat),  alpha = 0.7,
             shape = 0, size = 4,
             data = filter(sumflux, multi > quantile(sumflux$multi, 0.95, na.rm = TRUE))) +
  geom_point(aes(x = lon, y = lat),  alpha = 0.7,
             shape = 1, size = 4,
             data = filter(sumflux, Fn_r > quantile(sumflux$Fn_r, 0.95, na.rm = TRUE))) +
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
  ggplot(sumflux) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = get_scale(Fp_r), 
                 size = get_scores(Fp_r)),  alpha = 0.8,
             position = position_jitter(0,0), shape = 16) +
  geom_point(aes(x = lon, y = lat),  alpha = 0.7,
             shape = 0, size = 4,
             data = filter(sumflux, multi > quantile(sumflux$multi, 0.95, na.rm = TRUE))) +
  geom_point(aes(x = lon, y = lat),  alpha = 0.7,
             shape = 1, size = 4,
             data = filter(sumflux, Fp_r > quantile(sumflux$Fp_r, 0.95, na.rm = TRUE))) +
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
  ggplot(sumflux) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = get_scale(Gc_r), 
                 size = get_scores(Gc_r)),  alpha = 0.8,
             position = position_jitter(0,0), shape = 16) +
  geom_point(aes(x = lon, y = lat),  alpha = 0.7,
             shape = 0, size = 4,
             data = filter(sumflux, multi > quantile(sumflux$multi, 0.95, na.rm = TRUE))) +
  geom_point(aes(x = lon, y = lat),  alpha = 0.7,
             shape = 1, size = 4,
             data = filter(sumflux, Gc_r > quantile(sumflux$Gc_r, 0.95, na.rm = TRUE))) +
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
  ggplot(sumflux) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = get_scale(I_herb_r), 
                 size = get_scores(I_herb_r)),  alpha = 0.8,
             position = position_jitter(0,0), shape = 16) +
  geom_point(aes(x = lon, y = lat),  alpha = 0.7,
             shape = 0, size = 4,
             data = filter(sumflux, multi > quantile(sumflux$multi, 0.95, na.rm = TRUE))) +
  geom_point(aes(x = lon, y = lat),  alpha = 0.7,
             shape = 1, size = 4,
             data = filter(sumflux, I_herb_r > quantile(sumflux$I_herb_r, 0.95, na.rm = TRUE))) +
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
  ggplot(sumflux) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = get_scale(I_pisc_r), 
                 size = get_scores(I_pisc_r)),  alpha = 0.8,
             position = position_jitter(0,0), shape = 16) +
  geom_point(aes(x = lon, y = lat),  alpha = 0.7,
             shape = 0, size = 4,
             data = filter(sumflux, multi > quantile(sumflux$multi, 0.95, na.rm = TRUE))) +
  geom_point(aes(x = lon, y = lat),  alpha = 0.7,
             shape = 1, size = 4,
             data = filter(sumflux, I_pisc_r > quantile(sumflux$I_pisc_r, 0.95, na.rm = TRUE))) +
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
ggsave("output/plots/multimap2.png", multimap, width = 8, height = 8)

################# new analysis for fig 1 ####################

standard <- function(x){
  m <- mean(x, na.rm = TRUE)
  sd <- sd(x, na.rm = TRUE)
  st <- sapply(x, FUN = function(i){
    st <- (i - m)/sd
    return(st)
  })
  return(st)
}
## load data 
flux <- readd(summary_transect)
hp <- readd(herb_pisc)$summary_herb_pisc
sst <- read.csv("data/avSst.csv")

flux <- left_join(flux, hp) %>% left_join(sst)

fitN <- brm(standard(log(Fn)) ~ log(biomass_tot) + mean + (1|locality) + (1|sites), 
            data = flux, cores = 4)

fitP <- brm(standard(log(Fp)) ~ log(biomass_tot) + mean + (1|locality) + (1|sites), 
            data = flux, cores = 4)


summary(fitN)

res_Fn <- 
  fitN %>%
  spread_draws(r_locality[locality, Intercept]) %>%
  median_qi() %>%
  select(locality, Fn_r = r_locality)

res_Fp <- 
  fitP %>%
  spread_draws(r_locality[locality, Intercept]) %>%
  median_qi() %>%
  select(locality, Fp_r = r_locality)

res <- left_join(res_Fn, res_Fp)

ggplot(res) +
  geom_point(aes(x = Fn_r, y = Fp_r)) +
  geom_text_repel(aes(x = Fn_r, y = Fp_r, label = locality))


dt <- data.frame(
  x = 1:200
)
dt$y <- dt$x^0.6
plot(dt)
plot(dt$y, dt$y/dt$x)

##### Supplemental analysis remote reefs #####
loadd(summary_transect_complete)

hum <- read_csv("data/humanPopulation.csv") 

remote <- read_csv("data/humanPopulation.csv") %>%
  mutate(remote = case_when(linear_dist > 200 & pop_100 == 0 ~ "yes",
         TRUE ~ "no")) %>%
  select(sites, remote, grav_markets) %>%
  inner_join(summary_transect_complete) %>%
  left_join(residuals)

nd <- remote %>%
  select(mean) %>%
  unique() %>%
  mutate(biomass_tot = 100) 

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
  mutate(multif = as.numeric(Fn > Fn_ref & Fp > Fp_ref& Gc > Gc_ref& I_herb > I_herb_ref& I_pisc > I_pisc_ref)) #%>%
  #group_by(sites) #%>%
  #summarize_all(mean)

# ggplot(ndp) +
#   geom_boxplot(aes(x = as.character(multif), y = log(biomass_tot))) +
#   geom_hline(yintercept = log(100))

ndp$logbiomass <- log(ndp$biomass_tot)

fit_mf<- brm(multif ~ standard(mean) + standard(logbiomass) + standard(nspec) + standard(size_m) + standard(troph_m) + standard(imm_m) +
               standard(size_q3) + standard(troph_q3) + standard(imm_q1) + standard(imm_q3) + standard(troph_q1) + standard(size_q1) ,
             data = ndp, chain = 1, cores = 1, family = "bernoulli")
summary(fit_mf)

fit_mf1 <- brm(multif ~ mean + logbiomass,
               data = ndp, chain = 1, cores = 1, family = "bernoulli")
summary(fit_mf1)

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
result <- 
  extract_b(fit_mf, "mf")

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
  geom_hline(yintercept = 0, size = 1, color = "black") +
  geom_linerange(aes(x = variable, ymin = 0, ymax = mean), color = "grey40",
                 position = position_dodge(.7), size = 0.5, linetype = 2) + 
  geom_linerange(aes(x = variable, ymin = lq, ymax = uq), color = "grey40",
                 position = position_dodge(.7), size = 1, linetype = 1) + 
  geom_point(aes(x = variable, y = mean), color = "grey40",
             position = position_dodge(.7), size = 3) +
  coord_flip() +
  theme_minimal() +
  labs(y = "effect size", x = "") +
  scale_x_discrete(labels = c("trophic level (median)",
                              "trophic level (2.5%)",
                              "size (median)",
                              "immaturity (median)",
                              "richness",
                              "size (2.5%)",
                              "trophic level (97.5%)",
                              "sst",
                              "biomass")) +
  theme(legend.position = "bottom", axis.text = element_text(size = 12),
        axis.title = element_text(size = 12), legend.text = element_text(size = 12), 
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()
  )   
slopes
ggsave("output/plots/slopes_mf_mod.png")

fit_mf2 <- brm(multif ~ mean + logbiomass ,
                 data = ndp, chain = 1, cores = 1, family = "bernoulli")
summary(fit_mf2)

marginal_effects(fit_mf1, method = "fitted", points = TRUE)
me <- marginal_effects(fit_mf1, method = "fitted")
me$logbiomass

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
ggsave("output/plots/mf_pp_logbiomass.png")
  
sum(ndp$biomass_tot > 100)/nrow(ndp)
sum(ndp$biomass_tot > 450 & ndp$multif>0)/nrow(ndp)
sum(ndp$biomass_tot > 1000)/nrow(ndp)
sum(ndp$multif)/nrow(ndp)
test <- ndp %>%
  filter(multif>0)

f <- fitted(fit_mf1)
sum(f[,1]>0.5)/nrow(f)

me_b <- me$logbiomass
ggplot(ndp) +
  geom_point(aes(x = logbiomass, y = f[,1]), alpha = 0.2) +
  geom_smooth(aes(x = logbiomass, y = f[,1]), method = "loess")

ggplot(ndp) +
  geom_point(aes(x = log(grav_markets), y = f[,1]), alpha = 0.2) +
  geom_smooth(aes(x = log(grav_markets), y = f[,1]), method = "gam")


ggplot(remote) +
  geom_boxplot(aes(x = remote, y = f[,1]))

fit <- brm(log(Fn) ~ log(biomass_tot), data = summary_transect_complete)
fit_Fp <- brm(log(Fp) ~ log(biomass_tot) + mean + remote, data = remote, cores = 4)
fit_Fn <- brm(log(Fn) ~ log(biomass_tot) + mean + remote, data = remote, cores = 4)
fit_Gc <- brm(log(Gc) ~ log(biomass_tot) + mean + remote, data = remote, cores = 4)
fit_I_herb <- brm(log(I_herb) ~ log(biomass_tot) + mean + remote, data = remote[remote$I_herb>0,], cores = 4)
fit_I_pisc <- brm(log(I_pisc) ~ log(biomass_tot) + mean + remote, data = remote[remote$I_herb>0,], cores = 4)

ggplot(summary_transect_complete) +
  geom_point(aes(x = log(biomass_tot), y = log(Fp)), size = 0.5, alpha = 0.5) +
  geom_abline(aes(intercept = -11.64, slope = 1.088), color = "green") +
  geom_abline(aes(intercept = -11.46, slope = 1.05), color = "black")


####### multivariate correlation model ######

# flux <- summary_transect_complete
# 
# # run model
# fit <- brm(
#   mvbind(log(Fn), log(Fp), log(Gc), log(I_herb), log(I_pisc)) ~ (log(biomass_tot)), 
#               data = flux[flux$I_herb>0 & flux$I_pisc>0,][1:100,], cores = 4)
# fails


##### vulnerability per community #####

loadd(vulnerability)
loadd(contributions)
loadd(herb_pisc)
loadd(summary_transect)


con <- left_join(contributions, vulnerability) %>%
  left_join(herb_pisc$contributions_herb_pisc) %>%
  group_by(transect_id) %>%
  dplyr::summarize(
    vuln_Fn_fi = sum(Fn_p * vuln_fi, na.rm = TRUE),#/sum(vuln_fi, na.rm = TRUE),
    vuln_Fp_fi = sum(Fp_p * vuln_fi, na.rm = TRUE),#/sum(vuln_fi, na.rm = TRUE),
    vuln_Gc_fi = sum(Gc_p * vuln_fi, na.rm = TRUE),#/sum(vuln_fi, na.rm = TRUE),
    vuln_Iherb_fi = sum(I_herb_p[I_herb_p > 0] * vuln_fi[I_herb_p > 0], na.rm = TRUE),#/sum(vuln_fi[I_herb_p > 0], na.rm = TRUE),
    vuln_Ipisc_fi = sum(I_pisc_p[I_pisc_p > 0] * vuln_fi[I_pisc_p > 0], na.rm = TRUE),#/sum(vuln_fi[I_pisc_p > 0], na.rm = TRUE),
    vuln_Fn_cl = sum(Fn_p * vuln_climate, na.rm = TRUE),#/sum(vuln_climate, na.rm = TRUE),
    vuln_Fp_cl = sum(Fp_p * vuln_climate, na.rm = TRUE),#/sum(vuln_climate, na.rm = TRUE),
    vuln_Gc_cl = sum(Gc_p * vuln_climate, na.rm = TRUE),#/sum(vuln_climate, na.rm = TRUE),
    vuln_Iherb_cl = sum(I_herb_p[I_herb_p > 0] * vuln_climate[I_herb_p > 0], na.rm = TRUE),#/sum(vuln_climate[I_herb_p > 0], na.rm = TRUE),
    vuln_Ipisc_cl = sum(I_pisc_p[I_pisc_p > 0] * vuln_climate[I_pisc_p > 0], na.rm = TRUE))

con_long <- con %>%#/sum(vuln_climate[I_pisc_p > 0], na.rm = TRUE)) %>%
  pivot_longer(cols = 2:11) %>%
  separate(name, into = c("v", "fun", "impact"), sep = "_", remove = TRUE) %>%
  filter(value>0) %>%
  pivot_wider(names_from = impact, values_from = value)

ggplot(con_long) +
  geom_boxplot(aes(x = fun, y = (value))) +
  facet_wrap(~impact)

con <- left_join(con, unique(select(ungroup(contributions), transect_id, sites, locality, lat, lon)))

test <- filter(con) %>%
  group_by(sites, locality) %>%
  summarise_if(is.numeric, median) %>%
  group_by(locality) %>%
  summarise_if(is.numeric, median)

ggplot(test) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = (value), 
                 size = get_scores(value)),  alpha = 0.8,
             position = position_jitter(0,0), shape = 16) +
  # geom_point(aes(x = lon, y = lat),  alpha = 0.7,
  #            shape = 1, size = 4,
  #            data = filter(location_effect, r_loc_Fp > quantile(location_effect$r_loc_Fp, 0.95, na.rm = TRUE))) +
  scale_color_gradient(low = "green", high = "red", name = "") +
  coord_sf(ylim = c(-35, 35), expand = FALSE) +
  geom_text(aes(x = -175, y = 30, label = "P excretion (g P/m²day)"), size = 3, hjust = 0) +
  scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
  theme_worldmap() + 
  theme(legend.position = "bottom", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(vjust = -8, hjust = 0.005),
        plot.margin = unit(c(0.0001,0.0001,0.0001,0.0001), units = , "cm"))


res_long <- pivot_longer(residuals, 2:6, values_to = "residual") %>%
  mutate(fun = case_when(name == "Fn_r" ~ "Fn",
                         name == "Fp_r" ~ "Fp",
                         name == "Gc_r" ~ "Gc",
                         name == "I_herb_r" ~ "Iherb",
                         name == "I_pisc_r" ~ "Ipisc")) %>%
  select(-name)

con_long <- left_join(con_long, res_long)


ggplot(con_long) +
  geom_point(aes(x = fi, y = residual)) +
  facet_wrap(~fun, ncol = 1, scales = "free_y")
ggplot(con_long) +
  geom_point(aes(x = cl, y = residual)) +
  facet_wrap(~fun, ncol = 1, scales = "free_y")

vu <- data.frame(
  fun = c("Fn", "Fp", "Gc", "Iherb", "Ipisc"),
  m_fi = c(
    median(con$vuln_Fn_fi),
    median(con$vuln_Fp_fi),
    median(con$vuln_Gc_fi),
    median(con[con$vuln_Iherb_fi>0,]$vuln_Iherb_fi),
    median(con[con$vuln_Ipisc_fi>0,]$vuln_Ipisc_fi)
  ),
  m_cl = c(
    median(con$vuln_Fn_cl),
    median(con$vuln_Fp_cl),
    median(con$vuln_Gc_cl),
    median(con[con$vuln_Iherb_cl>0,]$vuln_Iherb_cl),
    median(con[con$vuln_Ipisc_cl>0,]$vuln_Ipisc_cl)
  )
) %>%
  mutate(m_fi = mean(m_fi),
         m_cl = mean(m_cl))

summary(c(
  (con$vuln_Fn_cl),
  (con$vuln_Fp_cl),
  (con$vuln_Gc_cl),
  (con[con$vuln_Iherb_cl>0,]$vuln_Iherb_cl),
  (con[con$vuln_Ipisc_cl>0,]$vuln_Ipisc_cl)
))

q_cl <- quantile(c(
  (con$vuln_Fn_cl),
  (con$vuln_Fp_cl),
  (con$vuln_Gc_cl),
  (con[con$vuln_Iherb_cl>0,]$vuln_Iherb_cl),
  (con[con$vuln_Ipisc_cl>0,]$vuln_Ipisc_cl)
), c(0.2, 0.4, 0.6, 0.8))
q_fi <- quantile(c(
  (con$vuln_Fn_fi),
  (con$vuln_Fp_fi),
  (con$vuln_Gc_fi),
  (con[con$vuln_Iherb_fi>0,]$vuln_Iherb_fi),
  (con[con$vuln_Ipisc_fi>0,]$vuln_Ipisc_fi)
), c(0.2, 0.4, 0.6, 0.8))

test <- con_long %>%
  left_join(res_long) %>%
  left_join(vu) %>%
  mutate(cat = case_when(fi > m_fi &
                           cl > m_cl ~ "High vulnerability to both",
                         cl > m_cl ~ "High vulnerability to climate change",
                         fi > m_fi ~ "High vulnerability to fishing",
                         TRUE ~ "Low vulnerability to both")) %>%
  mutate(high = residual > 0) %>%
  group_by(cat, fun) %>%
  dplyr::summarise(n_high = sum(high), n = length(unique(transect_id))) %>%
  ungroup() %>% dplyr::group_by(fun) %>%
  mutate(ntot = sum(n)) %>%
  mutate(prop_high = n_high / ntot, prop = n/ntot)

ggplot(test) +
  geom_bar(aes(x = fun, y = prop, fill = fun), 
           alpha = 0.2, stat = "identity") +
  geom_bar(aes(x = fun, y = prop_high, fill = fun), 
           alpha = 1, stat = "identity") +
  #geom_point(aes(x = fun, y = prop/2), size = 1) +
  facet_wrap(~cat) +
  scale_fill_fish_d(option = "Callanthias_australis",
                    labels = c("N excretion", 
                               "P excretion", 
                               "Production",
                               "Herbivory",
                               "Piscivory")) +
  labs(y = "Proportion of communities", x = "", fill = "") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) 

test2_fi <- con_long %>%
  left_join(res_long) %>%
  left_join(vu) %>%
  mutate(cat_fi = case_when(fi > q_fi[4] ~ "Very high",
                            fi > q_fi[3] ~ "High",
                            fi > q_fi[2] ~ "Medium",
                            fi > q_fi[1] ~ "Low",
                            TRUE ~ "Very low")) %>%
  mutate(high = residual > 0) %>%
  group_by(cat_fi, fun) %>%
  dplyr::summarise(n_high_fi = sum(high), n_fi = length(unique(transect_id))) %>%
  ungroup() %>% dplyr::group_by(fun) %>%
  mutate(ntot_fi = sum(n_fi)) %>%
  mutate(prop_high_fi = n_high_fi / ntot_fi, prop_fi = n_fi/ntot_fi) %>%
  mutate(cat = factor(cat_fi, levels = c("Very low", "Low", "Medium", "High", "Very high")))
test2_cl <- con_long %>%
  left_join(res_long) %>%
  left_join(vu) %>%
  mutate(cat_cl = case_when(cl > q_cl[4] ~ "Very high",
                            cl > q_cl[3] ~ "High",
                            cl > q_cl[2] ~ "Medium",
                            cl > q_cl[1] ~ "Low",
                            TRUE ~ "Very low")) %>%
  mutate(high = residual > 0) %>%
  group_by(cat_cl, fun) %>%
  dplyr::summarise(n_high_cl = sum(high), n_cl = length(unique(transect_id))) %>%
  ungroup() %>% dplyr::group_by(fun) %>%
  mutate(ntot_cl = sum(n_cl)) %>%
  mutate(prop_high_cl = n_high_cl / ntot_cl, prop_cl = n_cl/ntot_cl) %>%
  mutate(cat = factor(cat_cl, levels = c("Very low", "Low", "Medium", "High", "Very high")))

p1 <- ggplot(test2_fi) +
  # geom_bar(aes(x = cat_fi, y = -(prop - prop_high), fill = fun), 
  #          alpha = 0.5, stat = "identity", position = "dodge") +
  geom_bar(aes(x = cat, y = prop_fi, fill = fun), 
           alpha = 1, stat = "identity", position = "dodge") +
  scale_fill_fish_d(option = "Callanthias_australis",
                    labels = c("N excretion", 
                               "P excretion", 
                               "Production",
                               "Herbivory",
                               "Piscivory")) +
  labs(y = "Proportion of communities", x = "Vulnerability to fishing", fill = "") +
  #scale_y_continuous(breaks = c(-0.1, 0, 0.1, 0.2, 0.3), labels = c("0.1", "0", "0.1", "0.2", "0.3")) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(), legend.position = "top")  

p2 <- ggplot(test2_cl) +
  # geom_bar(aes(x = cat_cl, y = -(prop - prop_high), fill = fun), 
  #          alpha = 0.5, stat = "identity", position = "dodge") +
  geom_bar(aes(x = cat, y = prop_cl, fill = fun), 
           alpha = 1, stat = "identity", position = "dodge") +
  scale_fill_fish_d(option = "Callanthias_australis",
                    labels = c("N excretion", 
                               "P excretion", 
                               "Production",
                               "Herbivory",
                               "Piscivory")) +
  #scale_y_continuous(breaks = c(-0.1, 0, 0.1, 0.2, 0.3), labels = c("0.1", "0", "0.1", "0.2", "0.3")) +
  labs(y = "Proportion of communities", x = "Vulnerability to climate change", fill = "") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(), legend.position = "none")  

p1 + p2 + plot_layout(nrow = 2)

double <- 
  con_long %>%
  left_join(res_long) %>%
  left_join(vu) %>%
  mutate(cat_double = case_when(fi > q_fi[3] & cl > q_cl[3] ~ "High",
                                fi < q_fi[2] & cl < q_cl[2] ~ "Low",
                                TRUE ~ "rest")) %>%
  group_by(cat_double, fun) %>%
  dplyr::summarise( n = length(unique(transect_id))) %>%
  ungroup() %>% dplyr::group_by(fun) %>%
  mutate(ntot = sum(n)) %>%
  mutate(prop = n/ntot) %>%
  filter(!cat_double == "rest") %>%
  mutate(cat_double = factor(cat_double, levels = c("Low", "High")))


p3 <- ggplot(double) +
  geom_bar(aes(x = cat_double, y = prop, fill = fun), 
           alpha = 1, stat = "identity", position = "dodge", width = 0.5) +
  scale_fill_fish_d(option = "Callanthias_australis",
                    labels = c("N excretion", 
                               "P excretion", 
                               "Production",
                               "Herbivory",
                               "Piscivory")) +
  #scale_y_continuous(breaks = c(-0.1, 0, 0.1, 0.2, 0.3), labels = c("0.1", "0", "0.1", "0.2", "0.3")) +
  labs(y = "Proportion of communities", x = "Vulnerability to both", fill = "") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor = element_blank(), 
        
        legend.position = "none")  

p3

p1 + p2 + p3 + plot_layout(ncol = 1)
ggsave("output/plots/figure_4_vuln_com.png", width = 8, height = 10)


test2 <- test %>%
  mutate(cat = case_when(fi > m_fi &
                           vl > m_cl ~ "both",
                         fi < m_fi &
                           cl > m_cl ~ "climate",
                         fi > m_fi &
                           cl < m_cl ~ "fishing",
                         TRUE ~ "none")) %>%
  mutate(high = Gc_r > 0) %>%
  group_by(cat) %>%
  dplyr::summarise(n_high = sum(Gc_r > 0), n = n()) %>%
  mutate(prop_high = n_high / n, prop = n/sum(n))
test2 <- test %>%
  mutate(cat = case_when(vuln_Fp_fi > m_fi &
                           vuln_Fp_cl > m_cl ~ "both",
                         vuln_Fp_fi < m_fi &
                           vuln_Fp_cl > m_cl ~ "climate",
                         vuln_Fp_fi > m_fi &
                           vuln_Fp_cl < m_cl ~ "fishing",
                         TRUE ~ "none")) %>%
  mutate(high = Fp_r >0) %>%
  group_by(cat) %>%
  dplyr::summarise(n_high = sum(Fp_r > 0), n = n()) %>%
  mutate(prop_high = n_high / n, prop = n/sum(n))
test2 <- test %>%
  mutate(cat = case_when(vuln_Fn_fi > m_fi &
                           vuln_Fn_cl > m_cl ~ "both",
                         vuln_Fn_fi < m_fi &
                           vuln_Fp_cl > m_cl ~ "climate",
                         vuln_Fn_fi > m_fi &
                           vuln_Fn_cl < m_cl ~ "fishing",
                         TRUE ~ "none")) %>%
  mutate(high = Fn_r >0) %>%
  group_by(cat) %>%
  dplyr::summarise(n_high = sum(Fn_r > 0), n = n()) %>%
  mutate(prop_high = n_high / n, prop = n/sum(n))

test <- filter(con, fun == "Ipisc") %>%
  left_join(residuals) %>%
  group_by(sites, locality) %>%
  summarise_if(is.numeric, median) %>%
  group_by(locality) %>%
  summarise_if(is.numeric, median)


ggplot(test) + 
  geom_sf(data = world, color = NA, fill = "lightgrey") +
  geom_point(aes(x = lon, y = lat, color = (value), 
                 size = get_scores(value)),  alpha = 0.8,
             position = position_jitter(0,0), shape = 16) +
  # geom_point(aes(x = lon, y = lat),  alpha = 0.7,
  #            shape = 1, size = 4,
  #            data = filter(location_effect, r_loc_Fp > quantile(location_effect$r_loc_Fp, 0.95, na.rm = TRUE))) +
  scale_color_gradient(low = "green", high = "red", name = "", trans = "log10") +
  coord_sf(ylim = c(-35, 35), expand = FALSE) +
  geom_text(aes(x = -175, y = 30, label = "P excretion (g P/m²day)"), size = 3, hjust = 0) +
  scale_size_continuous(range = c(0.5, 3), guide = FALSE) +
  theme_worldmap() + 
  theme(legend.position = "bottom", 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.title = element_text(vjust = -8, hjust = 0.005),
        plot.margin = unit(c(0.0001,0.0001,0.0001,0.0001), units = , "cm"))


##### residuals #####
loadd(residuals)
human <- read_csv("data/humanPopulation.csv") %>%
  select(-lat, -lon, - locality, -region )%>%
  as.data.frame()

res2 <- residuals %>%
  drop_na() %>%
  mutate(hot = 
           Fn_r > mean(Fn_r) & 
           Fp_r > mean(Fp_r) &
           Gc_r > mean(Gc_r) &
           I_herb_r > mean(I_herb_r) &
           I_pisc_r > mean(I_pisc_r),
         cold = 
           Fn_r < median(Fn_r) & 
           Fp_r < median(Fp_r) &
           Gc_r < median(Gc_r) &
           I_herb_r < median(I_herb_r) &
           I_pisc_r < median(I_pisc_r),
         all_q = 
           Fn_r > uquar(Fn_r) & 
           Fp_r > uquar(Fp_r) &
           Gc_r > uquar(Gc_r) &
           I_herb_r > uquar(I_herb_r) &
           I_pisc_r > uquar(I_pisc_r))
  #        ) %>%
  # left_join(summary_transect)

uquar <- function(x) {
  quantile(x, 0.75)
}

sum(res2$hot) /nrow(res2)
sum(res2$cold) /nrow(res2)
sum(res2$all_q) /nrow(res2)

idhot <- res2 %>%
  filter(hot == TRUE) %>%
  left_join(summary_transect) %>%
  left_join(human)

idcold <- res2 %>%
  filter(cold == TRUE) %>%
  left_join(summary_transect) %>%
  left_join(human)

summary(log1p(idhot$pop_50))
summary(log1p(idcold$pop_50))

ggplot() +
  geom_density(aes(x = ((idhot$size_m))), fill = "red", alpha = 0.5) +
  geom_density(aes(x = ((idcold$size_m))), fill = "blue", alpha = 0.5) +
  geom_density(aes(x = ((res2$size_m))), fill = "grey", alpha = 0.5) 
  
ggplot() +
  geom_density(aes(x = ((idhot$troph_m))), fill = "red", alpha = 0.5) +
  geom_density(aes(x = ((idcold$troph_m))), fill = "blue", alpha = 0.5) +
  geom_density(aes(x = ((res2$troph_m))), fill = "grey", alpha = 0.5) 
ggplot() +
  geom_density(aes(x = ((idhot$troph_q3))), fill = "red", alpha = 0.5) +
  geom_density(aes(x = ((idcold$troph_q3))), fill = "blue", alpha = 0.5) +
  geom_density(aes(x = ((res2$troph_q3))), fill = "grey", alpha = 0.5) 

ggplot() +
  geom_density(aes(x = ((idhot$imm_m))), fill = "red", alpha = 0.5) +
  geom_density(aes(x = ((idcold$imm_m))), fill = "blue", alpha = 0.5) +
  geom_density(aes(x = ((res2$imm_m))), fill = "grey", alpha = 0.5) 

ggplot() +
  geom_density(aes(x = ((idhot$troph_q1))), fill = "red", alpha = 0.5) +
  geom_density(aes(x = ((idcold$troph_q1))), fill = "blue", alpha = 0.5) +
  geom_density(aes(x = ((res2$troph_q1))), fill = "grey", alpha = 0.5) 









loadd(vulnerability)
loadd(contributions)
loadd(herb_pisc)
loadd(summary_transect_imp)

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
  left_join(summary_transect_imp)

long <- con %>% 
  pivot_longer(cols = starts_with('vuln'), 
               names_to = c("null", "fun","type"), names_sep = "_") %>%
  pivot_wider(names_from = type, values_from = value)


rescor <- sum$rescor_pars
sum$random$locality

test <- left_join(degree_dominance, summary_transect_imp)

library(hitandrun)
te <- (simplex.sample(20, 100, sort=FALSE))[[1]]
rowSums(te)


ks_area <- function(data){
  
  if (length(data) == 0){
    return(NA)
  }else if (length(data) == 1){
    return(1)
  }else{
    data <- data.frame(data = data)
    drank <- data %>% 
      mutate(rank = dense_rank(desc(data))) %>%
      dplyr::arrange(rank) %>%
      mutate(cumsum = cumsum(data))
    
    A <- rank_area(drank$cumsum)
    
    r <- nrow(data) # species richness
    A_max <- 1 * (r - 1) # maximum if first species performs 100% (rectangle)
    A_min <- (r^2-1)/(2*r) # minimum if all sp contribute equally (trapezium)
    
    dd <- (A - A_min)/(A_max - A_min)
    
    return(dd)
  }
}


rank_area <- function(y) {
  n <- length(y)
  a <- y[2:n]
  b <- y[1:(n-1)]
  area <- sum(a + b)/2
  return(area)
}

list <- parallel::mclapply(21, function(j){
  te <- (simplex.sample(j, 1000, sort=FALSE))[[1]]
   dd <- (apply(te, 1, ks_area))
   data.frame(nspec = j, dd = dd)
}, mc.cores = 50)

vec <- plyr::ldply(list)

summary(vec$dd)

plot(5:100, vec)

hist(vec$dd)

ggplot(test) +
  geom_point(aes(x = dd_Fn, y = dd_Gc))

fit <- brm(dd_Fn ~ log(nspec),
    data = test, cores = 4, family = "beta",
    backend = "cmdstanr",
    threads = threading(10))


fit_sim <- update(fit, formula. = dd ~ log(nspec), newdata = vec)

fit_ddFp <- update(fit, formula. = dd_Fp ~ log(nspec), newdata = test)
fit_ddGc <- update(fit, formula. = dd_Gc ~ log(nspec), newdata = test)
fit_ddH <- update(fit, formula. = dd_I_herb ~ log(nspec), newdata = test, family = "zero_one_inflated_beta")
fit_ddP <- update(fit_ddH, formula. = dd_I_pisc ~ log(nspec), newdata = test)


summary(fit)

a <- conditional_effects(fit)
b <- conditional_effects(fit_sim)
c <- conditional_effects(fit_ddFp)
d <- conditional_effects(fit_ddGc)
e <- conditional_effects(fit_ddH)
f <- conditional_effects(fit_ddP)

ggplot() +
  geom_ribbon(data = a$nspec, aes(x = nspec, ymin = lower__, ymax = upper__), alpha = 0.6) +
  geom_ribbon(data = b$nspec, aes(x = nspec, , ymin = lower__, ymax = upper__), alpha = 0.6) +
  geom_line(data = a$nspec, aes(x = nspec, y = estimate__)) +
  geom_line(data = b$nspec, aes(x = nspec, y = estimate__), color = "red") +
  geom_ribbon(data =c$nspec, aes(x = nspec, ymin = lower__, ymax = upper__), alpha = 0.6) +
  geom_ribbon(data = d$nspec, aes(x = nspec, , ymin = lower__, ymax = upper__), alpha = 0.6) +
  geom_line(data = c$nspec, aes(x = nspec, y = estimate__)) +
  geom_line(data = d$nspec, aes(x = nspec, y = estimate__)) +
  geom_ribbon(data = e$nspec, aes(x = nspec, ymin = lower__, ymax = upper__), alpha = 0.6) +
  geom_ribbon(data = f$nspec, aes(x = nspec, , ymin = lower__, ymax = upper__), alpha = 0.6) +
  geom_line(data = e$nspec, aes(x = nspec, y = estimate__)) +
  geom_line(data = f$nspec, aes(x = nspec, y = estimate__)) +
  scale_x_continuous(lim = c(10,90)) +
  labs(x = "Richness", y = "Degree of dominance")


nd = data.frame(nspec = 50)

pred1 <- fitted(fit, newdata = nd, nsamples = 1000, summary = FALSE, method = "predict")
predsim <- fitted(fit_sim, newdata = nd, nsamples = 1000, summary = FALSE, method = "predict")
pred2 <- fitted(fit_ddFp, newdata = nd, nsamples = 1000, summary = FALSE, method = "predict")
pred3 <- fitted(fit_ddGc, newdata = nd, nsamples = 1000, summary = FALSE, method = "predict")
pred4 <- fitted(fit_ddH, newdata = nd, nsamples = 1000, summary = FALSE, method = "predict")
pred5 <- fitted(fit_ddP, newdata = nd, nsamples = 1000, summary = FALSE, method = "predict")

pred1 <- fitted(fit)
pred2 <- fitted(fit_ddFp)
pred3 <- fitted(fit_ddGc)
pred4 <- fitted(fit_ddH)
pred5 <- fitted(fit_ddP)
predsim <- fitted(fit_sim)



###### degree_dominance ######
  
  hist(vec$dd)
  
  fitdd1 <- brm(dd_Fn ~ 1 + (1|sites) + (1|locality), 
                cores = 4, family = "beta",
                backend = "cmdstanr",
                threads = threading(10),
                data = test)
  
  fitdd2 <- update(fitdd1, dd_Fp ~ 1 + (1|sites) + (1|locality), 
                   newdata = test)
  fitdd3 <- update(fitdd1, dd_Gc ~ 1 + (1|sites) + (1|locality), 
                   newdata = test)
  fitdd4 <- update(fitdd1, dd_I_herb ~ 1 + (1|sites) + (1|locality), 
                   newdata = test, family = "zero_one_inflated_beta")
  fitdd5 <- update(fitdd1, dd_I_pisc ~ 1 + (1|sites) + (1|locality), 
                   newdata = test, family = "zero_one_inflated_beta")
  


  inverse_logit <- function(x){
    exp(x)/(1+exp(x))
  }
  

  
  pr1 <- fitted(fitdd1, newdata = data.frame(locality = NA, sites = NA))
  pr2 <- fitted(fitdd2, newdata = data.frame(locality = NA, sites = NA))
  pr3 <- fitted(fitdd3, newdata = data.frame(locality = NA, sites = NA))
  pr4 <- fitted(fitdd4, newdata = data.frame(locality = NA, sites = NA))
  pr5 <- fitted(fitdd5, newdata = data.frame(locality = NA, sites = NA))
  
  cols <- fish(option = "Callanthias_australis", n = 5)
  
  plot1 <-
  ggplot() + 
    geom_jitter(aes(y = "e", x = dd_Fn), color = cols[1],
                data = test, alpha = 0.1, size = 0.1) +
    geom_jitter(aes(y = "d", x = dd_Fp), color = cols[2],
                data = test, alpha = 0.1, size = 0.1) +
    geom_jitter(aes(y = "c", x = dd_Gc), color = cols[3],
                data = test, alpha = 0.1, size = 0.1) +
    geom_jitter(aes(y = "b", x = dd_I_herb), color = cols[4],
                data = test, alpha = 0.1, size = 0.1, width = 0.01) +
    geom_jitter(aes(y = "a", x = dd_I_pisc), color = cols[5],
                data = test, alpha = 0.1, size = 0.1, width = 0.01) +
    geom_pointrange(aes(y = "e", x = pr1[1], xmin = pr1[3], xmax = pr1[4])) +
    geom_pointrange(aes(y = "d", x = pr2[1], xmin = pr2[3], xmax = pr2[4])) + 
    geom_pointrange(aes(y = "c", x = pr3[1], xmin = pr3[3], xmax = pr3[4])) +
    geom_pointrange(aes(y = "b", x = pr4[1], xmin = pr4[3], xmax = pr4[4])) +
    geom_pointrange(aes(y = "a", x = pr5[1], xmin = pr5[3], xmax = pr5[4])) +
    geom_vline(xintercept = mean(vec$dd), linetype = 2) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_discrete(labels = c("Piscivory", 
                                "Herbivory", 
                                "Production", 
                                "P excretion", 
                                "N excretion")) +
    labs(x = "Degree of dominance (communities)", y = "") +
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
          legend.position = "none", 
          panel.border = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line())
  
  spi_sp_long <- freq_dominance %>%
    pivot_longer(c(Fn_i, Fp_i, Gc_i, I_herb_i, I_pisc_i)) %>%
    drop_na() %>%
    filter(occ>4, value > 0) 
  
  a <- spi_sp_long %>% group_by(name) %>%
    summarize(n = n())
  
  b <- freq_dominance %>%
    pivot_longer(c(Fn_i, Fp_i, Gc_i, I_herb_i, I_pisc_i)) %>%
    drop_na() %>%
    filter(occ>4) %>% 
    group_by(name) %>%
    summarize(ntot = n())
  
  ab <- left_join(a,b) %>%
    mutate(prop = n/ntot)
  

  plot2 <- 
    ggplot(ab) +
    geom_col(aes(x = prop, y = reorder(name, desc(name)), fill = name),
             alpha = 0.9) +
    geom_text(aes(x = prop - 0.1, y = reorder(name, desc(name)), 
                  label = paste0(round(100*prop), " %")),
              size = 4) +
    scale_fill_fish_d(option = "Callanthias_australis") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_discrete(labels = c("Piscivory", 
                                "Herbivory", 
                                "Production", 
                                "P excretion", 
                                "N excretion")) +
    scale_color_fish_d(option = "Callanthias_australis") +
    theme(legend.position = "none", axis.title.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(x = "Proportion of species being dominant") +
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
          legend.position = "none", 
          panel.border = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line())
  
  fit_occ <- brm(value ~ name,
                 cores = 4, family = "beta",
                 backend = "cmdstanr",
                 threads = threading(12),
                 data = spi_sp_long)
  
  pr <- fitted(fit_occ, newdata = data.frame(name = unique(spi_sp_long$name)))
  
  plot3 <- 
  ggplot() + 
    geom_jitter(aes(y = reorder(name, desc(name)), x = value, color = name),
                data = spi_sp_long, alpha = 0.5, size = 0.2) +
    geom_pointrange(aes(y = 5, x = pr[1,1], xmin = pr[1,3], xmax = pr[1,4])) +
    geom_pointrange(aes(y = 4, x = pr[2,1], xmin = pr[2,3], xmax = pr[2,4])) +
    geom_pointrange(aes(y = 3, x = pr[3,1], xmin = pr[3,3], xmax = pr[3,4])) +
    geom_pointrange(aes(y = 2, x = pr[4,1], xmin = pr[4,3], xmax = pr[4,4])) +
    geom_pointrange(aes(y = 1, x = pr[5,1], xmin = pr[5,3], xmax = pr[5,4])) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_discrete(labels = c("Piscivory", 
                                "Herbivory", 
                                "Production", 
                                "P excretion", 
                                "N excretion")) +
    scale_color_fish_d(option = "Callanthias_australis") +
    theme(legend.position = "none", axis.title.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    labs(x = "Frequency of dominance per species") +
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
          legend.position = "none", 
          panel.border = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.grid.major.y = element_blank(),
          axis.line.x = element_line())
  
  plot <- plot1 + plot2 + plot3 + plot_annotation(tag_levels = "a")
  ggsave("output/plots/fig4_dd.png", plot, width = 9, height = 4)
 
  
   ####### vuln ######
  
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
    filter(vuln_Ipisc_cl > 0, vuln_Iherb_cl > 0) 
  
  con_long <- con %>%
    pivot_longer(cols = 2:13) %>%
    separate(name, into = c("v", "fun", "impact"), sep = "_", remove = TRUE) %>%
    pivot_wider(names_from = impact, values_from = value) %>%
    filter(!fun == "bm")
  
  ggplot(con) +
    geom_point(aes(x = vuln_Fp_cl, y = vuln_Fp_fi)) +
    geom_smooth(aes(x = vuln_Fp_cl, y = vuln_Fp_fi))
  
  con <- left_join(con, summary_transect_imp)
  ggplot(con) +
    geom_point(aes(x = log(nspec), y = vuln_Fp_cl)) +
    geom_smooth(aes(x = log(nspec), y = vuln_Fp_cl))
  
  con_long <- left_join(con_long, select(summary_transect_imp, sites, locality, transect_id, nspec, biomass_tot))
  
  
  fit_v <- brm(mvbind(fi, cl) ~ fun + fun:log(nspec) + fun:log(biomass_tot),#+ (1|sites) + (1|locality),
  data = con_long, cores = 4,
  backend = "cmdstanr",
  threads = threading(15))
  
  fe <- conditional_effects(fit_v)
  fe1 <- fe[[4]]
  fe2 <- fe[[3]]
  fe3 <- fe[[5]]
  fe4 <- fe[[6]]
  
  v1 <-
  ggplot(fe1) +
    geom_line(aes(x = log(effect1__), y = estimate__, color = effect2__)) +
    labs(x = "log(Richness)", y = "Vulnerability to climate change") +
    geom_ribbon(aes(x = log(effect1__), ymin = lower__, ymax = upper__, fill = effect2__), alpha = 0.4) +
    scale_fill_fish_d(name = "Function",
                       option = "Callanthias_australis",
                       labels = c("N excretion",
                                  "P excretion",
                                  "Production", 
                                  "Herbivory",
                                  "Piscivory" )) +
    scale_color_fish_d(name = "Function",
                       option = "Callanthias_australis",
                       labels = c("N excretion",
                                  "P excretion",
                                  "Production", 
                                  "Herbivory",
                                  "Piscivory" )) +
    theme_bw() +
    theme(panel.grid = element_blank())
    
  v2 <- 
    ggplot(fe2) +
    geom_line(aes(x = log(effect1__), y = estimate__, color = effect2__)) +
    geom_ribbon(aes(x = log(effect1__), ymin = lower__, ymax = upper__, fill = effect2__), alpha = 0.4) +
    scale_fill_fish_d(name = "Function",
                      option = "Callanthias_australis",
                      labels = c("N excretion",
                                 "P excretion",
                                 "Production", 
                                 "Herbivory",
                                 "Piscivory" )) +
    labs(x = "log(Richness)", y = "Vulnerability to fishing") +
    scale_color_fish_d(name = "Function",
                       option = "Callanthias_australis",
                       labels = c("N excretion",
                                  "P excretion",
                                  "Production", 
                                  "Herbivory",
                                  "Piscivory" )) +
    theme_bw() +
    theme(panel.grid = element_blank())
  v2
  v3 <- 
    ggplot(fe3) +
    geom_line(aes(x = log(effect1__), y = estimate__, color = effect2__)) +
    geom_ribbon(aes(x = log(effect1__), ymin = lower__, ymax = upper__, fill = effect2__), alpha = 0.4) +
    scale_fill_fish_d(name = "Function",
                      option = "Callanthias_australis",
                      labels = c("N excretion",
                                 "P excretion",
                                 "Production", 
                                 "Herbivory",
                                 "Piscivory" )) +
    labs(x = "log(Biomass)", y = "Vulnerability to fishing") +
    scale_color_fish_d(name = "Function",
                       option = "Callanthias_australis",
                       labels = c("N excretion",
                                  "P excretion",
                                  "Production", 
                                  "Herbivory",
                                  "Piscivory" )) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  v4 <- 
    ggplot(fe4) +
    geom_line(aes(x = log(effect1__), y = estimate__, color = effect2__)) +
    labs(x = "log(Biomass)", y = "Vulnerability to climate change") +
    geom_ribbon(aes(x = log(effect1__), ymin = lower__, ymax = upper__, fill = effect2__), alpha = 0.4) +
    scale_fill_fish_d(name = "Function",
                      option = "Callanthias_australis",
                      labels = c("N excretion",
                                 "P excretion",
                                 "Production", 
                                 "Herbivory",
                                 "Piscivory" )) +
    scale_color_fish_d(name = "Function",
                       option = "Callanthias_australis",
                       labels = c("N excretion",
                                  "P excretion",
                                  "Production", 
                                  "Herbivory",
                                  "Piscivory" )) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  plot <- 
  v2 + v1 + v3 + v4 + plot_annotation(tag_levels = "a") +
    plot_layout(guides = 'collect') &
    theme(legend.position = "bottom")
  
  ggsave("output/plots/fig5.png", width = 8, height = 6)
                          

  summary(fit_v)
  
fe[[1]]

pv1 <- 
  ggplot(fe[[1]]) +
  geom_jitter(aes(x = fi, y = reorder(fun, desc(fun)), color = fun), data = con_long, height = 0.4, alpha = 0.1, size = 0.1) +
  geom_pointrange(aes(y = (effect1__), x = estimate__, xmin = lower__, xmax = upper__)) +
  labs(y = "", x = "Vulnerability to fishing") +
  scale_color_fish_d(name = "Function",
                     option = "Callanthias_australis",
                     labels = c("N excretion",
                                "P excretion",
                                "Production", 
                                "Herbivory",
                                "Piscivory" )) +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))  +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.position = "none", 
        panel.border = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line())
pv1

pv2 <- 
  ggplot(fe[[2]]) +
  geom_jitter(aes(x = cl, y = reorder(fun, desc(fun)), color = fun), data = con_long, height = 0.4, alpha = 0.1, size = 0.1) +
  geom_pointrange(aes(y = (effect1__), x = estimate__, xmin = lower__, xmax = upper__)) +
  labs(y = "", x = "Vulnerability to climate change") +
  guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))  +
  scale_fill_fish_d(name = "Function",
                    option = "Callanthias_australis",
                    labels = c("N excretion",
                               "P excretion",
                               "Production", 
                               "Herbivory",
                               "Piscivory" )) +
  scale_color_fish_d(name = "Function",
                     option = "Callanthias_australis",
                     labels = c("N excretion",
                                "P excretion",
                                "Production", 
                                "Herbivory",
                                "Piscivory" )) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
        legend.position = "none", 
        panel.border = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(),
        axis.line.x = element_line())
pv2

plot <-
pv1 + pv2 + plot_layout( nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "right")
ggsave("output/plots/vuln_jitter.png", plot, width = 5, height = 4)
ggsave("output/plots/vuln_jitter.pdf", plot, width = 5, height = 4)

###### imputation ######
fit_imp <- brm(bf(mvbind(I_herb, I_pisc)|mi() ~ log(nspec) + log(biomass_tot)) +
                set_rescor(TRUE),
             backend = "cmdstanr",
             threads = threading(12),
             cores = 4,
             data = summary_transect_complete)

prior_summary(fit_imp)
stancode(fit_imp)
summary(fit_imp)

re <- ranef(mod_mv_siteloc)
rloc <- re$locality[,,1]
rsite <- re$sites[,,1]

get_variables(mod_mf_siteloc)

ran <- 
mod_mf_siteloc %>%
  spread_draws(n = 1000, r_sites[site, Intercept], r_locality[loc, Intercept]) %>%
  filter(site == "red_sea_2", loc == "red_sea") %>%
  mutate(r_tot = r_sites + r_locality) %>%
  group_by(site, loc) %>%
  summarize(rloc_m = mean(r_locality), 
            rloc_lb = quantile(r_locality, 0.025),
            rloc_ub = quantile(r_locality, 0.975),
            rsite_m = mean(r_sites), 
            rsite_lb = quantile(r_sites, 0.025),
            rsite_ub = quantile(r_sites, 0.975),
            rtot_m = mean(r_tot), 
            rtot_lb = quantile(r_tot, 0.025),
            rtot_ub = quantile(r_tot, 0.975),
            )
ran



loadd(mod_mvfun_bm)

nd <- data.frame(biomass_tot = seq(50,500, 4), mean = 26, locality = NA, sites = NA)

fit <- predict(mod_mvfun_bm, newdata = nd, summary = F, nsamples = 1000)
ef <- fit

ef[,,1] <- normalize(exp(ef[,,1]))
ef[,,2] <- normalize(exp(ef[,,2])) 
ef[,,3] <- normalize(exp(ef[,,3])) 
ef[,,4] <- normalize(exp(ef[,,4])) 
ef[,,5] <- normalize(exp(ef[,,5])) 

mf <- apply(ef, 1:2, geomean)
mf[,] <- normalize(mf[,])

f1_m <- apply(ef[,,1], 2, mean) 
mf_m <- apply(mf, 2, mean) 
mf_lb <- apply(mf, 2, quantile, 0.025) 
mf_ub <- apply(mf, 2, quantile, 0.975) 

ggplot() +
  geom_ribbon(aes(x = nd$biomass_tot, ymin = mf_lb, ymax = mf_ub), alpha = 0.5)+
  geom_line(aes(x = nd$biomass_tot, y= (mf_m))) 

lm(log(mf_m)~log(nd$biomass_tot))
lm(log(f1_m)~log(nd$biomass_tot))

sub <- ef[,1,]
hist(sub[,1])

test <- apply(sub, 1, geomean)

summary(test)

sub <- ef[,100,]
hist(sub[,1])

test <- apply(sub, 1, geomean)

summary(test)




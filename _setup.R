packages <- c(
  "purrr",
  "tidybayes",
  "ggplot2",
  "fishflux",
  "Rphylopars",
  "dplyr",
  "tidyr",
  "drake",
  "e1071",
  "brms",
  "utils",
  "patchwork",
  "fishualize",
  "ggnewscale",
  "tidytree",
  "ggtree",
  "forcats",
  "ggrepel",
  "tidyr",
  "readr",
  "Rcpp",
  "broom",
  "flextable",
  "officer",
  "drake")

to_install <- packages[!packages %in% installed.packages()]

install.packages(to_install)


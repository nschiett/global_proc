make_table_mod_mvfun_bm <- function(brmsfit) {
  
  sum <- summary(brmsfit)
  tab_fixed <- sum$fixed %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") %>%
    tidyr::separate(name, into = c("response", "term", "extra"), sep = "_") %>%
    dplyr::mutate(term = paste0(term, extra)) %>%
    dplyr::mutate(term = str_replace_all(term, "NA", "")) %>%
    dplyr::mutate(response = case_when(response == "logFn" ~ "log(N excretion)",
                                       response == "logFp" ~ "log(P excretion)",
                                       response == "logGc" ~ "log(production)",
                                       response == "logIherb" ~ "log(herbivory)",
                                       response == "logIpisc" ~ "log(piscivory)"),
                  term = case_when(term == "logbiomasstot" ~ "log(biomass)",
                                   term == "mean" ~ "sst",
                                   TRUE ~ term)) %>%
    dplyr::select(-extra) %>%
    dplyr::select(-7, -8, -9) %>%
    dplyr::mutate(type = "fixed")
  
  tab_sigma <- sum$spec_pars %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") %>%
    tidyr::separate(name, into = c("term", "response"), sep = "_") %>%
    dplyr::mutate(response = case_when(response == "logFn" ~ "log(N excretion)",
                                       response == "logFp" ~ "log(P excretion)",
                                       response == "logGc" ~ "log(production)",
                                       response == "logIherb" ~ "log(herbivory)",
                                       response == "logIpisc" ~ "log(piscivory)")) %>%
    dplyr::select(-7, -8, -9)
  
  tab_rescor<- sum$rescor_pars %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) 
  
  tab_cor_loc <- sum$random$locality %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) %>%
    dplyr::mutate(type = "random locality")
  
  tab_cor_site <- sum$random$sites %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) %>%
    dplyr::mutate(type = "random sites")
  
  tab <- full_join(tab_fixed, tab_sigma) %>%
    full_join(tab_rescor) %>%
    full_join(tab_cor_loc) %>%
    full_join(tab_cor_site) %>%
    mutate_all(tidyr::replace_na, "_") %>%
    dplyr::select(-Est.Error) %>%
    dplyr::mutate_if(is.numeric, round, 4) 
  
  write_csv(tab, here::here("output", "data", "model_output_mod_mvfun_bm.csv"))
  tab
}

make_table_mod_mvfun_com <- function(brmsfit) {
  
  sum <- summary(brmsfit)
  tab_fixed <- sum$fixed %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") %>%
    tidyr::separate(name, into = c("response", "term", "extra"), sep = "_") %>%
    dplyr::mutate(term = paste0(term,"_", extra)) %>%
    dplyr::mutate(term = str_replace_all(term, "_NA", "")) %>%
    dplyr::mutate(term = str_replace_all(term, "standard", "")) %>%
    dplyr::mutate(response = case_when(response == "Fnst" ~ "log(N excretion)",
                                       response == "Fpst" ~ "log(P excretion)",
                                       response == "Gcst" ~ "log(production)",
                                       response == "Iherbst" ~ "log(herbivory)",
                                       response == "Ipiscst" ~ "log(piscivory)"),
                  term = case_when(term == "logbiomass_tot" ~ "log(biomass)",
                                   term == "mean" ~ "sst",
                                   TRUE ~ term)) %>%
    dplyr::mutate(term = str_replace_all(term, "_q1", "_lb")) %>%
    dplyr::mutate(term = str_replace_all(term, "_q3", "_ub")) %>%
    dplyr::mutate(term = str_replace_all(term, "_m", "_median")) %>%
    dplyr::select(-extra) %>%
    dplyr::select(-7, -8, -9) %>%
    dplyr::mutate(type = "fixed")
  
  tab_fixed
  
  tab_sigma <- sum$spec_pars %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") %>%
    tidyr::separate(name, into = c("term", "response"), sep = "_") %>%
    dplyr::mutate(response = case_when(response == "Fnst" ~ "log(N excretion)",
                                       response == "Fpst" ~ "log(P excretion)",
                                       response == "Gcst" ~ "log(production)",
                                       response == "Iherbst" ~ "log(herbivory)",
                                       response == "Ipiscst" ~ "log(piscivory)")) %>%
    dplyr::select(-7, -8, -9)
  
  tab_rescor<- sum$rescor_pars %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) 
  
  tab_cor_loc <- sum$random$locality %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) %>%
    dplyr::mutate(type = "random locality")
  
  tab_cor_site <- sum$random$sites %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) %>%
    dplyr::mutate(type = "random sites")
  
  tab <- full_join(tab_fixed, tab_sigma) %>%
    full_join(tab_rescor) %>%
    full_join(tab_cor_loc) %>%
    full_join(tab_cor_site) %>%
    mutate_all(tidyr::replace_na, "_") %>%
    dplyr::select(-Est.Error) %>%
    dplyr::mutate(across(c(3,4,5), as.numeric)) %>%
    dplyr::mutate_if(is.numeric, round, 4) 
  
  write_csv(tab, here::here("output", "data", "model_output_mod_mvfun_com.csv"))
  
  tab
  
}

make_table_mod_mvfun_com2 <- function(brmsfit) {
  
  sum <- summary(brmsfit)
  tab_fixed <- sum$fixed %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") %>%
    tidyr::separate(name, into = c("response", "term", "extra"), sep = "_") %>%
    dplyr::mutate(term = paste0(term,"_", extra)) %>%
    dplyr::mutate(term = str_replace_all(term, "_NA", "")) %>%
    dplyr::mutate(term = str_replace_all(term, "standard", "")) %>%
    dplyr::mutate(response = case_when(response == "logFn" ~ "log(N excretion)",
                                       response == "logFp" ~ "log(P excretion)",
                                       response == "logGc" ~ "log(production)",
                                       response == "logIherb" ~ "log(herbivory)",
                                       response == "logIpisc" ~ "log(piscivory)"),
                  term = case_when(term == "logbiomass_tot" ~ "log(biomass)",
                                   term == "mean" ~ "sst",
                                   TRUE ~ term)) %>%
    dplyr::mutate(term = str_replace_all(term, "_q1", "_lb")) %>%
    dplyr::mutate(term = str_replace_all(term, "_q3", "_ub")) %>%
    dplyr::mutate(term = str_replace_all(term, "_m", "_median")) %>%
    dplyr::select(-extra) %>%
    dplyr::select(-7, -8, -9) %>%
    dplyr::mutate(type = "fixed")
  
  tab_fixed
  
  tab_sigma <- sum$spec_pars %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") %>%
    tidyr::separate(name, into = c("term", "response"), sep = "_") %>%
    dplyr::mutate(response = case_when(response == "logFn" ~ "log(N excretion)",
                                       response == "logFp" ~ "log(P excretion)",
                                       response == "logGc" ~ "log(production)",
                                       response == "logIherb" ~ "log(herbivory)",
                                       response == "logIpisc" ~ "log(piscivory)")) %>%
    dplyr::select(-7, -8, -9)
  
  tab_rescor<- sum$rescor_pars %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) 
  
  tab_cor_loc <- sum$random$locality %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) %>%
    dplyr::mutate(type = "random locality")
  
  tab_cor_site <- sum$random$sites %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) %>%
    dplyr::mutate(type = "random sites")
  
  tab <- full_join(tab_fixed, tab_sigma) %>%
    full_join(tab_rescor) %>%
    full_join(tab_cor_loc) %>%
    full_join(tab_cor_site) %>%
    mutate_all(tidyr::replace_na, "_") %>%
    dplyr::select(-Est.Error) %>%
    dplyr::mutate(across(c(3,4,5), as.numeric)) %>%
    dplyr::mutate_if(is.numeric, round, 4) 
  
  write_csv(tab, here::here("output", "data", "model_output_mod_mvfun_com_absolute.csv"))
  
  tab
  
}

# make_table_mod_mf_com <- function(brmsfit) {
#   
#   sum <- summary(brmsfit)
#   tab_fixed <- sum$fixed %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("term") %>%
#     dplyr::mutate(term = str_replace_all(term, "_NA", "")) %>%
#     dplyr::mutate(term = str_replace_all(term, "standard", "")) %>%
#     dplyr::mutate(
#                   term = case_when(term == "logbiomass_tot" ~ "log(biomass)",
#                                    term == "mean" ~ "sst",
#                                    TRUE ~ term)) %>%
#     dplyr::mutate(term = str_replace_all(term, "_q1", "_lb")) %>%
#     dplyr::mutate(term = str_replace_all(term, "_q3", "_ub")) %>%
#     dplyr::mutate(term = str_replace_all(term, "_m", "_median")) %>%
#     dplyr::select(-6, -7, -8) %>%
#     dplyr::mutate(type = "fixed")
#   
#   tab_sigma <- sum$spec_pars %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("term") %>%
#     dplyr::select(-6, -7, -8)
#   
#   tab_cor_loc <- sum$random$locality %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("term") %>%
#     dplyr::select(-6, -7, -8) %>%
#     dplyr::mutate(type = "random locality")
#   
#   tab_cor_site <- sum$random$sites %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("term") %>%
#     dplyr::select(-6, -7, -8) %>%
#     dplyr::mutate(type = "random sites")
#   
#   
#   tab <- full_join(tab_fixed, tab_sigma) %>%
#     full_join(tab_cor_loc) %>%
#     full_join(tab_cor_site) %>%
#     mutate_all(tidyr::replace_na, "_") %>%
#     dplyr::select(-Est.Error) %>%
#     dplyr::mutate(across(c(2,3,4), as.numeric)) %>%
#     dplyr::mutate_if(is.numeric, round, 4) 
#   
#   tab
#   
# }


# make_table_mod_mf_siteloc <- function(brmsfit) {
#   
#   sum <- summary(brmsfit)
#  
#   tab_sigma <- sum$spec_pars %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("term") %>%
#     dplyr::select(-6, -7, -8)
#   
#   tab_cor_loc <- sum$random$locality %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("term") %>%
#     dplyr::select(-6, -7, -8) %>%
#     dplyr::mutate(type = "random locality")
#   
#   tab_cor_site <- sum$random$sites %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column("term") %>%
#     dplyr::select(-6, -7, -8) %>%
#     dplyr::mutate(type = "random sites")
#   
#   tab <- 
#     full_join(tab_sigma, tab_cor_loc) %>%
#     full_join(tab_cor_site) %>%
#     mutate_all(tidyr::replace_na, "_") %>%
#     dplyr::select(-Est.Error) %>%
#     dplyr::mutate(across(c(2,3,4), as.numeric)) %>%
#     dplyr::mutate_if(is.numeric, round, 4) 
#   
#   tab
#   
# }

make_table_mod_mv_siteloc <- function(brmsfit) {
  
  sum <- summary(brmsfit)
  
  tab_sigma <- sum$spec_pars %>%
    as.data.frame() %>%
    tibble::rownames_to_column("name") %>%
    tidyr::separate(name, into = c("term", "response"), sep = "_") %>%
    dplyr::mutate(response = case_when(response == "logFn" ~ "log(N excretion)",
                                       response == "logFp" ~ "log(P excretion)",
                                       response == "logGc" ~ "log(production)",
                                       response == "logIherb" ~ "log(herbivory)",
                                       response == "logIpisc" ~ "log(piscivory)")) %>%
    dplyr::select(-7, -8, -9)
  
  tab_rescor<- sum$rescor_pars %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) 
  
  tab_cor_loc <- sum$random$locality %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) %>%
    dplyr::mutate(type = "random locality")
  
  tab_cor_site <- sum$random$sites %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    dplyr::select(-6, -7, -8) %>%
    dplyr::mutate(type = "random sites")
  
  tab <- full_join(tab_sigma, tab_rescor) %>%
    full_join(tab_cor_loc) %>%
    full_join(tab_cor_site) %>%
    mutate_all(tidyr::replace_na, "_") %>%
    dplyr::select(-Est.Error) %>%
    dplyr::mutate(across(c(3,4,5), as.numeric)) %>%
    dplyr::mutate_if(is.numeric, round, 4) 
  
  write_csv(tab, here::here("output", "data", "model_output_mod_mvfun_siteloc.csv"))
  
  tab
}



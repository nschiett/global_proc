
plan <- drake_plan(
  
  ##### set up output #####
  out_folder = dir.create("output"),
  out_folder_data = dir.create("output/data"),
  out_folder_plots = dir.create("output/plots"),
  out_folder_text = dir.create("output/text"),
  
  ##### Bioenergetic model #####
  
  # load all parameters
  sptl = get_sptl(),
  kmax = get_kmax(),
  lw = get_lw(sptl),
  cnp = get_cnp(),
  cnpdiet = get_cnpdiet(),
  metpar = get_metpar(kmax, sptl),
  artr = get_artr(sptl),
  
  # combine parameters
  params = combine_params(sptl, kmax, lw, cnp, cnpdiet, metpar, artr),
  output_params = write.csv(params, "output/data/params_sst_glob.csv", row.names = FALSE),
  tfish_unique = get_tfish_unique(),
  
  # Run fishflux bioenergetic model
  cnpflux = run_fishflux(tfish_unique, params),
  output_cnpflux = write.csv(cnpflux, "output/data/flux_results.csv", row.names = FALSE),
  
  ##### Summarize transect level #####
  
  # load all transect data
  tfish = read.csv("data/tfish_corr.csv"),
  
  # summarize results per transect
  summary_transect = summarise_pertransect(tfish, cnpflux, params),
  output_fluxsum = write.csv(summary_transect, "output/data/flux_summary_transect.csv", row.names = FALSE),
  
  contributions = get_contributions(tfish, cnpflux, params, summary_transect),
  contributions_site = get_contributions_site(tfish, cnpflux, params, summary_transect),
  sp_loc = unique(select(ungroup(contributions), bioregion, species)),
  
  # herbivory, piscivory 
  herb_pisc = get_herb_pisc(tfish, cnpflux, params, summary_transect),
  
  # combine summary data, sst 
  summary_transect_complete = combine_summary(summary_transect, herb_pisc),
  
  ##### models #####
  
  # site loc
  mod_mv_siteloc = fit_mvfun_siteloc(summary_transect_complete),
  #mod_mf_siteloc = fit_mf_siteloc(summary_transect_imp),
  
  # add biomass sst
  mod_mvfun_bm = fit_mvfun_bm(summary_transect_complete),

  # community analysis
  mod_mvfun_com = fit_mvfun_com(summary_transect_complete),
  mod_mvfun_com2 = fit_mvfun_com2(summary_transect_complete), #get absolute slopes

  # output tables 
  tab_mod_mv_siteloc = make_table_mod_mv_siteloc(mod_mv_siteloc),
  tab_mod_mvfun_bm = make_table_mod_mvfun_bm(mod_mvfun_bm),
  tab_mod_mvfun_com = make_table_mod_mvfun_com(mod_mvfun_com),
  #tab_mod_mvfun_com2 = make_table_mod_mvfun_com(mod_mvfun_com2),
  
  ##### contribution analysis #####
  #contr_family = get_cf(contributions, herb_pisc),
  degree_dominance = get_dd_site(contributions_site),
  sp_importance = get_importance_site(contributions_site),
  freq_dominance = get_fd_site(sp_importance),
  
  mod_dd = fit_mod_dd(degree_dominance),
  mod_fd = fit_mod_fd(freq_dominance),

  ##### Species vulnerability #####
  #vulnerability = get_vuln(sptl),
  #spi_vuln = get_spi_vuln(sp_importance, vulnerability)#,
  # 
  # ##### FIGURES #####
  coords = get_coords_siteloc(summary_transect_complete)
  # # main
  # fig1 = make_fig1(location_effect),
  # fig2 = make_fig2(commodels),
  # fig3 = make_fig3(contributions, herb_pisc, degree_dominance, freq_dominance),
  # fig4 = make_fig4(contributions, vulnerability, herb_pisc, residuals),
  # 
  # # Supplemental data figures
  # SI_fig1 = make_annex_fig1(summary_transect_complete, residuals, bmmodels),
  # SI_fig2 = make_annex_fig2(summary_transect_complete, residuals),
  # 
  # SI_rank_plots = make_rank_plots(contributions, herb_pisc),
  # SI_pp_plots = make_pp_plots(bmmodels, procmodels, commodels),
  # 
  # extra_diet_analysis = alt_diet(tfish, cnpflux),
  # 
  # ##### TEXT #####
  # main_text_doc = rmarkdown::render(knitr_in("text/main.Rmd"),
  #                                   output_format = "word_document",
  #                                   output_dir = "./output/text/",
  #                                   output_file = "Schiettekatte_global_functions_main.docx"),
  # # methods_text_doc = rmarkdown::render(knitr_in("text/methods.Rmd"),
  # #                                   output_format = "word_document",
  # #                                   output_dir = "./output/text/",
  # #                                   output_file = "Schiettekatte_global_functions_methods.docx"),
  # # suppl_methods_text_doc = rmarkdown::render(knitr_in("text/suppl_methods.Rmd"),
  # #                                   output_format = "word_document",
  # #                                   output_dir = "./output/text/",
  # #                                   output_file = "Schiettekatte_global_functions_suppl_methods.docx"),
  # suppl_mat_doc = rmarkdown::render(knitr_in("text/supplementary_materials.Rmd"),
  #                                             output_format = "word_document",
  #                                             output_dir = "./output/text/",
  #                                             output_file = "Schiettekatte_global_functions_suppl_materials.docx")
  #  
)

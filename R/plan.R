
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
  sp_loc = unique(select(ungroup(contributions), bioregion, species)),
  
  # herbivory, piscivory 
  herb_pisc = get_herb_pisc(tfish, cnpflux, params, summary_transect),
  
  # combine summary data, sst 
  summary_transect_complete = combine_summary(summary_transect, herb_pisc),
  
  ##### models #####
  # run models only with biomass ans sst
  bmmodels = run_bmmodels(summary_transect_complete),
  residuals = get_residuals(summary_transect_complete, procmodels),
  
  # run models proc ~ biomass et al
  procmodels = run_procmodels(summary_transect_complete),
  location_effect = get_location_effect(procmodels, summary_transect),
  
  # run models proc ~ community vars
  commodels = run_commodels(summary_transect_complete),
  commodels_real = run_commodels_abs(summary_transect_complete),
  
  ##### contribution analysis #####
  contr_family = get_cf(contributions, herb_pisc),
  degree_dominance = get_dd(contributions, herb_pisc),
  sp_importance = get_importance(contributions, herb_pisc),
  freq_dominance = get_fd(sp_importance),
  
  ##### Species vulnerability #####
  vulnerability = get_vuln(sptl),
  spi_vuln = get_spi_vuln(sp_importance, vulnerability),
  
  ##### FIGURES #####
  # main
  fig1 = make_fig1(location_effect),
  fig2 = make_fig2(commodels),
  fig3 = make_fig3(contributions, herb_pisc, degree_dominance, freq_dominance),
  fig4 = make_fig4(contributions, vulnerability, herb_pisc, residuals),
  
  # Supplemental data figures
  SI_fig1 = make_annex_fig1(summary_transect_complete, residuals, bmmodels),
  SI_fig2 = make_annex_fig2(summary_transect_complete, residuals),

  SI_rank_plots = make_rank_plots(contributions, herb_pisc),
  SI_pp_plots = make_pp_plots(bmmodels, procmodels, commodels),
  
  extra_diet_analysis = alt_diet(tfish, cnpflux),
  
  ##### TEXT #####
  main_text_doc = rmarkdown::render(knitr_in("text/main.Rmd"),
                                    output_format = "word_document",
                                    output_dir = "./output/text/",
                                    output_file = "Schiettekatte_global_functions_main.docx"),
  methods_text_doc = rmarkdown::render(knitr_in("text/methods.Rmd"),
                                    output_format = "word_document",
                                    output_dir = "./output/text/",
                                    output_file = "Schiettekatte_global_functions_methods.docx"),
  suppl_methods_text_doc = rmarkdown::render(knitr_in("text/suppl_methods.Rmd"),
                                    output_format = "word_document",
                                    output_dir = "./output/text/",
                                    output_file = "Schiettekatte_global_functions_suppl_methods.docx"),
  suppl_figs_doc = rmarkdown::render(knitr_in("text/suppl_tables_figures.Rmd"),
                                             output_format = "word_document",
                                             output_dir = "./output/text/",
                                             output_file = "Schiettekatte_global_functions_suppl_tables_figures.docx")
  
)

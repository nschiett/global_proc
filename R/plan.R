
plan <- drake_plan(
  
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
  
  # herbivory, piscivory --> in g dry mass
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
  
  # run models to get relative contribution of each ecosystem process
  #models_contributions_Fp = fit_contribution_Fp(contributions, sp_loc),
  #models_contributions_Fn = fit_contribution_Fn(contributions, sp_loc),
  #models_contributions_Fc = fit_contribution_Fc(contributions, sp_loc),
  #models_contributions_Ic = fit_contribution_Ic(contributions, sp_loc),
  #models_contributions_Gc = fit_contribution_Gc(contributions, sp_loc),
  
  # run models to get relative contribution herbivory piscivory
  #models_contributions_I_herb = fit_contribution_I_herb(herb_pisc, sp_loc),
  #models_contributions_I_pisc = fit_contribution_I_pisc(herb_pisc, sp_loc),
  
  # extract contributions from model
  # contributions_Fn = extract_contr(models_contributions_Fn),
  # contributions_Fp = extract_contr(models_contributions_Fp),
  # contributions_Gc = extract_contr(models_contributions_Gc),
  # contributions_I_herb = extract_contr(models_contributions_I_herb),
  # contributions_I_pisc = extract_contr(models_contributions_I_pisc),

  # combine contributions
  # contributions_sp_loc = combine_contributions(
  #   sp_loc, contributions_Fn, contributions_Fp,
  #   contributions_Gc, contributions_I_herb, contributions_I_pisc, herb_pisc),

  #contributions_sp_loc_occ = add_occurence(contributions, contributions_sp_loc, herb_pisc),
  
  ##### contribution analysis #####
  contr_family = get_cf(contributions, herb_pisc),
  degree_dominance = get_dd(contributions, herb_pisc),
  sp_importance = get_importance(contributions, herb_pisc),
  freq_dominance = get_fd(sp_importance),
  
  ##### Species vulnerability #####
  vulnerability = get_vuln(sptl),
  spi_vuln = get_spi_vuln(sp_importance, vulnerability),
  
  ##### FIGURES #####
  fig1 = make_fig1(location_effect),
  fig2 = make_fig2(commodels),
  fig3 = make_fig3(contributions, herb_pisc, degree_dominance, freq_dominance),
  fig4 = make_fig4(spi_vuln),
  #plot_covu = plot_con_vuln_all(contributions_sp_loc_occ, vulnerability),
  
  annex_fig1 = make_annex_fig1(summary_transect_complete, bmmodels),
  annex_fig2 = make_annex_fig2(summary_transect_complete, residuals),
  annex_fig3 = make_annex_fig3(summary_transect_complete, residuals),
  
  ##### TEXT #####
  main_text_doc = rmarkdown::render(knitr_in("text/main.Rmd"), 
                                    output_format = "word_document", 
                                    output_dir = "./output/text/",
                                    output_file = "Schiettekatte_global_functions_main.docx")
)

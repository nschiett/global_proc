
plan <- drake_plan(
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
  
  # load all transect data
  tfish = read.csv("data/tfish_corr.csv"),
  
  # summarize results per transect
  summary_transect = summarise_pertransect(tfish, cnpflux, params),
  output_fluxsum = write.csv(summary_transect, "output/data/flux_summary_transect.csv", row.names = FALSE),
  
  contributions = get_contributions(tfish, cnpflux, params, summary_transect),
  sp_loc = unique(select(ungroup(contributions), bioregion, species)),
  
  # run models to get relative contribution of each ecosystem process
  models_contributions_Fp = fit_contribution_Fp(contributions, sp_loc),
  models_contributions_Fn = fit_contribution_Fn(contributions, sp_loc),
  models_contributions_Fc = fit_contribution_Fc(contributions, sp_loc),
  models_contributions_Ic = fit_contribution_Ic(contributions, sp_loc),
  models_contributions_In = fit_contribution_In(contributions, sp_loc),
  models_contributions_Ip = fit_contribution_Ip(contributions, sp_loc),
  models_contributions_Gc = fit_contribution_Gc(contributions, sp_loc),
  models_contributions_Gn = fit_contribution_Gn(contributions, sp_loc),
  models_contributions_Gp = fit_contribution_Gp(contributions, sp_loc),
  models_contributions_Wp = fit_contribution_Wp(contributions, sp_loc),
  models_contributions_Wn = fit_contribution_Wn(contributions, sp_loc),
  models_contributions_Wc = fit_contribution_Wc(contributions, sp_loc),
  
  # herbivory, piscivory --> in g dry mass
  herb_pisc = get_herb_pisc(tfish, cnpflux, params, summary_transect)
  
  # run models to get relative contribution herbivory piscivory
  
)


#hello

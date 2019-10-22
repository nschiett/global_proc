
plan <- drake_plan(
  sptl = get_sptl(),
  kmax = get_kmax(),
  lw = get_lw(sptl),
  cnp = get_cnp(),
  cnpdiet = get_cnpdiet(),
  metpar = get_metpar(kmax, sptl),
  artr = get_artr(),
  params = combine_params(sptl, kmax, lw, cnp, cnpdiet, metpar, artr),
  output_params = write.csv(params, "params_sst_glob.csv", row.names = FALSE),
  tfish_unique = get_tfish_unique(),
  cnpflux = run_fishflux(tfish_unique, params),
  output_cnpflux = write.csv(cnpflux, "flux_results.csv", row.names = FALSE)
)

#make(plan)

small_config <- drake_config(plan)
vis_drake_graph(small_config, targets_only = TRUE)

make(plan)


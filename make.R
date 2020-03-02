source("R/packages.R")  # loads packages
source("R/functions_wrangling.R") # defines the create_plot() function
source("R/functions_analysis.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan
# options(clustermq.scheduler = "multicore") # optional parallel computing. Also needs parallelism = "clustermq"

make(
  plan, # defined in R/plan.R
  verbose = 2, lock_envir = FALSE
)

#make(plan, targets = c("models_contributions_Wc", "models_contributions_Wn", "models_contributions_Wp"))

config <- drake_config(plan)
drake::r_vis_drake_graph(targets_only = TRUE)

drake::r_make()

# drake::clean("models_contributions_Wn")
# make(plan, targets = "herb_pisc", lock_envir = FALSE)

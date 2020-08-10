source("R/packages.R")  # loads packages
source("R/functions_wrangling.R") # defines the create_plot() function
source("R/functions_analysis.R") # defines the create_plot() function
source("R/functions_plots.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan

config <- drake_config(plan, lock_envir = FALSE)
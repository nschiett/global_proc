source("R/packages.R")  # loads packages
source("R/functions_wrangling.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan
# options(clustermq.scheduler = "multicore") # optional parallel computing. Also needs parallelism = "clustermq"

make(
  plan, # defined in R/plan.R
  verbose = 2
)

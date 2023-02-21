all_studies = c(10, 16, 24, 48)
alpha = 3; beta = 3
#alpha = 8.5; beta = 8.5
#alpha = 20; beta = 20
#rate = 0.04
rate = 0.015

all_dims <- c()
for(studies in all_studies){

  file_id <- paste0(
    "~/Desktop/Comparison Simulation/Raw Results/Balanced/",
    #"F:/School/MSc/STA4000/RareEventsMeta Simulations/Comparisons/Final Results 4 (ACTUAL FINAL RESULTS)/Balanced/",
    "CP",
    "_K", studies,
    "_A", alpha,
    "_B", beta,
    "_R", rate,
    ".rds"
  )
  my_data <- readRDS(file_id)
  my_comp <- my_data$all_comp
  my_res <- my_data$all_res

  all_dims <- c(all_dims, dim(my_comp)[2])

}

all_dims/3

# Shouldn't these all have 6000 columns?
# There should be 2000 simulations for all settings, even if sometimes methods do not work.
# Shouldn't all these have # of columns that are multiples of 3?
# These seems to be a problem for all cases of alpha, beta, and rate above and balanced vs. unbalanced.


# You also have a line like:
# Why aren't the simulations working? Please check.
# my_comp <- my_comp[, colSums(is.na(my_comp))==0]




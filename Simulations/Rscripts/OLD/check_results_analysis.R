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


# Check the Type I error calculation for one setting.
studies = 48
alpha = 3; beta = 3
rate = 0.015


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

all_comp <- my_data$all_comp
all_rest <- my_data$all_res

# Remove NA colums.
all_comp_rm_NA <- all_comp[, colSums(is.na(all_comp))==0]
dim(all_comp_rm_NA)
comp_include_H0 <- rowMeans(all_comp_rm_NA[, seq(3, ncol(all_comp_rm_NA), by = 3)])
type_1_error <- 1 - comp_include_H0
type_1_error

# This seems to match the plot here: https://drive.google.com/file/d/1gSg0oOVTGcjat4brc20rKb4z4j-KXE13/view?usp=sharing
# The current results can be found in this folder: https://drive.google.com/drive/folders/1ftgW4TjavtV8uDEwfKWllmu67XIn95R_?usp=sharing


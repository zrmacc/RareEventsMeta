# Library.
library(RareEventsMeta)
library(optparse)

# Library for comparison methods.
library(meta)

# Don't drop double zero studies from meta-analysis - use this updated
# data generation function.

source("~/Documents/GitHub/RareEventsMeta/RareEventsMeta/R/data_gen3.R")


setwd("/Users/jgrons/Documents/GitHub/RareEventsMeta/Simulations/")

studies <- 48
study_sizes <- data.table::fread(file = "Configs/study_sizes.txt")

if(studies > 48){

  n1 <- rep(study_sizes$n1, studies/48)
  n2 <- rep(study_sizes$n2,  studies/48)

}else{

  n1 <- study_sizes$n1[1:studies]
  n2 <- study_sizes$n2[1:studies]


}

get_eff_N <- function(studies, n1, n2, alpha, beta, psi){

  res <- GenData(
    total_studies = studies,
    n1 = n1,
    n2 = n2,
    alpha = alpha,
    beta = beta,
    psi = psi
  )

  return(sum((res[, "events_1"] == 0) & (res[, "events_2"] == 0)) / studies)

}


alpha <- 1.1
beta <- 1.65
psi <- alpha/0.04
n_sims <- 2000

total_res <- c()
for(i in n_sims){
  set.seed(i)
  total_res <- c(total_res, get_eff_N(studies, n1, n2, alpha, beta, psi))
}

mean(total_res)








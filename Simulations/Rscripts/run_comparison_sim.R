###########################
###### Type 1 error #######
###########################
setwd("/Users/jgrons/Documents/GitHub/RareEventsMeta/Simulations/Rscripts")

## Rare base rate, high heterogeneity, small n

my_settings <- c(10, 3, 3, 0.002)
source("comparison_sim.R")

## Rare base rate, moderate heterogeneity, small n

my_settings <- c(10, 8.5, 8.5, 0.002)
source("comparison_sim.R")

## Rare base rate, low heterogeneity, small n

my_settings <- c(10, 20, 20, 0.002)
source("comparison_sim.R")

## Moderate base rate, high heterogeneity, small n

my_settings <- c(10, 3, 3, 0.02)
source("comparison_sim.R")

## Moderate base rate, moderate heterogeneity, small n

my_settings <- c(10, 8.5, 8.5, 0.02)
source("comparison_sim.R")

## Moderate base rate, low heterogeneity, small n

my_settings <- c(10, 20, 20, 0.02)
source("comparison_sim.R")


####################
###### Power #######
####################

## Rare base rate, high heterogeneity, small n

my_settings <- c(10, 6, 3, 0.002)
source("comparison_sim.R")

## Rare base rate, moderate heterogeneity, small n

my_settings <- c(10, 17, 8.5, 0.002)
source("comparison_sim.R")

## Rare base rate, low heterogeneity, small n

my_settings <- c(10, 40, 20, 0.002)
source("comparison_sim.R")

## Moderate base rate, high heterogeneity, small n

my_settings <- c(10, 6, 3, 0.02)
source("comparison_sim.R")

## Moderate base rate, moderate heterogeneity, small n

my_settings <- c(10, 17, 8.5, 0.02)
source("comparison_sim.R")

## Moderate base rate, low heterogeneity, small n

my_settings <- c(10, 40, 20, 0.02)
source("comparison_sim.R")

# Library.
library(RareEventsMeta)
library(optparse)

# Library for comparison methods.
library(meta)

# Don't drop double zero studies from meta-analysis - use this updated
# data generation function.

source("~/Documents/GitHub/RareEventsMeta/RareEventsMeta/R/data_gen3.R")


setwd("/Users/jgrons/Documents/GitHub/RareEventsMeta/Simulations/")

# -----------------------------------------------------------------------------
# Unpack simulation settings.
# -----------------------------------------------------------------------------

# Command line options.
opt_list <- list()

# Sample size.
opt <- make_option(c("--studies"), type = "integer", help = "Studies",  default = 48)
opt_list <- c(opt_list, opt)

# Alpha.
opt <- make_option(c("--alpha"), type = "numeric", help = "Alpha", default = 1.1)
opt_list <- c(opt_list, opt)

# Beta.
opt <- make_option(c("--beta"), type = "numeric", help = "Beta", default = 1.1)
opt_list <- c(opt_list, opt)

# Psi.
opt <- make_option(c("--psi"), type = "numeric", help = "Psi", default = 1.1 / 0.04)
opt_list <- c(opt_list, opt)

# Simulation replicates.
opt <- make_option(c("--reps"), type = "integer", help = "Replicates", default = 200)
opt_list <- c(opt_list, opt)

# Iterations.
opt <- make_option(c("--mc"), type = "integer", help = "MC iterations", default = 200)
opt_list <- c(opt_list, opt)

# Output directory.
opt <- make_option(c("--out"), type = "character", help = "Output stem", default = "Results/")
opt_list <- c(opt_list, opt)

# Option parsing.
t0 <- proc.time()
parsed_opts <- OptionParser(option_list = opt_list)
params <- parse_args(object = parsed_opts)

# Output stem.
file_id <- paste0(
  "CP",
  "_K", params$studies,
  "_A", params$alpha,
  "_B", params$beta,
  "_P", params$psi,
  ".rds"
)


# -----------------------------------------------------------------------------
# Simulation parameters.
# -----------------------------------------------------------------------------
#setwd('~/Documents/GitHub/RareEventsMeta/Simulations')

# Data Generation.
studies <- params$studies
alpha <- params$alpha
beta <- params$beta
psi <- params$psi
t1e <- 0.05

study_sizes <- data.table::fread(file = "Configs/study_sizes.txt")

if(studies > 48){
  
  n1 <- rep(study_sizes$n1, studies/48)
  n2 <- rep(study_sizes$n2,  studies/48)
  
}else{
  
  n1 <- study_sizes$n1[1:studies]
  n2 <- study_sizes$n2[1:studies]
  
  
}


# Simulations.
reps <- params$reps
mc <- params$mc
num_nu_vals <- 15

# -----------------------------------------------------------------------------
# Functions.
# -----------------------------------------------------------------------------

#' Data Generating Process
#'
#' Wraps data generation.
#'
#' @return Simulated data.

DGP <- function() {
  
  # Data.
  data <- GenData(
    total_studies = studies,
    n1 = n1,
    n2 = n2,
    alpha = alpha,
    beta = beta,
    psi = psi
  )
  
  #print(warning())
  # Remove study if events exceeds study size.
  sub <- subset(
    x = data,
    (events_1 < size_1) & (events_2 < size_2)
  )
  
  removed <- nrow(data) - nrow(sub)
  if (removed > 0) {
    msg <- paste0(removed, " studies removed due to excess events.\n")
    warning(msg)
  }
  return(sub)
}
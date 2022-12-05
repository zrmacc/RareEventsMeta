# Library.
library(RareEventsMeta)
library(optparse)

# -----------------------------------------------------------------------------
# Unpack simulation settings.
# -----------------------------------------------------------------------------

# Command line options.
opt_list <- list()

# Sample size.
opt <- make_option(c("--studies"), type = "integer", help = "Studies", default = 10)
opt_list <- c(opt_list, opt)

# Alpha.
opt <- make_option(c("--alpha"), type = "numeric", help = "Alpha", default = 8.5)
opt_list <- c(opt_list, opt)

# Beta.
opt <- make_option(c("--beta"), type = "numeric", help = "Beta", default = 8.5)
opt_list <- c(opt_list, opt)

# Base rate.
opt <- make_option(c("--rate"), type = "numeric", help = "Base rate", default = 0.006)
opt_list <- c(opt_list, opt)

# Simulation replicates.
opt <- make_option(c("--reps"), type = "integer", help = "Replicates", default = 500)
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
  "_R", params$rate,
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
rate <- params$rate
t1e <- 0.05

study_sizes <- data.table::fread(file = "Configs/study_sizes.txt")
n1 <- study_sizes$n1[1:studies]
n2 <- study_sizes$n2[1:studies]

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
    alpha2 = alpha,
    beta2 = beta,
    rate1 = rate
  )
  
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

eff_n <- c()
for(i in 1:5000){
  
  data <- DGP()
  eff_n <- c(eff_n, nrow(data))
  
}

summary(eff_n)

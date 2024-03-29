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
opt <- make_option(c("--alpha"), type = "numeric", help = "Alpha", default = 5)
opt_list <- c(opt_list, opt)

# Beta.
opt <- make_option(c("--beta"), type = "numeric", help = "Beta", default = 5)
opt_list <- c(opt_list, opt)

# Base rate.
opt <- make_option(c("--rate"), type = "numeric", help = "Base rate", default = 0.0038)
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


# -----------------------------------------------------------------------------
# Alpha, beta pairs corresponding to nu search sequence.
# These do no change across simulation replicates.
ab_vals <- NuSeq(
  alpha = alpha,
  beta = beta,
  num_nu_vals = num_nu_vals
)


# -----------------------------------------------------------------------------

#' Check Coverage.
#'
#' @param data Data.frame returned by `DGP`.
#' @return Vector of p-values of length `num_nu_vals`.

CheckCoverage <- function(data) {

  aux <- function(i) {
    out <- try(
      RunMC(
        size_1 = data$size_1,
        events_1 = data$events_1,
        size_2 = data$size_2,
        events_2 = data$events_2,
        reps = mc,
        alpha = ab_vals$alpha[i],
        beta = ab_vals$beta[i],
        p_only = TRUE
      )
    )
    if (class(out) == "try-error") {
      out <- NA
    }
    return(out)
  }

  pvals <- sapply(seq_len(num_nu_vals), aux)
  return(pvals)
}

# -----------------------------------------------------------------------------

#' Simulation loop.
Sim <- function(i) {
  data <- DGP()
  pvals <- CheckCoverage(data = data)
  pvals_all <- c(nrow(data), pvals, any(pvals >= 0.05))
  return(pvals_all)
}

# # For Jesses Testing
# set.seed(92047)
# all_res <- c()
# for(i in 51:600){
#   print(i)
#   all_res <- rbind(all_res, Sim(i))
# }
#
# t1 <- proc.time()
# elapsed <- t1-t0
# cat("Time elapsed: ", elapsed["elapsed"], "sec.\n")
#
# dim(all_res)
# colMeans(all_res)

results <- lapply(seq_len(reps), Sim)
results <- do.call(rbind, results)

# -----------------------------------------------------------------------------

maxNA <- function(x) {return(max(x, na.rm = TRUE))}
out <- data.frame(
  "studies" = studies,
  "rate" = rate,
  "alpha" = alpha,
  "beta" = beta,
  "reps" = reps,
  "mc" = mc,
  "na" = sum(is.na(results)),
  "coverage" = mean(apply(results, 1, maxNA) > t1e, na.rm = TRUE)
)

out_stem <- params$out
if (!dir.exists(out_stem)) {
  dir.create(out_stem, recursive = TRUE)
}
out_file <- paste0(out_stem, file_id)
saveRDS(object = out, file = out_file)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------
t1 <- proc.time()
elapsed <- t1-t0
cat("Time elapsed: ", elapsed["elapsed"], "sec.\n")

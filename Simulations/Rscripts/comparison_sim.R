# Library.
library(RareEventsMeta)
library(optparse)

# Library for comparison methods.
library(meta)

# Don't drop double zero studies from meta-analysis.
source("~/Documents/GitHub/RareEventsMeta/RareEventsMeta/R/data_gen2.R")
setwd("/Users/jgrons/Documents/GitHub/RareEventsMeta/Simulations/")
# -----------------------------------------------------------------------------
# Unpack simulation settings.
# -----------------------------------------------------------------------------

# Command line options.
opt_list <- list()

# Sample size.
opt <- make_option(c("--studies"), type = "integer", help = "Studies", default = my_settings[1])
opt_list <- c(opt_list, opt)

# Alpha.
opt <- make_option(c("--alpha"), type = "numeric", help = "Alpha", default = my_settings[2])
opt_list <- c(opt_list, opt)

# Beta.
opt <- make_option(c("--beta"), type = "numeric", help = "Beta", default = my_settings[3])
opt_list <- c(opt_list, opt)

# Base rate.
opt <- make_option(c("--rate"), type = "numeric", help = "Base rate", default = my_settings[4])
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
  beta = alpha, # If under H0, alpha = beta.
  # If under H1, we want to check H0 value.
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
#' Comparison methods.
IncludeNull <- function(CI, null_val = log(1)){

  lower_less <- I(CI[1] <= null_val) * 1
  upper_more <- I(CI[2] >= null_val) * 1

  return(lower_less * upper_more)
}

#' Comparison methods.
CompMethods <- function(data, data_dz_removed){

  # With continuity correction.
  rr <- metabin(data[,"events_1"], data[, "size_1"],
                data[,"events_2"], data[, "size_2"])
  rr_fixed <- c(rr$lower.fixed, rr$upper.fixed)
  rr_random <- c(rr$lower.random, rr$upper.random)

  or <- metabin(data[,"events_1"], data[, "size_1"],
                data[,"events_2"], data[, "size_2"], sm = "OR")
  or_fixed <- c(or$lower.fixed, or$upper.fixed)
  or_random <- c(or$lower.random, or$upper.random)

  peto <-  metabin(data[,"events_1"], data[, "size_1"],
                   data[,"events_2"], data[, "size_2"], method = "Peto")
  peto_fixed <- c(peto$lower.fixed, peto$upper.fixed)
  peto_random <- c(peto$lower.random, peto$upper.random)

  # With double zero studies removed.
  rr <- metabin(data_dz_removed[,"events_1"], data_dz_removed[, "size_1"],
                data_dz_removed[,"events_2"], data_dz_removed[, "size_2"])
  rr_fixed_dzr <- c(rr$lower.fixed, rr$upper.fixed)
  rr_random_dzr  <- c(rr$lower.random, rr$upper.random)

  or <- metabin(data_dz_removed[,"events_1"], data_dz_removed[, "size_1"],
                data_dz_removed[,"events_2"], data_dz_removed[, "size_2"], sm = "OR")
  or_fixed_dzr  <- c(or$lower.fixed, or$upper.fixed)
  or_random_dzr  <- c(or$lower.random, or$upper.random)

  peto <-  metabin(data_dz_removed[,"events_1"], data_dz_removed[, "size_1"],
                   data_dz_removed[,"events_2"], data_dz_removed[, "size_2"], method = "Peto")
  peto_fixed_dzr  <- c(peto$lower.fixed, peto$upper.fixed)
  peto_random_dzr  <- c(peto$lower.random, peto$upper.random)

  all_CIs <- rbind(rr_fixed,
                   rr_random,
                   rr_fixed_dzr,
                   rr_random_dzr,
                   or_fixed,
                   or_random,
                   or_fixed_dzr,
                   or_random_dzr,
                   peto_fixed,
                   peto_random,
                   peto_fixed_dzr,
                   peto_random_dzr
  )

  all_CIs_e <- cbind(all_CIs,
                     sapply(1:nrow(all_CIs), function(xx)
                       IncludeNull(all_CIs[xx, ])))

  return(all_CIs_e)
}



# -----------------------------------------------------------------------------

#' Simulation loop.
Sim <- function(i) {

  data <- DGP()

  data_dz_removed <- subset(
    x = data,
    !((events_1 == 0) & (events_2) == 0)
  )

  pvals <- CheckCoverage(data = data_dz_removed)

  pvals_all <- c(nrow(data_dz_removed), pvals, any(pvals >= 0.05))

  comp <- CompMethods(data, data_dz_removed)

  return(list(pvals_all = pvals_all,
              comp = comp))
}



set.seed(92047)
all_res <- c()
all_comp <- c()
for(i in 1:300){

  print(i)
  res <- Sim(i)
  pvals <- res$pvals_all
  comps <- res$comp

  all_res <- rbind(all_res, pvals)
  all_comp <- cbind(all_comp, comps)
#
#   if(i > 1){
#     print(colMeans(all_res))
#
#     print(rowMeans(all_comp[, seq(3, ncol(all_comp), by = 3)]))
#   }
}

t1 <- proc.time()
elapsed <- t1-t0
cat("Time elapsed: ", elapsed["elapsed"], "sec.\n")

dim(all_res)
colMeans(all_res)

rowMeans(all_comp[, seq(3, ncol(all_comp), by = 3)])

# -----------------------------------------------------------------------------

out <- data.frame(
  "studies" = studies,
  "rate" = rate,
  "alpha" = alpha,
  "beta" = beta,
  "reps" = reps,
  "mc" = mc
)

out_stem <- params$out
if (!dir.exists(out_stem)) {
  dir.create(out_stem, recursive = TRUE)
}
out_file <- paste0(out_stem, file_id)
saveRDS(object = list(all_res = all_res, all_comp = all_comp), file = out_file)


setwd("/Users/jgrons/Documents/GitHub/RareEventsMeta/Simulations/Rscripts")


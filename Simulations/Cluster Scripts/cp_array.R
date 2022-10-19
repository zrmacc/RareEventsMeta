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

# Balanced
opt <- make_option(c("--bl"), type = "integer", help = "Balanced design", default = 1)
opt_list <- c(opt_list, opt)

# Output directory.
opt <- make_option(c("--out"), type = "character", help = "Output stem", default = "Results/")
opt_list <- c(opt_list, opt)

opt <- make_option(c("--job"), type = "integer", help = "Array job", default = 1)
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
  "_mc", params$mc,
  "_reps", params$reps,
  "_bl", params$bl,
  "_job", params$job,
  ".rds"
)

# -----------------------------------------------------------------------------
# Simulation parameters.
# -----------------------------------------------------------------------------

# Data Generation.
studies <- params$studies
alpha <- params$alpha
beta <- params$beta
rate <- params$rate
t1e <- 0.05

study_sizes <- data.table::fread(file = "Configs/study_sizes.txt")

# Balanced or not.
bal <- params$bl

if (bal == 1) {
  n1 <- study_sizes$n1[1:studies]
  n2 <- n1
} else {
  n1 <- study_sizes$n1[1:studies]
  n2 <- study_sizes$n2[1:studies]
}

# Simulations.
reps <- params$reps
mc <- params$mc
num_nu_vals <- 15
start_index <- (params$job - 1)*reps + 1

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
  set.seed(i + 2021)
  data <- DGP()
  pvals <- CheckCoverage(data = data)
  return(pvals)
}

end_index <- start_index + reps - 1

results <- lapply(c(start_index:end_index), Sim)
results <- do.call(rbind, results)

# -----------------------------------------------------------------------------
# Save the results to a file. 

out_stem <- params$out
if (!dir.exists(out_stem)) {
  dir.create(out_stem, recursive = TRUE)
}

out_file <- paste0(out_stem, file_id)
saveRDS(object = results, file = out_file)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------
t1 <- proc.time()
elapsed <- t1 - t0
cat("Time elapsed: ", elapsed["elapsed"], "sec.\n")

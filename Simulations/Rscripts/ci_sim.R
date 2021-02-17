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
opt <- make_option(c("--reps"), type = "integer", help = "Replicates", default = 2)
opt_list <- c(opt_list, opt)

# Iterations.
opt <- make_option(c("--mc"), type = "integer", help = "MC iterations", default = 50)
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
  "CI",
  "_K", params$studies, 
  "_A", params$alpha,
  "_B", params$beta,
  "_R", params$rate,
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

study_sizes <- data.table::fread(file = "Configs/study_sizes.txt")
n1 <- study_sizes$n1[1:studies]
n2 <- study_sizes$n2[1:studies]

# Simulations.
reps <- params$reps
mc <- params$mc
maxit <- 250

# -----------------------------------------------------------------------------
# Simulation function.
# -----------------------------------------------------------------------------

sim <- function(i) {
  
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
  data <- subset(
    x = data,
    (events_1 < size_1) & (events_2 < size_2)
  )
  
  # Confidence interval.
  lower <- try(LowerBound(
    size_1 = data$size_1,
    events_1 = data$events_1, 
    size_2 = data$size_2,
    events_2 = data$events_2,
    reps = mc,
    step_size = 0.002,
    maxit = maxit
  ))
  if (class(lower) != "try-error") {
    lower <- min(lower$mu)
  } else {
    lower <- NA
  }
  
  upper <- try(UpperBound(
    size_1 = data$size_1,
    events_1 = data$events_1, 
    size_2 = data$size_2,
    events_2 = data$events_2,
    reps = mc,
    step_size = 0.005,
    maxit = maxit
  ))
  if (class(upper) != "try-error") {
    upper <- max(upper$mu)
  } else {
    upper <- NA
  }
  
  # Output normalized p-value.
  ci <- data.frame(
    lower = lower,
    upper = upper
  )
  return(ci)
}

results <- lapply(seq_len(reps), sim)
results <- do.call(rbind, results)
results$delta <- results$upper - results$lower

out <- data.frame(
  "studies" = studies,
  "rate" = rate,
  "alpha" = alpha,
  "beta" = beta,
  "reps" = reps,
  "mc" = mc,
  "na_lower" = sum(is.na(results$lower)),
  "mean_lower" = mean(results$lower, na.rm = TRUE),
  "med_lower" = median(results$lower, na.rm = TRUE),
  "na_upper" = sum(is.na(results$upper)),
  "mean_upper" = mean(results$upper, na.rm = TRUE),
  "med_upper" = median(results$upper, na.rm = TRUE),
  "min_len" = min(results$delta, na.rm = TRUE),
  "mean_len" = mean(results$delta, na.rm = TRUE),
  "med_len" = median(results$delta, na.rm = TRUE),
  "max_len" = max(results$delta, na.rm = TRUE)
)

out_stem <- params$out
if (!dir.exists(out_stem)) {dir.create(out_stem, recursive = TRUE)}
out_file <- paste0(out_stem, file_id)
saveRDS(object = out, file = out_file)

# -----------------------------------------------------------------------------
# End
# -----------------------------------------------------------------------------
t1 <- proc.time()
cat(t1-t0, "\n")

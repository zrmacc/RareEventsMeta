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

# -----------------------------------------------------------------------------
# Simulation function.
# -----------------------------------------------------------------------------

sim <- function(i) {
  
  # Data.
  data <- RareEventsMeta::GenData(
    total_studies = studies,
    n1 = n1,
    n2 = n2,
    alpha2 = alpha,
    beta2 = beta,
    rate1 = rate
  )
  print(dim(data))
  
  # Remove study if events exceeds study size.
  data <- subset(
    x = data,
    (events_1 < size_1) & (events_2 < size_2)
  )
  print(dim(data))
  
  # Get sequence of variance values to consider.
  mu <- alpha / (alpha + beta)
  boundary_nu <- min(mu^2*(1-mu)/(1+mu), mu*(1-mu)^2/(2-mu))
  nu_vals <- seq(1e-6, boundary_nu, length.out = 15)
  null_vals <- BoundaryAB(0.5, nu_vals)[1:length(nu_vals)]
  
  out <- sapply(null_vals, function(jj) try(RunMC(
    size_1 = data$size_1,
    events_1 = data$events_1,
    size_2 = data$size_2,
    events_2 = data$events_2,
    reps = mc,
    alpha = jj,
    beta = jj
  ))[5])
  
  if (class(out) != "try-error") {
    out <- out
  } else {
    out <- NA
  }
  
  # Output normalized p-value.
  return(out)
}

results <- sapply(seq_len(reps), sim)

out <- data.frame(
  "studies" = studies,
  "rate" = rate,
  "alpha" = alpha,
  "beta" = beta,
  "reps" = reps,
  "mc" = mc,
  "na" = sum(is.na(results)),
  "coverage" = mean(apply(results, 2, max) > t1e, na.rm = TRUE)
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
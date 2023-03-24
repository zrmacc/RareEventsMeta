# Library.
library(RareEventsMeta)
library(optparse)

# Library for comparison methods.
library(meta)

# Don't drop double zero studies from meta-analysis - use this updated 
# data generation function.
source("data_gen2.R")

#setwd("/Users/jgrons/Documents/GitHub/RareEventsMeta/Simulations/")
# -----------------------------------------------------------------------------
# Unpack simulation settings.
# -----------------------------------------------------------------------------

# Command line options.
opt_list <- list()

# Sample size.
opt <- make_option(c("--studies"), type = "integer", help = "Studies", default = 48)
opt_list <- c(opt_list, opt)

# Alpha.
opt <- make_option(c("--alpha"), type = "numeric", help = "Alpha", default = 3)
opt_list <- c(opt_list, opt)

# Beta.
opt <- make_option(c("--beta"), type = "numeric", help = "Beta", default = 3)
opt_list <- c(opt_list, opt)

# Base rate.
opt <- make_option(c("--rate"), type = "numeric", help = "Base rate", default = 0.015)
opt_list <- c(opt_list, opt)

# Simulation replicates.
opt <- make_option(c("--reps"), type = "integer", help = "Replicates", default = 1000)
opt_list <- c(opt_list, opt)

# Iterations.
opt <- make_option(c("--mc"), type = "integer", help = "MC iterations", default = 500)
opt_list <- c(opt_list, opt)

# Output directory.
opt <- make_option(c("--out"), type = "character", help = "Output stem", default = "Local Results 3/")
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

study_sizes <- data.table::fread(file = "study_sizes.txt")
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

# data = DGP()
# ort = tryCatch(metabin(data[,"events_1"], data[, "size_1"],
#         data[,"events_2"], data[, "size_2"], 
#         sm = "OR", 
#         allstudies = TRUE),
#         error = function(e){
#           return(rep(NA, 2))
#         })
# tryCatch(c(ort$lower.fixed, ort$upper.fixed), 
#          error = function(e){
#            return(rep(NA, 2))
#          })

#' Comparison methods.
CompMethods <- function(data){
  
  # ------------------------------------------------ #
  # Comparison to existing fixed effects approaches. #
  # ------------------------------------------------ #
  
  # MH odds ratio with continuity correction (include DZ studies).
  or <- tryCatch(metabin(data[,"events_1"], data[, "size_1"],
                         data[,"events_2"], data[, "size_2"], 
                         sm = "OR", 
                         allstudies = TRUE),
                 error = function(e){
                   return(rep(NA, 2))
                 })
  or_MH_cc <- tryCatch(c(or$lower.fixed, or$upper.fixed),
                       error = function(e){
                         return(rep(NA, 2))
                       })
  
  # MH odds ratio without continuity correction.
  or <- tryCatch(metabin(data[,"events_1"], data[, "size_1"],
                         data[,"events_2"], data[, "size_2"], 
                         sm = "OR", 
                         MH.exact = TRUE),
                 error = function(e){
                   return(rep(NA, 2))
                 })
  or_MH <- tryCatch(c(or$lower.fixed, or$upper.fixed),
                    error = function(e){
                      return(rep(NA, 2))
                    })
  
  # Peto method for odds ratio, fixed effects.
  or <- tryCatch(metabin(data[,"events_1"], data[, "size_1"],
                         data[,"events_2"], data[, "size_2"],
                         sm= "OR",
                         method = "Peto"),
                 error = function(e){
                   return(rep(NA, 2))
                 })
  or_peto_fixed <- tryCatch(c(or$lower.fixed, or$upper.fixed),
                            error = function(e){
                              return(rep(NA, 2))
                            })
  
  # ------------------------------------------------- #
  # Comparison to existing random effects approaches. #
  # ------------------------------------------------- #
  
  # DL method for odds ratio with continuity correction.
  or <- tryCatch(metabin(data[,"events_1"], data[, "size_1"],
                         data[,"events_2"], data[, "size_2"], 
                         sm = "OR", 
                         allstudies = TRUE),
                 error = function(e){
                   return(rep(NA, 2))
                 })
  or_dl <- tryCatch(c(or$lower.random, or$upper.random),
                    error = function(e){
                      return(rep(NA, 2))
                    })
  
  # Peto method for odds ratio, random effects.
  or <- tryCatch(metabin(data[,"events_1"], data[, "size_1"],
                         data[,"events_2"], data[, "size_2"],
                         sm= "OR",
                         method = "Peto"),
                 error = function(e){
                   return(rep(NA, 2))
                 })
  or_peto_rand <- tryCatch(c(or$lower.random, or$upper.random),
                           error = function(e){
                             return(rep(NA, 2))
                           })
  
  all_CIs <- rbind(or_MH_cc,
                   or_MH,
                   or_peto_fixed,
                   or_peto_rand,
                   or_dl
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
  
  comp <- CompMethods(data)
  # tryCatch(CompMethods(data),
  #                error = function(e){
  #                  return(rep(NA, 12))
  #                })
  
  return(list(pvals_all = pvals_all,
              comp = comp))
}


all_res <- c()
all_comp <- c()
for(i in 1:reps){
  
  res <- Sim(i)
  pvals <- res$pvals_all
  comps <- res$comp
  
  all_res <- rbind(all_res, pvals)
  all_comp <- cbind(all_comp, comps)
  
}

t1 <- proc.time()
elapsed <- t1-t0
cat("Time elapsed: ", elapsed["elapsed"], "sec.\n")

#dim(all_res)
#colMeans(all_res)

#rowMeans(all_comp[, seq(3, ncol(all_comp), by = 3)])

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


#setwd("/Users/jgrons/Documents/GitHub/RareEventsMeta/Simulations/Rscripts")


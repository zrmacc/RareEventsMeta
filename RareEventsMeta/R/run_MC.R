#' Run Monte Carlo Simulation.
#' 
#' Estimates p-values assessing the null hypothesis that the supplied \deqn{\alpha} and
#' \deqn{\beta} are the generative parameters for the observed event counts.
#' 
#' @param size_1 Size of study 1.
#' @param events_1 Events in study 1.
#' @param size_2 Size of study 2.
#' @param events_2 Events in study 2.
#' @param study Optional study identifier.
#' @param reps Replications.
#' @param alpha First shape parameter for beta distribution,  
#'   lambda_1/(lambda_1+lambda_2) ~ B(alpha, beta).
#' @param beta Second shape parameter for beta distribution,  
#'   lambda_1/(lambda_1+lambda_2) ~ B(alpha, beta).
#' @importFrom stats rbeta rbinom
#' @export
#' @return Numeric vector containing:
#' \itemize{
#'   \item `alpha` and `beta` of interest.
#'   \item Mean `mu` and variance `nu` corresponding to `alpha` and `beta`.
#'   \item P-values based on the normalized and un-normalized test statistics,
#'     `p_val_norm` and `p_val_unnorm`.
#' }

RunMC <- function(
  size_1,
  events_1,
  size_2,
  events_2,
  study = NULL,
  reps, 
  alpha, 
  beta
) {

  # Observed data.
  studies <- length(size_1)
  total_events <- events_1 + events_2
  
  # Observed moment estimators.
  obs_est <- MomentEst(size_1, events_1, size_2, events_2, study)
  
  # True mean and variance.
  mu <- alpha / (alpha + beta)
  nu <- mu * (1 - mu) / (alpha + beta + 1)

  # Observed test statistics.
  t_stat_obs_unnorm <- (obs_est$mu - mu)^2
  t_stat_obs_norm <- t_stat_obs_unnorm / obs_est$mu_se2

  # Simulate test statistics.
  loop <- function(i) {
    
    # Simulate data.
    theta <- -Logit(rbeta(studies, alpha, beta))
    pi <- Expit(-theta-log(size_2 / size_1))
    sim_events <- rbinom(studies, total_events, pi)
    
    # Moment estimators and statistics.
    sim_est <- MomentEst(size_1, sim_events, size_2, total_events - sim_events)
    t_stat_sim_unnorm <- (sim_est$mu - mu)^2
    t_stat_sim_norm <- t_stat_sim_unnorm / (sim_est$mu_se2)
    
    # Output.
    out <- c("t_stat_unnorm" = t_stat_sim_unnorm, "t_stat_norm" = t_stat_sim_norm)
    return(out)
  }
  sim <- lapply(seq_len(reps), loop)
  sim <- do.call(rbind, sim)

  # Empirical p-values.
  p_val_unnorm <- mean(sim[, "t_stat_unnorm"] >= t_stat_obs_unnorm)
  p_val_norm <- mean(sim[, "t_stat_norm"] >= t_stat_obs_norm)

  # Output.
  out <- c(
    "alpha" = alpha,
    "beta" = beta,
    "mu" = mu,
    "nu" = nu,
    "p_val_norm" = p_val_norm,
    "p_val_unnorm" = p_val_unnorm
  )
  return(out)
}

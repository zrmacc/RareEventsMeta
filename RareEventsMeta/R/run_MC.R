# Purpose: Monte Carlo test of the hypothesis that a given (alpha, beta) pair
# are the generative parameters for the observed data.
# Updated: 2021-03-16

# -----------------------------------------------------------------------------

#' Chi Squared Statistic
#' 
#' @param mu_hat Estimated mu.
#' @param mu True mu.
#' @param se2 Squared standard error of mu_hat.
#' @return Numeric chi2 statistic.

CalcChi2 <- function(mu_hat, mu, se2) {
  out <- (mu_hat - mu)^2 / se2
  return(out)
}


# -----------------------------------------------------------------------------

#' Run Monte Carlo Simulation.
#' 
#' Estimates p-values assessing the null hypothesis that the supplied \eqn{\alpha} and
#' \eqn{\beta} are the generative parameters for the observed event counts.
#' 
#' @param size_1 Size of study 1.
#' @param events_1 Events in study 1.
#' @param size_2 Size of study 2.
#' @param events_2 Events in study 2.
#' @param reps Replications.
#' @param alpha First shape parameter for beta distribution,  
#'   lambda_1/(lambda_1+lambda_2) ~ B(alpha, beta).
#' @param beta Second shape parameter for beta distribution,  
#'   lambda_1/(lambda_1+lambda_2) ~ B(alpha, beta).
#' @param mu Optional mu supplied in place of alpha and beta. Note that if mu is
#'   supplied, the values of alpha and beta are overwritten.
#' @param study Optional study identifier.
#' @param p_only Return p-value only? 
#' @importFrom stats rbeta rbinom
#' @export
#' @return Numeric vector containing:
#' \itemize{
#'   \item `alpha` and `beta` of interest.
#'   \item Mean `mu` and variance `nu` corresponding to `alpha` and `beta`.
#'   \item P-value `p` assessing the null hypothesis that 
#'     \eqn{\mu = \alpha / (\alpha + \beta)}.
#' }

RunMC <- function(
  size_1,
  events_1,
  size_2,
  events_2,
  reps, 
  alpha, 
  beta,
  mu = NULL,
  study = NULL,
  p_only = FALSE
) {

  # Check for mu.
  if (!is.null(mu)) {
    ab <- as.numeric(BoundaryAB(mu))
    alpha <- ab[1]
    beta <- ab[2]
  }
  
  # Observed data.
  studies <- length(size_1)
  total_events <- events_1 + events_2
  
  # Observed moment estimators.
  obs_est <- MomentEst(size_1, events_1, size_2, events_2, study)
  
  # True mean and variance.
  mu <- alpha / (alpha + beta)
  nu <- mu * (1 - mu) / (alpha + beta + 1)

  # Observed test statistics.
  t_stat_obs <- CalcChi2(obs_est$mu, mu, obs_est$mu_se2)

  # Simulate test statistics.
  loop <- function(i) {
    
    # Simulate data.
    theta <- -Logit(rbeta(studies, alpha, beta))
    pi <- Expit(-theta-log(size_2 / size_1))
    sim_events <- rbinom(studies, total_events, pi)
    
    # Moment estimators and statistics.
    sim_est <- MomentEst(size_1, sim_events, size_2, total_events - sim_events)
    t_stat_sim <- CalcChi2(sim_est$mu, mu, sim_est$mu_se2)
    return(t_stat_sim)
  }
  sim <- lapply(seq_len(reps), loop)
  sim <- do.call(c, sim)

  # Empirical p-values.
  p_val <- mean(sim >= t_stat_obs)
  
  # Output.
  if (p_only) {
    out <- c("p" = p_val)
  } else {
    out <- data.frame(
      "alpha" = alpha,
      "beta" = beta,
      "mu" = mu,
      "nu" = nu,
      "p" = p_val
    )
  }
  return(out)
}

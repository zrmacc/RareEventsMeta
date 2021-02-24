# Updated: 2020-12-18

#' Upper Bound for Mu.
#' 
#' Finds the upper boundary for mu via grid search.
#' 
#' @param size_1 Size of study 1.
#' @param events_1 Events in study 1.
#' @param size_2 Size of study 2.
#' @param events_2 Events in study 2.
#' @param reps MC replicates for p-value estimation.
#' @param t1e Type I error level.
#' @param mu0 Optional initial value for mu.
#' @param step_size Distance between successive estimates of mu.
#' @param maxit Maximum number of iterations to perform.
#' @importFrom stats qnorm
#' @export
#' @return Data.frame containing the results of \code{\link{RunMC}}
#'   for each \eqn{\mu} in the search interval.

UpperBound <- function(
  size_1,
  events_1,
  size_2,
  events_2,
  reps = 500,
  t1e = 0.05,
  step_size = 0.0002,
  maxit = 100,
  mu0 = NULL
) {
  
  # Create a wrapper function that maps directly from a value of mu
  # to the RunMC output.
  WrapMC <- function(mu) {
    ab <- as.numeric(BoundaryAB(mu))
    out <- RunMC(
      size_1 = size_1,
      events_1 = events_1,
      size_2 = size_2,
      events_2 = events_2,
      reps = reps,
      alpha = ab[1],
      beta = ab[2]
    )
    return(out)
  }
  
  # -------------------------------------------------------
  # Initiation.
  # -------------------------------------------------------
  
  # Try initial value of mu.
  if (!is.null(mu0)) {
    current_bound <- mu0
    current_mc <- WrapMC(current_bound)
    current_p <- current_mc["p_val_norm"]
    if (current_p <= t1e) {
      stop("Search was initiated outside the boundary.\n")
    }
  } else {
    
    # Observed moment estimators.
    obs_est <- MomentEst(size_1, events_1, size_2, events_2)
    
    # Moment estimator CI.
    crit <- qnorm(p = 1 - t1e / 2)
    mom_upper <- obs_est$mu + crit * obs_est$mu_se2
    current_bound <- min(mom_upper, 1 - step_size)
    current_mc <- WrapMC(current_bound)
    current_p <- current_mc["p_val_norm"]
  }
  
  # Initialize output.
  out <- NULL
  if (current_p > t1e) {
    out <- rbind(out, current_mc)
  } else {
    current_bound <- obs_est$mu
    current_p <- 1
  }
  
  # -------------------------------------------------------
  # Extension.
  # -------------------------------------------------------
  
  # Extend upper bound.
  for(i in 1:maxit) {
    new_bound <- current_bound + step_size
    if (new_bound >= 1) {break}
    new_mc <- WrapMC(new_bound)
    new_p <- new_mc["p_val_norm"]
    out <- rbind(out, new_mc)
    if (new_p <= t1e) {break}
    current_bound <- new_bound
  }

  # Send warning if max iteration was reached.
  if (i == maxit & new_p > t1e) {
    warning("Maximum iterations performed 
            without reaching upper bound; increase step_size.\n")
  }
  
  rownames(out) <- NULL
  out <- data.frame(out)
  return(out)
}

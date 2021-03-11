# Updated: 2020-12-18

#' Lower Bound for Mu.
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
#' @param keep_history Keep search history?#'
#' @param mu_num_extra_steps Number of steps for mu to check beyond boundary.
#' @param nu_num_extra_steps Number of steps for nu to check beyond boundary.
#' @importFrom stats qnorm
#' @export
#' @return Data.frame containing the results of \code{\link{RunMC}}
#'   for each \eqn{\mu} in the search interval.
#' @examples 
#' set.seed(2013)
#' data <- GenData(
#'   total_studies = 10,
#'   n1 = 100,
#'   n2 = 100,
#'   alpha2 = 10,
#'   beta2 = 10
#' )
#' # Note: use more high `reps` and smaller `step_size` for more accurate results.
#' LowerBound(
#'   size_1 = data$size_1,
#'   events_1 = data$events_1,
#'   size_2 = data$size_2,
#'   events_2 = data$events_2,
#'   reps = 50,
#'   step_size = 0.02
#' )

LowerBound <- function(
  size_1,
  events_1,
  size_2,
  events_2,
  reps = 500,
  t1e = 0.05,
  mu0 = NULL,
  step_size = 0.0002,
  maxit = 100,
  keep_history = FALSE,
  mu_num_extra_steps = 10,
  nu_num_extra_steps = 10
) {
  
  # Create a wrapper function that maps directly from a value of mu
  # to the RunMC output.
  WrapMC <- function(mu = 0.5, a = NULL, b = NULL) {
    
    if(is.null(a) & is.null(b)){
      ab <- as.numeric(BoundaryAB(mu))
    }
    
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
    current_p <- current_mc$p
    if (current_p <= t1e) {
      stop("Search was initiated outside the boundary.\n")
    }
  } else {
    
    # Observed moment estimators.
    obs_est <- MomentEst(size_1, events_1, size_2, events_2)
    
    # Moment estimator CI.
    crit <- qnorm(p = 1 - t1e / 2)
    mom_lower <- obs_est$mu - crit * obs_est$mu_se2
    current_bound <- max(step_size, mom_lower)
    current_mc <- WrapMC(current_bound)
    current_p <- current_mc$p
  }
  
  # Initialize output.
  if (current_p > t1e) {
    out <- current_mc
  } else {
    current_bound <- obs_est$mu
    current_p <- 1
  }
  
  # -------------------------------------------------------
  # Extension.
  # -------------------------------------------------------
  
  # Extend lower bound.
  for(i in 1:maxit) {
    new_bound <- current_bound - step_size
    if (new_bound <= 0) {break}
    new_mc <- WrapMC(new_bound)
    new_p <- new_mc$p
    if (keep_history) {
      out <- rbind(out, new_mc)
    } else {
      out <- new_mc
    }
    if (new_p <= t1e) {break}
    current_bound <- new_bound
  }

  # Send warning if max iteration was reached.
  if (i == maxit & new_p > t1e) {
    warning("Maximum iterations performed without reaching lower bound.
            Consider increasing the step_size or initializing the search.\n")
  }
  
  # -------------------------------------------------------
  # Double check the values beyond the boundary.
  # ------------------------------------------------------- 
  check_values <- seq(current_bound - step_size*mu_num_extra_steps,
                           current_bound - step_size, 
                           by = step_size)
  
  check_out <- NULL
  for(i in check_values){
    ab_vals <- NuSeq(num_nu_vals = nu_num_extra_steps, mu = i)
    for(j in 1:mu_num_extra_steps){
      temp_result <- WrapMC(mu = NULL, a = ab_vals[j, 1], b = ab_vals[j, 2])
      check_out <- rbind(check_out, temp_result)
    }
  }
 
  check_values_in_CI <- check_out[check_out[,'p'] > t1e, ]
  
  rownames(out) <- NULL
  out <- data.frame(out)
  
  rownames(check_values_in_CI) <- NULL
  check_values_in_CI <- data.frame(check_values_in_CI)
  
  return(list(out = out, beyond_boundary = check_values_in_CI))
}

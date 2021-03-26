# Updated: 2021-03-26

#' Outputs the confidence interval.
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
#' @param mu_extra_steps Number of steps for mu to check beyond boundary.
#' @param nu_extra_steps Number of steps for nu to check beyond boundary.
#' @importFrom stats qnorm
#' @export
#' @return Data.frame containing the confidence interval `lower`
#'   and `upper` bound.
#' 
#' @examples 
#' set.seed(2013)
#' data <- GenData(
#'   total_studies = 10,
#'   n1 = 100,
#'   n2 = 100,
#'   alpha2 = 10,
#'   beta2 = 10
#' )
#' # Note: use more `reps` and smaller `step_size` for more accurate results.
#' ExactConfInt(
#'   size_1 = data$size_1,
#'   events_1 = data$events_1,
#'   size_2 = data$size_2,
#'   events_2 = data$events_2,
#'   reps = 50,
#'   step_size = 0.02
#' )

ExactConfInt <- function(
  size_1,
  events_1,
  size_2,
  events_2,
  reps = 500,
  t1e = 0.05,
  mu0 = NULL,
  step_size = 0.0002,
  maxit = 100,
  mu_extra_steps = 0,
  nu_extra_steps = 0
){

  # Confidence interval lower bound.
  lower <- try(
    FindBound(
      size_1 = size_1,
      events_1 = events_1, 
      size_2 = size_2,
      events_2 = events_2,
      lower_bound = TRUE,
      reps = reps,
      t1e = t1e,
      mu0 = mu0,
      step_size = step_size,
      maxit = maxit,
      keep_history = FALSE,
      mu_extra_steps = mu_extra_steps,
      nu_extra_steps = nu_extra_steps
    )
  )
  if (class(lower) != "try-error") {
    
    # Lower bound candidates.
    lower_candidates <- c(
      lower$search_results$bound, 
      lower$boundary_region_check$mu[lower$boundary_region_check$p >= t1e]  
    )
    lower <- min(lower_candidates)
    
  } else {
    lower <- NA
  }
  
  # Confidence interval upper bound.
  upper <- try(
    FindBound(
      size_1 = size_1,
      events_1 = events_1, 
      size_2 = size_2,
      events_2 = events_2,
      lower_bound = FALSE,
      reps = reps,
      t1e = t1e,
      mu0 = mu0,
      step_size = step_size,
      maxit = maxit,
      keep_history = FALSE,
      mu_extra_steps = mu_extra_steps,
      nu_extra_steps = nu_extra_steps
    )
  )
  if (class(upper) != "try-error") {
    
    # Upper bound candidates.
    upper_candidates <- c(
      upper$search_results$bound,
      upper$boundary_region_check$mu[upper$boundary_region_check$p >= t1e]  
    )
    upper <- max(upper_candidates)    

  } else {
    upper <- NA
  }
  
  # Output.
  ci <- data.frame(
    lower = lower,
    upper = upper
  )
  
  return(ci)
}
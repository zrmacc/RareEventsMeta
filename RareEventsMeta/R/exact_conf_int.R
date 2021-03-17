# TODO: How do we handle the results of the boundary region check? 
# Updated: 2021-03-01

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
  maxit = 100
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
      keep_history = FALSE
    )
  )
  if (class(lower) != "try-error") {
    lower <- lower$search_results$bound
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
      keep_history = FALSE
    )
  )
  if (class(upper) != "try-error") {
    upper <-upper$search_results$bound
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
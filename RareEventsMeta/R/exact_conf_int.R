# Updated: 2020-12-18

#' Lower Bound for Mu.
#' 
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
#' @return Data.frame containing the results of \code{\link{RunMC}}
#'   for each \eqn{\mu} in the search interval.

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

  # Confidence interval.
  lower <- try(LowerBound(
    size_1 = data$size_1,
    events_1 = data$events_1, 
    size_2 = data$size_2,
    events_2 = data$events_2,
    reps = mc,
    step_size = step_size,
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
    step_size = step_size,
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
# Updated: 2021-03-26

# Internal: find one confidence bound (lower or upper) via FindBound; returns scalar or NA.
find_one_bound <- function(
  size_1,
  events_1,
  size_2,
  events_2,
  lower_bound,
  reps,
  t1e,
  mu0,
  step_size,
  maxit,
  mu_extra_steps,
  nu_extra_steps
) {
  result <- try(
    FindBound(
      size_1 = size_1,
      events_1 = events_1,
      size_2 = size_2,
      events_2 = events_2,
      lower_bound = lower_bound,
      reps = reps,
      t1e = t1e,
      mu0 = mu0,
      step_size = step_size,
      maxit = maxit,
      keep_history = FALSE,
      mu_extra_steps = mu_extra_steps,
      nu_extra_steps = nu_extra_steps
    ),
    silent = TRUE
  )
  if (inherits(result, "try-error")) {
    return(NA_real_)
  }
  candidates <- c(
    result$search_results$bound,
    result$boundary_region_check$mu[result$boundary_region_check$p >= t1e]
  )
  if (lower_bound) {
    return(min(candidates))
  } else {
    return(max(candidates))
  }
}

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
#' # Use fixed data with non-zero events so the example runs reliably.
#' data <- data.frame(
#'   size_1 = c(100, 100),
#'   events_1 = c(10, 12),
#'   size_2 = c(100, 100),
#'   events_2 = c(6, 7)
#' )
#' # Note: use more `reps` and smaller `step_size` for accurate results.
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
) {
  lower <- find_one_bound(
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
    mu_extra_steps = mu_extra_steps,
    nu_extra_steps = nu_extra_steps
  )
  upper <- find_one_bound(
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
    mu_extra_steps = mu_extra_steps,
    nu_extra_steps = nu_extra_steps
  )
  return(data.frame(lower = lower, upper = upper))
}
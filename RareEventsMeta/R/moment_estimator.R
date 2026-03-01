# Purpose: Moment estimator.

# -----------------------------------------------------------------------------
# Data preparation.
# -----------------------------------------------------------------------------

# Internal: one study's (events, total_events, weight) for PrepData. Not exported.
prep_one_study <- function(n_draw, n_other, e_draw, e_other) {
  if (e_other == 0) {
    return(list(events = 0, total_events = e_draw, weight = 1))
  }
  weight <- dhyper(0:e_other, e_other, n_other - e_other, n_draw)
  events <- 0:e_other
  total_events <- e_draw + events
  return(list(events = events, total_events = total_events, weight = weight))
}

#' Prepare Data for Simulation
#'
#' @param size_1 Size of study 1.
#' @param events_1 Events in study 1.
#' @param size_2 Size of study 2.
#' @param events_2 Events in study 2.
#' @param study Optional study identifier.
#'
#' @importFrom stats dhyper
#' @return Data.frame containing
#' \itemize{
#'   \item `study` identifier.
#'   \item `events` down-sampled number of events.
#'   \item `total_events` total number of events.
#'   \item `weight` based on the hypergeometric distribution.
#'   \item `n1` and `n2`, original study sizes.
#' }
PrepData <- function(
  size_1,
  events_1,
  size_2,
  events_2,
  study
) {

  loop <- function(i) {
    if (size_1[i] >= size_2[i]) {
      # More patients in group 1: draw size_2 from (size_1, size_2); count group-1 events.
      x <- prep_one_study(size_2[i], size_1[i], events_2[i], events_1[i])
    } else {
      # More patients in group 2: draw size_1; total_events = events_1 + (0:events_2).
      x <- prep_one_study(size_1[i], size_2[i], events_1[i], events_2[i])
      x$events <- rep(events_1[i], length(x$total_events))
      x$total_events <- events_1[i] + (0:events_2[i])
    }

    out <- data.frame(
      study = study[i],
      events = x$events,
      total_events = x$total_events,
      weight = x$weight,
      n1 = size_1[i],
      n2 = size_2[i],
      stringsAsFactors = FALSE
    )

    is_zero <- (out$total_events == 0)
    if (any(is_zero)) {
      out <- out[!is_zero, ]
      out$weight <- out$weight / sum(out$weight)
    }
    out <- out[out$weight > 0, ]
    return(out)
  }

  out <- lapply(seq_len(length(size_1)), loop)
  out <- do.call(rbind, out)
  return(out)
}


# -----------------------------------------------------------------------------
# Moment estimation.
# -----------------------------------------------------------------------------

#' Moment Estimation
#'
#' @param size_1 Size of study 1.
#' @param events_1 Events in study 1.
#' @param size_2 Size of study 2.
#' @param events_2 Events in study 2.
#' @param study Optional study identifier.
#'
#' @export
#' @return List containing:
#' \itemize{
#'   \item `mu`, estimated first moment.
#'   \item `mu_se2`, squared standard error of estimated first moment.
#'   \item `mu_cc`, estimated first moment with continuity correction.
#'   \item `mu2`, estimated second moment.
#'   \item `nu`, estimate of nu.
#' }

MomentEst <- function(
  size_1,
  events_1,
  size_2,
  events_2,
  study = NULL
) {

  # Create study identifier if not provided.
  if (is.null(study)) {
    study <- seq_len(length(size_1))
  }
  studies <- length(study)

  # Prepare data.
  data <- PrepData(
    size_1 = size_1,
    events_1 = events_1,
    size_2 = size_2,
    events_2 = events_2,
    study = study
  )

  # Raw first moment.
  mu <- sum(data$weight * data$events / data$total_events) / studies

  # Continuity correction.
  data$total_events <- data$total_events + 1
  # 1/(data$n1 + data$n2)
  data$events <- data$events + 0.5
  #1/data$n1

  # Continuity corrected first moment.
  mu_cc <- sum(data$weight * data$events / data$total_events) / studies

  # Estimate nu using first and second moments.
  num <- (sum(data$weight * (data$events / data$total_events)^2) / studies -
    sum(data$weight * mu_cc / data$total_events) / studies)
  denom <- sum(data$weight * (1 - 1 / (data$total_events))) / studies
  mu2 <- num / denom
  nu <- max(0, mu2 - mu_cc^2)

  # Standard error of first moment estimator.
  se2 <- sum(
    data$weight * (mu_cc * (1 - mu_cc) / data$total_events) +
      data$weight * (1 - 1 / data$total_events) * nu
  ) / studies^2

  # Output.
  out <- list(
    "mu" = mu,
    "mu_se2" = se2,
    "mu_cc" = mu_cc,
    "mu2" = mu2,
    "nu" = nu
  )
  return(out)
}

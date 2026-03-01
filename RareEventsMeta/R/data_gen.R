# -----------------------------------------------------------------------------
# Data Generation
# -----------------------------------------------------------------------------

#' Generate Data
#'
#' Simulates meta-analysis data with gamma base rates and Poisson event counts.
#'
#' @param total_studies Total number of studies.
#' @param n1 Sample size in group 1.
#' @param n2 Sample size in group 2.
#' @param alpha Shape parameter for group 1 (gamma).
#' @param beta Shape parameter for group 2 (gamma).
#' @param psi Rate parameter for gamma distributions.
#'
#' @importFrom stats rgamma rpois
#' @export
#' @return Data.frame with columns `study`, `size_1`, `events_1`, `size_2`, `events_2`.
GenData <- function(
  total_studies,
  n1,
  n2,
  alpha,
  beta,
  psi
) {


  base_rate_2 <- rgamma(total_studies, beta, psi)
  y2 <- rpois(total_studies, n2 * base_rate_2)

  base_rate_1 <- rgamma(total_studies, alpha, psi)
  y1 <- rpois(total_studies, n1 * base_rate_1)

  # Data for analysis.
  data <- data.frame(
    "study" = seq_len(total_studies),
    "size_1" = n1,
    "events_1" = y1,
    "size_2" = n2,
    "events_2" = y2
  )

  return(data)
}


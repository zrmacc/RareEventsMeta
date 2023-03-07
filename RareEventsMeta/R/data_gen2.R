# -----------------------------------------------------------------------------
# Data Generation
# -----------------------------------------------------------------------------

#' Generate Data
#'
#' @param total_studies Total number of studies.
#' @param n1 Sample size in group 1.
#' @param n2 Sample size in group 2.
#' @param alpha1 Shape parameters for group 1.
#' @param rate1 Baseline event rate in group 1.
#' @param alpha2 First shape parameter for beta distribution,
#'   lambda_1/(lambda_1+lambda_2) ~ B(alpha2, beta2).
#' @param beta2 Second shape parameter for beta distribution,
#'   lambda_1/(lambda_1+lambda_2) ~ B(alpha2, beta2).
#'
#' @importFrom stats rbeta rgamma rpois
#' @export

GenData <- function(
  total_studies,
  n1,
  n2,
  alpha1 = 1.44,
  rate1 = 0.01,
  alpha2,
  beta2
) {

  # # Baseline event rate in group 1.
  # base_rate <- rgamma(total_studies, alpha2, alpha2 / rate1)
  # 
  # # Events in group 1.
  # y2 <- rpois(total_studies, n2 * base_rate)
  # 
  # # Events in group 2.
  # #rr <- rbeta(total_studies, alpha2, beta2)
  # base_rate_2 <- rgamma(total_studies, beta2, beta2 / rate1)
  # y1 <- rpois(total_studies, n1 * base_rate_2)
  #               #base_rate * rr / (1 - rr))
  
  # Baseline event rate in group 1.
  base_rate <- rgamma(total_studies, alpha2, alpha2 / rate1)
  
  # Events in group 1.
  y2 <- rpois(total_studies, n2 * base_rate)
  
  # Events in group 2.
  rr <- rbeta(total_studies, alpha2, beta2)
  y1 <- rpois(total_studies, n1 * base_rate * rr / (1 - rr))
  
  
  # Data for analysis.
  data <- data.frame(
    "study" = seq_len(total_studies),
    "size_1" = n1,
    "events_1" = y1,
    "size_2" = n2,
    "events_2" = y2
  )

  # Drop double zero studies.
  # events_1 <- events_2 <- NULL
  # data <- subset(
  #   x = data,
  #   !((events_1 == 0) & (events_2) == 0)
  # )
  #
  # Output.
  return(data)
}


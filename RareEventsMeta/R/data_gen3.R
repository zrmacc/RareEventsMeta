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
  alpha,
  beta,
  psi
) {


  base_rate_2 <- rgamma(total_studies, beta, psi)
  y2 <- rpois(total_studies, n2 * base_rate_2)

  #if(alpha == beta){

  #  y1 <- rpois(total_studies, n1 * base_rate_2)

 # }else{

    base_rate_1 <- rgamma(total_studies, alpha, psi)
    y1 <- rpois(total_studies, n1 * base_rate_1)

 # }

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


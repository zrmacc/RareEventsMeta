#' Find Boundary Alpha Beta
#' 
#' @param mu Mean of heterogeneity parameter pi.
#' @param mu Mean of heterogeneity parameter pi.
#' @param tol Small tolerance to ensure alpha, beta > 1.
#' @export
#' @return Numeric vector containing alpha and beta.
ComputeAB <- function(mu, v = NULL, tol = 1e-6) {

  # Find nu along the boundary for given mu.
  if(is.null(v)){
  v1 <- mu^2 * (1 - mu) / (1 + mu)
  v2 <- mu * (1 - mu)^2 / (2 - mu)
  v <- min(v1, v2)
  }

  # Corresponding alpha and beta + some tolerance such that alpha, beta > 1.
  alpha <- mu * ((mu * (1 - mu) - v) / v) * (1 + tol)
  beta <- (1 - mu) * ((mu * (1 - mu) - v) / v) * (1 + tol)

  # Output.
  out <- c("alpha" = alpha, "beta" = beta)
  return(out)
}

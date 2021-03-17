#' Nu Search Sequence
#' 
#' Obtain a sequence of alpha, beta pairs with the same mu, but varying nu.
#' 
#' @param alpha Generative alpha.
#' @param beta Generative beta. 
#' @param num_nu_vals Number of nu values.
#' @param mu Value of mu if alpha, beta not entered.
#' @param tol Boundary tolerance level.
#' @export
#' @return Data.frame containing alpha, beta pairs for nu search sequence.

NuSeq <- function(
  alpha, 
  beta,
  mu = NULL,
  num_nu_vals = 10,
  tol = 1e3
) {
   
  if(is.null(mu)){
    mu <- alpha / (alpha + beta)
  }
  
  boundary_nu <- min(mu^2 * (1 - mu) / (1 + mu), mu * (1 - mu)^2 / (2 - mu))
  
  nu_vals <- seq(
    from = boundary_nu / tol, 
    to = boundary_nu,
    length.out = num_nu_vals
  )
  
  ab_vals <- lapply(nu_vals, function(nu) {return(BoundaryAB(mu, nu))})
  ab_vals <- data.frame(do.call(rbind, ab_vals))
  
  return(ab_vals)
}

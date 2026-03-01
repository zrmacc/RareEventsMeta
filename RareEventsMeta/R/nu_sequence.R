#' Nu Search Sequence
#' 
#' Obtain a sequence of alpha, beta pairs with the same mu, but varying nu.
#' 
#' @param alpha Generative alpha.
#' @param beta Generative beta. 
#' @param num_nu_vals Number of nu values.
#' @param mu Value of mu if alpha and beta are not supplied.
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
   
  if (is.null(mu)) {
    mu <- alpha / (alpha + beta)
  }

  nu_vals <- seq(
    from = boundary_nu(mu) / tol,
    to = boundary_nu(mu),
    length.out = num_nu_vals
  )

  ab_vals <- lapply(nu_vals, function(nu) {
    return(BoundaryAB(mu, nu))
  })
  ab_vals <- data.frame(do.call(rbind, ab_vals))
  
  return(ab_vals)
}

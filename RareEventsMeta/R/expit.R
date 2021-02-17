#' Logit
#'
#' Calculates \eqn{log(x) - log(1-x)}.
#'
#' @param x Numeric vector.
#' @return Numeric vector.
Logit <- function(x) {
  return(log(x / (1 - x)))
}

#' Expit
#'
#' Calculates \eqn{1 / (1 + exp(-x))}.
#'
#' @param x Numeric vector.
#' @return Numeric vector.
Expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

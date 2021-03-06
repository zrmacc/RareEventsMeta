% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nu_sequence.R
\name{NuSeq}
\alias{NuSeq}
\title{Nu Search Sequence}
\usage{
NuSeq(alpha, beta, mu = NULL, num_nu_vals = 10, tol = 1000)
}
\arguments{
\item{alpha}{Generative alpha.}

\item{beta}{Generative beta.}

\item{mu}{Value of mu if alpha, beta not entered.}

\item{num_nu_vals}{Number of nu values.}

\item{tol}{Boundary tolerance level.}
}
\value{
Data.frame containing alpha, beta pairs for nu search sequence.
}
\description{
Obtain a sequence of alpha, beta pairs with the same mu, but varying nu.
}

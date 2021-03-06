% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cholesky_banded.R
\name{chol_banded}
\alias{chol_banded}
\title{Banded Cholesky}
\usage{
chol_banded(A, k, checks_off = FALSE)
}
\arguments{
\item{A}{The input matrix}

\item{k}{bandwidth parameter}

\item{checks_off}{Boolean variable to turn on / off checks}
}
\description{
\code{chol_banded} computes the Cholesky decomposition of a banded
positive definite matrix.
}

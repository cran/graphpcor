% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/param_prec.R
\name{Lprec}
\alias{Lprec}
\alias{fillLprec}
\alias{theta2Lprec2C}
\title{Precision matrix parametrization helper functions.}
\usage{
Lprec(theta, p, ilowerL)

fillLprec(L, lfi)

theta2Lprec2C(theta, p, ilowerL)
}
\arguments{
\item{theta}{numeric vector of length \code{m}.}

\item{p}{numeric giving the dimention of Q.
If missing, \code{p = (1+sqrt(1+8*length(theta)))}
and Q is assumed to be dense.}

\item{ilowerL}{numeric vector as index
to (lower) L to be filled with \code{theta}.
Default is missing and Q is assumed to be dense.}

\item{L}{matrix as the lower triangle
containing the Cholesky decomposition of
a precision matrix}

\item{lfi}{indicator of fill-in elements}
}
\value{
matrix as the Cholesky factor of a
precision matrix as the inverse of a correlation

a matrix whose elements at the lower
triangle are the filled in elements of the
Cholesky decomposition of a precision matrix
and diagonal elements as \code{1:p}.
}
\description{
Precision matrix parametrization helper functions.
}
\details{
The precision matrix definition consider
\code{m} parameters for the lower part of L.
If Q is dense, then \code{m = p(p-1)/2}, else
\code{m = length(ilowerL)}.
Then the precision is defined as
\eqn{Q(\theta) = L(\theta)L(\theta)^T}
}
\section{Functions}{
\itemize{
\item \code{fillLprec()}: Function to fill-in a Cholesky matrix

\item \code{theta2Lprec2C()}: Internal function to build C

}}
\examples{
theta1 <- c(1, -1, -2)
Lprec(theta1)
theta2 <- c(0.5, -0.5, -1, -1)
Lprec(theta2, 4, c(2,4,7,12))
}

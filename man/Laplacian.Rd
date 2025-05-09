% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa_methods.R, R/graphpcor.R
\name{Laplacian}
\alias{Laplacian}
\alias{Laplacian.default}
\alias{Laplacian.matrix}
\title{The Laplacian of a graph}
\usage{
Laplacian(graph)

\method{Laplacian}{default}(graph)

\method{Laplacian}{matrix}(graph)
}
\arguments{
\item{graph}{an object that inherits a matrix class}
}
\value{
matrix as the Laplacian of a graph
}
\description{
The (symmetric) Laplacian of a graph is a
square matrix with dimention
equal the number of nodes.
It is defined as
\deqn{L_{ij} = n_i \textrm{ if } i=j, -1 \textrm{ if } i\sim j, 0 \textrm{ otherwise}}{%
      Lij = ni if i=j, -1 if i~j or 0 otherwise}
where i~j means that there is an edge
between nodes i and j and
n_i is the number of edges including node i.
}
\section{Methods (by class)}{
\itemize{
\item \code{Laplacian(default)}: The Laplacian default method (none)

\item \code{Laplacian(matrix)}: The Laplacian of a matrix

}}

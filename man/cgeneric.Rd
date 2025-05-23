% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa_classes.R, R/cgeneric.R,
%   R/cgeneric_utils.R, R/graphpcor.R, R/kronecker.R, R/treepcor.R
\docType{class}
\name{inla.cgeneric-class}
\alias{inla.cgeneric-class}
\alias{cgeneric}
\alias{cgeneric.default}
\alias{cgeneric.character}
\alias{cgeneric_get}
\alias{initial}
\alias{initial.inla.cgeneric}
\alias{mu}
\alias{mu.inla.cgeneric}
\alias{prior}
\alias{prior.inla.cgeneric}
\alias{graph}
\alias{graph.inla.cgeneric}
\alias{Q}
\alias{prec.inla.cgeneric}
\alias{cgeneric.graphpcor}
\alias{kronecker,inla.cgeneric,inla.cgeneric-method}
\alias{cgeneric.treepcor}
\title{\code{inla.cgeneric} class, short \code{cgeneric},
to define a \code{\link[INLA:cgeneric]{INLA::cgeneric()}} latent model}
\usage{
cgeneric(model, ...)

\method{cgeneric}{default}(model, debug = FALSE, useINLAprecomp = TRUE, libpath = NULL, ...)

\method{cgeneric}{character}(model, ...)

cgeneric_get(
  model,
  cmd = c("graph", "Q", "initial", "mu", "log_prior"),
  theta,
  optimize = TRUE
)

initial(model)

\method{initial}{inla.cgeneric}(model)

mu(model, theta)

\method{mu}{inla.cgeneric}(model, theta)

prior(model, theta)

\method{prior}{inla.cgeneric}(model, theta)

graph(model, ...)

\method{graph}{inla.cgeneric}(model, ...)

Q(model, ...)

\method{prec}{inla.cgeneric}(model, ...)

\method{cgeneric}{graphpcor}(...)

\S4method{kronecker}{inla.cgeneric,inla.cgeneric}(X, Y, FUN = "*", make.dimnames = FALSE, ...)

\method{cgeneric}{treepcor}(...)
}
\arguments{
\item{model}{an object \code{inla.cgeneric} object.}

\item{...}{arguments passed on.}

\item{debug}{integer, default is zero, indicating the verbose level.
Will be used as logical by INLA.}

\item{useINLAprecomp}{logical, default is TRUE, indicating if it is to
be used the shared object pre-compiled by INLA.
This is not considered if 'libpath' is provided.}

\item{libpath}{string, default is NULL, with the path to the shared object.}

\item{cmd}{an string to specify which model element to get}

\item{theta}{numeric vector with the model parameters.
If missing, the \code{\link[=initial]{initial()}} will be used.}

\item{optimize}{logical. If missing or FALSE,
the graph and precision are as a sparse matrix.
If TRUE, graph only return the row/col indexes and
precision return only the elements as a vector.}

\item{X}{\code{inla.cgeneric} or \code{inla.rgeneric}}

\item{Y}{\code{inla.cgeneric} or \code{inla.rgeneric}}

\item{FUN}{see \code{kronecker}}

\item{make.dimnames}{see \code{kronecker}}
}
\value{
a \code{inla.cgeneric}, \code{\link[=cgeneric]{cgeneric()}} object.

depends on \code{cmd}
}
\description{
This organize data needed on the C interface
for building latent models, which are characterized
from a given model parameters \eqn{\theta} and the
the following model elements.
\itemize{
\item \code{graph} to define the non-zero precision matrix pattern.
only the upper triangle including the diagonal is needed.
The order should be by line.
\item \code{Q} vector where the
\itemize{
\item first element (N) is the size of the matrix,
\item second element (M) is the number of non-zero
elements in the upper part (including) diagonal
\item the remaining (M) elements are the actual
precision (upper triangle plus diagonal) elements
whose order shall follow the graph definition.
}
\item \code{mu} the mean vector,
\item \code{initial} vector with
\itemize{
\item first element as the number of the parameters in the model
\item remaining elements should be the initials for the model parameters.
}
\item \code{log.norm.const} log of the normalizing constant.
\item \code{log.prior} log of the prior for the model parameters.
}

See details in \code{\link[INLA:cgeneric]{INLA::cgeneric()}}
}
\section{Methods (by class)}{
\itemize{
\item \code{cgeneric(default)}: This calls \code{\link[INLA:cgeneric]{INLA::inla.cgeneric.define()}}

\item \code{cgeneric(character)}: Method for when \code{model} is a character.
E.g. cgeneric(model = "generic0")
calls \link{cgeneric_generic0}

\item \code{cgeneric(graphpcor)}: The \code{cgeneric} method for \code{graphpcor} uses \code{\link[=cgeneric_graphpcor]{cgeneric_graphpcor()}}

\item \code{cgeneric(treepcor)}: The \code{cgeneric} method for \code{treepcor}, uses \code{\link[=cgeneric_treepcor]{cgeneric_treepcor()}}

}}
\section{Functions}{
\itemize{
\item \code{cgeneric_get()}: \code{cgeneric_get} is an internal function used by
\code{graph}, \code{prec}, \code{initial}, \code{mu} or \code{prior}
methods for \code{inla.cgeneric}

\item \code{initial()}: Retrieve the initial model parameter(s).

\item \code{initial(inla.cgeneric)}: Retrive the initial parameter(s) of an \code{inla.cgeneric} model.

\item \code{mu()}: Evaluate the mean.

\item \code{mu(inla.cgeneric)}: Evaluate the mean for an \code{inla.cgeneric} model.

\item \code{prior()}: Evaluate the log-prior.

\item \code{prior(inla.cgeneric)}: Evaluate the prior for an \code{inla.cgeneric} model

\item \code{graph()}: Retrieve the graph

\item \code{graph(inla.cgeneric)}: Retrieve the graph of an \code{inla.cgeneric} object

\item \code{Q()}: Evaluate \code{\link[=prec]{prec()}} on a model

\item \code{prec(inla.cgeneric)}: Evaluate \code{\link[=prec]{prec()}} on an \code{inla.cgeneric} object

\item \code{kronecker(X = inla.cgeneric, Y = inla.cgeneric)}: Kronecker (product) between two \code{inla.cgeneric} models as
a method for kronecker()

}}
\seealso{
\code{\link[INLA:cgeneric]{INLA::cgeneric()}}
}

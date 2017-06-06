\name{marg_param}
\alias{marg_param}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute marginal parametrization}
\description{Function that creates matrices C and M for the marginal parametrization of the probability vector for a vector of categorical variables (internal use)}
\usage{
marg_param(lev, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lev}{vector containing the number of levels of each variable}
  \item{type}{vector with elements "l", "g", "c", "r" indicating the type of logit}
}
\value{
  \item{C}{matrix of constrats (the first sum(lev)-length(r) elements are referred to univariate logits) }
  \item{M}{marginalization matrix with elements 0 and 1}
  \item{G}{corresponding design matrix for the corresponding log-linear model}
 }

\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), http://www.stat.unipg.it/bartolucci}


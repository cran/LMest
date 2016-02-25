\name{lk_sta}

\alias{lk_sta}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute the stationary log-likelihood}

\description{Function that computes the stationary log-likelihood  
(internal use).}

\usage{
lk_sta(tau, u, V, G2, outl = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau}{vector of parameters}
  \item{u}{vector of initial posterior probabilities}
  \item{V}{sum of the posterior probabilities}
  \item{G2}{design matrix}
  \item{outl}{to return additional output (la,PI)}
}

\value{
  \item{flk}{log-likelihood}
 
}

\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT) }

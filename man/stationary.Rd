\name{stationary}
\alias{stationary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stationary}

\description{Function that computes the derivatives of the log-likelihood of a stationary model (internal function).}
\usage{stationary(tau, k, G2, IPI)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau}{regression parameters}
  \item{k}{number of latent states}
  \item{G2}{internal argument}
  \item{IPI}{index of the transition matrix elements}

}
\value{
  \item{d0}{first derivatives of the log-likelihood}
  \item{d1}{second derivatives of the log-likelihood}
  }


\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)}

\name{lk_obs}

\alias{lk_obs}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute the observable log-likelihood of the basic LM model}

\description{Function that computes the observable log-likelihood of the basic LM model (internal use).}

\usage{
lk_obs(th, m, Bm, Cm, bv, k, S, R, yv, TT, r, mod)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{th}{vector of parameters}
  \item{m}{list of design matrices for the logits}
  \item{Bm}{design matrix for the logits}
  \item{Cm}{design matrix for the logits}
  \item{bv}{number of response categories}
  \item{k}{number of states}
  \item{S}{matrix of distinct response configurations}
  \item{R}{matrix of missing response configurations}
  \item{yv}{corresponding vector of frequencies}
  \item{TT}{number of times occasions}
  \item{r}{number of response variabels}
  \item{mod}{type of model}
}

\value{
  \item{lk}{log-likelihood}
  \item{sc}{score vector}
}

\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), http://www.stat.unipg.it/bartolucci}

\name{lk_obs_mixed}

\alias{lk_obs_mixed}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute the observable log-likelihood of the mixed LM model}

\description{Function that computes the observable log-likelihood of the mixed LM model with discrete random effect in the latent model (internal use).}

\usage{
lk_obs_mixed(th, nla, nPiv, nPi, nPsi, S, yv, r, k1, k2)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{th}{vector of parameters}
  \item{nla}{number of logits for the mass probability vector}
  \item{nPiv}{number of logits for the intial probabilities}
  \item{nPi}{number of logits for the transition probabilities}
  \item{nPsi}{number of logits for the conditional response probabilities}
  \item{S}{matrix of distinct response configurations}
  \item{yv}{corresponding vector of frequencies}
  \item{r}{number of response variables}
  \item{k1}{number of latent classes}
  \item{k2}{number of latent states}
}

\value{
  \item{lk}{log-likelihood}
  \item{sc}{score vector}
}

\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)}

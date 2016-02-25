\name{lk_obs_latent}

\alias{lk_obs_latent}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute the observable log-likelihood of the LM model with covariates in the latent model}

\description{Function that computes the observable log-likelihood of the LM model with covariates in the latent model (internal use).}

\usage{
lk_obs_latent(th, S, R, b, yv, Am, XXdis, Xlab, ZZdis, Zlab, param, fort = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{th}{vector of parameters}
  \item{S}{matrix of distinct response configurations}
  \item{R}{matrix of missing response configurations}
  \item{b}{number of response categories}
  \item{yv}{corresponding vector of frequencies}
  \item{Am}{design matrix for the logits}
  \item{XXdis}{design matrix used for estimation of the initial probabilities}
  \item{Xlab}{list of labels used for estimation of the initial probabilities}
  \item{ZZdis}{design matrix used for estimation of the transition probabilities}
  \item{Zlab}{list of labels used for estimation of the transition probabilities}
  \item{param}{type of parametrization for the transition probabilities ("multilogit" = standard multinomial logit for every row of the transition matrix, "difflogit" = multinomial logit based on the difference between two sets of parameters)}
  \item{fort}{to use fortran routine when possible}
}

\value{
  \item{lk}{log-likelihood}
  \item{sc}{score vector}
  \item{Psi}{conditional response probabilities}
  \item{be}{parameters on initial probabilities}
  \item{Ga}{parameters on transition probabilities}
  \item{Piv}{initial probabilities}
  \item{PI}{transition probabilities}
  \item{dlPsi}{matrix of derivatives of the logarithm of the conditional response probabilities}
  \item{dlPiv}{matrix of derivatives of the logarithm of the intial probabilities}
  \item{dlPI}{matrix of derivatives of the logarithm of the transition probabilities}
}

\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)}

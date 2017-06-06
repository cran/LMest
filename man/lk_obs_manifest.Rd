\name{lk_obs_manifest}

\alias{lk_obs_manifest}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute the observable log-likelihood of the LM model with covariates in the measurement model}

\description{Function that computes the observable log-likelihood of the LM model with covariates in the measurement model (internal use).}

\usage{
lk_obs_manifest(par, Y, Xd, yv, indn, lev, k, sup, G2, IPI, mod, outp = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{vector of regression parameters}
  \item{Y}{matrix of response variables}
  \item{Xd}{matrix of covariates}
  \item{yv}{vector of frequencies of response configurations}
  \item{indn}{index of the design matrix for each unit}
  \item{lev}{vector containing the number of levels of each variable}
  \item{k}{number of latent states}
  \item{sup}{support points}
  \item{G2}{design matrix }
  \item{IPI}{index of the transition matrix elements}
  \item{mod}{model (0 = LM with stationary transition, 1 = finite mixture)}
  \item{outp}{to compute the score of the observable log-likelihood}
}

\value{
  \item{lk}{log-likelihood}
  \item{U}{array containing the posterior distribution of the latent states for each response configuration and each pair of consecutive time occasions}
  \item{s}{score of the observable log-likelihood}
  }

\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), http://www.stat.unipg.it/bartolucci}

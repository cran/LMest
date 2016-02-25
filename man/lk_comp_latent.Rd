\name{lk_comp_latent}

\alias{lk_comp_latent}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Complete log-likelihood of the latent Markov model with covariates }

\description{Function that computes complete log-likelihood of the latent Markov model with covariates
in the distribution of the latent process (internal use).}

\usage{lk_comp_latent(S, R, yv, Piv, PI, Psi, k, fort = TRUE, der = FALSE,
       dlPsi = NULL, dlPiv = NULL, dlPI = NULL)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{matrix of distinct response configurations}
  \item{R}{matrix of missing response configurations}
  \item{yv}{corresponding vector of frequencies}
  \item{Piv}{initial probability matrix}
  \item{PI}{transition probability matrices}
  \item{Psi}{conditional response probability matrix}
  \item{k}{number of latent classes}
  \item{fort}{to use fortran routine when possible}
  \item{der}{to compute derivatives}
  \item{dlPsi}{matrix of derivatives of the logarithm of the conditional response probabilities}
  \item{dlPiv}{matrix of derivatives of the logarithm of the intial probabilities}
  \item{dlPI}{matrix of derivatives of the logarithm of the transition probabilities}
}

\value{
  \item{lk}{log-likelihood}
  \item{Phi}{matrix of the conditional probabilities of the observed response configurations}
  \item{L}{matrix of the forward probabilities}
  \item{pv}{vector of marginal probabilities}
  \item{dlk}{derivatives of the log-likelihood}
  \item{dlPhi}{matrix of derivatives of the log-conditional probabilities of the observed response configurations}
  \item{dlL}{matrix of derivatives of the log-forward probabilities}
  \item{dlL2}{matrix of second derivatives of the log-forward probabilities}
  \item{dlpv}{matrix of derivatives of the log-marginal probabilities}
}

\references{Baum, L. E., Petrie, T., Soules, G., and Weiss, N. (1970). 
A maximization technique occurring in the statistical analysis of probabilistic functions 
of Markov chains. Annals of Mathematical Statistics, 41, 164-171.}


\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), http://www.stat.unipg.it/bartolucci}
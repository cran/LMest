\name{prob_post_cov}

\alias{prob_post_cov}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute posterior probabilities.}

\description{Function that uses backward recursion to compute posterior probabilities (internal funtion).}

\usage{prob_post_cov(S, yv, Psi, Piv, PI, Phi, L, pv, der = FALSE, fort = TRUE,
       dlPhi = NULL, dlPiv = NULL, dlPI = NULL, dlL = NULL, dlL2 = NULL, dlpv = NULL)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{matrix of distinct response configurations}
  \item{yv}{corresponding vector of frequencies}
  \item{Psi}{conditional response probabilities}
  \item{Piv}{initial probability matrix}
  \item{PI}{transition probability matrices}
  \item{Phi}{matrix of the conditional probabilities of the observed response configurations}
  \item{L}{matrix of the forward probabilities}
  \item{pv}{vector of marginal probabilities}
  \item{der}{to compute derivatives}
  \item{fort}{to use fortran routine when possible}
  \item{dlPhi}{matrix of derivatives of the log-conditional probabilities of the observed response configurations}
  \item{dlPiv}{matrix of derivatives of the logarithm of the intial probabilities}
  \item{dlPI}{matrix of derivatives of the logarithm of the transition probabilities}
  \item{dlL}{matrix of derivatives of the log-forward probabilities}
  \item{dlL2}{matrix of second derivatives of the log-forward probabilities}
  \item{dlpv}{matrix of derivatives of the log-marginal probabilities}
}

\value{
  \item{U }{array containing the posterior distribution of the latent states for each response configuration and each pair of consecutive time occasions}
  \item{V }{array containing the posterior distribution of the latent states for each response configuration and time occasion}
  \item{dlU}{matrix of derivatives of the logarithm of U}
  \item{dlV}{matrix of derivatives of the logarithm of V}
}

\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), http://www.stat.unipg.it/bartolucci}

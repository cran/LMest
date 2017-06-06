\name{complk}

\alias{complk}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Complete log-likelihood of the basic latent Markov model}

\description{Function that computes complete log-likelihood of the basic latent Markov model
(internal use).}

\usage{complk(S, R, yv, piv, Pi, Psi, k)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{matrix of distinct response configurations (n x TT x r)}
  \item{R}{matrix of missing response configurations}
  \item{yv}{corresponding vector of frequencies}
  \item{piv}{vector of initial probabilities}
  \item{Pi}{transition probability matrix}
  \item{Psi}{conditional response probabilities }
  \item{k}{number of latent classes}
}

\value{
  \item{lk}{log-likelihood}
  \item{Phi}{matrix of the conditional probabilities of the observed response configurations}
  \item{L}{matrix of the forward probabilities}
  \item{pv}{vector of marginal probabilities}
}

\references{Baum, L. E., Petrie, T., Soules, G., and Weiss, N. (1970). 
A maximization technique occurring in the statistical analysis of probabilistic functions 
of Markov chains. Annals of Mathematical Statistics, 41, 164-171.}


\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), http://www.stat.unipg.it/bartolucci}
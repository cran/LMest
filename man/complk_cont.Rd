\name{complk_cont}

\alias{complk_cont}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Complete log-likelihood of the basic latent Markov model for continuous outcomes}

\description{Function that computes complete log-likelihood of the basic latent Markov model for continuous outcomes (internal use).}

\usage{complk_cont(Y, piv, Pi, Mu, Si, k)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{array of continuous outcomes (n x TT x r)}
  \item{piv}{vector of initial probabilities}
  \item{Pi}{transition probability matrix}
  \item{Mu}{matrix of conditional means for the response variables (r x k) }
  \item{Si}{var-cov matrix common to all states (r x r) }
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
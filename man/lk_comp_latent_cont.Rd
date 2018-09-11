\name{lk_comp_latent_cont}

\alias{lk_comp_latent_cont}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Complete log-likelihood of the latent Markov model for continuous outcomes with covariates }

\description{Function that computes complete log-likelihood of the latent Markov model for continuous outcomes with covariates in the distribution of the latent process (internal use).}

\usage{lk_comp_latent_cont(Y, yv, Piv, PI, Mu, Si, k)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{array of continuous response outcomes}
  \item{yv}{corresponding vector of frequencies}
  \item{Piv}{initial probability matrix}
  \item{PI}{transition probability matrices}
  \item{Mu}{matrix of conditional means of the response variables}
  \item{Si}{var-cov matrix common to all states}
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
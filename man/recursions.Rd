\name{recursions}
\alias{recursions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Recursions used by est_lm_basic}
\description{Implementation of the recursions to compute manifest probability of the responses and 
	posterior distribution of the latent states (internal function).}
\usage{recursions(S, R, yv, Psi, piv, Pi, k, lth, matr, Bm, Cm, bv, mod)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{matrix of distinct response configurations}
  \item{R}{matrix of missing response configurations}
  \item{yv}{corresponding vector of frequencies}
  \item{Psi}{matrix of conditional response probabilties}
  \item{piv}{initial probability vector}
  \item{Pi}{transition probability matrices}
  \item{k}{number of latent states}
  \item{lth}{internal argument}
  \item{matr}{internal argument}
  \item{Bm}{internal argument}
  \item{Cm}{internal argument}
  \item{bv}{internal argument}
  \item{mod}{model on the transition probabilities (0 for time-heter., 1 for time-homog., from 2 to (TT-1) partial homog. of that order)}
}
\value{
  \item{lk}{log-likelihood}
  \item{sc}{score vector}
  \item{F1}{internal argument}
  \item{F2}{internal argument}
  \item{F1d}{internal argument}
  \item{F2d}{internal argument}
}
\references{Baum, L. E., Petrie, T., Soules, G., and Weiss, N. (1970). 
A maximization technique occurring in the statistical analysis of probabilistic functions 
of Markov chains. Annals of Mathematical Statistics, 41, 164-171.}

\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), http://www.stat.unipg.it/bartolucci}

\name{rec1}
\alias{rec1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Recursions used by est_lm_cov_manifest}
\description{Implementation of the forward recursions to compute  
	posterior distribution of the latent states (internal function).}
\usage{rec1(Pio, las, PI)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pio}{matrix of manifest probabilities }
  \item{las}{stationary initial probability}
  \item{PI}{transition probability matrices}
}
\value{
  \item{Q}{matrices of posterior probabilities}
 }
 
\references{Baum, L. E., Petrie, T., Soules, G., and Weiss, N. (1970). 
A maximization technique occurring in the statistical analysis of probabilistic functions 
of Markov chains. Annals of Mathematical Statistics, 41, 164-171.}

\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), http://www.stat.unipg.it/bartolucci}

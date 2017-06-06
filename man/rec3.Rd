\name{rec3}
\alias{rec3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Recursions used by est_lm_cov_manifest}
\description{Implementation of the backward recursions to compute  
	posterior distribution of the latent states (internal function).}
\usage{rec3(Q, yv, PI, Pio, pim)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{posterior probabilities from the forward recursion} 
  \item{yv}{vector of frequencies of response configurations} 
  \item{PI}{transition probability matrices}
  \item{Pio}{matrix of manifest probabilities }
  \item{pim}{vector of marginal probabilities}
 
}
\value{
  \item{U}{array containing the posterior distribution of the latent states for each response configuration and time occasion}
  \item{V}{matrix containing the posterior distribution of the latent states}
 }
 
\references{Baum, L. E., Petrie, T., Soules, G., and Weiss, N. (1970). 
A maximization technique occurring in the statistical analysis of probabilistic functions 
of Markov chains. Annals of Mathematical Statistics, 41, 164-171.}

\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)}

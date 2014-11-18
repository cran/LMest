\name{draw_lm_mixed}

\alias{draw_lm_mixed}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Draws samples from the mixed LM model}

\description{Function that draws samples from the mixed LM model with specific parameters.}

\usage{
draw_lm_mixed(la, Piv, Pi, Psi, n, TT)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{la }{vector of mass probabilities for the first latent variable}
  \item{Piv }{matrix of initial probabilities of the latent Markov chain (k2 x k1)}
  \item{Pi }{set of transition matrices (k2 x k2 x k1)}
  \item{Psi }{array of conditional response probabitlies (l x k2 x r)}
  \item{n }{sample size}
  \item{TT }{number of time occasions}
}

\value{
  \item{Y }{matrix of response configurations unit by unit}
  \item{S }{matrix of distinct response configurations}
  \item{yv }{corresponding vector of frequencies}
}

\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia, http://www.stat.unipg.it/bartolucci}

\examples{
# draw a sample for 1000 units and only one response variable and 5 time occasions
k1 = 2; k2 = 3
la = rep(1/k1,k1)
Piv = matrix(1/k2,k2,k1)
Pi = array(0,c(k2,k2,k1))
Pi[,,1] = diag(k2)
Pi[,,2] = 1/k2
Psi = cbind(c(0.6,0.3,0.1),c(0.1,0.3,0.6))
out = draw_lm_mixed(la,Piv,Pi,Psi,n=1000,TT=5)
}
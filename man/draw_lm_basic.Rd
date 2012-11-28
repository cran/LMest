\name{draw_lm_basic}

\alias{draw_lm_basic}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Draws samples from LM model}

\description{Function that draws samples from the LM model with specific parameters.}

\usage{
draw_lm_basic(piv, Pi, Psi, n)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{piv }{vector of initial probabilities of the latent Markov chain}
  \item{Pi }{set of transition matrices}
  \item{Psi }{matrix of conditional response probabitlies}
  \item{n }{sample size}
}

\value{
  \item{Y }{matrix of response configurations unit by unit}
  \item{S }{matrix of distinct response configurations}
  \item{yv }{corresponding vector of frequencies}
}

\author{Francesco Bartolucci, University of Perugia, http://www.stat.unipg.it/bartolucci}

\examples{
# draw a sample for 1000 units and only one response variable
piv = c(0.7,0.3)
Pi = matrix(c(0.9,0.1,0.1,0.9),2,2)
Pi = array(Pi,c(2,2,6))
Pi[,,1] = 0
Psi = matrix(c(0.7,0.2,0.1,0.5,0.4,0.1),3,2)
Psi = array(Psi,c(3,2,1))
out = draw_lm_basic(piv, Pi, Psi, n=1000)
}
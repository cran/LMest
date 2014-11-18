\name{trans_par}
\alias{trans_par}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert matrix parametrization}
\description{Function that converts matrix parametrization (internal function).}
\usage{trans_par(par,lev,k,sup,G2,IPI,mod)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{regression parameters}
  \item{lev}{vector containing the number of levels of each variable}
  \item{k}{number of latent states}
  \item{sup}{vectors of support points}
  \item{G2}{internal argument}
  \item{IPI}{internal argument}
  \item{mod}{model (0 = LM with stationary transition, 1 = finite mixture)}

}
\value{
  \item{la}{vector of initial probabilities}
  \item{PI}{transition matrix}
  \item{rho}{parameter vector for AR}
  \item{si}{sigma of the AR process}	
  \item{par}{regression parameters}
  \item{lrho}{Fisher transformation of the autocorrelation coefficient}
  \item{tau}{regression parameters}
  }


\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)}

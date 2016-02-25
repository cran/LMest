\name{lk_ar_rho}

\alias{lk_ar_rho}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute complete log-likelihood for AR(1) latent process}

\description{Compute the complete log-likelihood for the transtion and intial probabilties under the constraint of AR(1) latent process
(internal use).}

\usage{lk_ar_rho(lrho, SUP, V, outp = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lrho}{Fisher transformation of the autocorrelation coefficient}
  \item{SUP}{matrix of the support points}
  \item{V}{array containing the posterior distribution of the latent states for each response configuration and time occasion}
  \item{outp}{to return addtional output (Wei, rho)}
}

\value{
  \item{flk}{Minus the complete log-likelihood}
}

\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)}
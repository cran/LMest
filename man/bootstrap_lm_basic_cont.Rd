\name{bootstrap_lm_basic_cont}

\alias{bootstrap_lm_basic_cont}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Parametric bootstrap for the basic LM model for continuous outcomes}

\description{Function that performs bootstrap parametric resampling to compute standard errors for the parameter estimates.}

\usage{bootstrap_lm_basic_cont(piv, Pi, Mu, Si, n, B = 100, start = 0, mod = 0, tol = 10^-6)}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{piv }{initial probability vector}
  \item{Pi }{probability transition matrices (k x k x TT)}
  \item{Mu }{matrix of conditional means for the response variables (r x k)}
  \item{Si }{var-cov matrix common to all states (r x r)}
  \item{n }{sample size}
  \item{B }{number of bootstrap samples}
  \item{start}{type of starting values (0 = deterministic, 1 = random)}
  \item{mod}{model on the transition probabilities (0 for time-heter., 1 for time-homog., from 2 to (TT-1) partial homog. of that order)}
  \item{tol}{tolerance level for convergence}
}

\value{
  \item{mMu }{average of bootstrap estimates of the conditional means of the response variables}
  \item{mSi }{average of bootstrap estimates of the var-cov matrix}
  \item{mpiv }{average of bootstrap estimates of the initial probability vector}
  \item{mPi }{average of bootstrap estimates of the transition probability matrices}
  \item{seMu }{standard errors for the conditional means of the response variables}
  \item{seSi }{standard errors for the var-cov matrix}
  \item{sepiv }{standard errors for the initial probability vector}
  \item{sePi }{standard errors for the transition probability matrices}
}

\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), http://www.stat.unipg.it/bartolucci}

\examples{
\dontrun{
# Example based on multivariate longitudinal continuous data 
# load data
require(mmm)
data(multiLongGaussian)
res = long2matrices(multiLongGaussian$ID,X=cbind(multiLongGaussian$X,multiLongGaussian$time),
      Y=cbind(multiLongGaussian$resp1, multiLongGaussian$resp2))
Y = res$YY
n = dim(Y)[1]
# fit of the Basic LM model for continuous outcomes
k = 3
out1 = est_lm_basic_cont(Y,k,mod=1)
out2 = bootstrap_lm_basic_cont(out1$piv,out1$Pi,out1$Mu,out1$Si,n,mod=1,B = 1000)
}
}
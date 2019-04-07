\name{est_lm_basic_cont}

\alias{est_lm_basic_cont}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate basic LM model for continuous outcomes}

\description{Main function for estimating the basic LM model for continuous outcomes.}

\usage{est_lm_basic_cont(Y, k, start = 0, mod = 0, tol = 10^-8, maxit = 1000,
                  out_se=FALSE, piv = NULL, Pi = NULL, Mu = NULL, Si = NULL)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{array of continuous outcomes (n x TT x r) - with missing responses, data are imputed using mix package}
  \item{k}{number of latent states}
  \item{start}{type of starting values (0 = deterministic, 1 = random, 2 = initial values in input)}
  \item{mod}{model on the transition probabilities (0 for time-heter., 1 for time-homog., from 2 to (TT-1) partial homog. of that order)}
  \item{tol}{tolerance level for convergence}
  \item{maxit}{maximum number of iterations of the algorithm}
  \item{out_se}{to compute the information matrix and standard errors}
  \item{piv }{initial value of the initial probability vector (if start=2)}
  \item{Pi }{initial value of the transition probability matrices (k x k x TT) (if start=2) }
  \item{Mu }{initial value of the conditional means (r x k) (if start=2)}
  \item{Si }{initial value of the var-cov matrix common to all states (r x r) (if start=2)}

}

\value{
  \item{lk }{maximum log-likelihood}
  \item{piv }{estimate of initial probability vector}
  \item{Pi }{estimate of transition probability matrices}
  \item{Mu }{estimate of conditional means of the response variables}
  \item{Si }{estimate of var-cov matrix common to all states}
  \item{np }{number of free parameters}
  \item{aic }{value of AIC for model selection}
  \item{bic }{value of BIC for model selection}
  \item{lkv }{log-likelihood trace at every step}
  \item{V }{array containing the posterior distribution of the latent states for each units and time occasion}
  \item{Y }{optional completed array with imputed (with missing responses)}
  \item{call}{command used to call the function}
}

\references{Bartolucci, F., Farcomeni, A. and Pennoni, F. (2013), Latent Markov Models for Longitudinal Data,
	Chapman and Hall/CRC press.}
\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), http://www.stat.unipg.it/bartolucci}

\examples{
# Example based on multivariate longitudinal continuous data 
# load data
require(mmm)
data(multiLongGaussian)
res = long2matrices(multiLongGaussian$ID,X=cbind(multiLongGaussian$X,multiLongGaussian$time),
      Y=cbind(multiLongGaussian$resp1, multiLongGaussian$resp2))
Y = res$YY

# fit of the Basic LM model for continuous outcomes
k = 3
out = est_lm_basic_cont(Y,k,mod=1,tol=10^-5)
summary(out)

}
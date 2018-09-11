\name{bootstrap_lm_cov_latent_cont}

\alias{bootstrap_lm_cov_latent_cont}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Parametric bootstrap for LM models for continuous outcomes with individual covariates in the latent model}

\description{Function that performs bootstrap parametric resampling to compute standard errors for the parameter estimates.}

\usage{bootstrap_lm_cov_latent_cont(X1, X2, param = "multilogit", Mu, Si, Be, Ga, B = 100)}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{X1 }{matrix of covariates affecting the initial probabilities (n x nc1)}
  \item{X2 }{array of covariates affecting the transition probabilities (n x TT-1 x nc2)}
  \item{param }{type of parametrization for the transition probabilities ("multilogit" = standard multinomial logit for every row of the transition matrix, "difflogit" = multinomial logit based on the difference between two sets of parameters)}
  \item{Mu }{matrix of conditional means for the response variables (r x k)}
  \item{Si }{var-cov matrix common to all states (r x r)}
  \item{Be }{parameters affecting the logit for the initial probabilities}
  \item{Ga }{parametes affecting the logit for the transition probabilities}
  \item{B}{number of bootstrap samples}
}

\value{
  \item{mMu }{average of bootstrap estimates of the conditional means for the response variables}
  \item{mSi }{average of bootstrap estimates of the var-cov matrix}
  \item{mBe }{average of bootstrap estimates of the parameters affecting the logit for the initial probabilities}
  \item{mGa }{average of bootstrap estimates of the parameters affecting the logit for the transition probabilities}
  \item{seMu }{standard errors for the conditional means }
  \item{seSi}{standard errors for the var-cov matrix}
  \item{seBe }{standard errors for the parameters in Be}
  \item{seGa }{standard errors for the parameters in Ga }
}

\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)}

\examples{
\dontrun{
# Example based on multivariate longitudinal continuous data 
# load data
require(mmm)
data(multiLongGaussian)
TT=4
res = long2matrices(multiLongGaussian$ID,X=cbind(multiLongGaussian$X,multiLongGaussian$time),
                    Y=cbind(multiLongGaussian$resp1, multiLongGaussian$resp2))
Y = res$YY
X1 =res$XX[,1,]
X2 =res$XX[,2:TT,]

# estimate the model
est = est_lm_cov_latent_cont(Y,X1,X2,k=3,output=TRUE)
out = bootstrap_lm_cov_latent_cont(X1,X2,Mu=est$Mu,Si=est$Si,Be=est$Be,Ga=est$Ga,B=1000)

}
}
\name{bootstrap_lm_cov_latent}

\alias{bootstrap_lm_cov_latent}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Parametric bootstrap for LM models with individual covariates in the latent model}

\description{Function that performs bootstrap parametric resampling to compute standard errors for the parameter estimates.}

\usage{bootstrap_lm_cov_latent(X1,X2,param="multilogit",Psi,Be,Ga,B=100,fort=TRUE)}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{X1 }{matrix of the initial probabilities (n x nc1)}
  \item{X2 }{array of the transition probabilities (n x TT-1 x nc2)}
  \item{param }{type of parametrization for the transition probabilities ("multilogit" = standard multinomial logit for every row of the transition matrix, "difflogit" = multinomial logit based on the difference between two sets of parameters)}
  \item{Psi }{matrix of conditional response probabilities}
  \item{Be }{parameters affecting the logit for the initial probabilities}
  \item{Ga }{parametes affecting the logit for the transition probabilities}
  \item{B}{number of bootstrap samples}
  \item{fort}{to use fortran routine when possible (FALSE for not use fortran)}
}

\value{
  \item{mPsi }{average of bootstrap estimates of the conditional response matrix}
  \item{mBe }{average of bootstrap estimates of the parameters affecting the logit for the initial probabilities}
  \item{mGa }{average of bootstrap estimates of the parameters affecting the logit for the transition probabilities}
  \item{sePsi }{standard errors for the conditional response matrix}
  \item{seBe }{standard errors for the parameters in Be}
  \item{seGa }{standard errors for the parameters in Ga }
}

\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)}

\examples{
\dontrun{
# Example based on self-rated health status (SRHS) data
# load SRHS data
data(data_SRHS_long)
data = data_SRHS_long

TT = 8 
head(data)
res = long2matrices(data$id,cbind(data$gender-1,data$race==2|data$race==3,
data$education==4,data$education==5,data$age-50,(data$age-50)^2/100),data$srhs)

# matrix of responses (with ordered categories from 0 to 4)
S = 5-res$YY

# matrix of covariates (for the first and the following occasions)
# colums are: gender,race,educational level (2 columns),age,age^2)
X1 =res$XX[,1,]
X2 =res$XX[,2:TT,]

# estimate the model
out1 = est_lm_cov_latent(S,X1,X2,k=2,output=TRUE,out_se=TRUE)

out2 = bootstrap_lm_cov_latent(X1,X2,Psi=out1$Psi,Be=out1$Be,Ga=out1$Ga,B=1000)
}
}
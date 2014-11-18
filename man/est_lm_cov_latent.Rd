\name{est_lm_cov_latent}

\alias{est_lm_cov_latent}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate LM model with covariates in the latent model}

\description{Main function for estimating LM model with covariates in the latent model.}

\usage{est_lm_cov_latent(S,X1,X2,yv=rep(1,ns),k,start=0,tol=10^-8,maxit=1000,param="multilogit",
                       Psi,Be,Ga,fort=TRUE,output=FALSE,out_se=FALSE,fixPsi=FALSE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{array of available configurations (n x T x r)}
  \item{X1}{design matrix for the initial probabilities (n x nc1)}
  \item{X2}{design matrix for the transition probabilities (n x T-1 x nc2)}
  \item{yv}{vector of frequencies of the available configurations}
  \item{k}{number of latent states}
  \item{start}{type of starting values (0 = deterministic, 1 = random, 2 = initial values in input)}
  \item{tol}{tolerance level for checking convergence of the algorithm}
  \item{maxit}{maximum number of iterations for convergence of the algorithm}
  \item{param}{type of parametrization for the transition probabilities ("multilogit" = standard multinomial logit for every row of the transition matrix, "difflogit" = multinomial logit based on the difference between two sets of parameters)}
 \item{Psi}{intial value of the matrix of the conditional response probabilities}
 \item{Be}{intial value of the parameters affecting the logit for the initial probabilities (if start=2)}
 \item{Ga}{intial value of the parametes affecting the logit for the transition probabilities (if start=2)}
 \item{fort}{to use fortran routine when possible (FALSE for not use fortran)} 
 \item{output}{to return additional output (V,PI,Piv,Ul)}
 \item{out_se}{to compute the information matrix and standard errors}
 \item{fixPsi}{TRUE if Psi is given in input and is not updated anymore}
}

\value{
  \item{lk }{maximum log-likelihood}
  \item{Be }{estimated array of the parameters affecting the logit for the initial probabilities}
  \item{Ga }{estimated array of the parameters affecting the logit for the transition probabilities}
  \item{Piv }{estimate of initial probability matrix}
  \item{PI }{estimate of transition probability matrices}
  \item{Psi }{estimate of conditional response probabilities}
  \item{np }{number of free parameters}
  \item{aic }{value of AIC for model selection}
  \item{bic }{value of BIC for model selection}
  \item{lkv }{log-likelihood trace at every step}
  \item{V }{array containing the posterior distribution of the latent states for each response configuration and time occasion}
  \item{Ul }{matrix containing the predicted sequence of latent states by the local deconding method}
  \item{sePsi}{standard errors for the conditional response matrix}
  \item{seBe}{standard errors for Be}
  \item{seGa}{standard errors for Ga}
  \item{call}{command used to call the function}
}

\references{Bartolucci, F., Farcomeni, A. and Pennoni, F. (2013), Latent Markov Models for Longitudinal Data,
	Chapman and Hall/CRC press.}
\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia, http://www.stat.unipg.it/bartolucci}

\examples{

\dontrun{
# Example based on self-rated health status (SRHS) data
# load SRHS data
data(data_SRHS)
data = as.matrix(data_SRHS)

n = nrow(data); T = 8 
# matrix of responses (with ordered categories from 0 to 4)
S = 5-data[,1:T]
# matrix of covariates (for the first and the following occasions)
# colums are: gender,race,educational level (2 columns),age,age^2)
t = 1
X1 = cbind(data[,9]-1,data[,10]==2|data[,10]==3,data[,11]==4,
data[,11]==5,data[,11+t]-50,(data[,11+t]-50)^2/100)
X2 = array(0,c(n,6,T-1))
for(t in 2:T){
  X2[,,t-1] = cbind(data[,9]-1,data[,10]==2|data[,10]==3,data[,11]==4,
  data[,11]==5,data[,11+t]-50,(data[,11+t]-50)^2/100)
}
X2 = aperm(X2,c(1,3,2))
# estimate the model
est2f = est_lm_cov_latent(S,X1,X2,k=2,output=TRUE,out_se=TRUE)
# matrice probabilita' media
PI = round(apply(est2f$PI[,,,2:8],c(1,2),mean),4)

# Transition probability matrix for white females with high educational level
ind1 = (data[,9]==2 & data[,10]==1 & data[,11]==5)
PI1 = round(apply(est2f$PI[,,ind1,2:8],c(1,2),mean),4)

# Transition probability matrix for non-white male, low educational level
ind2 = (data[,9]==1 & data[,10]>1 & data[,11]<4)
PI2 = round(apply(est2f$PI[,,ind2,2:8],c(1,2),mean),4)
}
}
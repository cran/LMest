\name{est_lm_cov_manifest}

\alias{est_lm_cov_manifest}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate LM model with covariates in the measurement model}

\description{Main function for estimating LM model with covariates in the measurement model based on a cumulative logit parameterization.}

\usage{est_lm_cov_manifest(S,X,lev,k,q,mod,tol=10^-8,maxit=1000,start=0,mu=NULL,al=NULL,
                     be=NULL,si=NULL,rho=NULL,la=NULL,PI=NULL,output=FALSE,out_se=FALSE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{array of available configurations (n x T x r)}
  \item{X}{array (np x npar x n x T) of covariates with eventually includes lagged response (np = number of marginal logits, npar = number of covariates)}
  \item{lev}{vector containing the number of levels of each variable}
  \item{k}{number of latent states}
  \item{q}{number of support points for AR}
  \item{mod}{model (0 = LM with stationary transition, 1 = finite mixture)}
  \item{tol}{tolerance for the convergence (optional) and tolerance of conditional probability if tol>1 then return}
  \item{maxit}{maximum number of iterations for convergence of the algorithm}
  \item{start}{equal to 1 for random starting values (optional)}
 \item{mu}{starting value for mu (optional)}
 \item{al}{starting value for al (optional)}
 \item{be}{starting value for be (optional)}
 \item{si}{starting value for si (optional)} 
 \item{rho}{starting value for rho (optional)}
 \item{la}{starting value for la (optional)}
 \item{PI}{starting value for PI (optional)}
 \item{output}{to return additional output (PRED0, PRED1)}
 \item{out_se}{TRUE for computing information matrix and standard error}
}

\value{
	\item{mu }{vector of cuptpoints}
	\item{al }{support points for the latent states}
	\item{be }{estimate of the vector of regression parameters}
	\item{si }{sigma of the AR process}
	\item{rho }{parameter vector for AR}
	\item{la }{vector of initial probabilities}
	\item{PI }{transition matrix}
 	\item{lk }{maximum log-likelihood}
 	\item{npa }{number of parameters}
 	\item{aic }{value of AIC index}
 	\item{bic }{value of BIC index}
 	\item{PRED0 }{prediction of latent state}
 	\item{PRED1 }{prediction of the overall latent effect}
 	\item{sebe }{standard errors for the regression parameters be}
 	\item{selrho }{standard errors for logit type transformation of rho}
 	\item{J1 }{information matrix}
	\item{call}{command used to call the function}
}

\references{Bartolucci, F., Farcomeni, A. and Pennoni, F. (2013), Latent Markov Models for Longitudinal Data,
	Chapman and Hall/CRC press.
	Bartolucci, F., Bacci, S. and Pennoni, F. (2013), Longitudinal analysis of the self-reported health status by mixture 
	latent autoregressive models, Journal of the Royal Statistical Society - series C, in press}
\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)}

\examples{

\dontrun{
# Example based on self-rated health status (SRHS) data

# load SRHS data
data(data_SRHS)
data = as.matrix(data_SRHS)

n = nrow(data); T = 8
r = 1	 		# number of response variables
lev = 5 		# number of response categories for each response variable
np = sum(lev)-r	# number of marginal logits
npar = 6      # number of covariates 
# *** design matrices X
X = array(0,c(np,npar,n,T))
for(i in 1:n) for(t in 1:T){
   X[,,i,t] = rep(1,4)\%o\%c(data[i,9]-1,(data[i,10]==2||data[i,10]==3),
   data[i,11]==4,data[i,11]==5,data[i,11+t]-50,(data[i,11+t]-50)^2)
}

S = 5-data[,1:T]

# *** fit stationary LM model
res0 = vector("list",10); tol = 10^-5;
for(k in 1:10){
  res0[[k]] = est_lm_cov_manifest(S,X,lev,k,1,mod=0,tol)
   save.image("example_SRHS3_temp.RData")
}

# *** fit the mixture latent auto-regressive model
tol = 0.005
res = vector("list",4)
k=1
q = 51
res[[k]]=est_lm_cov_manifest(S,X,lev,k,q,mod=1,tol,outp=T)
for(k in 2:4) res[[k]]=est_lm_cov_manifest(S,X,lev,k,q=61,mod=1,tol,output=T)
}
}
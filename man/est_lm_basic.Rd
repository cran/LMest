\name{est_lm_basic}

\alias{est_lm_basic}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate basic LM model}

\description{Main function for estimating the basic LM model.}

\usage{est_lm_basic(S, yv, k, start = 0, mod = 0, tol = 10^-6)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{matrix of available configurations (n by T by r)}
  \item{yv}{vector of frequencies of the available configurations}
  \item{k}{number of latent states}
  \item{start}{type of starting values (0 = deterministic, 1 = random)}
  \item{mod}{model on the transition probabilities (0 for time-heter., 1 for time-homog., 2-(T-1) partial homog. of that order)}
  \item{tol}{tolerance level for convergence}
}

\value{
  \item{lk }{maximum log-likelihood}
  \item{piv }{estimate of initial probability vector}
  \item{Pi }{estimate of transition probability matrices}
  \item{Psi }{estimate of conditional response probabilities}
  \item{np }{number of free parameters}
  \item{aic }{value of AIC for model selection}
  \item{bic }{value of BIC for model selection}
  \item{lkv }{log-likelihood trace at every step}
  \item{V }{array containing the posterior distribution of the latent states for each response configuration and time occasion}
}

\references{Bartolucci, F., Farcomeni, A. and Pennoni, F. (2010), An overview of latent Markov models for longitudinal categorical data, http://arxiv.org/abs/1003.2804.}
\author{Francesco Bartolucci, University of Perugia, http://www.stat.unipg.it/bartolucci}

\examples{
# Example of drug consumption data
# load data
data(data_drug)
data_drug = as.matrix(data_drug)
S = data_drug[,1:5]-1
yv = data_drug[,6]
n = sum(yv)
# fit of the Basic LC model
k = 3
out = est_lm_basic(S,yv,k,mod=1)

\dontrun{
# Example based on criminal data
# load criminal data
data = data(data_criminal)
out = long2wide(data_criminal,"id","time","sex",c("y1","y2","y3","y4","y5","y6","y7","y8","y9","y10"),aggr=T,full=0)	
XX = out$XX; for(j in 1:dim(XX)[1]) XX[j,] = max(XX[j,])
YY = out$YY
freq = out$freq
n1 = sum(freq[XX[,1]==1])
n2 = sum(freq[XX[,1]==2])
nc1 = floor(n1/0.34); nc2 = floor(n2/0.085);
dd = dim(YY)
YY1 = YY; YY = array(0,c(dd[1]+2,dd[2],dd[3]))	
YY[3:(2+dd[1]),,] = YY1 
XX = rbind(rep(1,dd[2]),rep(2,dd[2]),XX)
freq = c(nc1-n1,nc2-n2,freq)
n = sum(freq)
# fit basic LM model with increasing number of states to select the most suitable
Res0 = vector("list",7)
for(k in 1:7){
    Res0[[k]] = est_lm_basic(YY,freq,k,mod=1,tol=10^-4)
    save(list = ls(),file="example_criminal_temp.RData")
}
out1 = Res0[[6]]
}
}
\name{decoding}

\alias{decoding}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perform local and global deconding}

\description{Function that performs local deconding and global decoding (Viterbi) from the output for est_lm_basic, est_lm_cov_latent, est_lm_cov_manifest, and est_lm_mixed.}

\usage{
decoding(est,Y,X1=NULL,X2=NULL,fort=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est }{output from est_lm_basic, est_lm_cov_latent, est_lm_cov_manifest, or est_lm_mixed}
  \item{Y }{single vector of matrix of responses}
  \item{X1 }{matrix of covariates on the initial probabilities (est_lm_cov_latent) or on the responses (est_lm_cov_manifest)}
  \item{X2 }{array of covariates on the transition probabilites}
  \item{fort }{to use Fortran routines}
}

\value{
  \item{Ul }{matrix of local deconded states corresponding to each row of Y}
  \item{Ug }{matrix of global deconded states corresponding to each row of Y}
}

\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia, http://www.stat.unipg.it/bartolucci}

\examples{
# example for the output from est_lm_basic
data(data_drug)
data_drug = as.matrix(data_drug)
S = data_drug[,1:5]-1
yv = data_drug[,6]
n = sum(yv)
# fit of the Basic LC model
k = 3
est = est_lm_basic(S,yv,k,mod=1)
# deconding for a single sequence
out1 = decoding(est,S[1,])
# deconding for all sequences
out2 = decoding(est,S)


# example for the output from est_lm_cov_latent with difflogit parametrization
data(data_SRHS)
data = as.matrix(data_SRHS[1:200,])
n = nrow(data); T = 8 
# matrix of responses (with ordered categories from 0 to 4)
S = 5-data[,1:T]
# matrix of covariates (for the first and the following occasions)
# colums are: gender,race,educational level (2 columns),age,age^2)
t = 1
X1 = cbind(data[,9]-1,data[,10]==2|data[,10]==3,data[,11]==4,
data[,11]==5,data[,11+t]-50,(data[,11+t]-50)^2/100)
X2 = array(0,c(n,6,T-1))
for(t in 2:T) X2[,,t-1] = cbind(data[,9]-1,data[,10]==2|data[,10]==3,
data[,11]==4,data[,11]==5,data[,11+t]-50,(data[,11+t]-50)^2/100)
X2 = aperm(X2,c(1,3,2))
# estimate the model
est = est_lm_cov_latent(S,X1,X2,k=2,output=TRUE,param="difflogit")
# deconding for a single sequence
out1 = decoding(est,S[1,],X1[1,],X2[1,,])
# deconding for all sequences
out2 = decoding(est,S,X1,X2)
}

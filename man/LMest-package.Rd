\name{LMest-package}
\alias{LMest-package}
\alias{LMest}
\docType{package}
\title{Fit latent Markov models}
\description{Set of functions to fit latent Markov models in the basic version.}
\details{
\tabular{ll}{
Package: \tab LMest\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-11-27\cr
License: \tab Unlimited\cr
}
The most important functions is \code{est_lm_basic} which estimates the basic LM model 
}
\author{
Francesco Bartolucci, University of Perugia(IT), http://www.stat.unipg.it/bartolucci\cr
Maintainer: Francesco Bartolucci <bart@stat.unipg.it>
}
\references{Bartolucci, F., Farcomeni, A. and Pennoni, F. (2010), An overview of latent Markov models for longitudinal categorical data, http://arxiv.org/abs/1003.2804.}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
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
out1 = est_lm_basic(S,yv,k,mod=1)
\dontrun{
out2 = bootstrap_lm_basic(out1$piv,out1$Pi,out1$Psi,n,mod=1,B=1000)
}
}

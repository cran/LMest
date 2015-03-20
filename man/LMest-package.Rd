\name{LMest-package}
\alias{LMest-package}
\alias{LMest}
\docType{package}
\title{Fit latent Markov models}
\description{Set of functions to fit latent Markov models in the basic version and in the extended version with individual covariates.}
\details{
\tabular{ll}{
Package: \tab LMest\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2015-03-19\cr
License: \tab Unlimited\cr
}
The most important functions are \code{est_lm_basic}, \code{est_lm_cov_manifest}, and\code{est_lm_cov_latent} which estimate the basic LM model and its extensions including covariates
}
\author{
Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)
Maintainer: Francesco Bartolucci <bart@stat.unipg.it>
}
\references{Bartolucci, F., Farcomeni, A. and Pennoni, F. (2013), Latent Markov Models for Longitudinal Data,
	Chapman and Hall/CRC press.}
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
out1 = est_lm_basic(S,yv,k,mod=1,out_se=TRUE)
\dontrun{
out2 = bootstrap_lm_basic(out1$piv,out1$Pi,out1$Psi,n,mod=1,B=1000)
}
}

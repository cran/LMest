\name{long2matrices}
\alias{long2matrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{From data in the long format to data in array format}
\description{Function that transforms data in the long format to data in array format.}
\usage{
long2matrices(id,X,Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{vector of subjects id}
  \item{X}{matrix of covariates in long format}
  \item{Y}{matrix of responses in long format}
}
\value{
  \item{XX}{array of covariates (n x TT x nc)}
  \item{YY}{array of responses (n x TT x r)}
}

\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), http://www.stat.unipg.it/bartolucci}

\examples{
# Example based on SRHS data
# load SRHS data
data(data_SRHS_long)
data = data_SRHS_long
head(data)
X = cbind(data$gender-1,data$race==2|data$race==3,data$education==4,data$education==5,
data$age-50,(data$age-50)^2/100)
Y = data$srhs
res = long2matrices(data$id,X,Y)
}


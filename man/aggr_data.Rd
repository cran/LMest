\name{aggr_data}

\alias{aggr_data}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Function that aggregates data}

\description{Given a matrix of configurations (covariates and responses) unit by unit, this
	function finds the corresponding matrix
	of distinct configurations and the corresponding vector of frequencies.}

\usage{
aggr_data(data)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{matrix of covariate and response configurations unit by unit}
}

\value{
  \item{data_dis}{matrix of distinct configurations}
  \item{freq}{vector of corresponding frequencies}
  \item{label}{the index of each provided response among the distinct ones}
}


\author{Francesco Bartolucci, University of Perugia, http://www.stat.unipg.it/bartolucci}

\examples{
# draw a matrix of random responses and find distinct responses
X = matrix(sample(5,100,replace=TRUE),50,2)
out = aggr_data(X)
}
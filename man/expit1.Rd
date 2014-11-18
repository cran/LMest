\name{expit1}

\alias{expit1}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute the expit function with respect to a reference category.}

\description{Function to compute the inverse of the logit function with respect to a reference category.}

\usage{expit1(lp,ref=1)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lp}{a vector which contains the values to be transformed}
  \item{ref}{the reference category}
}

\value{
  \item{p}{The transformed vector}
  \item{Der}{The matrix of derivatives}
}

\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)}

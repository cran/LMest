\name{prob_multilogit}
\alias{prob_multilogit}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute multinomial probabilities}

\description{The function computes multinomial probabilities (internal function).}

\usage{prob_multilogit(Xdis, be, label, fort = TRUE,der = FALSE)}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{Xdis}{array of all covariates (two- or three-dimensional)}
  \item{be}{initial value of the parameter vector}
  \item{label}{label associated to every covariate configuration}
  \item{fort}{to use fortran routines when possible}
  \item{der}{to compute derivatives}
}

%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}

\value{
\item{P}{estimated matrix of probabilities}
\item{Pdis}{estimated matrix of distinct probabilities}
\item{dP}{derivatives of the estimated matrix of probabilities}
\item{dPdis}{derivatives of the estimated matrix of distinct probabilities}

}

\author{Francesco Bartolucci, Silvia Pandolfi - University of Perugia (IT)}

%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
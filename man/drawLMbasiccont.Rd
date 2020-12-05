\name{drawLMbasiccont}
\alias{drawLMbasiccont}

\title{Draw samples from the basic LM model for continuous outcomes}

\description{Function that draws samples from the basic LM model for continuous outcomes with specific parameters.}

\usage{
drawLMbasiccont(piv, Pi, Mu, Si, n, est = NULL,
                format = c("long","matrices"), seed = NULL)
}

\arguments{
  \item{piv}{vector of initial probabilities of the latent Markov chain}
  \item{Pi}{set of transition probabilities matrices (k x k x TT)}
  \item{Mu}{matrix of conditional means for the response variables (r x k)}
  \item{Si}{var-cov matrix common to all states (r x r)}
  \item{n}{sample size}
  \item{est}{object of class \code{LMbasiccont} (\code{\link{LMbasiccont-class}})}
  \item{format}{character string indicating the format of final responses matrix}
  \item{seed}{an integer value with the random number generator state}
}

\value{
  \item{Y}{array of continuous outcomes (n x TT x r)}
  \item{piv}{vector of initial probabilities of the latent Markov chain}
  \item{Pi}{set of transition probabilities matrices (k x k x TT)}
  \item{Mu}{matrix of conditional means for the response variables (r x k)}
  \item{Si}{var-cov matrix common to all states (r x r)}
  \item{n}{sample size}
  \item{est}{object of class \code{LMbasiccont} (\code{\link{LMbasiccont-class}})}
}

\author{Francesco Bartolucci,
  Silvia Pandolfi,
  Fulvia Pennoni,
  Alessio Farcomeni,
  Alessio Serafini
}

\examples{
\dontrun{
# draw a sample for 1000 units and 3 response variable
n <- 1000
TT <- 5
k <- 2
r <- 3 #number of response variables

piv <- c(0.7,0.3)
Pi <- matrix(c(0.9,0.1,0.1,0.9), k, k)
Pi <- array(Pi, c(k, k, TT))
Pi[,,1] <- 0
Mu <- matrix(c(-2,-2,0,0,2,2), r, k)
Si <- diag(r)
out <- drawLMbasiccont(piv, Pi, Mu, Si, n)

data(data_long_cont)

est <- lmestCont(responsesFormula = Y1 + Y2 + Y3 ~ NULL,
                 index = c("id", "time"),
                 data = data_long_cont,
                 k = 3,
                 modBasic = 1,
                 tol = 10^-5)

out2 <- drawLMbasiccont(est = est, n = 100, format = "long", seed = 4321)

}


}

\name{drawLMbasic}
\alias{drawLMbasic}

\title{Draw samples from the basic LM model}

\description{Function that draws samples from the basic LM model.}

\usage{
drawLMbasic(piv, Pi, Psi, n, est = NULL, format = c("long","matrices"), seed = NULL)
}

\arguments{
  \item{piv}{vector of initial probabilities of the latent Markov chain}
  \item{Pi}{set of transition probabilities matrices (k x k x TT)}
  \item{Psi}{array of conditional response probabitlies (mb x k x r)}
  \item{n}{sample size}
  \item{est}{object of class \code{LMbasic} (\code{\link{LMbasic-class}})}
  \item{format}{character string indicating the format of final responses matrix}
  \item{seed}{an integer value with the random number generator state}
}

\value{
  \item{Y}{matrix of response configurations unit by unit}
  \item{S}{matrix of distinct response configurations}
  \item{yv}{corresponding vector of frequencies}
   \item{piv}{vector of initial probabilities of the latent Markov chain}
  \item{Pi}{set of transition probabilities matrices (k x k x TT)}
  \item{Psi}{array of conditional response probabitlies (mb x k x r)}
  \item{n}{sample size}
  \item{est}{object of class \code{LMbasic} (\code{\link{LMbasic-class}})}
}

\author{Francesco Bartolucci,
  Silvia Pandolfi,
  Fulvia Pennoni,
  Alessio Farcomeni,
  Alessio Serafini
}


\examples{
\dontrun{

# draw a sample for 1000 units and only one response variable
n <- 1000
TT <- 6
k <- 2
r <- 1 #number of response variables
mb <- 3 #maximum number of response categories

piv <- c(0.7,0.3)
Pi <- matrix(c(0.9,0.1,0.1,0.9), k, k)
Pi <- array(Pi, c(k, k, TT))
Pi[,,1] <- 0
Psi <- matrix(c(0.7,0.2,0.1,0.5,0.4,0.1), mb, k)
Psi <- array(Psi, c(mb, k, r))

out <- drawLMbasic(piv, Pi, Psi, n = 1000)


data("data_SRHS_long")
SRHS <- data_SRHS_long[1:2400,]

SRHS$srhs <- 5 - SRHS$srhs

est <- lmest(responsesFormula = srhs ~ NULL,
             index = c("id","t"),
             data = SRHS,
             k = 3)

out1 <- drawLMbasic(est = est, format = "matrices", seed = 4321, n = 100)
}
}

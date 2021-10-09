\name{drawLMmixed}
\alias{drawLMmixed}

\title{Draws samples from the mixed LM model}

\description{Function that draws samples from the mixed LM model with specific parameters.}

\usage{
drawLMmixed(la, Piv, Pi, Psi, n, TT, est = NULL,
            format = c("long","matrices"), seed = NULL)
}

\arguments{
  \item{la }{vector of mass probabilities for the first latent variable}
  \item{Piv }{matrix of initial probabilities of the latent Markov chain (k2 x k1)}
  \item{Pi }{set of transition matrices (k2 x k2 x k1)}
  \item{Psi }{array of conditional response probabitlies (mb x k2 x r)}
  \item{n }{sample size}
  \item{TT }{number of time occasions}
  \item{est}{object of class \code{LMmixed} (\code{\link{LMmixed-class}})}
  \item{format}{character string indicating the format of final responses matrix}
  \item{seed}{an integer value with the random number generator state}
}

\value{
  \item{Y}{matrix of response configurations unit by unit}
  \item{S}{matrix of distinct response configurations}
  \item{yv}{corresponding vector of frequencies}
  \item{la }{vector of mass probabilities for the first latent variable}
  \item{Piv }{matrix of initial probabilities of the latent Markov chain (k2 x k1)}
  \item{Pi }{set of transition matrices (k2 x k2 x k1)}
  \item{Psi }{array of conditional response probabitlies (mb x k2 x r)}
  \item{n }{sample size}
  \item{TT }{number of time occasions}
  \item{est}{object of class \code{LMmixed} (\code{\link{LMmixed-class}})}
}

\author{Francesco Bartolucci,
  Silvia Pandolfi,
  Fulvia Pennoni,
  Alessio Farcomeni,
  Alessio Serafini
}

\examples{
# draw a sample for 1000 units and only one response variable and 5 time occasions
k1 <- 2
k2 <- 3
la <- rep(1/k1, k1)
Piv <- matrix(1/k2, k2, k1)
Pi <- array(0, c(k2, k2, k1))
Pi[,,1] <- diag(k2)
Pi[,,2] <- 1/k2
Psi <- cbind(c(0.6,0.3,0.1), c(0.1,0.3,0.6), c(0.3,0.6,0.1))
out <- drawLMmixed(la, Piv, Pi, Psi, n = 1000, TT = 5)


\dontrun{
# Example based on criminal data

data(data_criminal_sim)
data_criminal_sim = data.frame(data_criminal_sim)
# Estimate mixed LM model for females

responsesFormula <- lmestFormula(data = data_criminal_sim,
                                  response = "y")$responsesFormula
est <- lmestMixed(responsesFormula = responsesFormula,
                  index = c("id","time"),
                  k1 = 2,
                  k2 = 2,
                  data = data_criminal_sim[data_criminal_sim$sex == 2,])

out <- drawLMmixed(est = est, n = 100,  seed = 4321)

}
}

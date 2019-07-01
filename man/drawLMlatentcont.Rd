\name{drawLMlatentcont}
\alias{drawLMlatentcont}

\title{Draw samples from LM model for continuous outcomes with covariates in the latent model}

\description{Function that draws samples from the LM model for continuous outcomes with individual covariates with specific parameters.}

\usage{
drawLMlatentcont(Mu, Si, Be, Ga, latentFormula, data, index,
                 paramLatent = c("multilogit","difflogit"), est = NULL,
                 format = c("long","matrices"), fort = TRUE, seed = NULL)
}

\arguments{
  \item{Mu}{array of conditional means for the response variables (r x k)}
  \item{Si}{var-cov matrix common to all states (r x r)}
  \item{Be}{parameters affecting the logit for the initial probabilities}
  \item{Ga}{parametes affecting the logit for the transition probabilities}
  \item{latentFormula}{a symbolic description of the model to be fitted. A detailed description is given in \code{\link{lmestCont}}}
  \item{data}{a data frame in long format, with rows corresponding to observations and columns corresponding to variables, a column corresponding to time occasions and a column containing the unit identifier}
  \item{index}{a character vector with two elements indicating the name of the "id" column as first element and the "time" column as second element}
  \item{paramLatent}{type of parametrization for the transition probabilities ("multilogit" = standard multinomial logit for every row of the transition matrix, "difflogit" = multinomial logit based on the difference between two sets of parameters)}
   \item{est}{object of class \code{LMlatentcont} (\code{\link{LMlatentcont-class}})}
  \item{format}{character string indicating the format of final responses matrix}
  \item{fort}{to use fortran routine when possible (FALSE for not use fortran)}
  \item{seed}{an integer value with the random number generator state}
  }

\value{
  \item{Y}{array of continuous outcomes (n x TT x r)}
  \item{U}{matrix containing the sequence of latent states (n x TT)}
  \item{Mu}{array of conditional means for the response variables (r x k)}
  \item{Si}{var-cov matrix common to all states (r x r)}
  \item{Be}{parameters affecting the logit for the initial probabilities}
  \item{Ga}{parametes affecting the logit for the transition probabilities}
  \item{latentFormula}{a symbolic description of the model to be fitted. A detailed description is given in \code{\link{lmestCont}}}
  \item{data}{a data frame in long format, with rows corresponding to observations and columns corresponding to variables, a column corresponding to time occasions and a column containing the unit identifier}
  \item{est}{object of class \code{LMlatentcont} (\code{\link{LMlatentcont-class}})}

 }

\author{Francesco Bartolucci,
  Silvia Pandolfi,
  Fulvia Pennoni,
  Alessio Farcomeni,
  Alessio Serafini
}

\examples{

\dontrun{
require(mmm)
data(multiLongGaussian)

t <- rep(1:4, times  = max(multiLongGaussian$ID))
multiLongGaussian <- data.frame(t = t, multiLongGaussian)

est <- lmestCont(responsesFormula = resp1 + resp2 ~ NULL,
                 latentFormula = ~ X + time,
                 index = c("ID", "t"),
                 data = multiLongGaussian,
                 k = 3)

out <- drawLMlatentcont(est = est,format = "matrices", seed = 4321)
out1 <- drawLMlatentcont(latentFormula = ~ X + time,data = multiLongGaussian,
                         index = c("ID", "t"),
                         Mu = est$Mu,Si = est$Si,
                         Be = est$Be,Ga = est$Ga, fort=TRUE, seed = 4321, format = "matrices")

}
}

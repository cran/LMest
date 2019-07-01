\name{drawLMlatent}

\alias{drawLMlatent}

\title{Draw samples from LM model with covariates in the latent model}

\description{Function that draws samples from the LM model with individual covariates with specific parameters.}

\usage{
drawLMlatent(Psi, Be, Ga, latentFormula, data, index,
             paramLatent = c("multilogit","difflogit"), est = NULL,
             format = c("long","matrices"), fort = TRUE, seed = NULL)
}

\arguments{
 \item{Psi}{array of conditional response probabilities (mb x k x r)}
  \item{Be}{parameters affecting the logit for the initial probabilities}
  \item{Ga}{parametes affecting the logit for the transition probabilities}
  \item{latentFormula}{a symbolic description of the model to be fitted. Detailed description is given in \code{\link{lmest}}}
  \item{data}{a data frame in long format, with rows corresponding to observations and columns corresponding to variables, a column corresponding to time occasions and a column containing the unit identifier}
  \item{index}{a character vector with two elements indicating the name of the "id" column as first element and the "time" column as second element}
  \item{paramLatent}{type of parametrization for the transition probabilities ("multilogit" = standard multinomial logit for every row of the transition matrix, "difflogit" = multinomial logit based on the difference between two sets of parameters)}
  \item{est}{object of class \code{LMlatent} (\code{\link{LMlatent-class}})}
  \item{format}{character string indicating the format of final responses matrix}
  \item{fort}{to use fortran routine when possible (FALSE for not use fortran)}
  \item{seed}{an integer value with the random number generator state}
}

\value{
  \item{Y}{matrix of response configurations}
  \item{U}{matrix containing the sequence of latent states (n x TT)}
  \item{Psi}{array of conditional response probabilities (mb x k x r)}
  \item{Be}{parameters affecting the logit for the initial probabilities}
  \item{Ga}{parametes affecting the logit for the transition probabilities}
  \item{latentFormula}{a symbolic description of the model to be fitted. Detailed description is given in \code{\link{lmest}}}
  \item{data}{a data frame in long format, with rows corresponding to observations and columns corresponding to variables, a column corresponding to time occasions and a column containing the unit identifier}
  \item{est}{object of class \code{LMlatent} (\code{\link{LMlatent-class}})}

 }

\author{Francesco Bartolucci,
  Silvia Pandolfi,
  Fulvia Pennoni,
  Alessio Farcomeni,
  Alessio Serafini
}


\examples{
\dontrun{
data(data_SRHS_long)
dataSRHS <- data_SRHS_long
data_SRHS_long$srhs <- 5 - data_SRHS_long$srhs
est  <- lmest(responsesFormula = srhs ~ NULL,
              latentFormula =  ~
              I(gender - 1) +
              I( 0 + (race == 2) + (race == 3)) +
              I(0 + (education == 4)) +
              I(0 + (education == 5)) +
              I(age - 50) + I((age-50)^2/100),
              index = c("id","t"),
              data = data_SRHS_long,
              k = 2,
              paramLatent = "multilogit",
              start = 0)

out <- drawLMlatent(est = est,format = "matrices",seed = 4321)

out1 <- drawLMlatent(Psi = est$Psi,Be = est$Be,Ga = est$Ga,
                     data = data_SRHS_long, index = c("id","t"),
                     latentFormula =  ~
                     I(gender - 1) +
                     I( 0 + (race == 2) + (race == 3)) +
                     I(0 + (education == 4)) +
                     I(0 + (education == 5)) +
                     I(age - 50) + I((age-50)^2/100),
                     paramLatent = "multilogit", format = "matrice",
                     seed = 4321)

}
}

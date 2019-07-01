\name{bootstrap}
\alias{bootstrap}
\alias{bootstrap.LMbasic}
\alias{bootstrap.LMbasiccont}
\alias{bootstrap.LMlatent}
\alias{bootstrap.LMlatentcont}

\title{Parametric bootstrap}
\description{
Function that performs bootstrap parametric resampling to compute standard errors for the parameter estimates.
}
\usage{
bootstrap(est, ...)
\method{bootstrap}{LMbasic}(est, n = 1000, B = 100, seed = NULL, ...)
\method{bootstrap}{LMbasiccont}(est, n = 1000, B=100, seed = NULL, ...)
\method{bootstrap}{LMlatent}(est, B = 100, seed = NULL, ...)
\method{bootstrap}{LMlatentcont}(est, B = 100, seed = NULL, ...)
}
\arguments{
  \item{est}{an object obtained from a call to \code{\link{lmest}} and \code{\link{lmestCont}}}
  \item{n}{sample size}
  \item{B}{number of bootstrap samples}
  \item{seed}{an integer value with the random number generator state}
  \item{\dots}{further arguments}
}

\value{
Average of bootstrap estimates and standard errors for the model parameters in \code{est} object.
}

\author{Francesco Bartolucci,
  Silvia Pandolfi,
  Fulvia Pennoni,
  Alessio Farcomeni,
  Alessio Serafini
}

\examples{
\dontrun{
data("data_SRHS_long")
SRHS <- data_SRHS_long[1:2400,]
SRHS$srhs <- 5 - SRHS$srhs

out1 <- lmest(responsesFormula = srhs ~ NULL,
            index = c("id","t"),
            data = SRHS,
            k = 3,
            tol = 1e-8,
            start = 1,
            modBasic = 1,
            out_se = TRUE,
            seed = 123)

boot1 <- bootstrap(out1)

out2 <- lmest(responsesFormula = srhs ~ NULL,
             latentFormula =  ~
             I(gender - 1) +
             I( 0 + (race == 2) + (race == 3)) +
             I(0 + (education == 4)) +
             I(0 + (education == 5)) +
             I(age - 50) + I((age-50)^2/100),
             index = c("id","t"),
             data = SRHS,
             k = 2,
             paramLatent = "multilogit",
             start = 0)

boot2 <- bootstrap(out2)

require(mmm)
data(multiLongGaussian)

t <- rep(1:4, times  = max(multiLongGaussian$ID))
multiLongGaussian <- data.frame(t = t, multiLongGaussian)

out3 <- lmestCont(responsesFormula = resp1 + resp2 ~ NULL,
                  index = c("ID", "t"),
                  data = multiLongGaussian,
                  k =3,
                  modBasic=1,
                  tol=10^-5)

boot3 <- bootstrap(out3)

out4 <- lmestCont(responsesFormula = resp1 + resp2 ~ NULL,
                   latentFormula = ~ X + time,
                   index = c("ID", "t"),
                   data = multiLongGaussian,
                   k = 3,
                   output=TRUE)
boot4 <- bootstrap(out4)
}
}

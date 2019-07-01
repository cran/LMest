\name{se}
\alias{se}
\alias{se.LMbasic}
\alias{se.LMbasiccont}
\alias{se.LMlatent}
\alias{se.LMlatentcont}

\title{Standard errors}
\description{
  Function to compute standard errors for the parameter estimates.
}
\usage{
se(est, ...)
\method{se}{LMbasic}(est, ...)
\method{se}{LMbasiccont}(est, ...)
\method{se}{LMlatent}(est, ...)
\method{se}{LMlatentcont}(est, ...)
}
\arguments{
  \item{est}{an object obtained from a call to \code{\link{lmest}} and \code{\link{lmestCont}}}
  \item{\dots}{further arguments}
}

\value{
  Standard errors for estimates in \code{est} object.
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

# Categories rescaled to vary from 1 (“poor”) to 5 (“excellent”)

SRHS$srhs <- 5 - SRHS$srhs

out <- lmest(responsesFormula = srhs ~ NULL,
            index = c("id","t"),
            data = SRHS,
            k = 3,
            modBasic = 1,
            out_se = FALSE)
out.se <- se(out)

out1 <- lmest(responsesFormula = srhs ~ NULL,
            index = c("id","t"),
            data = SRHS,
            k = 3,
            modBasic = 1,
            out_se = TRUE)
out1.se <- se(out1)

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
out2.se <- se(out2)

require(mmm)
data(multiLongGaussian)

t <- rep(1:4, times  = max(multiLongGaussian$ID))
multiLongGaussian <- data.frame(t = t, multiLongGaussian)

out3 <- lmestCont(responsesFormula = resp1 + resp2 ~ NULL,
                  index = c("ID", "t"),
                  data = multiLongGaussian,
                  k =3,
                  modBasic = 1,
                  tol = 10^-5)
out3.se <- se(out3)

out4 <- lmestCont(responsesFormula = resp1 + resp2 ~ NULL,
                   latentFormula = ~ X + time | X + time,
                   index = c("ID", "t"),
                   data = multiLongGaussian,
                   k = 3,
                   output = TRUE)
out4.se <- se(out4)
}
}

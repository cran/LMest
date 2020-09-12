\name{plot}
\docType{methods}
\alias{plot.LMbasic}
\alias{plot.LMbasiccont}
\alias{plot.LMlatent}
\alias{plot.LMlatentcont}

\title{Plots for Generalized Latent Markov Models}

\description{Plots for outputs of LMest objects: \code{LMbasic}, \code{LMbasiccont}, \code{LMlatent}, and \code{LMlatentcont}}

\usage{

\method{plot}{LMbasic}(x,
                            what = c("modSel", "CondProb", "transitions","marginal"),
                            verbose=interactive(),...)
\method{plot}{LMlatent}(x,
                            what = c("modSel", "CondProb", "transitions","marginal"),
                            verbose=interactive(),...)
\method{plot}{LMbasiccont}(x,
                                what = c("modSel", "density", "transitions","marginal"),
                                components,verbose=interactive(),...)
\method{plot}{LMlatentcont}(x,
                                 what = c("modSel", "density", "transitions","marginal"),
                                 components, verbose=interactive(),...)
}



\arguments{
  \item{x}{an object of class \code{LMbasic}, \code{LMlatent}, \code{LMbasiccont} or \code{LMlatentcont}}
  \item{what}{a string indicating the type of plot. A detailed description is provided in the ‘Details’ section.}
  \item{components}{An integer or a vector of integers specifying the components (latent states) to be selected for the "density" plot.}
  \item{verbose}{A logical controlling if a text progress bar is displayed during the
    fitting procedure. By default is \code{TRUE} if the session is 
    interactive, and \code{FALSE} otherwise.}
    \item{...}{Unused argument.}
}
\details{The type of plots are the following:
\tabular{ll}{
\code{"modSel"} \tab plot of values of the Bayesian Information Criterion and of the Akaike Information Criterion for model selection\cr
\code{"CondProb"} \tab plot of the estimated  conditional response probabilities\cr
\code{"density"} \tab plot of the overall estimated density for continuous responses, with weights given by the estimated marginal distribution of the latent variable. For multivariate continuous responses a  contour plot is provided. If the argument \code{components} is specified, the density plot for the selected components results \cr
\code{"transitions"} \tab path diagram of the estimated transition probabilities \cr
\code{"marginal"} \tab plot of the estimated marginal distribution of the latent variable \cr
 }
If argument \code{what} is not specified, a menu of choices is proposed in an interactive session.
}


\author{Francesco Bartolucci, Silvia Pandolfi, Fulvia Pennoni, Alessio Farcomeni, Alessio Serafini}


\examples{
\dontrun{
### Plot of basic LM model

data("data_SRHS_long")
SRHS <- data_SRHS_long[1:2400,]

# Categories rescaled to vary from 0 (“poor”) to 4 (“excellent”)

SRHS$srhs <- 5 - SRHS$srhs

out <- lmest(responsesFormula = srhs ~ NULL,
            index = c("id","t"),
            data = SRHS,
            k = 1:3,
            start = 1,
            modBasic = 1,
            seed = 123)
out
summary(out)
plot(out)

### Plot of basic LM model for continuous responses

require(mmm)
data(multiLongGaussian)

t <- rep(1:4, times  = max(multiLongGaussian$ID))
multiLongGaussian <- data.frame(t = t, multiLongGaussian)

out1 <- lmestCont(responsesFormula = resp1 + resp2 ~ NULL,
                  index = c("ID", "t"),
                  data = multiLongGaussian,
                  k = 1:5,
                  modBasic=1,
                  tol=10^-5)

plot(out1,what="density")
plot(out1,what="density",components=c(1,3))

}
}

\name{mc_multinom}

\alias{mc_multinom}

\title{Estimate MC model with constraints}

\description{Function that estimates Markov chain model with suitable constraints on the transition matrix.}

\usage{
mc_multinom(formula, data, baseline = c("initial", "central"),
           model_int = c("all", "const", "dist1", "dist2", "dist3", "two", 
                         "symm", "rsymm", "diff"),
           model_cov = c("all", "const", "dist1", "dist2", "dist3", "two", 
                         "symm", "rsymm", "diff"),
           formula_init = NULL, reltol = 10^-10, output = FALSE)
}

\arguments{
  \item{formula }{A model formula specifying the response variable and the covariates affecting the transition probabilities. The response variable must be numeric, coded consecutively from \code{1} to \code{k}, where \code{k} is the number of observed states. See the \sQuote{Details} section for the specification of the covariates}
  \item{data }{A balanced data frame in long format. The first two required variables are \code{id}, identifying subjects, and \code{t}, identifying time occasions. Subject identifiers must be consecutive integers from \code{1} to the number of subjects; for each subject the time occasions must be the consecutive integers from \code{1} to the number of occasions available for that subject, which must be at least two; the rows must be ordered by subject and, within each subject, by time occasion. Missing values are not allowed in the variables entering \code{formula} and \code{formula_init}}
  \item{baseline }{Character string specifying the reference category in the multinomial logit parameterization for the transition probabilities. With \code{"initial"}, the transition towards the first state is used as the reference category for each row of the transition matrix. With \code{"central"}, the self-transition is used as the reference category for each row of the transition matrix.}
  \item{model_int }{Character string specifying the constraint structure for the intercepts of the transition logits. See \code{\link{design_matrices}} for the available options and their interpretation.}
  \item{model_cov }{Character string specifying the constraint structure for the regression coefficients associated with the transition covariates. See \code{\link{design_matrices}} for the available options and their interpretation}
  \item{formula_init }{Optional formula specifying the covariates affecting the initial probabilities. If \code{NULL}, the covariates specified in \code{formula} are also used for the initial probabilities}
  \item{reltol }{Positive numeric value specifying the relative convergence tolerance used in the Newton--Raphson iterations for the initial and transition components of the log-likelihood}
  \item{output }{Logical value. If \code{TRUE}, subject-specific initial and transition probability arrays are returned in addition to the parameter estimates}
}

\value{
  An object of class \code{"mc_multinom"} containing the following components:
  \item{la }{Estimated regression coefficients for the initial probabilities}
  \item{psi }{Estimated reduced vector of regression coefficients for the constrained transition logits}
  \item{eta }{Expanded vector of transition-logit parameters obtained by applying the design-matrix constraints to \code{psi}}
  \item{lk }{Total log-likelihood, equal to the sum of the initial-probability and transition-probability log-likelihoods}
  \item{lk1 }{Log-likelihood contribution of the initial probabilities}
  \item{lk2 }{Log-likelihood contribution of the transition probabilities}
  \item{np }{Number of free model parameters}
  \item{aic }{Akaike Information Criterion}
  \item{bic }{Bayesian Information Criterion}
  \item{sela }{Standard errors of the estimated coefficients \code{la}}
  \item{sepsi }{Standard errors of the estimated coefficients \code{psi}}
  \item{call }{Matched call used to fit the model}
  \item{Piv }{Matrix of subject-specific initial probabilities, with dimensions \code{n x k}. Returned only when \code{output = TRUE}}
  \item{PI }{Array of subject-specific transition probabilities, with dimensions \code{k x k x n x TT}, where \code{TT} is the maximum number of time occasions. Returned only when \code{output = TRUE}}
}

\details{
  The covariates are included in the model through \code{\link{model.matrix}}, so that factors are automatically coded by means of dummy variables and transformations of the covariates may be specified directly in the formula, as in \code{I(age - 50)} or \code{I((age - 50)^2/100)}. The model always includes an intercept and only one response variable is allowed. Covariates affecting the initial probabilities are specified through the separate argument \code{formula_init}; when this is \code{NULL}, the covariates in \code{formula} are used for both the initial and the transition probabilities, whereas \code{formula_init = y ~ 1} gives initial probabilities that do not depend on covariates.
}

\references{
  Bartolucci, F., Pandolfi, S., and Pennoni, F. (2026). Parsimonious parametrizations of transition matrices of Markov chain and hidden Markov models,  \emph{Annals of Operations Research}, 363, 233–266.
}

\author{Francesco Bartolucci and Silvia Pandolfi, University of Perugia (IT), Luca Brusa and Fulvia Pennoni, University of Milano-Bicocca (IT)}

\examples{
\dontrun{
data("data_SRHS_long")
srhs <- data_SRHS_long

# Covariates
srhs$gender <- factor(
  srhs$gender,
  levels = c(1, 2),
  labels = c("male", "female")
)

srhs$race <- factor(
  srhs$race,
  levels = c(1, 2, 3),
  labels = c("white", "black", "others")
)

srhs$education <- factor(
  srhs$education,
  levels = 1:5,
  labels = c(
    "high school",
    "general educational diploma",
    "high school graduate",
    "some college",
    "college and above"
  )
)

srhs$age_c <- srhs$age - mean(srhs$age)

# 1 = poor, ..., 5 = excellent

srhs$srhs1 <- 6L - as.integer(srhs$srhs)

# Consecutive numeric identifiers and ordering by subject and time

srhs$id <- as.integer(factor(srhs$id))
srhs <- srhs[order(srhs$id, srhs$t), ]

# Select a small reproducible balanced sample for a fast example

set.seed(2026)

n_example <- 2000L
example_ids <- sample(unique(srhs$id), n_example)

srhs_example <- srhs[srhs$id \%in\% example_ids, , drop = FALSE]

srhs_example$id <- match(srhs_example$id, example_ids)
srhs_example <- srhs_example[order(srhs_example$id, srhs_example$t), ]


# Example 1: gender and age affect both initial and transition probabilities

mod_cov <- mc_multinom(
  srhs1 ~ gender + age_c,
  data = srhs_example,
  baseline = "central",
  model_int = "dist1",
  model_cov = "dist1",
  output = TRUE)

print(mod_cov)
summary(mod_cov)


# Average initial and transition probabilities

colMeans(mod_cov$Piv)
round(apply(mod_cov$PI[, , , 2:max(srhs_example$t)], c(1, 2), mean),3)

# Example 2: covariates affect transitions only

mod_trans <- mc_multinom(
  srhs1 ~ gender + age_c,
  data = srhs_example,
  baseline = "central",
  model_int = "dist1",
  model_cov = "dist1",
  formula_init = srhs1 ~ 1)

print(mod_trans)
summary(mod_trans)
}
}

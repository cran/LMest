\name{design_matrices}

\alias{design_matrices}

\title{Build design matrices for \code{mc_multinom}}

\description{Function that builds the contrast and design matrices required by \code{\link{mc_multinom}} to parameterize, through multinomial logits, the initial and the transition probabilities of a Markov chain model.}

\usage{
design_matrices(k, nc, baseline, model_int, model_cov)
}

\arguments{
  \item{k }{Number of response categories (or observed states) in the Markov chain, coded as \code{1, ..., k}. Note that the structures \code{"dist1"}, \code{"dist2"}, \code{"dist3"}, \code{"two"} and \code{"diff"} (see \code{model_int}) rely on the distance or on the direction between categories, and are therefore meaningful only when the categories are ordered; the remaining structures do not require any ordering}
  \item{nc }{Number of covariates affecting the transition probabilities, excluding the intercept. Set \code{nc = 0} when no transition covariates are included}
  \item{baseline }{Character string specifying the reference category in the multinomial logit parameterization. With \code{"initial"}, category (1) is the reference category for both initial and transition probabilities. With \code{"central"}, category (1) remains the reference category for the initial probabilities, whereas the current state is the reference category for each row of the transition matrix}
  \item{model_int }{Character string specifying the constraint structure for the intercepts of the transition logits. Possible values are:
  \describe{
    \item{\code{"all"}}{unconstrained transition-specific intercepts, one for each of the \eqn{k(k-1)}{k(k-1)} transitions}
    \item{\code{"const"}}{one common intercept for all transitions}
    \item{\code{"dist1"}}{one parameter for each signed distance \eqn{v-u}{v-u} between the destination category \eqn{v} and the origin category \eqn{u} (\eqn{2(k-1)}{2(k-1)} parameters)}
    \item{\code{"dist2"}}{one parameter for each absolute distance \eqn{|v-u|}{|v-u|}, shared by the two directions (\eqn{k-1}{k-1} parameters)}
    \item{\code{"dist3"}}{as \code{"dist2"}, but the parameter enters with a positive sign for upward transitions and with a negative sign for downward ones (\eqn{k-1}{k-1} parameters)}
    \item{\code{"two"}}{separate parameters for upward and downward transitions}
    \item{\code{"symm"}}{one parameter for each pair of categories, shared by the two directions (\eqn{k(k-1)/2}{k(k-1)/2} parameters)}
    \item{\code{"rsymm"}}{as \code{"symm"}, but the parameter enters with opposite signs in the two directions}
    \item{\code{"diff"}}{effects expressed as differences between destination- and origin-specific parameters, with the parameter of the first category set to zero (\eqn{k-1}{k-1} parameters)}
  }}
  \item{model_cov }{Character string specifying the constraint structure for the regression coefficients of the transition covariates. The admissible values and their interpretation are the same as for \code{model_int}. This argument is used only when \code{nc > 0}}
}

\value{
  \item{G }{Contrast matrix used to parameterize the initial probabilities in the \code{\link{mc_multinom}} function}
  \item{Z }{Design matrix mapping the reduced parameter vector to the transition-logit parameters under the selected intercept and covariate constraints}
  \item{GG }{Array of contrast matrices used to parameterize the transition probabilities for each origin state}
  \item{Z1 }{Design matrix defining the constraint structure for the intercepts of the transition logits}
  \item{Z2 }{Design matrix defining the constraint structure for the covariate effects on the transition logits. This component is returned only when \code{nc > 0}}
}


\references{
  Bartolucci, F., Pandolfi, S., and Pennoni, F. (2026). Parsimonious parametrizations of transition matrices of Markov chain and hidden Markov models,  \emph{Annals of Operations Research}, 363, 233–266.
}

\author{Francesco Bartolucci, Silvia Pandolfi, University of Perugia (IT), Luca Brusa and Fulvia Pennoni, University of Milano-Bicocca (IT)}

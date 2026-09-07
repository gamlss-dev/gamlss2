\name{ALD}
\alias{ALD}
\alias{fit_ALD}

\title{Asymmetric Laplace Quantile Regression}

\description{
  Asymmetric Laplace family for quantile regression and a wrapper for
  estimating one or more conditional quantile models.
}

\usage{
ALD(tau = 0.5)

fit_ALD(formula, data, tau, anchor = 0.5, max_step = 0.05,
  restarts = c("always", "on_failure", "never"), start = NULL,
  monotone = FALSE, verbose = interactive(), ...)
}

\arguments{
  \item{tau}{For \code{ALD}, one finite quantile probability strictly between
    zero and one. For \code{fit_ALD}, a numeric vector of such probabilities.}
  \item{formula}{A model formula accepted by \code{\link{gamlss2}}.}
  \item{data}{A data frame containing the variables in \code{formula}.}
  \item{anchor}{Quantile probability from which the sequence of fits starts.}
  \item{max_step}{Maximum distance between successive quantiles along the
    fitting path. Intermediate bridge models are fitted when necessary.
    Use \code{Inf} to disable bridge models.}
  \item{restarts}{Character. With \code{"always"}, each warm-started fit is
    compared with a cold fit; \code{"on_failure"} uses a cold fit only when
    the warm fit fails; and \code{"never"} uses continuation only.}
  \item{start}{Optional starting values for the anchor model, as accepted by
    \code{\link{gamlss2}}.}
  \item{monotone}{Logical. If \code{TRUE}, returned fitted values and
    predictions are rearranged pointwise to be nondecreasing in \code{tau}.}
  \item{verbose}{Logical. Should progress information be printed?}
  \item{\dots}{Further arguments passed to \code{\link{gamlss2}}.}
}

\details{
  The asymmetric Laplace density is proportional to
  \eqn{\exp\{-\rho_\tau((y-\mu)/\sigma)\}}, where
  \eqn{\rho_\tau(u) = u\{\tau-I(u<0)\}}. Thus, \eqn{\mu} is the conditional
  \eqn{\tau}-quantile and \eqn{\sigma>0} is a scale parameter.

  \code{fit_ALD} first estimates the anchor quantile and then proceeds toward
  the lower and upper requested quantiles. Coefficients and smoothing
  parameters are transferred between successive fits. Candidate fits are
  compared using their asymmetric Laplace log-likelihood.

  Monotone rearrangement affects only the values returned by
  \code{fitted()} and \code{predict()}; the individual fitted models are
  retained unchanged.
}

\value{
  \code{ALD} returns an object of class \code{"gamlss2.family"}.

  \code{fit_ALD} returns an object of class \code{"ALDfit"}, containing the
  individual models, fitted values, requested quantiles, and fitting
  diagnostics.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link[stats]{predict}}.
}

\examples{\dontrun{## load data
data("SpirometryUS", package = "gamlss2")

## subset for females
d <- subset(SpirometryUS, gender == "female")

## fit one conditional median model
m50 <- gamlss2(fev1 ~ s(age), data = d, family = ALD(tau = 0.5))

## estimate several conditional quantiles
qu <- c(0.025, 0.5, 0.975)
m <- fit_ALD(fev1 ~ s(age), data = d, tau = qu,
  max_step = 0.25, monotone = TRUE, verbose = FALSE)

## plot fitted quantile curves
fit <- fitted(m)
i <- order(d$age)
plot(fev1 ~ age, data = d)
matlines(d$age[i], fit[i, ], lty = 1, lwd = 2,
  col = c(2, 4, 2))
}}

\keyword{models}
\keyword{regression}


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

\references{
  Yu K, Moyeed RA (2001). \dQuote{Bayesian Quantile Regression.}
  \emph{Statistics & Probability Letters}, \bold{54}(4), 437--447.
  \doi{10.1016/S0167-7152(01)00124-9}

  Chernozhukov V, Fernandez-Val I, Galichon A (2010).
  \dQuote{Quantile and Probability Curves Without Crossing.}
  \emph{Econometrica}, \bold{78}(3), 1093--1125.
  \doi{10.3982/ECTA7880}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link[stats]{predict}}, \code{\link{mcmc}}.
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
  max_step = 0.1, monotone = TRUE, verbose = FALSE)

## plot fitted quantile curves
fit <- fitted(m)
i <- order(d$age)
plot(fev1 ~ age, data = d)
matlines(d$age[i], fit[i, ], lty = 1, lwd = 2,
  col = c(2, 4, 3))

## estimate lower limits model using full MCMC
set.seed(1328)
m2.5 <- mcmc(m$models[["tau=0.025"]],
  n.iter = 12000, burnin = 2000, thin = 10)

summary(m2.5)

## predict credible intervals
nd <- data.frame(age = seq(6, 80, length = 300))
nd$p <- predict(m2.5, newdata = nd,
  FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

## plot
plot(fev1 ~ age, data = d)
matlines(nd$age, nd$p, lty = 1, lwd = 2, col = 2)
}}

\keyword{models}
\keyword{regression}


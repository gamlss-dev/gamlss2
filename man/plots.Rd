\name{plot.gamlss2}
\alias{plot.gamlss2}
\alias{plot.gamlss2.list}

\title{
  Plotting GAMLSS
}

\description{
  Plotting methods for objects of class \code{"gamlss2"} and \code{"gamlss2.list"}, which can be
  used for effect plots of model terms or residual plots.
  Note that effect plots of model terms with more than two covariates are not supported, for this
  purpose use the \code{\link[stats]{predict}} method.
}

\usage{
\method{plot}{gamlss2}(x, parameter = NULL,
  which = "effects", terms = NULL,
  scale = TRUE, spar = TRUE, ...)

\method{plot}{gamlss2.list}(x, parameter = NULL, which = "effects",
  terms = NULL, spar = TRUE, legend = TRUE, ...)
}

\arguments{
  \item{x}{An object of class \code{"gamlss2"} or \code{"gamlss2.list"}, which can be
    created by using the \code{c()} method combining \code{"gamlss2"} objects.
    See th examples.}
  \item{parameter}{Character or integer. For which parameter/model/what should the plots be created?
    Note that instead of argument \code{parameter} plots can also be specified passing
    argument \code{model} and \code{what} to \dots.}
  \item{which}{Character or integer, selects the type of plot: \code{"effects"} produces
    effect plots of (special) model terms,
    \code{"hist-resid"} shows a histogram of residuals, \code{"qq-resid"} shows
    a quantile-quantile plot of residuals, \code{"scatter-resid"} shows a scatter plot of
    residuals with fitted values for the distribution mean
    (or median, if available in the family object).}
  \item{terms}{Character or integer. For which model term should the plot(s) be created?}
  \item{scale}{If set to \code{1}, effect plots all have the same scale on the y-axis. If set to
    \code{0} each effect plot has its own scale for the y-axis.}
  \item{spar}{Should graphical parameters be set?}
  \item{legend}{Should a legend be added using multiple model plots?}
  \item{\dots}{Arguments such as
    \code{lwd}, \code{lty}, \code{col}, \code{legend = TRUE} (for multiple model plots),
    a.o., depending on the type of plot. See the examples.}
}

\seealso{
  \code{\link{gamlss2}}.
}

\examples{
\dontrun{## Load data.
data("film90", package = "gamlss.data")

## Model formula.
f <-  ~ s(lboopen) + s(lnosc)
f <- rep(list(f), 4)
f[[1]] <- update(f[[1]], lborev1 ~ .)

## Estimate model.
b1 <- gamlss2(f, data = film90, family = BCCG)

## Plot effects (default).
plot(b1)

## Plot specific effect.
plot(b1, parameter = "sigma")
plot(b1, model = "sigma")
plot(b1, model = "nu", term = 1)
plot(b1, model = "nu", term = 2)
plot(b1, model = "nu", term = "lnosc")
plot(b1, term = "lnosc")

## Plot all residual diagnostics.
plot(b1, which = "resid")

## Single diagnostic plots.
plot(b1, which = "hist-resid")
plot(b1, which = "qq-resid")
plot(b1, which = "wp-resid")
plot(b1, which = "scatter-resid")

## Estimate another model.
b2 <- gamlss2(f, data = film90, family = BCPE)

## Compare estimated effects.
plot(c(b1, b2))
plot(c(b1, b2), term = "lboopen",
  col = c(1, 4), lwd = 3, lty = 1,
  pos = c("topleft", "topright", "bottomleft", "bottomright"))
plot(c(b1, b2), model = "sigma")
plot(c(b1, b2), model = "sigma", term = 2)
plot(c(b1, b2), model = c("mu", "nu"))
}
}

\keyword{aplot}

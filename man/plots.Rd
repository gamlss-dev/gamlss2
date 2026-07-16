\name{plot.gamlss2}
\alias{plot.gamlss2}
\alias{plot.gamlss2.list}

\title{Plotting GAMLSS}

\description{
  Plotting methods for objects of class \code{"gamlss2"} and
  \code{"gamlss2.list"}. These methods support effect plots for model
  terms and residual diagnostics. Effect plots for model terms with more
  than two covariates are not supported; in such cases, use
  \code{\link[stats]{predict}} and plot the results manually.
}

\usage{
\method{plot}{gamlss2}(x, parameter = NULL,
  which = "effects", terms = NULL,
  scale = TRUE, spar = TRUE, ...)

\method{plot}{gamlss2.list}(x, parameter = NULL, which = "effects",
  terms = NULL, spar = TRUE, legend = TRUE, ...)
}

\arguments{
  \item{x}{An object of class \code{"gamlss2"} or \code{"gamlss2.list"}. Objects of class
    \code{"gamlss2.list"} can be created by combining \code{"gamlss2"} objects with
    \code{c()}.
    See the examples.}
  \item{parameter}{Character or integer specifying for which distribution parameter the plots should be created.
    The arguments \code{model} and \code{what}, if supplied through \code{\dots}, are treated as
    aliases for \code{parameter} to support alternative calling styles.}
  \item{which}{Character or integer selecting the type of plot.
    \code{"effects"} produces effect plots of model terms,
    \code{"resid"} shows the available residual diagnostics,
    \code{"hist-resid"} produces a histogram of residuals,
    \code{"qq-resid"} a quantile-quantile plot of residuals,
    \code{"wp-resid"} a worm plot of residuals, and
    \code{"scatter-resid"} a scatter plot of residuals against fitted values
    for the distribution mean (or median, if available in the family object).
    In addition, \code{"selection"} plots the GAIC path stored for a fitted
    stepwise or selection object, and \code{"samples"} visualizes posterior or
    simulation draws for objects that contain \code{samples}.}
  \item{terms}{Character or integer. Which model term(s) should be plotted?}
  \item{scale}{Logical. If \code{TRUE}, effect plots all have the same scale on the y-axis.
    If \code{FALSE}, each effect plot has its own y-axis scale.}
  \item{spar}{Logical. Should graphical parameters be set automatically?}
  \item{legend}{Logical. Should a legend be added for plots with multiple models?}
  \item{\dots}{Additional graphical arguments such as
    \code{lwd}, \code{lty}, \code{col}, \code{legend = TRUE} (for multiple model plots),
    among others, depending on the type of plot. Arguments \code{model} and
    \code{what} are treated as aliases for \code{parameter}.}
}

\details{
  For objects of class \code{"gamlss2"}, \code{which = "effects"} is the default
  and plots the partial effects stored in \code{results(object)}. If no effect
  information is available, residual diagnostics are shown instead.

  For objects of class \code{"gamlss2.list"}, plots compare effect estimates
  across several fitted models. Residual diagnostics are only available for
  single-model objects.
}

\seealso{
  \code{\link{gamlss2}}.
}

\examples{\dontrun{## model formula
f <- Ozone ~ s(Temp) + s(Wind) + s(Solar.R)| . | 1 | 1

## estimate model
m1 <- gamlss2(f, data = airquality, family = BCCG)

## plot effects (default)
plot(m1)

## plot specific effect
plot(m1, parameter = "sigma")
plot(m1, model = "sigma")
plot(m1, model = "mu", terms = 1)
plot(m1, model = "mu", terms = 2)
plot(m1, model = "mu", terms = "Solar")
plot(m1, terms = "Solar.R")

## plot all residual diagnostics
plot(m1, which = "resid")

## single diagnostic plots
plot(m1, which = "hist-resid")
plot(m1, which = "qq-resid")
plot(m1, which = "wp-resid")
plot(m1, which = "scatter-resid")

## if available, plot stored stepwise-selection path
## plot(m1, which = "selection")

## estimate another model
m2 <- gamlss2(f, data = airquality, family = BCPE)

## compare estimated effects
plot(c(m1, m2))
plot(c(m1, m2), terms = "Wind",
  col = c(1, 4), lwd = 3, lty = 1,
  pos = c("topright", "bottomleft"))
plot(c(m1, m2), model = "sigma")
plot(c(m1, m2), model = "sigma", terms = 2)
plot(c(m1, m2), model = c("mu", "sigma"))
}
}

\keyword{aplot}

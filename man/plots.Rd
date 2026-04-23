\name{plot.gamlss2}
\alias{plot.gamlss2}
\alias{plot.gamlss2.list}

\title{Plotting GAMLSS}

\description{
  Plotting methods for objects of class \code{"gamlss2"} and
  \code{"gamlss2.list"}. These methods can be used for effect plots of model
  terms and for residual diagnostics. Effect plots for model terms with more
  than two covariates are not supported; in that case, use
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
  \item{x}{An object of class \code{"gamlss2"} or \code{"gamlss2.list"}, which can be
    created by using the \code{c()} method combining \code{"gamlss2"} objects.
    See the examples.}
  \item{parameter}{Character or integer. For which distribution parameter should the plots be created?
    The arguments \code{model} and \code{what}, if supplied through \code{\dots}, are treated as
    aliases for \code{parameter} to support alternative calling styles.}
  \item{which}{Character or integer selecting the type of plot.
    \code{"effects"} produces effect plots of model terms,
    \code{"hist-resid"} a histogram of residuals,
    \code{"qq-resid"} a quantile-quantile plot of residuals, and
    \code{"scatter-resid"} a scatter plot of residuals against fitted values
    for the distribution mean (or median, if available in the family object).}
  \item{terms}{Character or integer. Which model term(s) should be plotted?}
  \item{scale}{Logical. If \code{TRUE}, effect plots all have the same scale on the y-axis.
    If \code{FALSE}, each effect plot has its own y-axis scale.}
  \item{spar}{Logical. Should graphical parameters be set automatically?}
  \item{legend}{Logical. Should a legend be added for plots with multiple models?}
  \item{\dots}{Arguments such as
    \code{lwd}, \code{lty}, \code{col}, \code{legend = TRUE} (for multiple model plots),
    a.o., depending on the type of plot. See the examples.}
}

\seealso{
  \code{\link{gamlss2}}.
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
\dontrun{## load data
data("film90", package = "gamlss.data")

## model formula
f <-  ~ s(lboopen) + s(lnosc)
f <- rep(list(f), 4)
f[[1]] <- update(f[[1]], lborev1 ~ .)

## estimate model
b1 <- gamlss2(f, data = film90, family = BCCG)

## plot effects (default)
plot(b1)

## plot specific effect
plot(b1, parameter = "sigma")
plot(b1, model = "sigma")
plot(b1, model = "nu", term = 1)
plot(b1, model = "nu", term = 2)
plot(b1, model = "nu", term = "lnosc")
plot(b1, term = "lnosc")

## plot all residual diagnostics
plot(b1, which = "resid")

## single diagnostic plots
plot(b1, which = "hist-resid")
plot(b1, which = "qq-resid")
plot(b1, which = "wp-resid")
plot(b1, which = "scatter-resid")

## estimate another model
b2 <- gamlss2(f, data = film90, family = BCPE)

## compare estimated effects
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

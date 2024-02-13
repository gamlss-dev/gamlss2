\name{quantile.gamlss2}
\alias{quantile.gamlss2}

\title{
  Quantiles for GAMLSS
}

\description{
  The function computes estimated quantiles and optionally produces a plot.
}

\usage{
\method{quantile}{gamlss2}(x, probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
  variable = NULL, newdata = NULL,
  plot = FALSE, data = TRUE,
  n = 100L, ...)
}

\details{
  The function applies the predict method to determine the parameters of the response distribution.
  It then computes the quantiles as specified in the argument \code{probs}.
}

\value{
  A data frame of the estimated quantiles.
}

\arguments{
  \item{x}{An object of class \code{"gamlss2"}.}
  \item{probs}{Numeric vector of probabilities with values in [0,1].}
  \item{variable}{Logical or integer, should quantiles be plotted using
    the covariate data? Note that the variable option is only possible for
    single covariate models.}
  \item{newdata}{Data frame that should be used for computing the quantiles.}
  \item{plot}{Logical, should a plot be shown?}
  \item{data}{Logical, should the raw data be added to the plot?}
  \item{n}{Integer, number of observations that should be used to compute
    an equidistant grid for the selected \code{variable}.}
  \item{\dots}{Arguments such as
    \code{col}, \code{legend = TRUE/FALSE}. See the examples.}
}

\seealso{
  \code{\link{gamlss2}}.
}

\examples{
\dontrun{## Load data.
data("film90", package = "gamlss.data")

## Model formula.
f <-  ~ s(lboopen)
f <- rep(list(f), 4)
f[[1]] <- update(f[[1]], lborev1 ~ .)

## Estimate model.
b <- gamlss2(f, data = film90, family = BCPE)

## Compute quantiles using "newdata".
nd <- film90[1:10, ]
print(quantile(b, newdata = nd))

## Plot sorted quantiles.
quantile(b, plot = TRUE)

## Quantile plot using covariate data.
quantile(b, plot = TRUE, variable = TRUE)

## Plot without raw data.
quantile(b, plot = TRUE, variable = TRUE, data = FALSE)
}
}

\keyword{aplot}

\name{quantile.gamlss2}
\alias{quantile.gamlss2}

\title{Quantiles for GAMLSS}

\description{
  The function computes estimated quantiles and optionally produces a plot.
}

\usage{
\method{quantile}{gamlss2}(x, probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
  variable = NULL, newdata = NULL,
  plot = FALSE, data = TRUE,
  n = 100L, ...)
}

\value{
  A data frame with one row per prediction and one column per probability.
  With \code{variable = NULL}, a single probability returns a numeric vector.
  When a covariate is selected, it is included as the first column.
  With \code{plot = TRUE}, the result is returned invisibly.
}

\details{
  Quantiles are computed from the predicted response distributions.

  With \code{plot = TRUE} and no selected \code{variable}, the x-axis shows
  observation ranks as percentages, ordered by fitted median; the y-axis is
  on the response scale. Lines show the fitted conditional quantiles at
  \code{probs}, and points show observed responses when \code{data = TRUE}.
  The same ordering is used for the responses and all quantiles, so each
  position refers to one observation. If \code{newdata} does not contain the
  response, only the fitted quantile curves are shown.
}

\arguments{
  \item{x}{An object of class \code{"gamlss2"}.}
  \item{probs}{Numeric vector of probabilities with values in [0,1].}
  \item{variable}{Covariate selector for single-covariate models: \code{TRUE}
    uses the sole covariate; an integer or character selects by column position
    or name. \code{NULL} (default) or \code{FALSE} computes quantiles at the
    supplied or fitted observations without a covariate grid.}
  \item{newdata}{Data frame that should be used for computing the quantiles.}
  \item{plot}{Logical, should an ordered quantile plot or a covariate plot be shown?}
  \item{data}{Logical, should observed responses be added to the plot when available?}
  \item{n}{Integer, number of observations that should be used to compute
    an equidistant grid for the selected \code{variable}.}
  \item{\dots}{Arguments such as
    \code{col}, \code{legend = TRUE/FALSE}. See the examples.}
}

\seealso{
  \code{\link{gamlss2}}.
}

\examples{\dontrun{## load data
data("SpirometryUS", package = "gamlss2")

## subset for female
d <- subset(SpirometryUS, gender == "female")

## model formula
f <- fev1 ~ s(age) | . | . | .

## estimate model
m <- gamlss2(f, data = d, family = BCPE)

## compute quantiles using "newdata"
nd <- d[1:10, ]
print(quantile(m, newdata = nd))

## plot quantiles ordered by the fitted median
quantile(m, plot = TRUE)

## quantile plot using covariate data
quantile(m, plot = TRUE,
  variable = TRUE, pos = "topright")

## plot without raw data
quantile(m, plot = TRUE,
  variable = TRUE, pos = "topright", data = FALSE)
}
}

\keyword{aplot}
\keyword{regression}


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

\details{
  The function applies the predict method to determine the parameters of the response distribution.
  It then computes the quantiles as specified in the argument \code{probs}.
}

\value{
  The estimated quantiles. For multiple probabilities, a data frame is returned.
  For a single probability and \code{variable = NULL}, the result is simplified
  to a numeric vector.
}

\arguments{
  \item{x}{An object of class \code{"gamlss2"}.}
  \item{probs}{Numeric vector of probabilities with values in [0,1].}
  \item{variable}{Logical, integer, or character. Should quantiles be plotted
    against a covariate? If \code{TRUE}, the single model covariate is used.
    Alternatively, a covariate can be selected by position or name. This option
    is only available for single-covariate models.}
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

## plot sorted quantiles
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

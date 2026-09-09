\name{marginal_predict}
\alias{marginal_predict}

\title{Marginal Prediction Profiles for gamlss2 Models}

\description{
Compute one-variable-at-a-time prediction profiles from a fitted
\pkg{gamlss2} model. Each selected covariate is varied over a grid while the
other model covariates are held at representative values.
}

\usage{
marginal_predict(object, newdata = NULL, variables = NULL,
  n = 100L, continuous = median, at = NULL, values = NULL, \dots)
}

\arguments{
  \item{object}{A fitted object inheriting from class \code{"gamlss2"} or
    \code{"bamlss2"}.}
  \item{newdata}{An optional data frame defining the observed ranges, levels,
    and representative values of the model covariates. By default, the model
    frame of \code{object} is used.}
  \item{variables}{A character vector naming the covariates for which profiles
    should be computed. By default, all covariates used by the model are
    included.}
  \item{n}{A single integer greater than one giving the grid size for numeric,
    \code{Date}, and \code{POSIXct} covariates. A covariate that is constant
    has a one-point grid.}
  \item{continuous}{A function, or the name of a function, used to compute the
    representative value of covariates not currently being varied. The default
    is \code{median}; \code{mean} is another common choice. The function must
    return one value.}
  \item{at}{An optional named list specifying values at which covariates should
    be held fixed. Each element must contain one value.}
  \item{values}{An optional named list specifying the prediction grid for one
    or more of the selected \code{variables}.}
  \item{\dots}{Further arguments passed to \code{\link{predict.gamlss2}}, such
    as \code{type}, \code{model}, \code{terms}, or \code{se.fit}.}
}

\details{
For each selected covariate, numeric and date-time covariates are varied over
an equally spaced grid spanning their observed range. Factors are varied over
their levels, and character and logical covariates over their observed values.
The \code{values} argument can be used to replace any of these grids.

Covariates not being varied are fixed at \code{at}, when supplied. Otherwise,
numeric and date-time covariates are summarized by \code{continuous}, factors
are set to their first level, character variables to their first observed
value, and logical variables to \code{FALSE} when that value is observed.

Here, \emph{marginal prediction} means a conditional one-variable profile:
predictions are evaluated at one representative combination of the remaining
covariates. They are not averages over the empirical joint distribution of
those covariates. In models with interactions,
the profile is conditional on the interacting covariates' fixed values; use
\code{at} to choose scientifically meaningful values.
}

\value{
A named list with one data frame per selected covariate. Each data frame
contains the focal covariate followed by the output from
\code{predict(object, newdata = ..., \dots)}. Vector predictions, such as those
from \code{type = "response"}, are stored in a column named \code{.fit}.
The representative values computed for the covariates are available in the
\code{"fixed"} attribute of the returned list. Each is used whenever that
covariate is not the focal variable.
}

\seealso{
\code{\link{predict.gamlss2}}, \code{\link{prodist.gamlss2}},
\code{\link{quantile.gamlss2}}
}

\examples{
## Same air-quality model as in the first gamlss2() example.
air <- subset(airquality, !is.na(Ozone))
m <- gamlss2(
  Ozone ~ s(Temp) + s(Wind) | s(Temp),
  data = air,
  family = GA
)

## Conditional mean profiles for temperature and wind.
mp <- marginal_predict(
  m,
  variables = c("Temp", "Wind"),
  n = 50,
  type = "response"
)
head(mp$Temp)
attr(mp, "fixed")

plot(.fit ~ Temp, data = mp$Temp, type = "l",
  xlab = "Temperature (degrees F)", ylab = "Fitted mean ozone")

## Supply a custom grid and condition on a particular wind speed.
mp2 <- marginal_predict(
  m,
  variables = "Temp",
  values = list(Temp = seq(60, 95, by = 5)),
  at = list(Wind = 5),
  type = "response"
)
mp2$Temp
}

\keyword{regression}
\keyword{models}

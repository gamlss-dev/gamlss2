\name{MCS}
\alias{MCS}

\title{Motorcycle Counts at Sonnenberg}

\description{
  This dataset contains daily motorcycle counts at Sonnenberg, located in the Harz region
  of Germany. It covers a period of nearly three years, from 2021-01-01 to 2023-11-30.
}

\usage{data("MCS")}

\format{
The \code{MCS} data contains the following variables:
\describe{
  \item{date:}{Date, the date of the recording.}
  \item{counts:}{Integer, the number of motorcycles on that day.}
  \item{tempmin:}{Numeric, minimum temperature in \eqn{^{\circ}C}.}
  \item{tempmax:}{Numeric, maximum temperature in \eqn{^{\circ}C}.}
  \item{temp:}{Numeric, mean temperature in \eqn{^{\circ}C}.}
  \item{humidity:}{Numeric, mean relative humidity in percent.}
  \item{tempdew:}{Numeric, average dewpoint temperature in \eqn{^{\circ}C}.}
  \item{cloudiness:}{Numeric, average cloud cover in percent.}
  \item{rain:}{Numeric, amount of precipitation in mm (snow and rain).}
  \item{sunshine:}{Numeric, sunshine duration in minutes.}
  \item{wind:}{Numeric, mean wind speed in m/s.}
  \item{windmax:}{Numeric, maximum wind speed in m/s.}
}
}

\examples{
data("MCS")

## Visualize seasonal trends.
## First, extract the day of the year.
MCS$yday <- as.integer(format(MCS$date, "\%j")) - 1L

## Plot the counts over the motorcycle season.
plot(counts ~ yday, data = MCS)

## Visualize count distribution.
barplot(table(MCS$counts))

## Estimation of a negative binomial seasonal model
## for motorcycle counts using cyclic splines.
## Model formula:
f <- counts ~ s(yday, bs = "cc") |
  s(yday, bs = "cc")

## Estimate count model.
b <- gamlss2(f, data = MCS, family = NBI)

## Plot estimated effects.
plot(b)

## Residual diagnostics.
plot(b, which = "resid")

## Visualize estimated quantiles.
## (1) Predict parameters.
nd <- data.frame(yday = 1:365)
par <- predict(b, newdata = nd)

## (2) Compute quantiles.
p <- sapply(c(0.05, 0.5, 0.95), function(q) family(b)$q(q, par))

## Plot.
plot(counts ~ yday, data = MCS,
  pch = 19, col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
matplot(nd$yday, p, type = "l", lty = c(2, 1, 2),
  lwd = 2, col = 4, add = TRUE)
}

\keyword{datasets}


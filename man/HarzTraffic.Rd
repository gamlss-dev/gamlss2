\name{HarzTraffic}
\alias{HarzTraffic}

\encoding{UTF-8}

\title{Traffic Counts at Sonnenberg in the Harz Region}

\description{
  This dataset contains daily tarrfic counts close to the Sonnenberg, located in the Harz region
  of Germany. It covers a period of nearly three years, from 2021-01-01 to 2023-11-30.
}

\usage{data("MCS")}

\format{
The \code{MCS} data contains the following variables:
\describe{
  \item{date:}{Date, the date of the recording.}
  \item{yday:}{Integer, the day of the year.}
  \item{bikes:}{Integer, the number of motorcycles on that day.}
  \item{cars:}{Integer, the number of cars on that day.}
  \item{trucks:}{Integer, the number of trucks on that day.}
  \item{others:}{Integer, the number of other vehicles on that day.}
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

\references{
  Weather Data:
  \describe{
    \item{Data Source:}{Deutscher Wetterdienst (DWD), Climate Data Center (CDC).}
    \item{Licence:}{CC BY 4.0}
    \item{URL:}{\url{https://opendata.dwd.de/climate_environment/CDC/}}
    \item{Station:}{Wernigerode (5490; Sachsen-Anhalt)}
    \item{Position:}{10.7686/51.8454/233 (lon, lat, alt, EPSG 4326)}
  }

  Traffic Data:
  \describe{
    \item{Data Source:}{Bundesanstalt für Strassenwesen (BASt)}
    \item{Licence:}{CC BY 4.0}
    \item{URL:}{\url{https://www.bast.de},
      \url{https://www.bast.de/DE/Verkehrstechnik/Fachthemen/v2-verkehrszaehlung/Verkehrszaehlung.html}}
  }
}

\examples{
data("HarzTraffic")

## Plot bike counts over the season.
plot(counts ~ yday, data = HarzTraffic)

## Visualize count distribution.
barplot(table(HarzTraffic$counts))

## Estimation of a negative binomial seasonal model
## for motorcycle counts using cyclic splines.
## Model formula:
f <- counts ~ s(yday, bs = "cc") |
  s(yday, bs = "cc")

\dontrun{## Estimate count model.
b <- gamlss2(f, data = HarzTraffic, family = NBI)

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
plot(counts ~ yday, data = HarzTraffic,
  pch = 19, col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
matplot(nd$yday, p, type = "l", lty = c(2, 1, 2),
  lwd = 2, col = 4, add = TRUE)
}}

\keyword{datasets}


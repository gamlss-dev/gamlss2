\name{HarzTraffic}
\alias{HarzTraffic}

\encoding{UTF-8}

\title{Traffic Counts at Sonnenberg in the Harz Region}

\description{
  This dataset contains daily traffic counts near Sonnenberg in the Harz region
  of Germany. It covers nearly three years, from 2021-01-01 to 2023-11-30, and
  combines traffic counts with weather covariates. The data are useful for
  illustrating seasonal count regression models, for example with cyclic
  smoothers over the day of the year.
}

\usage{data("HarzTraffic", package = "gamlss2")}

\format{
A data frame containing 1057 observations on 16 variables.
\describe{
  \item{date}{Date, the date of the record.}
  \item{yday}{Integer, the day of the year extracted from \code{date}. This is
    convenient for modeling within-year seasonality.}
  \item{bikes}{Integer, the number of motorcycles on that day.}
  \item{cars}{Integer, the number of cars on that day.}
  \item{trucks}{Integer, the number of trucks on that day.}
  \item{others}{Integer, the number of other vehicles on that day.}
  \item{tempmin}{Numeric, minimum temperature in \eqn{^{\circ}C}.}
  \item{tempmax}{Numeric, maximum temperature in \eqn{^{\circ}C}.}
  \item{temp}{Numeric, mean temperature in \eqn{^{\circ}C}.}
  \item{humidity}{Numeric, mean relative humidity in percent.}
  \item{tempdew}{Numeric, average dewpoint temperature in \eqn{^{\circ}C}.}
  \item{cloudiness}{Numeric, average cloud cover in percent.}
  \item{rain}{Numeric, amount of precipitation in mm (snow and rain).}
  \item{sunshine}{Numeric, sunshine duration in minutes.}
  \item{wind}{Numeric, mean wind speed in m/s.}
  \item{windmax}{Numeric, maximum wind speed in m/s.}
}
}

\details{
  The traffic variables record daily counts by vehicle type. The weather
  variables summarize daily conditions measured at a nearby station. In package
  examples, \code{bikes} is often modeled as a seasonal response because
  motorcycle traffic shows a pronounced annual cycle and substantial variation
  across weather conditions.
}

\source{
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
    \item{URL:}{\url{https://www.bast.de}}
  }
}

\examples{
## seasonal variation of motorcycle counts at Sonnenberg/Harz
data("HarzTraffic", package = "gamlss2")
plot(bikes ~ yday, data = HarzTraffic)

## count distribution
barplot(table(HarzTraffic$bikes))

## negative binomial seasonal model using cyclic splines
m <- gamlss2(bikes ~ s(yday, bs = "cc") | s(yday, bs = "cc"),
  data = HarzTraffic, family = NBI)

## visualize effects
plot(m)

## residual diagnostics
plot(m, which = "resid")

## fitted parameters for each day of the year
nd <- data.frame(yday = 1:365)

## corresponding quantiles
p <- quantile(m, newdata = nd, probs = c(0.05, 0.5, 0.95))

## visualization
plot(bikes ~ yday, data = HarzTraffic, pch = 19, col = gray(0.1, alpha = 0.3))
matlines(nd$yday, p, lty = c(2, 1, 2), lwd = 2, col = 4)
}

\keyword{datasets}


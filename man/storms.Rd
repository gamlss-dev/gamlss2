\name{storms}
\alias{storms}

\encoding{UTF-8}

\title{Severe Storms in Germany}

\description{
  According to the Beaufort scale, severe storms occur from a wind speed of 24.5-28.4 m/s.
  This dataset contains annual severe storm counts from weather stations in Germany
  from 1981 to 2021.
}

\usage{data("storms", package = "gamlss2")}

\format{
A data frame containing 3494 observations on 8 variables.
\describe{
  \item{id}{Factor, the weather station id.}
  \item{county}{Character, the county in Germany where the weather station is located.}
  \item{state}{Character, the state in Germany where the weather station is located.}
  \item{year}{Integer, the year the observation was measured.}
  \item{counts}{Integer, the number of severe storms in this year.}
  \item{alt}{Numeric, the altitude in meters above sea level of the weather station.}
  \item{lon}{Numeric, the longitude coordinate of the weather station.}
  \item{lat}{Numeric, the latitude coordinate of the weather station.}
}
}

\source{
  Severe Storms Data:
  \describe{
    \item{Data Source:}{Deutscher Wetterdienst (DWD), Climate Data Center (CDC).}
    \item{Licence:}{CC BY 4.0}
    \item{URL:}{\url{https://opendata.dwd.de/climate_environment/CDC/}}
    \item{Coordinate Reference System:}{EPSG 4326}
  }
}

\examples{
## load the data
data("storms", package = "gamlss2")

## yearly observations
plot(counts ~ year, data = storms)

## count distribution
barplot(table(storms$counts))

## NBI model
## model formula including spatial effect
f <- counts ~ s(year) + s(alt) + s(lon, lat) |
  s(year) + s(alt) + s(lon, lat)

## estimate model
b <- gamlss2(f, data = storms, family = NBI)

## estimated effects
plot(b)

## residual diagnostics
plot(b, which = "resid")
}

\keyword{datasets}


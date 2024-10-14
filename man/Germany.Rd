\name{Germany}
\alias{Germany}

\encoding{UTF-8}

\title{Map of Germany}

\description{
  The map contains the counties of Germany. The data was originally taken from
  GADM (\url{https://gadm.org/}) and slightly simplified to reduce disk space.
}

\usage{data("Germany", package = "gamlss2")}

\format{
A class \code{"sf"} data frame containing 403 counties of Germany.
\describe{
  \item{id}{Factor, a county id.}
  \item{county}{Character, the county in Germany where the weather station is located.}
  \item{state}{Character, the state in Germany where the weather station is located.}
  \item{geometry}{The polygon information.}
}
}

\source{
  Map of Germany:
  \describe{
    \item{Data Source:}{GADM}
    \item{Licence:}{CC BY}
    \item{URL:}{\url{https://gadm.org/}}
    \item{Coordinate Reference System:}{Longitude/latitude and the WGS84 datum.}
  }
}

\examples{
\dontshow{ if(!requireNamespace("sf")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
## load sf package for plotting
library("sf")

## load the data
data("Germany", package = "gamlss2")

## plot the map
plot(st_geometry(Germany))
}

\keyword{datasets}


\name{response.name}
\alias{response.name}

\title{
  Auxiliary Functions.
}

\description{
  Various auxiliary functions to facilitate the work of fitted models.
}

\usage{
response.name(formula)
}

\arguments{
  \item{formula}{A \code{\link[stats]{formula}}, \code{\link[Formula]{Formula}},
    or a fitted model object.}
}

\value{
  Function \code{response.name} extracts the response name as a character vector.
}

\seealso{
  \code{\link{gamlss2}}
}

\examples{
## Required package.
library("Formula")

## Simple formula
f <- y ~ x1 + x2 + log(x3)
rn <- response.name(f)
print(rn)

## Multiple responses.
f <- y1|y2|y3 ~ x1 + s(x2) + x3 + te(log(x3), x4) | x2 + ti(x5)
rn <- response.name(f)
print(rn)

f <- list(
  y1 ~ x1 + s(x2) + x3 + te(log(x3), x4),
  y2  ~ x2 + sqrt(x5),
  y3  ~ z2 + x1 + exp(x3) + s(x10)
)
rn <- response.name(f)
print(rn)
}

\keyword{utilities}
\keyword{regression}

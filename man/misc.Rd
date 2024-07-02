\name{response_name}
\alias{response_name}

\title{Auxiliary Functions for Formulas and Model Objects}

\description{
  Various auxiliary functions to facilitate the work with formulas and fitted
  model objects.
}

\usage{
response_name(formula)
}

\arguments{
  \item{formula}{A \code{\link[stats]{formula}}, \code{\link[Formula]{Formula}},
    or a fitted model object.}
}

\value{
  Function \code{response_name} extracts the response name as a character vector.
}

\seealso{
  \code{\link{gamlss2}}
}

\examples{
## basic formula
f <- y ~ x1 + x2 + log(x3)
response_name(f)

## formula with multiple responses
f <- y1 | y2 | y3 ~ x1 + s(x2) + x3 + te(log(x3), x4) | x2 + ti(x5)
response_name(f)

## list of formulas
f <- list(
  y1 ~ x1 + s(x2) + x3 + te(log(x3), x4),
  y2  ~ x2 + sqrt(x5),
  y3  ~ z2 + x1 + exp(x3) + s(x10)
)
response_name(f)
}

\keyword{utilities}
\keyword{regression}


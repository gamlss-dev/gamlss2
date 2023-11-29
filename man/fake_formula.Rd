\name{fake_formula}
\alias{fake_formula}

\title{
  Create a "Fake Formula".
}

\description{
  Create a "fake formula" from a \code{\link[stats]{formula}}, a \code{\link[Formula]{Formula}},
  or a \code{\link[base]{list}} of \code{\link[stats]{formula}}s. The function extracts
  all necessary variables (transformation of variables), to build a \code{\link[stats]{model.frame}}.
  The function also extracts all special model terms within the formulas, the information can
  be used to setup any special model term specification list.
}

\usage{
fake_formula(formula, specials = NULL,
  nospecials = FALSE, onlyspecials = FALSE)
}

\arguments{
  \item{formula}{A \code{\link[stats]{formula}}, \code{\link[Formula]{Formula}},
    or a \code{\link[base]{list}} of \code{\link[stats]{formula}}s.}
  \item{specials}{Character, vector of names of special functions in the formula,
    see \code{\link[stats]{terms.formula}}.}
  \item{nospecials}{Logical, should variables of special model terms be part
    of the "fake formula"?}
  \item{onlyspecials}{Logical, should only the special model terms be returned?}
}

\value{
  Depending on the input formula,
  the function returns a \code{\link[stats]{formula}} or \code{\link[Formula]{Formula}}.
  If \code{onlyspecials = TRUE} a vector or list of special model term names is returned.
}

\seealso{
  \code{\link{gamlss2}}
}

\examples{
## Required package.
library("Formula")

## Simple formula
f <- y ~ x1 + x2 + log(x3)
ff <- fake_formula(f)
print(ff)

## Including special model terms.
f <- y ~ x1 + s(x2) + x3 + te(log(x3), x4)
ff <- fake_formula(f)
print(ff)

## Multiple formulas.
f <- y ~ x1 + s(x2) + x3 + te(log(x3), x4) | x2 + sqrt(x5)
ff <- fake_formula(f)
print(ff)

## Can be collapsed.
print(formula(ff, collapse = TRUE))
print(formula(ff, collapse = TRUE, update = TRUE))

## A list of formulas.
f <- list(
  y ~ x1 + s(x2) + x3 + te(log(x3), x4),
    ~ x2 + sqrt(x5),
    ~ z2 + x1 + exp(x3)
)
ff <- fake_formula(f)
print(ff)

## Extract single parts of rhs.
formula(ff, rhs = 1)
formula(ff, rhs = 2)
formula(ff, rhs = 3)

## Multiple responses.
f <- y1|y2|y3 ~ x1 + s(x2) + x3 + te(log(x3), x4) | x2 + ti(x5)
ff <- fake_formula(f)
print(ff)

f <- list(
  y1 ~ x1 + s(x2) + x3 + te(log(x3), x4),
  y2  ~ x2 + sqrt(x5),
  y3  ~ z2 + x1 + exp(x3) + s(x10)
)
ff <- fake_formula(f)

## Extract only no special variables.
ff <- fake_formula(f, nospecials = TRUE)
print(ff)

## Extract only special terms.
ff <- fake_formula(f, onlyspecials = TRUE)
print(ff)
}

\keyword{utilities}
\keyword{regression}


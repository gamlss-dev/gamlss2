\name{fake_formula}
\alias{fake_formula}

\title{Extended Processing of "Fake" Formulas}

\description{
  Create a simplified formula from a \code{\link[stats]{formula}}, a
  \code{\link[Formula]{Formula}}, or a \code{\link[base]{list}} of
  \code{\link[stats]{formula}}s. The result retains the variables and
  transformations needed to build a \code{\link[stats]{model.frame}} while
  separating out special model terms. This is mainly used internally when
  \pkg{gamlss2} prepares formulas for fitting.
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

\details{
  The function is primarily an internal helper for formula preprocessing in
  \pkg{gamlss2}, but it can also be useful for inspecting how complex formulas
  are decomposed before fitting.

  Special model terms such as \code{s()}, \code{te()}, \code{ti()}, \code{re()},
  or package-specific specials are identified and handled separately from
  ordinary variables. Depending on the control arguments:
  \itemize{
    \item if \code{nospecials = FALSE}, variables used inside registered special
      terms are reintroduced into the returned formula so they remain available
      in the model frame;
    \item if \code{nospecials = TRUE}, the returned formula excludes the special
      terms and their variables;
    \item if \code{onlyspecials = TRUE}, the function returns the names of the
      detected special terms instead of a formula object.
  }

  For multi-part or multi-parameter formulas, the output preserves the overall
  formula structure.
}

\value{
  Depending on the input, the function returns a
  \code{\link[stats]{formula}} or \code{\link[Formula]{Formula}}. If
  \code{onlyspecials = TRUE}, a character vector or list of special model term
  names is returned.
}

\seealso{
  \code{\link{gamlss2}}
}

\examples{
## basic formula
f <- y ~ x1 + x2 + log(x3)
ff <- fake_formula(f)
print(ff)

## including special model terms
f <- y ~ x1 + s(x2) + x3 + te(log(x3), x4)
ff <- fake_formula(f)
print(ff)

## multiple parts on the right-hand side
f <- y ~ x1 + s(x2) + x3 + te(log(x3), x4) | x2 + sqrt(x5)
ff <- fake_formula(f)
print(ff)

## collapse all formula parts
print(formula(ff, collapse = TRUE))
print(formula(ff, collapse = TRUE, update = TRUE))

## list of formulas
f <- list(
  y ~ x1 + s(x2) + x3 + te(log(x3), x4),
    ~ x2 + sqrt(x5),
    ~ z2 + x1 + exp(x3)
)
ff <- fake_formula(f)
print(ff)

## extract separate parts on the right-hand side
formula(ff, rhs = 1)
formula(ff, rhs = 2)
formula(ff, rhs = 3)

## formula with multiple responses and multiple parts
f <- y1 | y2 | y3 ~ x1 + s(x2) + x3 + te(log(x3), x4) | x2 + ti(x5)
ff <- fake_formula(f)
print(ff)

## list of formulas with multiple responses
f <- list(
  y1 ~ x1 + s(x2) + x3 + te(log(x3), x4),
  y2 ~ x2 + sqrt(x5),
  y3 ~ z2 + x1 + exp(x3) + s(x10)
)
ff <- fake_formula(f)

## extract only without special terms
ff <- fake_formula(f, nospecials = TRUE)
print(ff)

## extract only special terms
ff <- fake_formula(f, onlyspecials = TRUE)
print(ff)
}

\keyword{utilities}
\keyword{regression}

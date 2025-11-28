\name{fake_formula}
\alias{fake_formula}

\title{Extended Processing of "Fake" Formulas}

\description{
  Create a simplified or "fake" version of a model formula for internal use in
  \pkg{gamlss2}. Many model formulas in \pkg{gamlss2} may contain complex terms
  such as smoothing functions, tensor-product interactions, distributional
  components, or multiple right-hand sides. These extended formulas cannot be
  passed directly to \code{\link[stats]{model.frame}}, which expects a
  traditional \code{\link[stats]{formula}} as used in \code{\link[stats]{lm}} or
  \code{\link[stats]{glm}}.

  \code{fake_formula()} therefore "fakes" a standard formula by rewriting the
  original model specification into a form that is acceptable to
  \code{model.frame}, while preserving all variable transformations
  (e.g., \code{log()}, \code{exp()}, arithmetic expressions) and extracting any
  special model terms. These special terms can then be handled separately during
  model setup and fitting.

  In summary, the function separates the parsing of model terms from the
  creation of the data frame: it provides a clean, simplified formula that
  ensures correct extraction and evaluation of variables.
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

\note{
  In some versions of the RStudio IDE, printing or inspecting a
  \code{fake_formula()} object may trigger a message of the form

    \code{'length = 2' in coercion to 'logical(1)'}.

  This is due to the way RStudio internally inspects objects of class
  \code{\link[Formula]{Formula}}, whose \code{length()} method returns a
  two-element vector (reflecting the number of left- and right-hand side
  components). The message does not indicate a problem with
  \code{fake_formula()}, and the returned object is valid and works correctly
  with \code{model.frame()} and subsequent processing. The message does not
  appear in a standard R session outside RStudio.
}

\seealso{
  \code{\link{gamlss2}}
}

\examples{
## basic formula, log(x3) should be kept
f <- y ~ x1 + x2 + log(x3)
ff <- fake_formula(f)
print(ff)

## including special model terms
## again, keep log(x3)
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
print(ff)

## extract only without special terms
ff <- fake_formula(f, nospecials = TRUE)
print(ff)

## extract only special terms
ff <- fake_formula(f, onlyspecials = TRUE)
print(ff)
}

\keyword{utilities}
\keyword{regression}


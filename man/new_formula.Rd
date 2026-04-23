\name{new_formula}
\alias{new_formula}

\title{Extracting a New Formula After Selection Algorithms}

\description{
   The generic function extracts the selected model terms from a fitted
   selection procedure and returns the corresponding reduced model formula.
}

\usage{
new_formula(object, ...)
}

\arguments{
  \item{object}{A fitted model object produced by a selection routine such as
    \code{\link{step_gamlss2}} or \code{\link{select_gamlss2}}.}
  \item{\dots}{Not used yet.}
}

\details{
  The function is intended for fitted models that store the outcome of an
  automatic term-selection step in \code{object$selection$formula}. It extracts
  the selected right-hand side terms and restores the original response on the
  left-hand side.

  For multi-parameter distributional models, the returned object is a
  \code{\link[Formula]{Formula}} with one right-hand side part for each
  distribution parameter. This makes it convenient to inspect the final model
  or to refit it without repeating the selection step.
}

\value{
  A \code{\link[Formula]{Formula}} containing the selected model terms.
}

\examples{
## no runnable example; see \code{\link{step_gamlss2}} and
## \code{\link{select_gamlss2}} for end-to-end usage
}

\seealso{
  \code{\link{step_gamlss2}}, \code{\link{select_gamlss2}}, \code{\link{gamlss2}}
}

\keyword{utilities}
\keyword{regression}

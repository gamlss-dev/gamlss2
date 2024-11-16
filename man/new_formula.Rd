\name{new_formula}
\alias{new_formula}

\title{Extracting a New Formula After Selection Algorithms}

\description{
   The generic function extracts the selected model terms after applying variable
   selection algorithms and constructs a new model formula based on those selected terms. 
}

\usage{
new_formula(object, ...)
}

\arguments{
  \item{object}{A fitted model.}
  \item{\dots}{Not used yet.}
}

\value{
  A \code{link[Formula]{Formula}} of selected model terms.
}

\seealso{
  \code{\link{step_gamlss2}}, \code{\link{select_gamlss2}}, \code{\link{gamlss2}}
}

\keyword{utilities}
\keyword{regression}


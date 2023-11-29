\name{gamlss2.control} 
\alias{gamlss2.control}

\title{Control Parameters for GAMLSS}

\description{
  Various parameters that control fitting of GAMLSS using \code{\link{gamlss2}}.
}

\usage{
gamlss2.control(optimizer = RS, step = 1, trace = TRUE,
  flush = TRUE, light = FALSE, expand = TRUE,
  model = TRUE, x = TRUE, y = TRUE, \dots)
}

\arguments{
  \item{optimizer}{Function, the optimizer function that should be used for fitting.}
  \item{step}{Numeric, a step length control parameter.}
  \item{trace}{Logical, should information be printed while the algorithm is running?}
  \item{flush}{Logical, use \code{\link[utils]{flush.console}} for displaying the
    current output in the console.}
  \item{light}{Logical, if set to \code{light = TRUE}, no model frame, response,
    model matrix and other design matrices will be part of the return value.}
  \item{expand}{Logical, if fewer formulas are supplied than there are parameters
    of the distribution, should formulas with intercept only formulas be added?}
  \item{model}{Logical, should the \emph{model frame} be included as component of the
    returned object.}
  \item{x}{Logical, indicating whether the model matrix should be included
    as component of the returned object.}
  \item{y}{Logical, should the response be included as component of the returned object.}
  \item{\dots}{Further control parameters to be part of the return value.}
}

\details{
  The control parameters in \code{gamlss2.control} can also be extended, e.g., if another
  optimization function is used, newly specified control parameters are automatically
  passed on to this function.
} 

\value{
  A list with the arguments specified.
}

\seealso{
  \code{\link{RS}}, \code{\link{gamlss2}}
}

\examples{## Load the abdominal circumference data
data("abdom", package = "gamlss.data")

## Specify the model Formula.
f <- y ~ s(x) | s(x)

## Estimate model with different step length control.
b1 <- gamlss2(f, data = abdom, family = BCT, step = 1)
b2 <- gamlss2(f, data = abdom, family = BCT, step = 0.9)
}

\keyword{models}
\keyword{regression}


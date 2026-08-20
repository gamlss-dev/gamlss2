\name{gamlss2_control} 
\alias{gamlss2_control}

\title{Control Parameters}

\description{
  Control parameters for fitting GAMLSS models with \code{\link{gamlss2}}.
}

\usage{
gamlss2_control(optimizer = RS, trace = TRUE,
  flush = TRUE, light = FALSE, expand = TRUE,
  model = TRUE, x = TRUE, y = TRUE,
  fixed = FALSE, \dots)
}

\arguments{
  \item{optimizer}{Function, the optimizer to be used for fitting.}
  \item{trace}{Logical, should information be printed while the algorithm is running?}
  \item{flush}{Logical, whether to use \code{\link[utils]{flush.console}} to display
    current output in the console.}
  \item{light}{Logical, if set to \code{light = TRUE}, no model frame, response,
    model matrix and other design matrices will be part of the return value.}
  \item{expand}{Logical, if fewer formulas are supplied than there are parameters
    of the distribution, should intercept-only formulas be added automatically?}
  \item{model}{Logical, should the \emph{model frame} be included as a component of the
    returned object.}
  \item{x}{Logical, indicating whether the model matrix should be included
    as a component of the returned object.}
  \item{y}{Logical, should the response be included as a component of the returned object.}
  \item{fixed}{Named logical vector indicating which parameters should be fixed during estimation.
    See \code{\link{gamlss2_start}} for examples.}
  \item{\dots}{Further control parameters to be included in the return value, for example
    parameters used by the optimizer function \code{\link{RS}}.}
}

\details{
  The set of control parameters can be extended. For example, if a different
  optimizer is used, newly specified control parameters are passed on to that
  optimizer automatically.

  For the \code{RS} and \code{CG} optimizers, argument
  \code{demmler.reinsch} controls the use of a weighted Demmler--Reinsch
  reparameterization when selecting the smoothing parameter of an eligible
  univariate, single-penalty smooth. The default, \code{"auto"}, starts with
  the direct solver and switches only when the current smoothing-parameter
  search reaches the edge of its one-decade search window. This indicates that
  another costly optimizer restart is needed and amortizes the setup of the
  reparameterization. Set \code{demmler.reinsch = TRUE} to force it for eligible
  updates, or \code{FALSE} to disable it. To preserve the numerical trajectory of
  nonlinear distributional models, the initial outer iteration continues to use the
  direct solver.

  Smooths with fixed smoothing parameters or multiple penalties, including
  adaptive and tensor-product smooths, always use the direct solver.
} 

\value{
  A list with the arguments specified.
}

\seealso{
  \code{\link{RS}}, \code{\link{gamlss2}}, \code{\link{gamlss2_start}}
}

\examples{\dontrun{## specify the model formula
f <- dist ~ s(speed) | s(speed)

## estimate model with different step length
## control in the RS algorithm
m1 <- gamlss2(f, data = cars, family = GA, step = 1)
m2 <- gamlss2(f, data = cars, family = GA, step = 0.1)
}}

\keyword{models}
\keyword{regression}

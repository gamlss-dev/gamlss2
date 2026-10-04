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
  \item{optimizer}{Function, the optimizer to be used for fitting. Use
    \code{\link{JR}} for joint REML estimation.}
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

  Optimizer-specific controls and algorithm details are documented in
  \code{\link{RS}} and \code{\link{CG}}, and in \code{\link{JR}} for
  joint REML estimation. Smoothing parameter selection is described in
  \code{\link{gamlss2_smoothing}}.
} 

\value{
  A list with the arguments specified.
}

\seealso{
  \code{\link{RS}}, \code{\link{CG}}, \code{\link{JR}},
  \code{\link{gamlss2}}, \code{\link{gamlss2_start}},
  \code{\link{gamlss2_smoothing}}
}

\examples{\dontrun{## specify the model formula
f <- dist ~ s(speed) | s(speed)

## specify fitting controls
ctrl <- gamlss2_control(trace = FALSE, maxit = c(50, 100))
m <- gamlss2(f, data = cars, family = GA, control = ctrl)
}}

\keyword{models}
\keyword{regression}

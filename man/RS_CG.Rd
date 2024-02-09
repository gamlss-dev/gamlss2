\name{RS} 
\alias{RS}
\alias{CG}

\title{Rigby and Stasinopoulos (RS) & Cole and Green (CG) Algorithm}

\description{
  The function \code{RS()} implements the algorithm of Rigby and Stasinopoulos,
  the function \code{CG()} implements the algorithm of Cole and Green for estimating
  a GAMLSS with \code{\link{gamlss2}}.
}

\usage{
## Rigby and Stasinopoulos algorithm.
RS(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)

## Cole and Green algorithm.
CG(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)
}

\arguments{
  \item{x}{The full model matrix to be used for fitting.}
  \item{y}{The response vector or matrix.}
  \item{specials}{A named list of special model terms, e.g., including design and
    penalty matrices for fitting smooth terms using \code{\link[mgcv]{smooth.construct}}.}
  \item{family}{A family object, see \code{\link{gamlss2.family}}.}
  \item{offsets}{If supplied, a list or data frame of possible model offset.}
  \item{weights}{If supplied, a numeric vector of weights.}
  \item{start}{Starting values, either for the parameters of the response distribution or,
    if specified as a named list in which each element of length one is named with
    \code{"(Intercept)"}, the respective intercepts are initialized. If starting values are
    specified as a named list, data frame or matrix, where each element/column is a
    vector with the same length as the number of observations in the data,
    the respective predictors are initialized with these. See the examples for
    \code{\link{gamlss2}}.}
  \item{xterms}{A named list specifying the linear model terms. Each named list
    element represents one parameter as specified in the family object.}
  \item{sterms}{A named list specifying the special model terms. Each named list
    element represents one parameter as specified in the family object.}
  \item{control}{Further control arguments as specified within the call of \code{\link{gamlss2}}.
    See the details.}
}

\details{
  Functions \code{RS()} and \code{CG()} are called within \code{\link{gamlss2}}. Both functions
  implement a backfitting algorithm for estimating GAMLSS. For algorithm details see
  Rigby and Stasinopoulos (2005).

  The functions use the following control arguments:
  \itemize{
    \item \code{eps}: Numeric vector of length 2, the stopping criterion.
      Default is \code{eps = c(1e-05, 1e-05)} for the outer and the inner backfitting loop.
    \item \code{maxit}: Integer vector of length 2, the maximum number of iterations of
      the outer and inner backfitting loop. Default is \code{maxit = c(100, 10)}.
    \item \code{step}: Numeric, the step length control parameter. Default is \code{step = 1}.
      Note that if \code{step} is set smaller than 1, it might be appropriate to lower
      the stopping criterion \code{eps}, too.
    \item \code{CG}: Integer, the number of iterations when to start the CG correction.
      Default is \code{CG = Inf}.
    \item \code{trace}: Logical, should information be printed while the algorithm is running?
    \item \code{flush}: Logical, use \code{\link[utils]{flush.console}} for displaying the
      current output in the console.
  }

  To facilitate the development of new algorithms for \pkg{gamlss2}, users can exchange
  them using the \code{optimizer} argument in \code{\link{gamlss2_control}}. Users
  developing new model fitting functions are advised to use these functions as templates
  and pass them to \code{\link{gamlss2_control}}. Alternatively, users can replace the
  optimizer function by adding a named list element, \code{"optimizer"}, to the family object.
  For instructions on setting up new families in \pkg{gamlss2}, see \code{\link{gamlss2.family}}.
} 

\value{
  Functions \code{RS()} and \code{CG()} return a named list of class \code{"gamlss2"}
  containing the following objects:
  \item{fitted.values}{A data frame of the fitted values of the modeled parameters
    of the selected distribution.}
  \item{fitted.specials}{A named list, one element for each parameter of the distribution,
    containing the fitted model object information of special model terms.}
  \item{fitted.linear}{A named list, one element for each parameter of the distribution,
    containing the information on fitted linear effects.}
  \item{coefficients}{A named list, one element for each parameter of the distribution,
    containing the estimated parameters of the linear effects.}
  \item{elapsed}{The elapsed runtime of the algorithm.}
  \item{iterations}{How many iterations the algorithm performed.}
  \item{logLik}{The final value of the log-likelihood of the model.}
  \item{control}{All control arguments used as supplied from function \code{\link{gamlss2_control}}.}
}

\references{
  Rigby, R.A. and Stasinopoulos, D.M. (2005). Generalized Additive Models for Location,
  Scale and Ahape (with Discussion). \emph{Journal of the Royal Statistical Society},
  Series C (Applied Statistics), \bold{54}: 507--554.
  \doi{10.1111/j.1467-9876.2005.00510.x}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2_control}}, \code{\link{gamlss2.family}}
}

\examples{\dontrun{## Load the abdominal circumference data
data("abdom", package = "gamlss.data")

## Specify the model Formula.
f <- y ~ s(x) | s(x) | s(x) | s(x)

## Estimate model using RS (default).
b <- gamlss2(f, data = abdom, family = BCT, optimizer = RS)

## Now with CG.
b <- gamlss2(f, data = abdom, family = BCT, optimizer = CG)

## First 2 RS iterations and afterwards switch to CG.
b <- gamlss2(f, data = abdom, family = BCT, CG = 2)
}}

\keyword{models}
\keyword{regression}


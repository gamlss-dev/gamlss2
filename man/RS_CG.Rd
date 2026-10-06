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
## rigby and stasinopoulos algorithm
RS(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)

## cole and green algorithm
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
    \item \code{eps}: Numeric vector of length 1 to 3, the stopping criterion
      for the global, per-parameter, and CG correction loops. The CG tolerance
      defaults to the per-parameter tolerance.
    \item \code{maxit}: Integer vector of length 1 to 3, the maximum number of global
      parameter sweeps and per-parameter IRLS/backfitting iterations. Default is
      \code{maxit = c(20, 50)}. For the CG algorithm, an optional third element
      limits the internal correction sweeps and defaults to the second element.
      RS checks global convergence after each parameter sweep and ignores the
      third element.
    \item \code{step}: Numeric, the step length control parameter. Default is
      \code{step = 1}. With fixed smoothing parameters, it scales complete RS
      parameter updates or CG correction sweeps after the first update.
      Smoothing-parameter changes and their coefficients are accepted or
      rejected together. Smaller steps may require a smaller \code{eps}.
    \item \code{autostep}: Logical, should an update that decreases its
      objective be halved automatically? Default is \code{TRUE}. Updates use
      the likelihood, penalized likelihood, or smoothing-selection criterion
      as appropriate; see \code{\link{gamlss2_smoothing}}.
    \item \code{CG}: Integer, the outer iteration in which to start the CG
      correction. Thus, \code{CG = 5} uses RS for the first four outer
      iterations and CG in the fifth. \code{CG = 0} starts with CG immediately;
      the default is \code{CG = Inf}.
    \item \code{trace}: Logical, should information be printed while the algorithm is running?
    \item \code{rs.cache}: Logical, reuse generated family evaluations and weighted
      matrix calculations with identical inputs within a fit? Default is
      \code{TRUE}; set \code{FALSE} to disable caching for comparisons.
    \item \code{flush}: Logical, use \code{\link[utils]{flush.console}} for displaying the
      current output in the console.
    \item \code{ridge}: Logical, should automatic ridge penalization be applied
      only to linear effects, without penalizing the intercept? For each parameter of
      the distribution the optimum ridge penalty is estimated using an information criterion.
      Possible options are \code{criterion = c("aic", "aicc", "bic", "gaic", "gcv")}. The
      default is \code{criterion = "gaic"} and argument \code{K = 2}, which can
      be set in \code{\link{gamlss2_control}}.
  }

  For smoothing parameter selection, numerical solvers, and examples, see
  \code{\link{gamlss2_smoothing}}.

  Use these functions as templates for custom fitting algorithms, supplied via
  \code{optimizer} in \code{\link{gamlss2_control}} or the family object's
  \code{"optimizer"} element. For new families, see \code{\link{gamlss2.family}}.
} 

\value{
  Functions \code{RS()} and \code{CG()} return a named list of class \code{"gamlss2"}
  containing the following objects:
  \item{fitted.values}{A data frame of the fitted values of the predictors
    of the selected distribution.}
  \item{fitted.specials}{A named list, one element for each parameter of the distribution,
    containing the fitted model object information of special model terms.}
  \item{fitted.linear}{A named list, one element for each parameter of the distribution,
    containing the information on fitted linear effects.}
  \item{coefficients}{A named list, one element for each parameter of the distribution,
    containing the estimated parameters of the linear effects.}
  \item{elapsed}{The elapsed runtime of the algorithm.}
  \item{iterations}{How many iterations the algorithm performed.}
  \item{converged}{Logical, did the global iteration satisfy the stopping criterion?}
  \item{logLik}{The final value of the log-likelihood of the model.}
  \item{control}{All effective control arguments used by the optimizer.}
}

\references{
  Rigby RA, Stasinopoulos DM (2005). \dQuote{Generalized Additive Models for Location,
  Scale and Shape (with Discussion).} \emph{Journal of the Royal Statistical Society,
  Series C (Applied Statistics)}, \bold{54}, 507--554.
  \doi{10.1111/j.1467-9876.2005.00510.x}
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2_control}},
  \code{\link{gamlss2.family}}, \code{\link{gamlss2_smoothing}}
}

\examples{\dontrun{## count response example using HarzTraffic data
data("HarzTraffic", package = "gamlss2")

## specify the model formula
f <- bikes ~ s(yday, bs = "cc") + s(temp) + s(rain) | .

## estimate model using RS (default)
m1 <- gamlss2(f, data = HarzTraffic, family = NBI, optimizer = RS)

## now with CG
m2 <- gamlss2(f, data = HarzTraffic, family = NBI, optimizer = CG)

## first 2 RS iterations and afterwards switch to CG
m3 <- gamlss2(f, data = HarzTraffic, family = NBI, CG = 3)
}}

\keyword{models}
\keyword{regression}


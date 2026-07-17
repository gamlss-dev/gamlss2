\name{stepwise} 
\alias{stepwise}
\alias{step_gamlss2}

\title{Stepwise Model Term Selection Using GAIC}

\description{
  The optimizer function \code{stepwise()} performs stepwise model term selection
  using a Generalized Akaike Information Criterion (GAIC). Estimation is based on the
  Rigby and Stasinopoulos (RS) and Cole and Green (CG) algorithms, as implemented in
  function \code{\link{RS}}.
}

\usage{
## wrapper function for stepwise GAMLSS estimation
step_gamlss2(formula, ..., K = 2,
  strategy = c("both.linear", "both"), keeporder = FALSE,
  cores = 1L)

## stepwise optimizer function
stepwise(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)
}

\arguments{
  \item{formula}{A model formula for \code{\link{gamlss2}}.}
  \item{\dots}{Arguments passed to \code{\link{gamlss2}}.}
  \item{K}{Numeric, the penalty for the \code{\link{GAIC}}.}
  \item{strategy}{Character, the strategy to be used by the stepwise algorithm.
    Possible options are \code{"forward.linear"}, \code{"forward"}, \code{"backward"},
    \code{"backward.linear"}, \code{"replace"}, \code{"replace.linear"},
    \code{"both"}, \code{"both.linear"}, \code{"all"} and
    \code{"all.linear"}. See the details.}
  \item{keeporder}{Logical. For the different stepwise strategies,
    should the updates be performed sequentially according to the order of the parameters
    of the response distribution as specified in the family (see \code{\link{gamlss2.family}}),
    or should the selection search be performed across all parameters?}
  \item{cores}{Integer. If \code{cores > 1L}, \code{\link[parallel]{mclapply}} is used
    to speed up computations within the selection steps.}
  \item{x}{The full model matrix to be used for fitting.}
  \item{y}{The response vector or matrix.}
  \item{specials}{A named list of special model terms, e.g., including design and
    penalty matrices for fitting smooth terms using \code{\link[mgcv]{smooth.construct}}.}
  \item{family}{A family object, see \code{\link{gamlss2.family}}.}
  \item{offsets}{If supplied, a list or data frame of model offsets.}
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
  \item{control}{Further control arguments as specified within the call
    of \code{\link{gamlss2}}.}
}

\details{
  The wrapper function \code{step_gamlss2()} calls \code{\link{gamlss2}}
  using the \code{stepwise()} optimizer function.

  The stepwise algorithm can apply the following strategies:
  \enumerate{
    \item Each predictor must include an intercept.
    \item In a forward selection step, model terms with the highest improvement
      on the GAIC are selected.
    \item In a replacement step, each model term is tested to see if an exchange with
      a term not yet selected improves the GAIC.
    \item In a backward step, model terms are removed if this further improves the GAIC.
    \item In a bidirectional step, model terms can be either added or removed.
    \item In addition, the forward, backward and replace selection step can be combined.
  }

  The selected strategies are iterated until no further improvement is achieved.

  The different strategies can be selected using argument \code{strategy}; see
  the examples. Possible values are \code{"both"}, \code{"forward"},
  \code{"backward"}, \code{"replace"}, and \code{"all"}.
  Here, \code{strategy = "all"} combines the forward, backward, and replace selection steps.

  In addition, each of the steps 2-4 can be applied to linear model terms only,
  prior to performing the steps for all model terms. This can be done by additionally
  setting \code{strategy} to \code{"both.linear"}, \code{"forward.linear"},
  \code{"backward.linear"}, \code{"replace.linear"}, or \code{"all.linear"}.

  The default is \code{strategy = c("both.linear", "both")} and \code{keeporder = FALSE}.

  Note that each of the steps 2--4 can be performed while maintaining the order of the
  parameters of the response distribution, i.e., if the \code{keeporder = TRUE}
  argument is set, then the parameters will be updated in the order specified in the
  \code{\link{gamlss2.family}}. Using backward elimination, the model terms are
  deselected in reverse order.
}

\value{
  The optimizer function \code{stepwise()} returns the final model as a named list
  of class \code{"gamlss2"}. See the return value of function \code{\link{RS}}.
  The wrapper function \code{step_gamlss2()} also returns the final model.
}

\seealso{
  \code{\link{new_formula}}, \code{\link{gamlss2}}, \code{\link{gamlss2_control}}, \code{\link{RS}}
}

\examples{\dontrun{## because of possible linear interactions,
## scale the covariates first
air <- subset(airquality, !is.na(Ozone))
air$Wind <- scale(air$Wind)
air$Temp <- scale(air$Temp)

## the formula defines the searching scope
f <- Ozone ~ Wind + Temp + Wind:Temp +
  ti(Wind) + ti(Temp) + ti(Wind, Temp) | .

## estimate a Gamma model using the stepwise algorithm
m1 <- step_gamlss2(f, data = air, family = GA, K = 2)

## same with
## m <- gamlss2(f, data = air, family = NO,
##   optimizer = stepwise, K = 2)

## show the new formula of selected model terms
new_formula(m1)

## final model summary
summary(m1)

## effect plots
plot(m1)

## diagnostic plots
plot(m1, which = "resid")

## plot GAIC
plot(m1, which = "selection")

## use forward linear, replace and backward strategy
m2 <- step_gamlss2(f, data = air, family = GA, K = 2,
  strategy = c("forward.linear", "replace", "backward"))

## plot GAIC
plot(m2, which = "selection")
}}

\keyword{models}
\keyword{regression}

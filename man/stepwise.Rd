\name{stepwise} 
\alias{stepwise}
\alias{step_gamlss2}

\title{Stepwise Model Term Selection Using GAIC}

\description{
  The optimizer function \code{stepwise()} performs stepwise model term selection
  using a Generalized Akaike Information Criterion (GAIC). Estimation is based on the
  Rigby and Stasinopoulos (RS) & Cole and Green (CG) algorithm as implemented in
  function \code{\link{RS}}.
}

\usage{
## Wrapper function for stepwise GAMLSS estimation.
step_gamlss2(formula, ..., K = 2,
  strategy = c("both.linear", "both"), keeporder = FALSE,
  cores = 1L)

## Stepwise optimizer function.
stepwise(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)
}

\arguments{
  \item{formula}{A model formula for \code{\link{gamlss2}}.}
  \item{\dots}{Arguments passed to \code{\link{gamlss2}}.}
  \item{K}{Numeric, the penalty for the \code{\link{GAIC}}.}
  \item{strategy}{Character, the strategy that should be applied for the stepwise algorithm.
    Possible options are \code{"forward.linear"}, \code{"forward"}, \code{"backward"},
    \code{"backward.linear"}, \code{"replace"}, \code{"replace.linear"},
    \code{"both"}, \code{"both.linear"}. See the details.}
  \item{keeporder}{Logical, For the different strategies of the stepwise algorithm,
    should the updates be performed sequentially according to the order of the parameters
    of the response distribution as specified in the family (see \code{\link{gamlss2.family}}),
    or should the selection search be performed across all parameters?}
  \item{cores}{Integer, if \code{cores > 1L}, function \code{\link[parallel]{mclapply}}
    function is used to speed up computations using multiple cores within the selection steps.}
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
      a model term not yet selected will improve the GAIC.
    \item In a backward step, model terms are deselected,
      if the GAIC can be further improved.
    \item In a bidirectional step, model terms can be either added or removed.
    \item In addition, the forward, backward and replace selection step can be combined.
  }

  The selected strategies are iterated until no further improvement is achieved.

  The different strategies can be selected using argument \code{strategy}. Please see 
  the examples. Possible values are \code{strategy = c("both", "forward", "backward", "replace", "all")}.
  Here, \code{strategy = "all"} combines the forward, backward and replace selection step.

  In addition, each of the steps 2-4 can be applied to linear model terms only,
  prior to performing the steps for all model terms. This can be done by additionally
  setting \code{strategy = c("both.linear", "forward.linear", "backward.linear", "replace.linear", "all.linear")}.

  The default is \code{strategy = c("both.linear", "both")} and \code{keeporder = FALSE}.

  Note that each of the steps 2-4 can be performed while maintaining the order of the 
  parameters of the response distribution, i.e., if the \code{keeporder = TRUE}
  argument is set, then the parameters will be updated in the order specified in the
  \code{\link{gamlss2.family}}. Using backward elimination, the model terms are
  deselected in reverse order.
}

\value{
  The optimizer function \code{stepwise()} returns the final model as named list
  of class \code{"gamlss2"}. See the return value of function \code{\link{RS}}.
  The wrapper function \code{step_gamlss2()} also returns the final model.
}

\seealso{
  \code{\link{new_formula}}, \code{\link{gamlss2}}, \code{\link{gamlss2_control}}, \code{\link{RS}}
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
\dontrun{## load the Munich rent data
data("rent", package = "gamlss.data")

## because of possible linear interactions,
## scale the covariates first
rent$Fl <- scale(rent$Fl)
rent$A <- scale(rent$A)

## the Formula defines the searching scope
f <- R ~ Fl + A + Fl:A + loc + s(Fl) + s(A) + te(Fl, A) |
  Fl + A + loc + Fl:A + s(Fl) + s(A) + te(Fl, A)

## estimate a Gamma model using the stepwise algorithm
b <- step_gamlss2(f, data = rent, family = GA, K = 2)

## same with
## b <- gamlss2(f, data = rent, family = GA, optimizer = stepwise, K = 2)

## show the new formula of selected model terms
new_formula(b)

## final model summary
summary(b)

## effect plots
plot(b)

## diagnostic plots
plot(b, which = "resid")

## plot GAIC
plot(b, which = "selection")

## use forward linear, replace and backward strategy
b <- step_gamlss2(f, data = rent, family = GA, K = 2,
  strategy = c("forward.linear", "replace", "backward"))

## more complex model
## note, the third parameter
## nu does not include any model terms
f <- R ~ Fl + A + Fl:A + loc + s(Fl) + s(A) + te(Fl, A) |
  Fl + A + loc + Fl:A + s(Fl) + s(A) + te(Fl, A) |
  1 |
  Fl + A + loc + Fl:A + s(Fl) + s(A) + te(Fl, A)

## model using the BCT family
b <- step_gamlss2(f, data = rent, family = BCT,
  K = 2, strategy = c("forward.linear", "both"),
  keeporder = TRUE)

## plot GAIC
plot(b, which = "selection")
}}

\keyword{models}
\keyword{regression}


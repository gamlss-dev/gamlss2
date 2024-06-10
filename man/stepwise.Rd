\name{stepwise} 
\alias{stepwise}
\alias{stepGAMLSS}

\title{Stepwise Model Term Selection Using GAIC}

\description{
  The optimizer function \code{stepwise()} performs stepwise model term selection
  using a Generalized Akaike Information Criterion (GAIC). Estimation is based on the
  Rigby and Stasinopoulos (RS) & Cole and Green (CG) algorithm as implemented in
  function \code{\link{RS}}.
}

\usage{
## Wrapper function for stepwise GAMLSS estimation.
stepGAMLSS(formula, ..., K = 2,
  strategy = c("both.linear", "both"), keeporder = FALSE)

## Stepwise optimizer function.
stepwise(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)
}

\arguments{
  \item{formula}{A model formula for \code{\link{gamlss2}}.}
  \item{\dots}{Arguments passed to \code{\link{gamlss2}}.}
  \item{K}{Numeric, the penalty for the \code{\link{GAIC}}.}
  \item{strategy}{Character, the strategy that should be applied for the stepwise algorithm. Possible options are \code{"forward.linear"}, \code{"forward"}, \code{"backward"},
    \code{"backward.linear"}, \code{"replace"}, \code{"replace.linear"},
    \code{"both"}, \code{"both.linear"}. See the details.}
    \item{keeporder}{Logical, For the different strategies of the stepwise algorithm,
    should the updates be performed sequentially according to the order of the parameters of the response distribution as specified in the family (see \code{\link{gamlss2.family}}),
    or should the selection search be performed across all parameters?}
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
  The wrapper functions simply calls \code{\link{gamlss2}} using the \code{stepwise()} optimizer
  function.

  The stepwise algorithm uses the following strategy:
  \enumerate{
    \item Each predictor must include an intercept.
    \item In a forward selection step, model terms with the highest improvement
      on the GAIC are selected.
    \item In a replacement step, each model term is tested to see if an exchange with
      a model term not yet selected will improve the GAIC.
    \item In a backward step, model terms are deselected,
      if the GAIC can be further improved.
    \item In a bidirectional step, model terms can be either added or removed.
  }

  Items 2-5 are iterated until no further improvement is achieved.

  The different strategies can be selected using argument \code{strategy}. Please see 
  the examples. Possible values are \code{strategy = c("both", "forward", "backward", "replace")}.

  In addition, each of the steps 2-4 can be applied to linear model terms only,
  prior to performing the steps for all model terms. This can be done by additionally
  setting \code{strategy = c("both.linear", "forward.linear", "backward.linear", "replace.linear")}.

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
  The wrapper function \code{stepGAMLSS()} also returns the final model.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{gamlss2_control}}, \code{\link{RS}}
}

\examples{\dontrun{## Load the Munich rent data
data("rent", package = "gamlss.data")

## Because of possible linear interactions,
## scale the covariates first.
rent$Fl <- scale(rent$Fl)
rent$A <- scale(rent$A)

## The formula defines the searching scope.
f <- R ~ Fl + A + Fl:A + loc + s(Fl) + s(A) + te(Fl, A) |
  Fl + A + loc + Fl:A + s(Fl) + s(A) + te(Fl, A)

## Estimate a Gamma model using the stepwise algorithm.
b <- stepGAMLSS(f, data = rent, family = GA, K = 2)

## Same with.
## b <- gamlss2(f, data = rent, family = GA, optimizer = stepwise, K = 2)

## Final model summary.
summary(b)

## Effect plots.
plot(b)

## Diagnostic plots.
plot(b, which = "resid")

## Plot GAIC.
plot(b, which = "selection")

## Use forward linear, replace and backward strategy.
b <- stepGAMLSS(f, data = rent, family = GA, K = 2,
  strategy = c("forward.linear", "replace", "backward"))

## More complex model. Note, the third parameter
## nu does not include any model terms.
f <- R ~ Fl + A + Fl:A + loc + s(Fl) + s(A) + te(Fl, A) |
  Fl + A + loc + Fl:A + s(Fl) + s(A) + te(Fl, A) |
  1 |
  Fl + A + loc + Fl:A + s(Fl) + s(A) + te(Fl, A)

## Model using the BCT family.
b <- stepGAMLSS(f, data = rent, family = BCT,
  K = 2, strategy = c("forward.linear", "both"),
  keeporder = TRUE)

## Plot GAIC.
plot(b, which = "selection")
}}

\keyword{models}
\keyword{regression}


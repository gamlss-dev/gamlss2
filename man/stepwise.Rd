\name{stepwise} 
\alias{stepwise}

\title{Stepwise Model Term Selection Using GAIC}

\description{
  The optimizer function \code{stepwise()} performs stepwise model term selection
  using a Generalized Akaike Information Criterion (GAIC). Estimation is based on the
  Rigby and Stasinopoulos (RS) & Cole and Green (CG) algorithm as implemented in
  function \code{\link{RS}}.
}

\usage{
stepwise(x, y, specials, family, offsets,
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
  Note that the GAIC penalty is controlled by the \code{"K"} argument, which can be easily
  specified in the \code{\link{gamlss2}} call (using \code{\link{gamlss2_control}}). The default
  is the BIC penalty. See also the examples.

  The stepwise algorithm uses the following strategy:
  \itemize{
    \item Each predictor must include an intercept.
    \item First, in a forward selection step, the linear terms with the highest improvement
      on the GAIC are selected.
    \item Second, in a forward selection step, all linear or
      special model terms (see \code{\link{specials}}) with the highest improvement
      on the GAIC are selected.
    \item Third, each model term is tested to see if an exchange with
      a model term not yet selected will improve the GAIC.
    \item Fourth, in a backward step, linear or nonlinear (special) terms are deselected,
      if the GAIC can be further improved.
  }

  The last three items are iterated until no further improvement is achieved.

  The different strategies can be selected using an argument named
  \code{strategy} in the \code{\link{gamlss2_control}} function. Please see 
  the examples. Possible values are
  \code{strategy = c("forward.linear", "forward", "backward", "replace")}.
} 

\value{
  The function returns the final model as named list of class \code{"gamlss2"}. See
  the return value of function \code{\link{RS}}.
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
b <- gamlss2(f, data = rent, family = GA, optimizer = stepwise, K = 2)

## Final model summary.
summary(b)

## Effect plots.
plot(b)

## Diagnostic plots.
plot(b, which = "resid")

## Plot GAIC.
plot(b, which = "selection")

## Use only forward linear and replace strategy.
b <- gamlss2(f, data = rent, family = NO,
  optimizer = stepwise, K = 2,
  strategy = c("forward.linear", "replace"))
}}

\keyword{models}
\keyword{regression}


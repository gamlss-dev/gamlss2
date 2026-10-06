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
step_gamlss2(formula, ..., fixed = NULL, K = 2,
  strategy = c("both.linear", "both"), keeporder = FALSE,
  cores = 1L)

## stepwise optimizer function
stepwise(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)
}

\arguments{
  \item{formula}{A model formula for \code{\link{gamlss2}}.}
  \item{\dots}{Arguments passed to \code{\link{gamlss2}}.}
  \item{fixed}{An optional one-sided formula, \code{Formula} object, or list of
    formulas specifying model terms that must always remain in the model. Terms
    must also occur in the search scope given by \code{formula}. For
    multi-parameter families, right-hand sides or list elements correspond to
    the distributional parameters in their usual order.}
  \item{K}{Numeric, the penalty for the \code{\link{GAIC}}.}
  \item{strategy}{Character vector of selection strategies; see the details.}
  \item{keeporder}{Logical. Update distributional parameters sequentially
    rather than search across all parameters? See the details for the order.}
  \item{cores}{Integer, number of workers for candidate fits. The default
    \code{cores = 1L} runs sequentially. Larger values use
    \code{\link[future.apply]{future_lapply}} with a temporary
    \code{\link[future]{multisession}} plan, including on Windows. The caller's
    future plan is restored on exit. This requires the suggested packages
    \pkg{future} and \pkg{future.apply}.}
  \item{x, y, specials, family, offsets, weights, start, xterms, sterms}{
    Model components and starting values passed to the optimizer, as documented
    for \code{\link{RS}}.}
  \item{control}{Further control arguments as specified within the call
    of \code{\link{gamlss2}}.}
}

\details{
  The wrapper function \code{step_gamlss2()} calls \code{\link{gamlss2}}
  using the \code{stepwise()} optimizer function.

  Argument \code{strategy} selects one or more strategies:
  \describe{
    \item{\code{"forward"}}{Add the term that most improves the GAIC.}
    \item{\code{"backward"}}{Remove the term that most improves the GAIC.}
    \item{\code{"replace"}}{Exchange a selected term for an unselected term
      if this improves the GAIC.}
    \item{\code{"both"}}{Consider both additions and removals.}
    \item{\code{"all"}}{Combine forward, backward, and replacement steps.}
  }

  Append \code{".linear"} to any strategy to restrict it to linear terms.
  These strategies run before those applied to all model terms. The default is
  \code{strategy = c("both.linear", "both")} and \code{keeporder = FALSE}.
  Selected strategies are iterated until the GAIC no longer improves.

  With \code{keeporder = TRUE}, parameters are updated in the order specified
  by \code{\link{gamlss2.family}}, or in reverse order for backward elimination.
  Otherwise, selection searches across all parameters.

  Terms supplied in \code{fixed} are included in the initial model and are never
  eligible for removal or replacement. The intercept is always retained and
  need not be included in \code{fixed}. Fixing a linear interaction also fixes
  all of its available lower-order marginal terms.

  Linear selection is performed on formula terms, not individual columns of the
  model matrix. Thus, all columns representing a factor or an interaction with
  a factor are added, removed, or replaced together. Linear interactions obey
  strong hierarchy: an interaction can enter only after its available
  lower-order terms have entered, and a lower-order term cannot be removed while
  an included interaction depends on it.
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

## always retain both linear main effects
m2 <- step_gamlss2(f, data = air, family = GA, K = 2,
  fixed = ~ Wind + Temp | .)

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

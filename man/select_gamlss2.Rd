\name{select_gamlss2}
\alias{select_gamlss2}
\alias{sRS}

\title{Smooth Model Term Selection with Additional Penalties}

\description{
  The function \code{select_gamlss2()} allows for penalizing all \pkg{mgcv}
  model terms with an additional shrinkage penalty, which can shrink some terms to
  zero, effectively selecting them out of the model. In addition to this penalty,
  model terms are selected based on two criteria: the estimated degrees of freedom
  of the term and the percentage of the predictor range covered by the model term.
  These two thresholds aim to mimic a natural selection process, similar to what
  one might do by inspecting summaries and effect plots, ensuring that only
  relevant terms are retained in the final model. 
}

\usage{
## Wrapper function.
select_gamlss2(formula, ..., criterion = "BIC", thres = c(0.9, 0.2))

## Modified RS optimizer function.
sRS(x, y, specials, family, offsets,
  weights, start, xterms, sterms, control)
}

\arguments{
  \item{formula}{A GAM-type \code{\link[stats]{formula}} or
    \code{\link[Formula]{Formula}}. All smooth terms of the \pkg{mgcv} package
    are supported, see also \code{\link[mgcv]{formula.gam}}.
    For \code{gamlss.list()} \code{formula} is a list of
    \code{\link[stats]{formula}}s.}
  \item{\dots}{Arguments passed to argument \code{control} in function \code{sRS()}.}
  \item{criterion}{The information criterion to be used for estimating the
    shrinkage parameters. This can also be a vector of length 2, where the first
    element specifies the criterion to be used during the selection step, and the
    second element specifies the criterion to be used during the refitting step.
    Possible options are \code{"BIC"}, \code{"GCV"}, \code{"AIC"}, \code{"AICc"} and
    \code{"GAIC"} with user defined penalty \code{K}, default is \code{K = 2}.
  \item{thres}{A vector of thresholds used for model term selection. The first
    element controls the minimum allowed estimated degrees of freedom for a model
    term to enter the model. The second element specifies the minimum percentage of
    the total estimated predictor range required for a term to be included in
    the model.}
  \item{x}{The full model matrix to be used for fitting.}
  \item{y}{The response vector or matrix.}
  \item{specials}{A named list of special model terms, e.g., including design and
    penalty matrices for fitting smooth terms using \code{\link[mgcv]{smooth.construct}}.}
  \item{family}{A family object, see \code{\link{gamlss2.family}}.}
  \item{offsets}{If supplied, a list or data frame of possible model offset.}
  \item{weights}{If supplied, a numeric vector of weights.}
  \item{start}{Starting values, either for the parameters of the response distribution or,
    if specified as a named list in which each element of length one is named with
    \code{"(Intercept)"}, the respective intercepts are initialized. If starting
    values are specified as a named list, data frame or matrix, where each
    element/column is a vector with the same length as the number of observations
    in the data, the respective predictors are initialized with these.
    See the examples for \code{\link{gamlss2}}.}
  \item{xterms}{A named list specifying the linear model terms. Each named list
    element represents one parameter as specified in the family object.}
  \item{sterms}{A named list specifying the special model terms. Each named list
    element represents one parameter as specified in the family object.}
  \item{control}{Further control arguments as specified within the call
    of \code{\link{gamlss2}}.}
}

\details{
  The function \code{select_gamlss2} selects model terms by identifying those
  with estimated degrees of freedom greater than a pre-specified threshold
  (e.g., \code{thres = 0.9}). In addition, model term selection can also be based
  on the percentage of the total predictor range covered by the model term to
  ensure that the term represents a substantial portion of the total effect.
  After the selection step, the model is refitted using only the selected terms,
  excluding the additional penalty applied during the initial fitting process.
  The additional penalty for automatic term selection is described in
  \code{?gam.selection}. Please note that this is experimental, and careful
  consideration should always be given to the modeling process.
}

\value{
  An object of class \code{"gamlss2"}.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{new_formula}}.
}

\examples{\dontrun{## simulate data
set.seed(123)

## number of observations
n <- 1000

## covariates
k <- 100
x <- matrix(runif(n * k, -3, 3), ncol = k)
colnames(x) <- paste0("x", 1:k)
d <- as.data.frame(x)

## true effects
d$f1 <- sin(d$x1)
d$f2 <- exp(d$x2)/15 - 1
d$f3 <- -d$x3 / 3
d$f4 <- d$x4^2/5 - 1 
d$f5 <- cos(d$x5)

## true parameters
mu <- with(d, f1 + f3 + f4)
sigma <- with(d, 1.5 + f2 + f4 + f5)

## simulate response
d$y <- rnorm(n, mean = mu, sd = exp(sigma))

## model formula
f <- paste("~", paste0("s(x", 1:k, ")", collapse = "+"))
f <- as.formula(f)
f <- list(f, f)
f[[1]] <- update(f[[1]], y ~ .)

## estimate model
b <- select_gamlss2(f, data = d, family = NO)

## plot selected estimated effects
plot(b)

## final model
new_formula(b)

## example taken from ?gam.selection
library("mgcv")
set.seed(3)
n <- 200

## simulate data
dat <- gamSim(1, n = n, scale = .15, dist = "poisson")

## spurious
dat$x4 <- runif(n, 0, 1)
dat$x5 <- runif(n, 0, 1)

## formula
f <- y ~ s(x0) + s(x1) + s(x2) + s(x3) + s(x4) + s(x5)

## estimate model
b1 <- gam(f, data = dat, family = poisson,
  select = TRUE, method = "REML")
summary(b1)
plot(b1, pages = 1)

## same with gamlss2
b2 <- select_gamlss2(f, data = dat, family = PO)

## plot selected effects
plot(b2)

## final model
new_formula(b2)
}}

\keyword{regression}
\keyword{models}


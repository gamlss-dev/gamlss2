\name{special_terms} 
\alias{special_terms}
\alias{special_fit}
\alias{special_predict}
\alias{specials}

\title{Special Model Terms for GAMLSS}

\description{
  The \pkg{gamlss2} package provides infrastructure to include special model terms
  for the optimizer functions \code{\link{RS}} and \code{\link{CG}},
  e.g., such as neural networks, trees and forests. The infrastructure assumes
  that such special model terms provide their own fitting and predict method.
}

\usage{
## Generic fitting method.
special_fit(x, \dots)

## Generic predict method.
special_predict(x, \dots)

## Extractor function for fitted special terms.
specials(object, model = NULL, terms = NULL, elements = NULL, ...)
}

\arguments{
  \item{x}{A model term object as supplied in the formula in the \code{\link{gamlss2}} call.}
  \item{object}{A fitted \code{\link{gamlss2}} object.}
  \item{model}{Character or integer, specifies the model for which fitted special terms should be extracted.}
  \item{terms}{Character or integer, specifies the special model terms that should be extracted.}
  \item{elements}{Character, specifies which elements of a fitted special term should be
    extracted. If \code{elements = "names"}, the corresponding element names are extracted.}
  \item{\dots}{Arguments needed for the \code{special_fit()} function to facilitate the
    fitting of the model term, see the details. Similarly, for the \code{special_predict()} function,
    the \code{\dots} argument encompasses the objects for computing predictions for the model term.}
}

\details{
  To implement a new special term, the first step is to write a formula constructor function for
  the new model term. For example, consider the implementation below, which demonstrates how to
  create a neural network model term. Additionally, the name of the new model term constructor
  must be passed to the \code{specials} argument of the function \code{\link{fake_formula}}.
  Please note that in the provided example, no new special name is passed because \code{"n"}
  is already registered in \code{\link{fake_formula}}.

  Afterwards, a fitting and a predict method for the new special model term needs to be implemented.
  Please also refer to the example below, implementing these functions for a neural network
  model term.

  The following describes the detailed arguments and return values.

  A method for \code{special_fit()} has the following arguments:
  \itemize{
    \item \code{x}: The special model term object, containing all the data for fitting.
    \item \code{z}: The current working response/residual from the backfitting step.
    \item \code{w}: The current working weights from the backfitting step.
    \item \code{y}: The response vector/matrix, e.g., used to evaluate the log-likelihood.
    \item \code{eta}: The current named list of predictors.
    \item \code{j}: Character, the parameter name for which the model term needs to be updated.
    \item \code{family}: The family object of the model, see \code{\link{gamlss2.family}}.
    \item \code{control}: A named list of control arguments, see \code{\link{gamlss2_control}}.
  }
  Note that for setting up a special model term only the first three arguments a mandatory,
  all other arguments are optional. The function must at least return a named list containing
  the \code{"fitted.values"} to work with \code{\link{RS}} and \code{\link{CG}}.

  A method for \code{special_predict()} has the following arguments:
  \itemize{
    \item \code{x}: Depending on the return value of function \code{special_fit()}, the
      fitted model term object, see the examples.
    \item{data}: The data for which predictions should be computed.
    \item{se.fit}: Logical, should standard errors of the predictions be computed.
  }
  Note that function \code{special_predict()} should return a data frame with named
  colums \code{"fit"}, \code{"lower"} and \code{"upper"}, \code{"lower"} and \code{"upper"}
  are optional.
}

\seealso{
  \code{\link{gamlss2}}, \code{\link{RS}}, \code{\link{gamlss2_control}}, \code{\link{gamlss2.family}}
}

\examples{\dontshow{ if(!requireNamespace("gamlss.data")) { if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) { stop("not all packages required for the example are installed") } else q() }}
## example special term for neural networks,
## the constructor function is used in the formula
## when calling gamlss2()
n <- function(formula, ...)
{
  stopifnot(requireNamespace("nnet"))

  ## list for setting up the special model term
  st <- list()

  ## list of control arguments
  ctr <- list(...)
  if(is.null(ctr$size))
    ctr$size <- 50
  if(is.null(ctr$maxit))
    ctr$maxit <- 1000
  if(is.null(ctr$decay))
    ctr$decay <- 0.1
  if(is.null(ctr$trace))
    ctr$trace <- FALSE
  if(is.null(ctr$MaxNWts))
    ctr$MaxNWts <- 10000
  if(is.null(ctr$scale))
    ctr$scale <- TRUE

  ## put all information together
  st$control <- ctr
  st$formula <- formula
  st$term <- all.vars(formula)
  st$label <- paste0("n(", paste0(gsub(" ", "", as.character(formula)), collapse = ""), ")")
  st$data <- model.frame(formula)

  ## scale per default!
  if(ctr$scale) {
    sx <- list()
    for(j in colnames(st$data)) {
      if(!is.factor(st$data[[j]])) {
        sx[[j]] <- range(st$data[[j]])
        st$data[[j]] <- (st$data[[j]] - sx[[j]][1]) / diff(sx[[j]])
      }
    }
    st$scalex <- sx
  }

  ## assign the "special" class and the new class "n"
  class(st) <- c("special", "n")

  return(st)
}

## set up the special "n" model term fitting function
special_fit.n <- function(x, z, w, control, ...)
{
  ## model formula needs to be updated
  .fnns <- update(x$formula, response_z ~ .)

  ## assign current working response
  x$data$response_z <- z
  x$data$weights_w <- w

  ## possible weights from last iteration
  Wts <- list(...)$transfer$Wts

  ## estimate model
  nnc <- parse(text = paste0('nnet::nnet(formula = .fnns, data = x$data, weights = weights_w,',
      'size = x$control$size, maxit = x$control$maxit, decay = x$control$decay,',
      'trace = x$control$trace, MaxNWts = x$control$MaxNWts, linout = TRUE',
      if(!is.null(Wts)) ', Wts = Wts)' else ')'))

  rval <- list("model" = eval(nnc))

  ## get the fitted.values
  rval$fitted.values <- predict(rval$model)

  ## transferring the weights for the next backfitting iteration
  ## note, "transfer" can be used to transfer anything from one
  ## iteration to the next
  rval$transfer <- list("Wts" = rval$model$wts)

  ## center fitted values
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift

  ## degrees of freedom
  rval$edf <- length(coef(rval$model))

  ## possible scaling
  rval$scalex <- x$scalex

  ## assign class for predict method
  class(rval) <- "n.fitted"

  return(rval)
}

## finally, the predict method
special_predict.n.fitted <- function(x, data, se.fit = FALSE, ...)
{
  if(!is.null(x$scalex)) {
    for(j in names(x$scalex)) {
      data[[j]] <- (data[[j]] - x$scalex[[j]][1]) / diff(x$scalex[[j]])
    }
  }
  p <- predict(x$model, newdata = data, type = "raw")
  p <- p - x$shift
  if(se.fit)
    p <- data.frame("fit" = p)
  return(p)
}

\dontrun{## example with data
data("abdom", package = "gamlss.data")

## specify the model Formula
f <- y ~ n(~x) | n(~x) | n(~x) | n(~x)

## estimate model,
## set the seed for reproducibility
## note, data should be scaled!
set.seed(123)
b <- gamlss2(f, data = abdom, family = BCT)

## visualize estimated effects
plot(b, which = "effects")

## plot diagnostics
plot(b, which = "resid")

## predict parameters
par <- predict(b)

## predict quantiles
pq <- sapply(c(0.05, 0.5, 0.95), function(q) family(b)$q(q, par))

## plot
plot(y ~ x, data = abdom, pch = 19,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
matplot(abdom$x, pq, type = "l", lwd = 2,
  lty = 1, col = 4, add = TRUE)

## another example using the Munich rent data
data("rent", package = "gamlss.data")

## model Formula
f <- R ~ n(~Fl+A,size=10,decay=0.7) | n(~Fl+A,size=10,decay=0.7)

## estimate model
set.seed(456)
b <- gamlss2(f, data = rent, family = GA)

## plot estimated effects
plot(b, which = "effects", persp = FALSE)

## diagnostics
plot(b, which = "resid")

## predict using new data
n <- 50
nd <- with(rent, expand.grid(
  "Fl" = seq(min(Fl), max(Fl), length = n),
  "A" = seq(min(A), max(A), length = n)
))

## predict parameters of the GA distribution
par <- predict(b, newdata = nd)

## compute median rent R estimate
nd$fit <- family(b)$q(0.5, par)

## visualize
library("lattice")

p1 <- wireframe(fit ~ Fl + A, data = nd,
  screen = list(z = 50, x = -70, y = -10),
  aspect = c(1, 0.9), drape = TRUE,
  main = "n(~Fl+A)",
  xlab = "Floor", ylab = "YoC",
  zlab = "Rent")

p2 <- levelplot(fit ~ Fl + A, data = nd,
  contour = TRUE,
  main = "n(~Fl+A)", xlab = "Floor", ylab = "YoC")

print(p1, split = c(1, 1, 2, 1), more = TRUE)
print(p2, split = c(2, 1, 2, 1), more = FALSE)

## extract fitted special terms,
## fitted NN for parameter mu
specials(b, model = "mu", elements = "model")

## same for sigma
specials(b, model = "sigma", elements = "model")

## return element names of fitted special term list
specials(b, model = "sigma", elements = "names")
}}

\keyword{models}
\keyword{regression}


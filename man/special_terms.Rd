\name{special_terms} 
\alias{special_terms}
\alias{special_fit}
\alias{special_predict}

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
}

\arguments{
  \item{x}{A model term object as supplied in the formula in the \code{\link{gamlss2}} call.}
  \item{\dots}{Arguments needed for the \code{special_fit()} function to facilitate the
    fitting of the model term, see the details. Similarly, for the \code{special_predict()} function,
    the \code{\dots} argument encompasses the objects for computing predictions for the model term.}
}

\details{
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

\examples{## Example special term for neural networks.
## The constructor function is used in the formula
## when calling gamlss2().
n <- function(formula, ...)
{
  stopifnot(requireNamespace("nnet"))

  ## List for setting up the special model term.
  st <- list()

  ## List of control arguments.
  ctr <- list(...)
  if(is.null(ctr$size))
    ctr$size <- 100
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

  ## Put all information together.
  st$control <- ctr
  st$formula <- formula
  st$term <- all.vars(formula)
  st$label <- paste0("n(", paste0(gsub(" ", "", as.character(formula)), collapse = ""), ")")
  st$data <- model.frame(formula)

  ## Scale per default!
  if(ctr$scale) {
    sx <- list()
    for(j in colnames(data)) {
      if(!is.factor(data[[j]])) {
        sx[[j]] <- range(data[[j]])
        data[[j]] <- (data[[j]] - sx[[j]][1]) / diff(sx[[j]])
      }
    }
    st$scalex <- sx
  }

  ## Assign the "special" class and the new class "n".
  class(st) <- c("special", "n")

  return(st)
}

## Set up the special "n" model term fitting function
special_fit.n <- function(x, z, w, control, ...)
{
  ## Model formula needs to be updated.
  .fnns <- update(x$formula, response_z ~ .)

  ## Assign current working response.
  x$data$response_z <- z
  x$data$weights_w <- w

  ## Possible weights from last iteration.
  Wts <- list(...)$transfer$Wts

  ## Estimate model.
  nnc <- parse(text = paste0('nnet::nnet(formula = .fnns, data = x$data, weights = weights_w,',
      'size = x$control$size, maxit = x$control$maxit, decay = x$control$decay,',
      'trace = x$control$trace, MaxNWts = x$control$MaxNWts, linout = TRUE',
      if(!is.null(Wts)) ', Wts = Wts)' else ')'))

  rval <- list("model" = eval(nnc))

  ## Get the fitted.values.
  rval$fitted.values <- predict(rval$model)

  ## Transferring the weights for the next backfitting iteration.
  ## Note, "transfer" can be used to transfer anything from one
  ## iteration to the next.
  rval$transfer <- list("Wts" = rval$model$wts)

  ## Center fitted values.
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift

  ## Degrees of freedom.
  rval$edf <- x$control$size

  ## Possible scaling.
  rval$scalex <- x$scalex

  ## Assign class for predict method.
  class(rval) <- "n.fitted"

  return(rval)
}

## Finally, the predict method.
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

\dontrun{## Example with data.
data("abdom", package = "gamlss.data")

## Specify the model Formula.
f <- y ~ n(~x) | n(~x) | n(~x) | n(~x)

## Estimate model.
## Set the seed for reproducibility.
## Note, data should be scaled!
set.seed(123)
b <- gamlss2(f, data = abdom, family = BCT)

## Plot estimated effects.
plot(b, which = "effects")

## Plot diagnostics.
plot(b, which = "resid")

## Predict parameters.
par <- predict(b)

## Predict quantiles.
pq <- NULL
for(q in c(0.05, 0.5, 0.95))
  pq <- cbind(pq, family(b)$q(q, par))

## Plot.
plot(y ~ x, data = abdom, pch = 19,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
matplot(abdom$x, pq, type = "l", lwd = 2,
  lty = 1, col = 4, add = TRUE)

## Another example using the Munich rent data.
data("rent", package = "gamlss.data")

## Model formula.
f <- R ~ n(~Fl+A,size=10,decay=0.7) | n(~Fl+A,size=10,decay=0.7)

## Estimate model.
set.seed(456)
b <- gamlss2(f, data = rent, family = GA)

## Plot estimated effects.
plot(b, which = "effects", persp = FALSE)

## Diagnostics.
plot(b, which = "resid")

## Predict using new data.
n <- 50
nd <- with(rent, expand.grid(
  "Fl" = seq(min(Fl), max(Fl), length = n),
  "A" = seq(min(A), max(A), length = n)
))

## Predict parameters of the GA distribution.
par <- predict(b, newdata = nd)

## Compute median rent R estimate.
nd$fit <- family(b)$q(0.5, par)

## Visualize.
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
}}

\keyword{models}
\keyword{regression}


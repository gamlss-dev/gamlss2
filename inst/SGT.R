## Install the gamlss2 package.
if(!require("gamlss2")) {
  devtools::install_github("gamlss-dev/gamlss2")
}

library("gamlss2")

## Skewed Generalized T Distribution.
SGT <- function(...) {
  stopifnot(requireNamespace("sgt"))

  fam <- list(
    "names" = c("mu", "sigma", "lambda", "p", "q"),
    "links" = c("mu" = "identity", "sigma" = "log",
       "lambda" = "rhogit", "p" = "log", "q" = "log"),
    "d" = function(y, par, log = FALSE, ...) {
      sgt::dsgt(y, mu = par$mu, sigma = par$sigma, lambda = par$lambda,
        p = par$p, q = par$q, log = log)
    },
    "p" = function(y, par, lower.tail = TRUE, log.p = FALSE) {
      sgt::psgt(y, mu = par$mu, sigma = par$sigma,
        lambda = par$lambda, p = par$p, q = par$q,
        lower.tail = lower.tail, log.p = log.p)
    },
    "q" = function(p, par, lower.tail = TRUE, log.p = FALSE) {
      sgt::qsgt(p, mu = par$mu, sigma = par$sigma,
        lambda = par$lambda, p = par$p, q = par$q,
        lower.tail = lower.tail, log.p = log.p)
    }
  )

  class(fam) <- "gamlss2.family"

  return(fam)
}

## Load the abdominal circumference data
data("abdom", package = "gamlss.data")

## Specify the model Formula.
f <- y ~ s(x) | s(x) | s(x) | s(x) | s(x)

## Estimate model.
b <- gamlss2(f, data = abdom, family = SGT,
  maxit = c(100, 100), eps = 1e-08)

## Plot estimated effects.
plot(b, which = "effects")

## Plot diagnostics.
plot(b, which = "resid")

## Predict parameters.
par <- predict(b, type = "parameter")

## Predict quantiles.
pq <- NULL
for(q in c(0.05, 0.5, 0.95))
  pq <- cbind(pq, family(b)$q(q, par))

## Plot.
plot(y ~ x, data = abdom, pch = 19,
  col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
matplot(abdom$x, pq, type = "l", lwd = 2,
  lty = 1, col = 4, add = TRUE)


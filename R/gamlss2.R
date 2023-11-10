## Generic gamlss method.
gamlss2 <- function(x, ...)
{
  UseMethod("gamlss2")
}

## Formula method.
gamlss2.formula <- function(formula, data, family = NO,
  subset, na.action, weights, offset,
  model = TRUE, x = TRUE, y = TRUE,
  control = gamlss2.control(...), ...)
{
  ## Process specific formulas.
  if(!is.null(control$sigma.formula) | !is.null(control$nu.formula) | !is.null(control$tau.formula)) {
    if(!inherits(formula, "list")) {
      formula <- list("mu" = formula)
    }
    formula[["sigma"]] <- if(is.null(control$sigma.formula)) ~1 else control$sigma.formula
    formula[["nu"]] <- if(is.null(control$nu.formula)) ~1 else control$nu.formula
    formula[["tau"]] <- if(is.null(control$tau.formula)) ~1 else control$tau.formula
    formula <- formula[c("mu", "sigma", "nu", "tau")]
    names(formula) <- NULL
  }

  ## Call.
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## Formula.
  if(is.list(formula)) {
    formula <- do.call("as.Formula", formula)
  }
  formula <- as.Formula(formula)
  if(is.list(formula)) {
    formula <- do.call("as.Formula", formula)
  }
  if(length(formula)[2L] < 2L) {
    formula <- as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  } else {
    simple_formula <- FALSE
  }

  mf$formula <- fake_formula(formula)

  ## Evaluate model.frame.
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  rval <- list("call" = call)

  ## Response and model.matrix.
  ff <- fake_formula(formula, nospecials = TRUE)
  mt <- terms(ff, data = data)
  Y <- model.response(mf)
  X <- model.matrix(mt, mf)

  if(any(is.na(X)))
    stop("detected 'NA' values in data!")
  if(any(!is.finite(X)))
    stop("detected 'Inf' values in data!")

  ## Model fitting here!

  ## Further model information.
  rval$formula <- formula
  rval$family <- family
  rval$terms <- mt
  rval$xlevels <- .getXlevels(mt, mf)

  if(!control$light) {
    if(model) {
      rval$model <- mf
    }
    if(y) {
      rval$y = Y
    }
    if(x) {
      rval$x <- X
    }
  }

  class(rval) <- "gamlss2"

  return(rval)
}

## List method.
gamlss2.list <- function(formula, ...)
{
  gamlss2.formula(formula, ...)
}

## Default method.
gamlss2.default <- function(x, y, specials, family = NO,
  variables = NULL, weights = NULL, offset = NULL,
  control = gamlss2.control(...), ...)
{

}

## Control parameters.
gamlss2.control <- function(c.crit = 0.001, n.cyc = 20,
  mu.step = 1, sigma.step = 1, nu.step = 1, 
  tau.step = 1, gd.tol = Inf, iter = 0,
  trace = TRUE, autostep = TRUE, 
  save = TRUE, light = FALSE, ...)
{
  ctr <- as.list(environment())
  ctr <- c(ctr, list(...))

  if(c.crit <= 0) {
    warning("the value of c.crit supplied is zero or negative the default value of 0.001 was used instead")
    c.crit <- 0.001
  }
   if(n.cyc < 1) {
    warning("the value of no cycles supplied is zero or negative the default value of 20 was used instead")
    n.cyc <- 20
  }
   if(iter < 0) {
    warning("the value of no iterations  supplied is  negative the default value of 0 was used instead")
    iter <- 0
  }
   if(mu.step > 1 | mu.step < 0) {
    warning("the value of mu.step supplied is less than zero or more than one the default value of 1 was used instead")
    mu.step <- 1
  }
   if(sigma.step > 1 | sigma.step < 0) {
    warning("the value of sigma.step supplied is less than zero or more than one the default value of 1 was used instead")
    sigma.step <- 1
  }
   if(nu.step > 1 | nu.step < 0) {
    warning("the value of nu.step supplied is less than zero or more than one the default value of 1 was used instead")
    nu.step <- 1
  }
   if(tau.step > 1 | tau.step < 0) {
    warning("the value of tau.step supplied is less than zero or more than one the default value of 1 was used instead")
    tau.step <- 1
  }
   if(gd.tol < 0) {
    warning("the value of gd.tol supplied is less than zero the default value of Inf was used instead")
    gd.tol <- Inf
  }

  return(ctr)
}

## Testing.
if(FALSE) {
  library("Formula")
  library("gamlss.dist")
  source("gamlss2.R")
  source("fake_formula.R")

  d <- bamlss::GAMart()

  f <- list(num ~ x1 + x2, ~ x3)  

  b1 <- gamlss2(f, data = d)
  b2 <- gamlss2(num ~ x1 + x2, sigma.formula = ~ x3, data = d)

  print(head(b1$model))
  print(head(b2$model))
}


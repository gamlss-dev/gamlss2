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
  m <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
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

  ## Process weights and offsets.
  weights <- model.weights(mf)
  if(!is.null(weights)) {
    if(length(weights) == 1) 
      weights <- rep.int(weights, nrow(mf))
    weights <- as.vector(weights)
    names(weights) <- rownames(mf)
  }

  expand_offset <- function(offset) {
    off <- NULL
    if(!is.null(offset)) {
      if(length(offset) == 1) 
        offset <- rep.int(offset, n)
      off <- as.vector(offset)
    }
    return(off)
  }

  ## Process variables and special term information.
  Xterms <- offsets <- list()
  for(i in 1:(length(ff)[2])) {
    Xterms[[i]] <- attr(terms(formula(ff, rhs = i, lhs = 0)), "term.labels")
    if(attr(terms(formula(ff, rhs = i, lhs = 0)), "intercept") > 0)
      Xterms[[i]] <- c("(Intercept)", Xterms[[i]])
    offi <- expand_offset(model.offset(model.part(ff, 
      data = mf, rhs = i, terms = TRUE)))
    offsets[[i]] <- if(length(offi)) offi else numeric(0)
  }
  Sterms <- fake_formula(formula, onlyspecials = TRUE)

  if(!missing(offset)) {
    anyoff <- any(sapply(offsets, function(x) length(x) > 0))
    if(anyoff)
      stop("multiple offsets supplied, either use argument offset or specify offsets in the formula!")
    cn <- NULL
    if(!is.list(offset) | !is.data.frame(offset)) {
      if(!is.matrix(offset)) {
        offset <- data.frame(offset)
      } else {
        offset <- as.data.frame(offset)
      }
    }
    offset <- as.data.frame(offset)
    if(nrow(offset) < 2)
      offset <- offset[rep(1L, nrow(mf)), , drop = FALSE]
    rownames(offset) <- rownames(mf)
  }

  ## Process special terms.
  Specials <- special_terms(Sterms, mf)

  ## Evaluate and complete family.
  ## Note, families structure is a bit different
  ## in order to support more than 4 parameter models.
  family <- complete_family(family)

  ## Set names.
  names(Xterms) <- names(Sterms) <- family$names
  Xterms0 <- Xterms
  if(length(offsets)) {
    names(offsets) <- family$names
    offsets <- do.call("cbind", offsets)
  }
  if(!missing(offset)) {
    if(!all(colnames(offset) %in% family$names))
      colnames(offset) <- family$names[1:ncol(offset)]
    offsets <- offset
  }

  ## Process factors and other linear model terms.
  xlev <- .getXlevels(mt, mf)
  for(i in names(Xterms)) {
    ## Factors.
    for(j in names(xlev)) {
      if(j %in% Xterms[[i]]) {
        xl <- paste0(j, xlev[[j]])
        xl <- xl[xl %in% colnames(X)]
        Xterms[[i]][Xterms[[i]] == j] <- NA
        Xterms[[i]] <- Xterms[[i]][!is.na(Xterms[[i]])]
        Xterms[[i]] <- c(Xterms[[i]], xl)
      }
    }
    ## Others.
    for(j in Xterms[[i]]) {
      if(!(j %in% colnames(X))) {
        if(any(ij <- grepl(j, colnames(X), fixed = TRUE))) {
          Xterms[[i]][Xterms[[i]] == j] <- NA
          Xterms[[i]] <- Xterms[[i]][!is.na(Xterms[[i]])]
          Xterms[[i]] <- c(Xterms[[i]], colnames(X)[ij])
        }
      }
    }
  }

  ## Optionally, use optimizer function provided from family
  optimizer <- if(is.null(family$optimizer)) {
    control$optimizer
  } else {
    family$optimizer
  }

  ## Estimation.
  rval <- optimizer(X, Y, Specials, family,
    offsets, weights, Xterms, Sterms, control)

  ## Further model information.
  rval$formula <- formula
  rval$family <- family
  rval$terms <- mt
  rval$xlevels <- xlev
  rval$xterms <- Xterms0
  rval$sterms <- Sterms
  rval$specials <- Specials

  ## Return model.frame, X and y.
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

  class(rval) <- c(class(rval), "gamlss2")

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
gamlss2.control <- function(optimizer = RS, step = 1,
  trace = TRUE, flush = TRUE, light = FALSE, ...)
{
  ctr <- as.list(environment())
  ctr <- c(ctr, list(...))

  return(ctr)
}

## Testing.
if(FALSE) {
  d <- bamlss::GAMart()

  f <- list(num ~ x1 + x2, ~ x3)  

  b1 <- gamlss2(f, data = d)
  b2 <- gamlss2(num ~ x1 + x2, sigma.formula = ~ x3, data = d)
  b3 <- gamlss2(num ~ x1 + x2 + s(x3) + te(sqrt(lon),lat), data = d)

  print(head(b1$model))
  print(head(b2$model))
  print(head(b3$model))

  f <- list(
    num ~ fac + x1 + pb(x2) + s(x3) + te(sqrt(lon),lat),
        ~ s(x1) + x2 + pb(sqrt(x3))
  )

  b <- gamlss2(f, data = d)
}


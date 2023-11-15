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
  mf <- call <- match.call(expand.dots = FALSE)
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
  rval$call <- call
  rval$formula <- formula
  rval$fake_formula <- fake_formula(formula)
  rval$terms <- terms(rval$fake_formula)
  rval$family <- family
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

  class(rval) <- unique(c(class(rval), "gamlss2"))

  return(rval)
}

## List method.
gamlss2.list <- function(formula, ...)
{
  gamlss2.formula(formula, ...)
}

## Control parameters.
gamlss2.control <- function(optimizer = RS, step = 1,
  trace = TRUE, flush = TRUE, light = FALSE, ...)
{
  ctr <- as.list(environment())
  ctr <- c(ctr, list(...))

  return(ctr)
}

## A model.frame method.
model.frame.gamlss2 <- function(formula, ...)
{
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  if(length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- quote(model.frame)
    fcall$xlev <- formula$xlevels
    fcall$formula <- terms(formula)
    fcall$formula <- update(fcall$formula, NULL ~ .)
    fcall[names(nargs)] <- nargs
    env <- if(is.null(environment(formula$terms))) {
      parent.frame()
    } else {
      environment(formula$terms)
    }
    return(eval(fcall, env))
  } else {
    return(formula$model)
  }
}

## The model.matrix.
model.matrix.gamlss2 <- function(object, data = NULL, ...)
{
  if(!is.null(data))
    object$x <- NULL
  if(n_match <- match("x", names(object), 0L)) {
    object[[n_match]]
  } else {
    object$terms <- terms(fake_formula(object$formula, nospecials = TRUE))
    data <- model.frame(object, xlev = object$xlevels, data = data, ...)
    if(exists(".GenericCallEnv", inherits = FALSE)) {
      NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
    } else {
      dots <- list(...)
      dots$data <- dots$contrasts.arg <- NULL
      do.call("model.matrix", c(list(object = object, 
        data = data, contrasts.arg = object$contrasts), dots))
    }
  }
}

## Family extractor.
family.gamlss2 <- function(object, ...) object$family

## Multiple grep.
grep2 <- function (pattern, x, ...) 
{
  i <- NULL
  for(p in pattern)
    i <- c(i, grep(p, x, ...))
  unique(i)
}

## Predict method.
predict.gamlss2 <- function(object, 
  parameter = NULL, newdata = NULL, type = c("link", "parameter", "response", "terms"), 
  terms = NULL, se.fit = FALSE, ...)
{
  type <- match.arg(type)

  if(!is.null(newdata)) {
    mf <- model.frame(object, data = newdata)
  } else {
    mf <- model.frame(object)
  }
  X <- model.matrix(fake_formula(object$formula, nospecials = TRUE), data = mf)

  family <- object$family

  if(is.null(parameter)) {
    parameter <- list(...)$what
    if(is.null(parameter))
    parameter <- list(...)$model
    if(is.null(parameter))
      parameter <- family$names
  }
  if(!is.character(parameter))
    parameter <- family$names[parameter]
  parameter <- family$names[pmatch(parameter, family$names)]

  p <- list()
  for(j in parameter) {
    p[[j]] <- rep(0, length.out = nrow(mf))
    tj <- if(is.null(terms)) {
      c(object$xterms[[j]], object$sterms[[j]])
    } else {
      grep2(terms, c(object$xterms[[j]], object$sterms[[j]]), fixed = TRUE, value = TRUE)
    }
    if(length(tj)) {
      if(length(object$xterms[[j]])) {
        xn <- NULL
        for(i in tj) {
          xn <- c(xn, grep(i, object$xterms[[j]], fixed = TRUE, value = TRUE))      
        }
        if(length(xn)) {
          p[[j]] <- p[[j]] + drop(X[, xn, drop = FALSE] %*% coef(object)[[j]][xn])
        }
      }
      if(length(object$sterms[[j]])) {
        xn <- NULL
        for(i in tj) {
          xn <- c(xn, grep(i, object$sterms[[j]], fixed = TRUE, value = TRUE))      
        }
        if(length(xn)) {
          for(i in xn) {
            if(inherits(object$specials[[i]], "mgcv.smooth")) {
              Xs <- PredictMat(object$specials[[i]], data = mf, n = nrow(mf))
              co <- object$fitted.specials[[j]][[i]]$coefficients
              p[[j]] <- p[[j]] + drop(Xs %*% co)
            } else {
              cs <- object$fitted.specials[[j]][[i]]$coefficients
              p[[j]] <- p[[j]] + cs$fun(mf[[cs$name]])
            }
          }
        }
      }
    }
  }

  if(type %in% c("parameter", "response")) {
    p <- family$map2par(p)
  }

  if(type == "response") {
    fm <- family$mean
    if(is.null(fm)) {
      if(!is.null(family$q)) {
        fm <- function(par) family$q(0.5, par)
      }
    }
    p <- if(is.null(fm)) p[[1L]] else fm(p)
  }

  if(is.list(p)) {
    if(length(p) < 2) {
      p <- p[[1L]]
    } else {
      p <- as.data.frame(p)
    }
  }

  return(p)
}

## Testing.
if(FALSE) {
  d <- bamlss::GAMart(n = 10000, sd = -1)

  f <- list(
    y ~ s(x1) + pb(x2) + s(x3) + te(lon,lat,k=10),
      ~ s(x1) + pb(x2)
  )

  b <- gamlss2(f, data = d, maxit = c(100, 100))

  fx2 <- predict(b, newdata = d[1:100, ], parameter = "sigma", term = "x2")
  plot(fx2 ~ d$x2[1:100])
}

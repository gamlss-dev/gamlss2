## Generic gamlss method.
gamlss2 <- function(x, ...)
{
  UseMethod("gamlss2")
}

## Formula method.
gamlss2.formula <- function(formula, data, family = NO,
  subset, na.action, weights, offset,
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

  ## Evaluate and complete family.
  ## Note, families structure is a bit different
  ## in order to support more than 4 parameter models.
  family <- complete_family(family)

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

  ## Expand formula.
  if((length(attr(formula, "rhs")) < length(family$names)) & control$expand) {
    k <- length(family$names) - length(attr(formula, "rhs"))
    attr(formula, "rhs") <- c(attr(formula, "rhs"), as.list(rep(1, k)))
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
  if(!("(Intercept)" %in% colnames(X)))
    X <- cbind("(Intercept)" = 1.0, X)

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
        offset <- rep.int(offset, nrow(X))
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

  ## Set names.
  names(Xterms) <- family$names[1:length(Xterms)]
  names(Sterms) <- family$names[1:length(Sterms)]
  Xterms0 <- Xterms
  if(length(offsets)) {
    names(offsets) <- family$names[1:length(offsets)]
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
    if(control$model) {
      rval$model <- mf
    }
    if(control$y) {
      rval$y = Y
    }
    if(control$x) {
      rval$x <- X
    }
    rval$results <- results(rval)
  }

  class(rval) <- unique(c(class(rval), "gamlss2"))

  return(rval)
}

## List method.
gamlss2.list <- function(x, ...)
{
  cl <- match.call()
  cl$formula <- do.call("as.Formula", x)
  cl$x <- FALSE
  cl[[1L]] <- as.name("gamlss2.formula")
  eval.parent(cl)
}

## Control parameters.
gamlss2.control <- function(optimizer = RS, step = 1,
  trace = TRUE, flush = TRUE, light = FALSE, expand = TRUE,
  model = TRUE, x = TRUE, y = TRUE, ...)
{
  ctr <- as.list(environment())
  ctr <- c(ctr, list(...))

  return(ctr)
}

## A model.frame method.
model.frame.gamlss2 <- function(formula, ...)
{
  dots <- list(...)
  if(is.null(dots$keepresponse))
    dots$keepresponse <- FALSE
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  if(length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- quote(model.frame)
    fcall$xlev <- formula$xlevels
    fcall$formula <- terms(formula)
    if(!dots$keepresponse)
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
    return(object[[n_match]])
  } else {
    object$terms <- terms(fake_formula(object$formula, nospecials = TRUE))
    data <- model.frame(object, xlev = object$xlevels, data = data, ...)
    if(exists(".GenericCallEnv", inherits = FALSE)) {
      X <- NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
    } else {
      dots <- list(...)
      dots$data <- dots$contrasts.arg <- NULL
      X <- do.call("model.matrix", c(list(object = object, 
        data = data, contrasts.arg = object$contrasts), dots))
    }
    if(!("(Intercept)" %in% colnames(X)))
      X <- cbind("(Intercept)" = 1.0, X)
    return(X)
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

  tt <- type == "terms"

  p <- list()
  for(j in parameter) {
    p[[j]] <- if(tt) NULL else rep(0, length.out = nrow(mf))
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
          for(i in seq_along(xn)) {
            if(!is.null(object$xlevels)) {
              if(xn[i] %in% names(object$xlevels)) {
                xnl <- paste0(xn[i], object$xlevels[[xn[i]]])
                xnl <- xnl[xnl %in% colnames(X)]
                xn <- c(xn[-i], xnl)
              }
            }
          }
          if(tt) {
            ft <- t(t(X[, xn, drop = FALSE]) * coef(object)[[j]][xn])
            p[[j]] <- cbind(p[[j]], ft)
          } else {
            p[[j]] <- p[[j]] + drop(X[, xn, drop = FALSE] %*% coef(object)[[j]][xn])
          }
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
              fit <- drop(Xs %*% co)
            } else {
              cs <- object$fitted.specials[[j]][[i]]$coefficients
              if(inherits(cs, "random")) {
                vn <- as.character(as.call(as.call(parse(text = i))[[1L]])[[2L]])
                xv <- mf[[vn]]
                fit <- cs$coef[as.character(xv)]
              } else {
                fit <- cs$fun(mf[[cs$name]])
              }
            }
            if(tt) {
              fit <- matrix(fit, ncol = 1L)
              colnames(fit) <- i
              rownames(fit) <- rownames(mf)
              p[[j]] <- cbind(p[[j]], fit)
            } else {
              p[[j]] <- p[[j]] + fit
            }
          }
        }
      }
    }
  }

  if(type %in% c("parameter", "response") & !tt) {
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
      if(!tt)
        p <- as.data.frame(p)
    }
  }

  return(p)
}

## Residuals.
residuals.gamlss2 <- function(object, type = c("quantile", "response", "parameter"), newdata = NULL, ...)
{
  family <- family(object)

  if(!is.null(family$residuals)) {
    res <- family$residuals(object, type = type, ...)
    if(length(class(res)) < 2) {
      if(inherits(res, "numeric"))
        class(res) <- c("gamlss2.residuals", class(res))
    }
  } else {
    type <- match.arg(type)

    object$model <- NULL

    object$y <- model.response(model.frame(object, keepresponse = TRUE, data = newdata))

    if(is.null(object$y))
      stop("response variable is missing, cannot compute residuals!")

    nobs <- nrow(object$y)
    y <- if(is.data.frame(object$y)) {
      if(ncol(object$y) < 2) {
        object$y[[1]]
      } else object$y
    } else {
      object$y
    }

    par <- predict(object, drop = FALSE, newdata = newdata, ...)
    for(j in family$names)
      par[[j]] <- make.link2(family$links[j])$linkinv(par[[j]])

    if(type == "quantile") {
      if(is.null(family$p)) {
        type <- "response"
        warning(paste("no $p() function in family '", family$family,
          "', cannot compute quantile residuals, computing response resdiuals instead!", sep = ""))
      } else {
        discrete <- FALSE
        if(!is.null(family$type)) {
          if(tolower(family$type) == "discrete")
            discrete <- TRUE
        }
        if(family$family == "binomial")
          discrete <- TRUE
        if(discrete) {
          ymin <- min(y, na.rm = TRUE)
          a <- family$p(ifelse(y == ymin, y, y - 1), par)
          a <- ifelse(y == ymin, 0, a)
          b <- family$p(y, par)
          u <- runif(length(y), a, b)
          u <- ifelse(u > 0.999999, u - 1e-16, u)
          u <- ifelse(u < 1e-06, u + 1e-16, u)
          res <- qnorm(u)
        } else {
          prob <- family$p(y, par)
          thres <- 0.999999999999999
          prob[prob > thres] <- thres
          prob[prob < (1 - thres)] <- 1 - thres
          res <- qnorm(prob)
          if(any(isnf <- !is.finite(res))) {
            warning("non finite quantiles from probabilities, set to NA!")
            res[isnf] <- NA
          }
        }
        attr(res, "type") <- "Quantile"
      }
    }

    if(type == "response") {
      mu <- if(is.null(family$mu)) {
        function(par, ...) { par[[1]] }
      } else family$mu
      res <- y - mu(par)
      attr(res, "type") <- "Response"
    }

    if(type == "parameter") {
      stop("Not implemented yet!")
    }
   
    class(res) <- c("gamlss2.residuals", class(res))
  }

  if(any(j <- !is.finite(res)))
    res[j] <- NA

  return(res)
}


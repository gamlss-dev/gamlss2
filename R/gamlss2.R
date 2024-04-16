## Generic gamlss method.
gamlss2 <- function(x, ...)
{
  UseMethod("gamlss2")
}

## Formula method.
gamlss2.formula <- function(formula, data, family = NO,
  subset, na.action, weights, offset, start = NULL,
  control = gamlss2_control(...), ...)
{
  ## Save environments.
  menv <- parent.frame()
  fenv <- environment(formula)

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
    environment(formula) <- fenv
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

  ## Check for "." in formula.
  for(i in 1:length(formula)[2L]) {
    rhs <- formula(formula, rhs = i)
    if(as.character(rhs[3]) == ".") {
      if(!inherits(data, "environment")) {
        yn <- NULL
        for(j in 1:length(formula)[1L])
          yn <- c(yn, as.character(formula(formula, lhs = j))[2L])
        vn <- names(data)
        vn <- vn[!(vn %in% yn)]
        attr(formula, "rhs")[[i]] <- as.call(str2lang(paste(vn, collapse = "+")))
      } else {
        stop('using "." in formula but no data argument supplied!')
      }
    }
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
  Specials <- special_terms(Sterms, mf, binning = control$binning,
    digits = control$digits, select = control$select)

  ## Process by variables using mgcv::smoothCon().
  olab <- sapply(Specials, function(x) if(is.list(x)) x$orig.label else "")
  nt <- names(olab)
  ulab <- unique(olab)
  ulab <- ulab[ulab != ""]
  for(j in seq_along(Sterms)) {
    if(length(Sterms[[j]])) {
      for(i in ulab) {
        ii <- which(Sterms[[j]] == i)
        if(length(ii)) {
          Sterms[[j]] <- as.list(Sterms[[j]])
          uti <- unique(as.character(nt[olab == i]))
          Sterms[[j]][[ii]] <- uti
          Sterms[[j]] <- unlist(Sterms[[j]])
        }
      }
    }
  }

  ## Set names.
  names(Xterms) <- family$names[1:length(Xterms)]
  names(Sterms) <- family$names[1:length(Sterms)]
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
  rval <- optimizer(x = X, y = Y, specials = Specials, family = family,
    offsets = offsets, weights = weights, start = start, xterms = Xterms, sterms = Sterms,
    control = control)

  ## Further model information.
  rval$call <- call
  rval$formula <- formula
  rval$fake_formula <- fake_formula(formula)
  rval$terms <- terms(merge_formula(formula(rval$fake_formula, collapse = TRUE), as.formula(mt)))
  environment(rval$terms) <- menv
  rval$family <- family
  rval$xlevels <- xlev
  rval$contrasts <- attr(X, "contrasts")
  rval$na.action <- attr(mf, "na.action")
  attr(Xterms, "terms") <- mt
  rval$xterms <- Xterms
  rval$sterms <- Sterms
  rval$specials <- Specials
  rval$df <- get_df(rval)

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
  } else {
    rval$fitted.values <- NULL
    if(!is.null(rval$fitted.linear)) {
      for(j in names(rval$fitted.linear))
        rval$fitted.linear[[j]]$fitted.values <- NULL
    }
    if(!is.null(rval$specials)) {
      for(j in names(rval$specials)) {
        if(!is.null(rval$specials[[j]][["X"]])) {
          rval$specials[[j]][["X"]] <- NULL
          rval$specials[[j]][["Xu"]] <- NULL
        }
      }
    }
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
gamlss2_control <- function(optimizer = RS,
  trace = TRUE, flush = TRUE, light = FALSE, expand = TRUE,
  model = TRUE, x = TRUE, y = TRUE, fixed = FALSE, ...)
{
  ctr <- as.list(environment())
  ctr <- c(ctr, list(...))

  if(is.null(ctr$binning))
    ctr$binning <- FALSE
  if(is.null(ctr$digits))
    ctr$digits <- Inf
  if(is.null(ctr$initialize))
    ctr$initialize <- FALSE
  if(is.null(ctr$nullmodel))
    ctr$nullmodel <- TRUE

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
    if(is.null(data)) {
      data <- model.frame(object, xlev = object$xlevels, ...)
    }
    dots <- list(...)
    dots$data <- dots$contrasts.arg <- NULL
    mt <- terms(update(terms(object), NULL ~ .))
    X <- do.call(stats::model.matrix.default, c(list(object = list("terms" = mt), 
      data = data, contrasts.arg = object$contrasts), dots))
    if(!("(Intercept)" %in% colnames(X)))
      X <- cbind("(Intercept)" = 1.0, X)
    return(X)
  }
}

## Family extractor.
family.gamlss2 <- function(object, ...) object$family

## A simple printing method.
print.gamlss2 <- function(x, ...)
{
  x$call[[1]] <- as.name("gamlss2")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("---\n")
  print(x$family, full = FALSE)
  cat("*--------\n")
  info1 <- c(
    paste("n =", x$nobs),
    paste("df = ", round(x$df, digits = 2)),
    paste("res.df = ", round(x$nobs - x$df, digits = 2))
  )
  info2 <- c(
    paste("logLik =", round(x$logLik, digits = 4)),
    paste("Deviance =", round(-2 * x$logLik, digits = 4)),
    paste("AIC =", round(-2 * x$logLik + 2*x$df, digits = 4))
  )
  cat(info1)
  cat("\n")
  cat(info2)
  cat("\n")
}

## Merging formulas.
is_formula <- function(x) inherits(x, "formula")

merge_formula <- function(x, y, ...)
{
  if(!is_formula(x) || length(x) != 3)
    stop("First argument is invalid")
  if(!is_formula(y)) stop("Second argument is invalid")
  if(length(list(...))) warning("extraneous arguments discarded")
  is.gEnv <- function(e) identical(e, .GlobalEnv)

  str <- paste(c(deparse(x[[2]]), "~",
    deparse(x[[3]]), "+",
    deparse(y[[length(y)]])), collapse = "")
  f <- as.formula(str)
  ex <- environment(x)
  ey <- environment(y)
  if(!is.gEnv(ex)) {
      environment(f) <- ex
      if(!is.gEnv(ey) && !identical(ex,ey)) {
          warning("`x' and `y' have different environments; x's is used")
      }
  } else if(!is.gEnv(ey))
      environment(f) <- ey
  f
}

## Combine method.
c.gamlss2 <- function(...)
{
  objects <- list(...)
  x <- NULL
  for(i in 1L:length(objects))
    x <- c(x, objects[i])
  Call <- match.call()
  names(x) <- as.character(Call[-1L])
  class(x) <- c("gamlss2.list")
  return(x)
}


## Generic gamlss method.
gamlss2 <- function(formula, ...)
{
  UseMethod("gamlss2")
}

## Formula method.
gamlss2.formula <- function(formula, data, family = NO,
  subset, na.action, weights, offset, start = NULL, knots = NULL,
  control = gamlss2_control(...), ...)
{
  ## Save environments.
  menv <- parent.frame()
  fenv <- environment(formula)

  ## Check more formula versions.
  for(j in c("sigma.f", "nu.f", "tau.f")) {
    if(length(i <- grep(j, names(control), value = TRUE, fixed = TRUE))) {
      k <- grep(i, c("sigma.formula", "nu.formula", "tau.formula"), value = TRUE, fixed = TRUE)
      control[[k]] <- control[[i]]
    }
  }

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
  family <- complete_family(family, .links = control$links)

  ## Use numeric hessian?
  if(isTRUE(control$numhessian)) {
    family$update <- NULL
    family$hessian <- NULL
    family <- complete_family(family, .links = control$links)
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
  if(length(formula)[2L] < 2L & FALSE) {
    formula <- as.Formula(formula(formula), ~ 1)
  }

  ## Expand formula.
  if((length(attr(formula, "rhs")) < length(family$names)) & control$expand) {
    k <- length(family$names) - length(attr(formula, "rhs"))
    attr(formula, "rhs") <- c(attr(formula, "rhs"), as.list(rep(1, k)))
  }

  ## Check for "." in formula.
  for(i in seq_len(length(formula)[2L])) {
    rhs <- formula(formula, rhs = i)
    if(as.character(rhs[3]) == ".") {
      if(!inherits(data, "environment")) {
        yn <- NULL
        for(j in seq_len(length(formula)[1L]))
          yn <- c(yn, as.character(formula(formula, lhs = j))[2L])
        if(i < 2) {
          vn <- names(data)
          vn <- vn[!(vn %in% yn)]
          attr(formula, "rhs")[[i]] <- as.call(str2lang(paste(vn, collapse = "+")))
        } else {
          attr(formula, "rhs")[[i]] <- attr(formula, "rhs")[[i - 1L]]
        }
      } else {
        stop('using "." in formula but no data argument supplied!')
      }
    }
  }

  ## Ordinary formulas need no special-term rewriting. This is the common
  ## path, and avoids three costly recursive Formula reconstructions.
  formula_names <- all.names(formula, functions = TRUE)
  has_specials <- any(formula_names %in% .gamlss2_special_names)
  has_list_term <- "list" %in% formula_names
  if(has_specials || has_list_term) {
    mf_formula <- fake_formula(formula)
    ff <- fake_formula(formula, nospecials = TRUE)
    Sterms <- fake_formula(formula, onlyspecials = TRUE)
  } else {
    ## A colon must be expanded only in the model-frame formula; it remains an
    ## interaction in the linear design formula.
    mf_formula <- if(":" %in% formula_names) fake_formula(formula) else formula
    ff <- formula
    Sterms <- rep.int(list(character(0L)), length(formula)[2L])
  }
  ## model.frame.Formula() immediately collapses all formula parts to an
  ## ordinary terms object. Supplying that equivalent plain formula directly
  ## avoids the extra S3 and Formula dispatch.
  mf$formula <- formula(mf_formula, collapse = TRUE)

  ## Evaluate model.frame.
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## Response and model.matrix.
  mt <- X <- list()
  for(j in seq_len(length(ff)[2L])) {
    ffj <- formula(ff, lhs = 0, rhs = j)
    mt[[j]] <- terms(ffj, data = data)
    X[[j]] <- model.matrix(mt[[j]], mf)
  }

  xnames <- unlist(lapply(X, colnames), use.names = FALSE)
  unique_xnames <- sort(unique(xnames))
  first_xnames <- colnames(X[[1L]])
  if(all(unique_xnames %in% first_xnames)) {
    ## With expanded formulas, later parameter predictors are commonly just
    ## an intercept already present in the first matrix. Reuse that matrix
    ## instead of allocating a duplicate cbind followed by another subset.
    X <- if(identical(unique_xnames, first_xnames)) {
      X[[1L]]
    } else {
      X[[1L]][, unique_xnames, drop = FALSE]
    }
  } else {
    X <- do.call("cbind", X)
    X <- X[, unique_xnames, drop = FALSE]
  }

  if(!("(Intercept)" %in% unique_xnames))
    X <- cbind("(Intercept)" = 1.0, X)
  xnames <- colnames(X)
  if(anyNA(X))
    stop("detected 'NA' values in data!")
  ## Avoid allocating a logical matrix as large as X in the usual all-finite
  ## case. anyNA() above preserves the distinct NA diagnostic.
  if(length(X) && (!is.finite(min(X)) || !is.finite(max(X))))
    stop("detected 'Inf' values in data!")

  Y <- model.response(mf)
  if(is.null(Y)) {
    rn <- response_name(formula)
    Y <- mf[, rn]
  }

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
  for(i in seq_len(length(ff)[2L])) {
    Xterms[[i]] <- attr(mt[[i]], "term.labels")
    if(attr(mt[[i]], "intercept") > 0L)
      Xterms[[i]] <- c("(Intercept)", Xterms[[i]])

    ## model.part() rebuilds the terms object and subsets the whole model
    ## frame. The already computed terms object identifies the same evaluated
    ## offset columns directly.
    oi <- attr(mt[[i]], "offset")
    offi <- NULL
    if(length(oi)) {
      variables <- as.list(attr(mt[[i]], "variables"))[-1L]
      offset_names <- vapply(variables[oi], deparse,
        character(1L), width.cutoff = 500L)
      offi <- mf[[offset_names[1L]]]
      if(length(offset_names) > 1L) {
        for(k in offset_names[-1L])
          offi <- offi + mf[[k]]
      }
    }
    offi <- expand_offset(offi)
    offsets[[i]] <- if(length(offi)) offi else numeric(0)
  }

  mt <- mt[seq_along(family$names)]
  names(mt) <- family$names

  if(!missing(offset)) {
    anyoff <- any(lengths(offsets) > 0L)
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
    if(nrow(offset) < 2)
      offset <- offset[rep(1L, nrow(mf)), , drop = FALSE]
    rownames(offset) <- rownames(mf)
  }

  ## Process special terms.
  Specials <- special_terms(Sterms, mf, binning = control$binning,
    digits = control$digits, select = control$select, knots = knots)

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
  if((length(family$names) < 2)) {
    Xterms <- list(unlist(Xterms))
    names(Xterms) <- family$names
  } else {
    names(Xterms) <- family$names[1:length(Xterms)]
  }
  if((length(family$names) < 2)) {
    Sterms <- list(unlist(Sterms))
    names(Sterms) <- family$names
  } else {
    names(Sterms) <- family$names[1:length(Sterms)]
  }
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
  xlev <- lapply(mt, function(x) .getXlevels(x, mf))

  for(i in names(Xterms)) {
    ## Factors.
    for(j in names(xlev[[i]])) {
      if(j %in% Xterms[[i]]) {
        xl <- xl0 <- paste0(j, xlev[[i]][[j]])
        xl <- xl[xl %in% xnames]
        if(attr(mt[[i]], "intercept") > 0L && all(xl %in% xnames) && length(xl0) == length(xl)) {
          xl <- xl[-1L]
        }
        if(length(xl)) {
          Xterms[[i]][Xterms[[i]] == j] <- NA
          Xterms[[i]] <- Xterms[[i]][!is.na(Xterms[[i]])]
          Xterms[[i]] <- c(Xterms[[i]], xl)
        }
      }
    }
    ## Others.
    for(j in Xterms[[i]]) {
      if(!(j %in% xnames)) {
        if(any(ij <- grepl(j, xnames, fixed = TRUE))) {
          Xterms[[i]][Xterms[[i]] == j] <- NA
          Xterms[[i]] <- Xterms[[i]][!is.na(Xterms[[i]])]
          Xterms[[i]] <- c(Xterms[[i]], xnames[ij])
        }
      }
      ## Interactions.
      if(grepl(":", j, fixed = TRUE)) {
        xl <- strsplit(j, ":", fixed = TRUE)[[1]][2]
        xl <- paste0(":", xl)
        xl <- grep(xl, xnames, fixed = TRUE, value = TRUE)
        Xterms[[i]] <- unique(c(Xterms[[i]], xl))
      }
    }
  }

  ## Drop.
  for(i in names(Xterms)) {
    Xterms[[i]] <- Xterms[[i]][Xterms[[i]] %in% xnames]
  }

  attr(Xterms, "xlevels") <- xlev

  ## Optionally, use optimizer function provided from family
  optimizer <- if(is.null(family$optimizer)) {
    control$optimizer
  } else {
    family$optimizer
  }

  ## Response sanity check.
  if(!is.null(family$valid.response)) {
    if(!family$valid.response(Y)) {
      stop(paste0("please check the response, family '",
        family$family, "' validity check not passed!"))
    }
  }

  ## Check for binomial families.
  if(family$family %in% .bi.list) {
    fenv <- environment(family[["d"]])
    ybd <- get_y_bd(Y)
    Y <- ybd$y
    attr(Y, "bd") <- ybd$bd
    fenv$bd <- ybd$bd
    environment(optimizer) <- fenv
  }

  ## Track runtime.
  tstart <- proc.time()

  ## Estimation.
  rval <- optimizer(x = X, y = Y, specials = Specials, family = family,
    offsets = offsets, weights = weights, start = start, xterms = Xterms, sterms = Sterms,
    control = control)

  ## Runtime.
  elapsed <- as.numeric((proc.time() - tstart)["elapsed"])

  ## Further model information.
  rval$call <- cl
  rval$formula <- formula
  rval$fake_formula <- mf_formula
  rval$terms <- mt ## terms(merge_formula(formula(rval$fake_formula, collapse = TRUE), as.formula(mt)))
  environment(rval$terms) <- menv
  rval$family <- family
  rval$xlevels <- xlev
  rval$contrasts <- attr(X, "contrasts")
  rval$na.action <- attr(mf, "na.action")
  attr(Xterms, "terms") <- mt
  if(is.null(rval$selection) | isTRUE(rval$selection$select)) {
    rval$xterms <- Xterms
    rval$sterms <- Sterms
    rval$specials <- Specials
  } else {
    attr(rval$xterms, "terms") <- mt
  }
  rval$df <- get_df(rval)
  rval$weights <- weights
  rval$elapsed <- elapsed

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
    rval$results <- results(rval, data = mf)
  } else {
    rval$fitted.values <- NULL
    rval$weights <- NULL
    if(!is.null(rval$fitted.linear)) {
      for(j in names(rval$fitted.linear))
        rval$fitted.linear[[j]]$fitted.values <- NULL
    }
    if(!is.null(rval$specials)) {
      for(j in names(rval$specials)) {
        if(!is.null(rval$specials[[j]][["X"]])) {
          rval$specials[[j]][["X"]] <- NULL
          ##rval$specials[[j]][["Xu"]] <- NULL
        }
      }
    }
  }

  class(rval) <- unique(c(class(rval), "gamlss2"))

  return(rval)
}

## List method.
gamlss2.list <- function(formula, ...)
{
  cl <- match.call()
  fl <- list()
  for(j in seq_along(formula)) {
    fl[[j]] <- as.character(formula[[j]])
    if(length(fl[[j]]) > 2L) { 
      fl[[j]] <- as.formula(paste(fl[[j]][2L], "~", fl[[j]][3L]))
    } else {
      fl[[j]] <- as.formula(paste("~", fl[[j]][2L]))
    }
  }
  cl$formula <- do.call("as.Formula", fl)
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
  if(is.null(ctr$demmler.reinsch))
    ctr$demmler.reinsch <- "auto"
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
  nargs <- intersect(c("data", "na.action", "subset"), names(dots))
  if(length(nargs) || is.null(formula$model)) {
    fcall <- formula$call
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    drop.unused.levels <- list(...)$drop.unused.levels
    if(is.null(drop.unused.levels))
      drop.unused.levels <- FALSE
    fcall$drop.unused.levels <- drop.unused.levels
    fcall[[1L]] <- quote(model.frame)
    xlev <- list()
    for(j in seq_along(formula$xlevels)) {
      for(i in names(formula$xlevels[[j]]))
        xlev[[i]] <- formula$xlevels[[j]][[i]]
    }
    xlev <- xlev[unique(names(xlev))]
    fcall$xlev <- xlev
    fcall$formula <- formula$fake_formula
    if(!dots$keepresponse) {
      fcall$formula <- update(fcall$formula, NULL ~ .)
      fcall$formula <- formula(as.Formula(fcall$formula), lhs = 0)
    }
    fcall$formula <- formula(as.Formula(fcall$formula), collapse = TRUE, update = TRUE)
    no_weights <- list(...)$no_weights
    if(isTRUE(no_weights)) {
      fcall["weights"] <- NULL
    }
    fcall[nargs] <- dots[nargs]
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
    if(is.list(data))
      data <- as.data.frame(data)
    dots <- list(...)
    dots$data <- dots$contrasts.arg <- NULL
    mt <- object$terms
    X <- list()
    for(j in names(mt)) {
      X[[j]] <- do.call(stats::model.matrix.default, c(list(object = list("terms" = mt[[j]]), 
        data = data, contrasts.arg = object$contrasts), dots))
    }
    X <- do.call("cbind", X)
    X <- X[, sort(unique(colnames(X))), drop = FALSE]
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
  x$call[[1]] <- if(inherits(x, "bamlss2")) as.name("bamlss2") else as.name("gamlss2")
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
  return(invisible(NULL))
}

## Merging formulas.
is_formula <- function(x) inherits(x, "formula")

merge_formula <- function(x, y, ...)
{
  if(!is_formula(x) || length(x) != 3)
    stop("first argument is invalid!")
  if(!is_formula(y)) stop("second argument is invalid!")
  if(length(list(...))) warning("extraneous arguments discarded!")
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
          warning("`x' and `y' have different environments; x's is used!")
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

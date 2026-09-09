## Predict method.
predict.gamlss2 <- function(object, 
  model = NULL, newdata = NULL,
  type = c("parameter", "link", "response", "terms"), 
  terms = NULL, se.fit = FALSE, drop = TRUE, ...,
  level = NULL, interval.cache = NULL)
{
  if(!is.null(level)) {
    return(predict_wald(object, model, newdata, match.arg(type), terms,
      drop, list(...), level, interval.cache))
  }

  ## FIXME: se.fit, terms ...
  samples <- NULL
  if(se.fit || !is.null(list(...)$FUN) || inherits(object, "bamlss2")) {
    if(is.null(object$samples)) {
      R <- list(...)$R
      if(is.null(R))
        R <- 200L
      seed <- list(...)$seed
      if(is.null(seed))
        seed <- 123
      if(is.logical(seed)) {
        seed <- if(!seed) NULL else seed <- 123
      }
      if(!is.null(seed))
        set.seed(seed)
      samples <- sampling(object, R = R, full = TRUE)
    } else {
      samples <- object$samples
    }
  }

  if(!is.null(samples)) {
    burnin <- list(...)$burnin
    if(!is.null(burnin)) {
      burnin <- as.integer(burnin)
      samples <- samples[-seq.int(burnin), , drop = FALSE]
    }
  }

  FUN <- list(...)$FUN
  if(is.null(FUN))
    FUN <- mean
  if(se.fit) {
    FUN <- function(x) {
       c("fit" = mean(x, na.rm = TRUE), "se" = sd(x, na.rm = TRUE))
    }
  }

  type <- match.arg(type)

  ## Extract the model frame.
  if(!is.null(newdata)) {
    mf <- try(model.frame(object, data = newdata,
      keepresponse = object$family$family %in% .bi.list, ...), silent = TRUE)
    if(inherits(mf, "try-error")) {
      mf <- model.frame(object, data = newdata, ...)
    }
  } else {
    mf <- if(is.null(object$model)) {
      model.frame(object)
    } else {
      object$model
    }
  }

  ## Check for binomial families.
  if(object$family$family %in% .bi.list) {
    fenv <- environment(object$family[["d"]])
    fenv$bd <- get_y_bd(model.response(mf))$bd
  }

  ## Linear effects design matrix.
  X <- model.matrix(object, data = mf)

  family <- object$family

  ## Which parameter model to predict?
  if(is.null(model)) {
    model <- list(...)$what
    if(is.null(model))
      model <- list(...)$parameter
    if(is.null(model))
      model <- family$names
  }
  if(!is.character(model))
    model <- family$names[model]
  model <- family$names[pmatch(model, family$names)]

  if((type == "response") && (length(family$names) > 1L)) {
    if(length(model) != length(family$names))
      stop('Predictions on the response scale require all distributional parameters. Please omit the "model" argument or specify all parameters.')
  }

  tt <- type == "terms"

  ## Predict all specified parameters.
  p <- list()
  for(j in model) {
    p[[j]] <- if(tt) {
      NULL
    } else {
      if(is.null(samples)) {
        rep(0.0, length.out = nrow(mf))
      } else {
        matrix(0.0, nrow = nrow(mf), ncol = nrow(samples))
      }
    }
    terms <- gsub(" ", "", terms)
    if(length(terms) < 1L)
      terms <- NULL
    tj <- if(is.null(terms)) {
      c(object$xterms[[j]], object$sterms[[j]])
    } else {
      if(isTRUE(list(...)$nogrep)) {
        terms
      } else {
        grep2(terms, c(object$xterms[[j]], object$sterms[[j]]), fixed = TRUE, value = TRUE)
      }
    }
    tj <- unique(tj)
    if(length(tj)) {
      ## Linear effects.
      if(length(object$xterms[[j]])) {
        xn <- NULL
        if(isTRUE(list(...)$nogrep)) {
          xn <- object$xterms[[j]][object$xterms[[j]] %in% tj]
        } else {
          for(i in tj) {
            xn <- c(xn, grep(i, object$xterms[[j]], fixed = TRUE, value = TRUE))      
          }
        }
        if(length(xn)) {
          xn2 <- list()
          for(i in seq_along(xn)) {
            if(!is.null(object$xlevels)) {
              if(!is.null(object$xlevels[[xn[i]]])) {
                xn2[[xn[i]]] <- paste0(xn[i], object$xlevels[[xn[i]]])
              }
            }
          }
          if(length(xn2)) {
            xnn <- xn
            xn <- as.list(xn)
            names(xn) <- xnn
            for(i in names(xn2))
              xn[[i]] <- xn2[[i]]
            xn <- unlist(xn)
            names(xn) <- NULL
            xn <- xn[xn %in% colnames(X)]
          }
          xn <- unique(xn)
          if(tt) {
            if(is.null(samples)) {
              ft <- t(t(X[, xn, drop = FALSE]) * object$coefficients[[j]][xn])
            } else {
              ij <- paste0(j, ".p.", xn)
              ps <- list()
              for(l in seq_along(ij)) {
                ps[[l]] <- apply(samples, 1, function(beta) {
                  X[, xn[l], drop = TRUE] * beta[ij[l]]
                })
                if(nrow(ps[[l]]) != nrow(mf))
                  ps[[l]] <- t(ps[[l]])
                ps[[l]] <- apply(ps[[l]], 1, FUN)
                if(!is.null(dim(ps[[l]]))) {
                  if(nrow(ps[[l]]) != nrow(mf))
                    ps[[l]] <- t(ps[[l]])
                } else {
                  ps[[l]] <- matrix(ps[[l]], ncol = 1L)
                }
                rownames(ps[[l]]) <- rownames(mf)
                colnames(ps[[l]]) <- paste0(xn[l], if(ncol(ps[[l]]) > 1) "." else "", colnames(ps[[l]]))
              }
              ft <- do.call("cbind", ps)
            }
            p[[j]] <- cbind(p[[j]], ft)
          } else {
            if(is.null(samples)) {
              p[[j]] <- p[[j]] + drop(X[, xn, drop = FALSE] %*% object$coefficients[[j]][xn])
            } else {
              ij <- paste0(j, ".p.", xn)
              ps <- apply(samples, 1, function(beta) {
                X[, xn, drop = FALSE] %*% beta[ij]
              })
              p[[j]] <- p[[j]] + ps
            }
          }
        }
      }
      ## Special effects.
      if(length(object$sterms[[j]])) {
        if(isTRUE(list(...)$nogrep)) {
          xn <- object$sterms[[j]][object$sterms[[j]] %in% tj]
        } else {
          xn <- NULL
          for(i in tj) {
            xn <- c(xn, grep(i, object$sterms[[j]], fixed = TRUE, value = TRUE))      
          }
        }
        xn <- unique(xn)
        if(length(xn)) {
          for(i in xn) {
            fit <- if(is.null(samples)) {
              rep(0.0, nrow(mf))
            } else {
              matrix(0.0, nrow = nrow(mf), ncol = nrow(samples))
            }
            if(inherits(object$specials[[i]], "mgcv.smooth")) {
              if(!is.null(object$fitted.specials[[j]][[i]]$selected)) {
                Xs <- PredictMat(object$specials[[i]], data = mf, n = nrow(mf))
                if(is.null(samples)) {
                  co <- object$fitted.specials[[j]][[i]]$coefficients
                  fit <- drop(Xs %*% co)
                } else {
                  ij <- paste0(paste0(j, ".s.", i), ".", 1:ncol(Xs))
                  fit <- apply(samples, 1, function(beta) {
                    Xs %*% beta[ij]
                  })
                }
              }
            } else {
              if(inherits(object$specials[[i]], "special")) {
                if(!is.null(object$fitted.specials[[j]][[i]])) {
                  fit <- special_predict(object$fitted.specials[[j]][[i]], data = mf)
                }
              } else {
                cs <- object$fitted.specials[[j]][[i]]$coefficients
                if(inherits(cs, "random")) {
                  vn <- as.character(as.call(as.call(parse(text = i))[[1L]])[[2L]])
                  xv <- mf[[vn]]
                  fit <- cs$coef[as.character(xv)]
                } else {
                  fit <- try(cs$fun(mf[[cs$name]]), silent = TRUE)
                  if(inherits(fit, "try-error")) {
                    fit <- try(predict(cs, newdata = mf), silent = TRUE)
                    if(inherits(fit, "try-error")) {
                      warning(paste0("cannot predict model term '", i, "'!"))
                      fit <- rep(0.0, nrow(mf))
                    }
                  }
                }
              }
            }
            if(tt) {
              if(is.null(samples)) {
                fit <- matrix(fit, ncol = 1L)
                colnames(fit) <- i
              } else {
                fit <- apply(fit, 1, FUN)
                if(!is.null(dim(fit))) {
                  if(nrow(fit) != nrow(mf)) {
                    fit <- t(fit)
                  }
                } else {
                  fit <- matrix(fit, ncol = 1)
                }
                colnames(fit) <- paste0(i, if(ncol(fit) > 1L) "." else "", colnames(fit))
              }
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

  ## Map to parameter scale.
  if(type %in% c("parameter", "response") & !tt) {
    p <- family$map2par(p)
  }

  ## Compute mean or median predictions.
  if(type == "response") {
    fm <- family$mean
    warn <- NULL
    if(is.null(fm)) {
      warn <- 'Prediction with type = "response", however, the mean function is missing in the family!'
      if(!is.null(family$q)) {
        warn <- paste(warn, 'Using the median instead!')
        fm <- function(par) family$quantile(par = par, 0.5)
      }
    }
    if(is.null(fm)) {
      warning('Prediction with type = "response", however, the mean function is missing in the family! Using the first parameter instead!')
    } else {
      if(!is.null(warn))
        warning(warn)
    }
    if(!is.null(samples)) {
      yp <- matrix(0.0, nrow(p[[1L]]), ncol(p[[1L]]))
      for(i in 1:ncol(p[[1L]])) {
        pi <- sapply(names(p), function(j) { p[[j]][, i] })
        pi <- as.data.frame(pi)
        yp[, i] <- fm(pi)
      }
      p <- t(apply(yp, 1, FUN))
    } else {
      p <- if(is.null(fm)) p[[1L]] else fm(p)
    }
  }

  if(!is.null(samples) & !tt) {
    for(j in names(p)) {
      p[[j]] <- apply(p[[j]], 1, FUN)
      if(!is.null(dim(p[[j]]))) {
        if((nrow(p[[j]]) != nrow(mf)) || any(rownames(p[[j]]) %in% c("fit", "se")))
          p[[j]] <- t(p[[j]])
      }
    }
  }

  ## Drop dimension if only one parameter is predicted.
  if(is.list(p)) {
    if((length(p) < 2 & drop)) {
      p <- p[[1L]]
    } else {
      if(!tt)
        p <- as.data.frame(p)
    }
  }
  if(is.matrix(p)) {
    if(!is.null(colnames(p))) {
      p <- as.data.frame(p)
    }
  }

  return(p)
}

## Multiple grep.
grep2 <- function (pattern, x, ...) 
{
  i <- NULL
  for(p in pattern)
    i <- c(i, grep(p, x, ...))
  unique(i)
}

## Extract fitted values.
fitted.gamlss2 <- function(object, newdata = NULL,
  type = c("parameter", "link"), model = NULL, ...)
{
  type <- match.arg(type)

  if(is.null(newdata) & !is.null(object$fitted.values)) {
    fit <- object$fitted.values
  } else {
    fit <- predict(object, newdata = newdata, type = "link")
  }

  if(type == "parameter")
    fit <- family(object)$map2par(fit)

  if(is.null(model)) {
    model <- list(...)$what
    if(is.null(model))
      model <- list(...)$parameter
    if(is.null(model))
      model <- object$family$names
  }
  if(!is.character(model))
    model <- object$family$names[model]
  model <- object$family$names[pmatch(model, object$family$names)]

  return(fit[, model])
}

## Compute marginal prediction profiles for each covariate.
marginal_predict <- function(object, newdata = NULL, variables = NULL,
  n = 100L, continuous = median, at = NULL, values = NULL, ...)
{
  if(!inherits(object, c("gamlss2", "bamlss2")))
    stop("'object' must inherit from class 'gamlss2' or 'bamlss2'.")

  if(is.null(newdata))
    newdata <- model.frame(object)
  newdata <- as.data.frame(newdata)
  if(!nrow(newdata))
    stop("'newdata' must contain at least one row.")

  ## Covariates used by the model, excluding the response.
  tt <- try(terms(object$fake_formula), silent = TRUE)
  vars <- if(inherits(tt, "try-error")) {
    all.vars(object$fake_formula)
  } else {
    all.vars(delete.response(tt))
  }
  vars <- unique(vars)

  missing_vars <- setdiff(vars, names(newdata))
  if(length(missing_vars)) {
    stop("variable(s) not found in 'newdata': ",
      paste(missing_vars, collapse = ", "))
  }
  if(!length(vars))
    stop("the model does not contain any covariates.")

  if(is.null(variables)) {
    variables <- vars
  } else {
    variables <- as.character(variables)
    if(!length(variables) || anyNA(variables) || any(!nzchar(variables)))
      stop("'variables' must contain at least one non-missing variable name.")
    variables <- unique(variables)
    bad <- setdiff(variables, vars)
    if(length(bad))
      stop("unknown variable(s): ", paste(bad, collapse = ", "))
  }

  if(length(n) != 1L || is.na(n) || !is.numeric(n) || !is.finite(n) ||
      n < 2L || n > .Machine$integer.max || n != floor(n))
    stop("'n' must be a single integer greater than 1.")
  n <- as.integer(n)

  if(is.character(continuous) && length(continuous) == 1L)
    continuous <- match.fun(continuous)
  if(!is.function(continuous))
    stop("'continuous' must be a function, e.g. median or mean.")

  check_named_list <- function(x, name, allowed) {
    if(is.null(x))
      return(list())
    if(!is.list(x) || is.null(names(x)) ||
        anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x)))
      stop("'", name, "' must be a named list with unique, non-empty names.")
    bad <- setdiff(names(x), allowed)
    if(length(bad))
      stop("unknown variable(s) in '", name, "': ",
        paste(bad, collapse = ", "))
    x
  }

  at <- check_named_list(at, "at", vars)
  values <- check_named_list(values, "values", variables)

  ## Keep factor levels/classes intact.
  coerce_value <- function(value, x, name) {
    if(is.factor(x)) {
      value <- as.character(value)
      if(anyNA(value) || any(!value %in% levels(x)))
        stop("invalid level for '", name, "'.")
      return(factor(value,
        levels = levels(x),
        ordered = is.ordered(x)))
    }

    if(is.character(x))
      return(as.character(value))

    if(inherits(x, "Date"))
      return(as.Date(value, origin = "1970-01-01"))

    if(inherits(x, "POSIXct"))
      return(as.POSIXct(value,
        origin = "1970-01-01",
        tz = attr(x, "tzone")))

    if(is.logical(x))
      return(as.logical(value))

    if(is.numeric(x) && !is.numeric(value))
      stop("invalid numeric value for '", name, "'.")

    value
  }

  check_value <- function(value, name, single = FALSE) {
    if(!length(value))
      stop("no values supplied for '", name, "'.")
    if(single && length(value) != 1L)
      stop("the value for '", name, "' must have length 1.")
    if(anyNA(value))
      stop("missing values are not allowed for '", name, "'.")
    if((is.numeric(value) || inherits(value, c("Date", "POSIXct"))) &&
        any(!is.finite(value)))
      stop("non-finite values are not allowed for '", name, "'.")
    value
  }

  ## Value used when a covariate is held fixed.
  fixed_value <- function(x, name) {
    if(name %in% names(at)) {
      z <- at[[name]]
      if(length(z) != 1L)
        stop("'at[[\"", name, "\"]]' must have length 1.")
      z <- coerce_value(z, x, name)
      return(check_value(z, name, single = TRUE))
    }

    ## First factor level by default.
    if(is.factor(x)) {
      z <- factor(levels(x)[1L],
        levels = levels(x),
        ordered = is.ordered(x))
      return(check_value(z, name, single = TRUE))
    }

    if(is.character(x)) {
      z <- unique(x[!is.na(x)])
      if(!length(z))
        stop("no non-missing values available for '", name, "'.")
      return(z[1L])
    }

    if(is.logical(x)) {
      z <- unique(x[!is.na(x)])
      if(!length(z))
        stop("no non-missing values available for '", name, "'.")
      return(if(FALSE %in% z) FALSE else z[1L])
    }

    z <- try(continuous(x, na.rm = TRUE), silent = TRUE)
    if(inherits(z, "try-error"))
      z <- continuous(x[!is.na(x)])

    if(length(z) != 1L || is.na(z))
      stop("could not compute a representative value for '",
        name, "'.")

    z <- coerce_value(z, x, name)
    check_value(z, name, single = TRUE)
  }

  ## Values over which the focal covariate is varied.
  focal_values <- function(x, name) {
    if(name %in% names(values)) {
      z <- coerce_value(values[[name]], x, name)
      return(check_value(z, name))
    }

    if(is.factor(x)) {
      z <- factor(levels(x),
        levels = levels(x),
        ordered = is.ordered(x))
      return(check_value(z, name))
    }

    if(is.character(x))
      return(check_value(unique(x[!is.na(x)]), name))

    if(is.logical(x))
      return(check_value(sort(unique(x[!is.na(x)])), name))

    if(inherits(x, c("Date", "POSIXct"))) {
      if(all(is.na(x)))
        stop("no non-missing values available for '", name, "'.")
      r <- range(x, na.rm = TRUE)
      z <- if(r[1L] == r[2L]) r[1L] else
        seq(r[1L], r[2L], length.out = n)
      return(check_value(z, name))
    }

    if(is.numeric(x)) {
      r <- range(x, na.rm = TRUE)

      if(any(!is.finite(r)))
        stop("no finite values available for '", name, "'.")

      if(r[1L] == r[2L])
        return(r[1L])

      return(seq(r[1L], r[2L], length.out = n))
    }

    check_value(unique(x[!is.na(x)]), name)
  }

  ## Construct one representative observation.
  base <- newdata[1L, , drop = FALSE]
  fixed <- list()

  for(j in vars) {
    fixed[[j]] <- fixed_value(newdata[[j]], j)
    base[[j]] <- fixed[[j]]
  }

  res <- vector("list", length(variables))
  names(res) <- variables

  for(i in variables) {
    xi <- focal_values(newdata[[i]], i)

    nd <- base[rep(1L, length(xi)), , drop = FALSE]
    nd[[i]] <- xi
    rownames(nd) <- NULL

    fit <- predict(object, newdata = nd, ...)

    if(is.atomic(fit) && is.null(dim(fit))) {
      fit <- data.frame(.fit = fit)
    } else {
      fit <- as.data.frame(fit)
    }

    if(nrow(fit) != nrow(nd))
      stop("prediction for '", i,
        "' returned an unexpected number of rows.")

    res[[i]] <- cbind(nd[i], fit)
    rownames(res[[i]]) <- NULL
  }

  ## Store the conditioning values for reference.
  attr(res, "fixed") <- fixed

  res
}

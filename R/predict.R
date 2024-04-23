## Predict method.
predict.gamlss2 <- function(object, 
  model = NULL, newdata = NULL, type = c("parameter", "link", "response", "terms"), 
  terms = NULL, se.fit = FALSE, drop = TRUE, ...)
{
  ## FIXME: se.fit, terms ...
  samples <- NULL
  if(se.fit) {
    R <- list(...)$R
    if(is.null(R))
      R <- 100L
    samples <- sampling(object, R = R, full = TRUE)
  }

  type <- match.arg(type)

  ## Extract the model frame.
  if(!is.null(newdata)) {
    mf <- model.frame(object, data = newdata)
  } else {
    mf <- if(is.null(object$model)) {
      model.frame(object)
    } else {
      object$model
    }
  }

  ## Linear effects design matrix.
  X <- model.matrix(object, data = mf)

  family <- object$family

  ## Which parameter model to predict?
  if(is.null(model)) {
    model <- list(...)$what
    if(is.null(model))
      model <- family$names
  }
  if(!is.character(model))
    model <- family$names[model]
  model <- family$names[pmatch(model, family$names)]

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
    tj <- if(is.null(terms)) {
      c(object$xterms[[j]], object$sterms[[j]])
    } else {
      grep2(terms, c(object$xterms[[j]], object$sterms[[j]]), fixed = TRUE, value = TRUE)
    }
    tj <- unique(tj)
    if(length(tj)) {
      ## Linear effects.
      if(length(object$xterms[[j]])) {
        xn <- NULL
        for(i in tj) {
          xn <- c(xn, grep(i, object$xterms[[j]], fixed = TRUE, value = TRUE))      
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
          if(tt) {
            ft <- t(t(X[, xn, drop = FALSE]) * object$coefficients[[j]][xn])
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
        xn <- NULL
        for(i in tj) {
          xn <- c(xn, grep(i, object$sterms[[j]], fixed = TRUE, value = TRUE))      
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

  ## Map to parameter scale.
  if(type %in% c("parameter", "response") & !tt) {
    p <- family$map2par(p)
  }

  ## Compute mean or median predictions.
  if(type == "response") {
    fm <- family$mean
    if(is.null(fm)) {
      if(!is.null(family$q)) {
        fm <- function(par) family$q(0.5, par)
      }
    }
    p <- if(is.null(fm)) p[[1L]] else fm(p)
  }

  if(!is.null(samples) & !tt) {
    FUN <- list(...)$FUN
    if(is.null(FUN))
      FUN <- mean
    if(se.fit) {
      FUN <- function(x) {
        c("fit" = mean(x, na.rm = TRUE), "se" = sd(x, na.rm = TRUE))
      }
    }
    for(j in names(p)) {
      p[[j]] <- apply(p[[j]], 1, FUN)
      if(!is.null(dim(p[[j]]))) {
        if(nrow(p[[j]]) != nrow(mf))
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
      model <- object$family$names
  }
  if(!is.character(model))
    model <- object$family$names[model]
  model <- object$family$names[pmatch(model, object$family$names)]

  return(fit[, model])
}


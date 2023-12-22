## Predict method.
predict.gamlss2 <- function(object, 
  parameter = NULL, newdata = NULL, type = c("link", "parameter", "response", "terms"), 
  terms = NULL, se.fit = FALSE, ...)
{
  ## FIXME: se.fit, terms ...

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

  ## Which parameter to predict?
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

  ## Predict all specified parameters.
  p <- list()
  for(j in parameter) {
    p[[j]] <- if(tt) NULL else rep(0, length.out = nrow(mf))
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
            p[[j]] <- p[[j]] + drop(X[, xn, drop = FALSE] %*% object$coefficients[[j]][xn])
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
            if(inherits(object$specials[[i]], "mgcv.smooth")) {
              if(object$fitted.specials[[j]][[i]]$selected) {
                Xs <- PredictMat(object$specials[[i]], data = mf, n = nrow(mf))
                co <- object$fitted.specials[[j]][[i]]$coefficients
                fit <- drop(Xs %*% co)
              }
            } else {
              if(inherits(object$specials[[i]], "special")) {
                fit <- special_predict(object$fitted.specials[[j]][[i]], data = mf)
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

  ## Drop dimension if only one parameter is predicted.
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
  type = c("link", "parameter"), ...)
{
  type <- match.arg(type)

  if(is.null(newdata) & !is.null(object$fitted.values)) {
    fit <- object$fitted.values
  } else {
    fit <- predict(object, newdata = newdata)
  }

  if(type == "parameter")
    fit <- family(object)$map2par(fit)

  return(fit)
}


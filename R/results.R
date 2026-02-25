## Compute results for linear effects
## and smooth terms for plotting and
## summary statistics.
results <- function(x, ...)
{
  UseMethod("results")
}

## Extract linear and special term information.
results.gamlss2 <- function(x, ...)
{
  res <- list()
  np <- x$family$names
  if(is.null(x$model))
    x$model <- model.frame(x)

  if(length(x$sterms)) {
    res$effects <- list()
    k <- 1L
    for(j in np) {
      if(length(x$sterms[[j]]) & (j %in% names(x$fitted.specials))) {
        for(i in x$sterms[[j]]) {
          ## For mgcv smooths.
          if("mgcv.smooth" %in% class(x$specials[[i]])) {
            dim <- x$specials[[i]]$dim
            by <- x$specials[[i]]$by
            if(dim < 3 & !is.null(x$fitted.specials[[j]][[i]])) {
              if(dim > 1) {
                xc <- unlist(lapply(x$model[, x$specials[[i]]$term, drop = FALSE], function(x) {
                  if(inherits(x, "matrix"))
                    return("numeric")
                  else
                    return(class(x))
                }))
                if(all(xc %in% c("numeric", "integer"))) {
                  nd <- expand.grid(seq(min(x$model[[x$specials[[i]]$term[1L]]]),
                    max(x$model[[x$specials[[i]]$term[1L]]]), length = 50),
                    seq(min(x$model[[x$specials[[i]]$term[2L]]]),
                    max(x$model[[x$specials[[i]]$term[2L]]]), length = 50))
                } else {
                  next
                }
              } else {
                if(!is.factor(x$model[[x$specials[[i]]$term]])) {
                  xr <- range(x$model[[x$specials[[i]]$term]])
                  nd <- data.frame(seq(xr[1L], xr[2L], length = 300L))
                } else {
                  xf <- sort(unique(x$model[[x$specials[[i]]$term]]))
                  nd <- data.frame(xf)
                }
              }
              names(nd) <- x$specials[[i]]$term
              if(by != "NA") {
                if(is.factor(x$model[[x$specials[[i]]$by]])) {
                  by.level <- x$specials[[i]]$by.level
                  xlevels <- levels(x$model[[x$specials[[i]]$by]])
                  nd[[by]] <- factor(by.level, levels = xlevels)
                } else {
                  nd[[by]] <- 1.0
                }
              }
              X <- PredictMat(x$specials[[i]], nd, n = nrow(nd))
              se <- rowSums((X %*% x$fitted.specials[[j]][[i]]$vcov) * X)
              se <- 2 * sqrt(se)
              nd$fit <- drop(X %*% coef(x$fitted.specials[[j]][[i]]))
              nd$lower <- nd$fit - se
              nd$upper <- nd$fit + se
              if(by == "NA") {
                lab <- strsplit(x$specials[[i]]$label, "")[[1L]]
                lab <- paste0(lab[-length(lab)], collapse = "")
                lab <- paste0(lab, ",", round(x$fitted.specials[[j]][[i]]$edf, 2L), ")")
              } else {
                lab <- x$specials[[i]]$label
                lab <- gsub("):",
                  paste0(",", round(x$fitted.specials[[j]][[i]]$edf, 2L), "):"), lab,
                  fixed = TRUE)
              }
              lab <- paste0(j, ".", lab)
              attr(nd, "label") <- lab
              res$effects[[lab]] <- nd
            }
          }
          ## gamlss smooth terms.
          if("smooth" %in% class(x$specials[[i]])) {
            ## FIXME: vcov?

            if(inherits(x$fitted.specials[[j]][[i]]$coefficients, "pb")) {
              xn <- attr(x$specials[[i]], "Name")
              xr <- range(x$model[[xn]])
              nd <- data.frame(seq(xr[1L], xr[2L], length = 300L))
              names(nd) <- xn
              nd$fit <- x$fitted.specials[[j]][[i]]$coefficients$fun(nd[[xn]])
              lab <- paste0(j, ".", i)
              attr(nd, "label") <- lab
              res$effects[[lab]] <- nd
            }
          }

          ## special terms.
          if(("special" %in% class(x$specials[[i]])) & (i %in% names(x$fitted.specials[[j]]))) {
            dim <- length(x$specials[[i]]$term)
            if(dim > 2)
              next

            nd <- list()

            if(dim > 1) {
              xc <- unlist(lapply(x$model[, x$specials[[i]]$term, drop = FALSE], class))
              if(all(xc %in% c("numeric", "matrix", "array"))) {
                nd <- expand.grid(seq(min(x$model[[x$specials[[i]]$term[1L]]]),
                  max(x$model[[x$specials[[i]]$term[1L]]]), length = 50),
                  seq(min(x$model[[x$specials[[i]]$term[2L]]]),
                  max(x$model[[x$specials[[i]]$term[2L]]]), length = 50))
              } else {
                next
              }
              nd <- as.data.frame(nd)
              names(nd) <- x$specials[[i]]$term
            } else {
              if(!is.null(dim(x$model[[x$specials[[i]]$term]]))) {
                if(ncol(x$model[[x$specials[[i]]$term]]) < 2L)
                  x$model[[x$specials[[i]]$term]] <- as.numeric(x$model[[x$specials[[i]]$term]])
              }
              if(!is.matrix(x$model[[x$specials[[i]]$term]])) {
                if(!is.factor(x$model[[x$specials[[i]]$term]])) {
                  xr <- range(x$model[[x$specials[[i]]$term]])
                  nd <- data.frame(seq(xr[1L], xr[2L], length = 300L))
                } else {
                  xf <- sort(unique(x$model[[x$specials[[i]]$term]]))
                  nd <- data.frame(xf)
                }
                nd <- as.data.frame(nd)
                names(nd) <- x$specials[[i]]$term
              } else {
                 nd <- list()
                 nd[[x$specials[[i]]$term]] <- x$model[[x$specials[[i]]$term]]
              }
            }

            p <- special_predict(x$fitted.specials[[j]][[i]], data = nd, se.fit = TRUE)

            if(is.null(dim(p))) {
              nd$fit <- as.numeric(p)
            } else {
              if(is.matrix(p))
                p <- as.data.frame(p)
              nd <- cbind(nd, p)
            }
            lab <- strsplit(x$specials[[i]]$label, "")[[1L]]
            lab <- paste0(lab[-length(lab)], collapse = "")
            lab <- paste0(lab, ",", round(x$fitted.specials[[j]][[i]]$edf, 2L), ")")
            lab <- paste0(j, ".", lab)
            attr(nd, "label") <- lab
            res$effects[[lab]] <- nd
          }
        }
      }
    }
  }
  if(length(x$xterms)) {
    if(length(res$effects) < 1L)
      res$effects <- list()
    res$effects <- c(res$effects, results_linear(x, NULL, x$model))
  }

  return(res)
}

'%||%' <- function(a, b) if (!is.null(a)) a else b

results_linear <- function(x, parameter = NULL, data, ...)
{
  if(is.null(parameter)) {
    parameter <- list(...)$what
    if(is.null(parameter)) parameter <- list(...)$model
    if(is.null(parameter)) parameter <- x$family$names
  }
  if(!is.character(parameter))
    parameter <- x$family$names[parameter]
  parameter <- x$family$names[pmatch(parameter, x$family$names)]
  parameter <- parameter[!is.na(parameter)]
  if(length(parameter) < 1L)
    stop("argument parameter is specified wrong!")

  ff <- fake_formula(formula(x), nospecials = TRUE)
  vn <- all.vars(ff)

  if(!length(vn)) {
    return(NULL)
  }

  env <- environment(formula(x))

  yn <- deparse(formula(as.Formula(formula(x)), rhs = 0)[[2L]])
  data <- data[, names(data) != yn, drop = FALSE]
  if(ncol(data) < 1L) {
    return(NULL)
  }

  nd <- list()
  for(j in names(data)) {
    if(is.numeric(data[[j]])) {
      nd[[j]] <- seq(min(data[[j]]), max(data[[j]]), length = 300L)
    } else {
      if(is.character(data[[j]])) {
        nd[[j]] <- factor(rep(unique(data[[j]]), length.out = 300L))
      } else {
        nd[[j]] <- factor(rep(levels(data[[j]]), length.out = 300L),
          levels = levels(data[[j]]))
      }
    }
  }
  rm(data)
  nd <- as.data.frame(nd)
  X <- model.matrix(x, data = nd)
  p <- list()
  for(j in seq_along(parameter)) {
    V <- x$fitted.linear[[parameter[j]]]$vcov
    if(is.null(V)) {
      next
    }
    cj <- x$fitted.linear[[parameter[j]]]$coefficients
    for(i in vn) {
      if(i %in% all.vars(formula(ff, lhs = 0, rhs = j))) {
        ii <- grep(i, colnames(X), value = TRUE)
        if(attr(terms(formula(ff, lhs = 0, rhs = j)), "intercept") > 0L)
          ii <- c("(Intercept)", ii)
        Xj   <- X[, ii, drop = FALSE]
        Vsub <- V[ii, ii, drop = FALSE]
        bsub <- cj[ii]
        if(is.factor(nd[[i]])) {
          Xjc <- Xj - matrix(colMeans(Xj), nrow(Xj), ncol(Xj), byrow = TRUE)
          fit <- as.vector(Xjc %*% bsub)
          vj  <- rowSums((Xjc %*% Vsub) * Xjc)
          sj  <- sqrt(pmax(vj, 0))
        } else {
          Xjc <- sweep(Xj, 2, colMeans(Xj), FUN = "-")
          fit <- as.vector(Xjc %*% bsub)
          vj  <- rowSums((Xjc %*% Vsub) * Xjc)
          sj  <- sqrt(pmax(vj, 0))
        }
        z <- qnorm(0.975)
        upper <- fit + z * sj
        lower <- fit - z * sj
        ji <- paste0(parameter[j], ".", i)
        p[[ji]] <- data.frame(nd[[i]], "fit" = fit, "lower" = lower, "upper" = upper)
        colnames(p[[ji]])[1L] <- i
        p[[ji]] <- p[[ji]][!duplicated(p[[ji]][[i]]), , drop = FALSE]
        p[[ji]] <- p[[ji]][order(p[[ji]][[i]]), , drop = FALSE]
        rownames(p[[ji]]) <- NULL
        attr(p[[ji]], "label") <- paste(ji, "effect")
        attr(p[[ji]], "linear") <- TRUE
      }
    }
  }

  p
}


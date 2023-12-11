## Compute results for linear effects
## and smooth terms for plotting and
## summary statistics.
results <- function(x, ...)
{
  UseMethod("results")
}

## Extract linear and special term information.
results.gamlss2 <- function(x)
{
  res <- list()
  np <- x$family$names
  if(is.null(x$model))
    x$model <- model.frame(x)

  if(length(x$sterms)) {
    res$effects <- list()
    k <- 1L
    for(j in np) {
      if(length(x$sterms[[j]])) {
        for(i in x$sterms[[j]]) {
          ## For mgcv smooths.
          if("mgcv.smooth" %in% class(x$specials[[i]])) {
            dim <- x$specials[[i]]$dim
            by <- x$specials[[i]]$by
            if(dim < 3) {
              if(dim > 1) {
                xc <- unlist(lapply(x$model[, x$specials[[i]]$term, drop = FALSE], class))
                if(all(xc == "numeric")) {
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
                lab <- paste0(j, ".", lab)
              } else {
                lab <- x$specials[[i]]$label
                lab <- gsub("):",
                  paste0(",", round(x$fitted.specials[[j]][[i]]$edf, 2L), "):"), lab,
                  fixed = TRUE)
              }
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
          if("special" %in% class(x$specials[[i]])) {
            dim <- length(x$specials[[i]]$term)

            if(dim > 2)
              next

            nd <- list()

            if(dim > 1) {
              xc <- unlist(lapply(x$model[, x$specials[[i]]$term, drop = FALSE], class))
              if(all(xc == "numeric")) {
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

            nd <- as.data.frame(nd)
            names(nd) <- x$specials[[i]]$term
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

  return(res)
}


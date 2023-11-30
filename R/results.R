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
    res$specials <- list()
    k <- 1L
    for(j in np) {
      if(length(x$sterms[[j]])) {
        for(i in x$sterms[[j]]) {
          if("mgcv.smooth" %in% class(x$specials[[i]])) {
            dim <- x$specials[[i]]$dim
            by <- x$specials[[i]]$by
            if(dim < 2) {
              if(!is.factor(x$model[[x$specials[[i]]$term]])) {
                xr <- range(x$model[[x$specials[[i]]$term]])
                nd <- data.frame(seq(xr[1L], xr[2L], length = 300L))
              } else {
                nd <- data.frame(unique(x$model[[x$specials[[i]]$term]]))
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
              se <- sqrt(se)
              cv <- qt(1 - 0.01/2, x$fitted.specials[[j]][[i]]$df)
              nd$fit <- drop(X %*% coef(x$fitted.specials[[j]][[i]]))
              nd$lower <- nd$fit - cv * se
              nd$upper <- nd$fit + cv * se
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
              res$specials[[lab]] <- nd
            }
          }
          if("smooth" %in% class(x$specials[[i]])) {
            ## FIXME: vcov?
          }
        }
      }
    }
  }

  return(res)
}


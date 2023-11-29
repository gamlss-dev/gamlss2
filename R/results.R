## Compute results for linear effects
## and smooth terms for plotting and
## summary statistics.
results <- function(x, ...)
{
  UseMethod("results")
}

results.gamlss2 <- function(x)
{
  res <- list()
  np <- x$family$names
  if(is.null(x$model))
    x$model <- model.frame(x)

  if(length(x$sterms)) {
    res$specials <- list()
    for(j in np) {
      if(length(x$sterms[[j]])) {
        res$specials[[j]] <- list()
        for(i in x$sterms[[j]]) {
          if("mgcv.smooth" %in% class(x$specials[[i]])) {
            dim <- x$specials[[i]]$dim
            if(dim < 2) {
              xr <- range(x$model[[x$specials[[i]]$term]])
              nd <- data.frame(seq(xr[1L], xr[2L], length = 300L))
              names(nd) <- x$specials[[i]]$term
              X <- PredictMat(x$specials[[i]], nd, n = 300L)
              nd$fit <- drop(X %*% coef(x$fitted.specials[[j]][[i]]))
              res$specials[[j]][[i]] <- nd
            }
          }
        }
      }
    }
  }

  return(res)
}

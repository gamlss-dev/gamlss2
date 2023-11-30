## A plotting method.
plot.gamlss2 <- function(x, ...)
{
  if(is.null(x$results))
    x$results <- results(x)
  par(mfrow = n2mfrow(length(x$results$specials)))
  for(j in names(x$results$specials)) {
    cn <- names(x$results$specials[[j]])
    if(!is.factor(x$results$specials[[j]][[1L]])) {
      matplot(x$results$specials[[j]][, 1L],
        as.matrix(x$results$specials[[j]][, c("fit", "lower", "upper")]),
        xlab = cn[1L], ylab = j, type = "l",
        lty = c(1, 2, 2), col = 1, ...)
    } else {
      z <- list("stats" = t(x$results$specials[[j]][, c("lower", "lower", "fit", "upper", "upper")]),
        "n" = as.numeric(table(x$results$specials[[j]][[1L]])),
        "names" = levels(x$results$specials[[j]][[1L]]))
      bxp(z, ...)
    }
  }
}

## Extract coefficients.
coef.gamlss2 <- function(object, full = FALSE, drop = TRUE, ...)
{
  co <- object$coefficients

  cos <- list()
  for(i in family(object)$names) {
    if(!is.null(co[[i]])) {
      cos[[i]] <- list("p" = co[[i]])
    }
    if(!is.null(object$fitted.specials[[i]]) & full) {
      if(is.list(cos[[i]]))
        cos[[i]]$s <- list()
      else
        cos[[i]] <- list("s" = list())
      for(j in names(object$fitted.specials[[i]])) {
        cij <- object$fitted.specials[[i]][[j]]$coefficients
        if(!is.null(cij)) {
          names(cij) <- as.character(1:length(cij))
          cos[[i]]$s[[j]] <- cij
        }
      }
      if(length(cos[[i]]$s) < 1L)
        cos[[i]]$s <- NULL
    }
  }
  co <- if(drop) unlist(cos) else cos

  return(co)
}

## Variance-covariance matrix.
vcov.gamlss2 <- function(object, full = FALSE, ...)
{
  y <- if(is.null(object$y)) {
    model.response(model.frame(object, keepresponse = TRUE))
  } else {
    object$y
  }
  x <- if(is.null(object[["x"]])) {
    model.matrix(object)
  } else {
    object[["x"]]
  }

  n <- if(is.null(dim(y))) length(y) else nrow(y)

  par <- coef(object, full = TRUE, drop = TRUE)
  lpar <- par2list(par)

  family <- object$family
  nx <- family$names
  eta <- list()
  for(i in nx) {
    eta[[i]] <- rep(0, n)
    if(!full) {
      if(!is.null(lpar[[i]]$s)) {
        for(j in names(lpar[[i]]$s)) {
          fit <- drop(object$specials[[j]]$X %*% lpar[[i]]$s[[j]])
          if(!is.null(object$specials[[j]]$binning)) {
            fit <- fit[object$specials[[j]]$binning$match.index]
          }
          eta[[i]] <- eta[[i]] + fit
        }
      }
    }
  }

  loglik <- function(par) {
    par <- par2list(par)
    for(i in nx) {
      if(!is.null(par[[i]]$p)) {
        eta[[i]] <- drop(x[, names(par[[i]]$p), drop = FALSE] %*% par[[i]]$p)
      }
      if(full) {
        if(!is.null(par[[i]]$s)) {
          for(j in names(par[[i]]$s)) {
            fit <- drop(object$specials[[j]]$X %*% par[[i]]$s[[j]])
            if(!is.null(object$specials[[j]]$binning)) {
              fit <- fit[object$specials[[j]]$binning$match.index]
            }
            eta[[i]] <- eta[[i]] + fit
          }
        }
      }
    }
    ll <- family$loglik(y, family$map2par(eta))
    return(ll)
  }

  par <- coef(object, full = full, drop = TRUE)

  v <- solve(as.matrix(optimHess(par, fn = loglik, control = list(fnscale = -1))))

  return(v)
}

## Little helper function.
par2list <- function(par)
{
  np <- names(par)
  pn <- strsplit(np, ".", fixed = TRUE)
  nx <- sapply(pn, function(x) x[1])
  ps <- sapply(pn, function(x) x[2])
  lab <- sapply(pn, function(x) x[3])
  lab[ps == "p"] <- ""
  pl <- list()
  for(i in unique(nx)) {
    if(any(j <- (ps == "p") & (nx == i))) {
      pj <- paste0(i, ".p.")
      pl[[i]]$p <- par[grep(pj, np, fixed = TRUE)]
      names(pl[[i]]$p) <- gsub(pj, "", names(pl[[i]]$p), fixed = TRUE)
    }
    if(any(j <- (ps == "s") & (nx == i))) {
      pl[[i]]$s <- list()
      for(k in unique(lab[j])) {
        pj <- paste0(i, ".s.", k, ".")
        pl[[i]]$s[[k]] <- par[grep(pj, np, fixed = TRUE)]
        names(pl[[i]]$s[[k]]) <- gsub(pj, "", names(pl[[i]]$s[[k]]), fixed = TRUE)
      }
    }
  }
  pl
}

## Summary extractor function.
summary.gamlss2 <- function(object, ...)
{
  df.res <- object$nobs - object$df
  v <- vcov(object, full = FALSE)
  par <- coef(object, full = FALSE)
  se <- sqrt(abs(diag(v)))
  tvalue <- par / se
  pvalue <- 2 * pt(-abs(tvalue), df.res)
  ct <- cbind(par, se, tvalue, pvalue)
  dimnames(ct) <- list(names(par),
    c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  nx <- unique(sapply(strsplit(names(par), ".", fixed = TRUE),
    function(x) x[1L]))
  ctl <- list()
  for(i in nx) {
    ctl[[i]] <- ct[grep(paste0(i, "."), rownames(ct), fixed = TRUE), , drop = FALSE]
    rownames(ctl[[i]]) <- gsub(paste0(i, ".p."), "", rownames(ctl[[i]]), fixed = TRUE)
  }
  sg <- object[c("call", "family", "df", "nobs", "logLik", "iterations", "elapsed")]
  sg$call[[1L]] <- as.name("gamlss2")
  sg$coefficients <- ctl
  if(!is.null(object$fitted.specials)) {
    sg$specials <- list()
    for(i in names(object$fitted.specials)) {
      for(j in names(object$fitted.specials[[i]])) {
        sg$specials[[i]] <- rbind(sg$specials[[i]], object$fitted.specials[[i]][[j]]$edf)
      }
      rownames(sg$specials[[i]]) <- names(object$fitted.specials[[i]])
      colnames(sg$specials[[i]]) <- "edf"
    }
  }
  sg$elapsed <- object$elapsed
  class(sg) <- "summary.gamlss2"
  return(sg)
}

## Summary printing method.
print.summary.gamlss2 <- function(x,
  digits = max(3L, getOption("digits") - 3L),
  symbolic.cor = x$symbolic.cor, 
  signif.stars = getOption("show.signif.stars"), ...)
{
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("---\n")
  print(x$family, full = FALSE)
  for(i in x$family$names) {
    cat("*--------\n")
    cat("Parameter:", i, "\n")
    if(length(x$coefficients[[i]])) {
      cat("---\nCoefficients:\n")
      printCoefmat(x$coefficients[[i]], digits = digits,
        signif.stars = signif.stars, na.print = "NA",
        signif.legend = i == tail(x$family$names, 1L), ...)
    }
    if(length(x$specials[[i]])) {
      cat("---\nSmooth terms:\n")
      printCoefmat(t(x$specials[[i]]))
    }
  }
  cat("*--------\n")
  info1 <- c(
    paste("n =", x$nobs),
    paste("df = ", round(x$df, digits = 2)),
    paste("res.df = ", round(x$nobs - x$df, digits = 2))
  )
  info2 <- c(
    ## paste("logLik =", round(x$logLik, digits = 4)),
    paste("Deviance =", round(-2 * x$logLik, digits = 4)),
    paste("AIC =", round(-2 * x$logLik + 2*x$df, digits = 4))
  )
  rt <- x$elapsed
  rt <- if(rt > 60) {
    paste(formatC(format(round(rt / 60, 2), nsmall = 2), width = 5), "min", sep = "")
  } else paste(formatC(format(round(rt, 2), nsmall = 2), width = 5), "sec", sep = "")
  info3 <- c(
    paste("elapsed =", rt)
  )
  cat(info1)
  cat("\n")
  cat(info2)
  cat("\n")
  cat(info3)
  cat("\n")
}

## Confint method.
confint.gamlss2 <- function(object, parm, level = 0.95, ...)
{
  co <- coef(object, full = FALSE, drop = TRUE)
  v <- vcov(object, full = FALSE)
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  se <- sqrt(abs(diag(v)))
  ci <- co + se %o% qnorm(a)
  colnames(ci) <- paste0(round(a * 100, 3), "%")
  rownames(ci) <- gsub(".p.", ".", rownames(ci), fixed = TRUE)
  if(!missing(parm)) {
    ci <- ci[unique(grep2(parm, rownames(ci), value = TRUE, fixed = TRUE)), ]
  }
  return(ci)
}


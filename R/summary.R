## Extract coefficients.
coef.gamlss2 <- function(object, full = FALSE, drop = TRUE, ...)
{
  co <- object$coefficients

  family <- family(object)
  if(is.null(family))
    family <- list("names" = names(co))

  cos <- list()
  for(i in family$names) {
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
        lambdas <- NULL
        if(isTRUE(list(...)$lambdas)) {
          lambdas <- object$fitted.specials[[i]][[j]]$lambdas
          names(lambdas) <- paste0("lambda", 1:length(lambdas))
        }
        if(!is.null(cij)) {
          if(is.null(names(cij))) {
            names(cij) <- as.character(1:length(cij))
          }
          cos[[i]]$s[[j]] <- cij
        }
        if(!is.null(lambdas))
          cos[[i]]$s[[j]] <- c(cos[[i]]$s[[j]], lambdas)
      }
      if(length(cos[[i]]$s) < 1L)
        cos[[i]]$s <- NULL
    }
  }

  model <- list(...)$model
  if(is.null(model)) {
    what <- list(...)$what
    if(!is.null(what))
      model <- what
  }
  if(!is.null(model)) {
    if(!is.character(model))
      model <- family$names[model]
    model <- family$names[pmatch(model, family$names)]
    cos <- cos[model]    
  }

  if(drop) {
    nc <- names(cos)
    cos <- unlist(cos)
    if((length(nc) < 2L) & !isFALSE(list(...)$dropall) & !is.null(cos)) {
      names(cos) <- gsub(paste0(nc, "."), "", names(cos), fixed = TRUE)
      for(j in c("p.", "s."))
        names(cos) <- gsub(j, "", names(cos), fixed = TRUE)
    }
  }

  if(is.null(cos))
    cos <- numeric(0)

  class(cos) <- "coef.gamlss2"

  return(cos)
}

print.coef.gamlss2 <- function(x, ...)
{
  print(unclass(x), ...)
}

## Variance-covariance matrix.
vcov.gamlss2 <- function(object, type = c("vcov", "cor", "se", "coef"), full = FALSE, ...)
{
  type <- match.arg(type)

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

  par <- coef(object, full = TRUE, drop = TRUE, dropall = FALSE)
  lpar <- par2list(par)

  family <- object$family
  nx <- family$names

  loglik <- function(par) {
    par <- par2list(par)
    eta <- list()
    for(i in nx) {
      if(!is.null(par[[i]]$p)) {
        eta[[i]] <- drop(x[, names(par[[i]]$p), drop = FALSE] %*% par[[i]]$p)
      }
      if(full | TRUE) { ## FIXME: always add?
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
    ll <- family$logLik(y, family$map2par(eta))
    return(ll)
  }

  gradient <- function(par) {
    npar <- names(par)
    par <- par2list(par)
    eta <- list()
    for(i in nx) {
      if(!is.null(par[[i]]$p)) {
        eta[[i]] <- drop(x[, names(par[[i]]$p), drop = FALSE] %*% par[[i]]$p)
      }
      if(full | TRUE) { ## FIXME: always add?
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
    par2 <- family$map2par(eta)
    g <- NULL
    for(i in nx) {
      score <- family$score[[i]](y, par2)
      if(!is.null(par[[i]]$p)) {
        g <- c(g, colSums(x[, names(par[[i]]$p), drop = FALSE] * score))
      }
      if(full | TRUE) { ## FIXME: always add?
        if(!is.null(par[[i]]$s)) {
          for(j in names(par[[i]]$s)) {
            if(!is.null(object$specials[[j]]$binning)) {
              fit <- fit[object$specials[[j]]$binning$match.index]
              g <- c(g, colSums(object$specials[[j]]$X[object$specials[[j]]$binning$match.index, ] * score))
            } else {
              g <- c(g, colSums(object$specials[[j]]$X * score))
            }
          }
        }
      }
    }
    names(g) <- npar
    return(g)
  }

  par <- coef(object, full = full, drop = TRUE, dropall = FALSE)

  ## Set NA coefficients to zero.
  p_na <- is.na(par)
  par[p_na] <- 0.0

  if(type == "coef")
    return(par)

  control <- list(...)
  control$fnscale <- -1

  H <- as.matrix(optimHess(par, fn = loglik, gr = gradient, control = control))
  if(ncol(H) > 1L) {
    v <- try(solve(H), silent = TRUE)
    if(inherits(v, "try-error")) {
      H <- H + diag(1e-05, ncol(H))
      v <- try(solve(H), silent = TRUE)
      if(inherits(v, "try-error")) {
        H <- H + diag(1e-03, ncol(H))
        v <- solve(H)
      }
    }
  } else {
    v <- 1 / H
  }

  v <- -v

  if(type == "cor") {
    dv <- sqrt(diag(v))
    v <- v / (dv %*% t(dv))
  }

  if(type == "se")
    v <- sqrt(diag(v))

  v[p_na, p_na] <- NA

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
  par <- coef(object, full = FALSE, dropall = FALSE)
  se <- sqrt(abs(diag(v)))
  tvalue <- par / se
  pvalue <- 2 * pt(-abs(tvalue), df.res)
  ct <- cbind(par, se, tvalue, pvalue)
  dimnames(ct) <- list(names(par),
    c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ctl <- NULL
  if(nrow(ct) > 0L) {
    nx <- unique(sapply(strsplit(names(par), ".", fixed = TRUE),
      function(x) x[1L]))
    ctl <- list()
    for(i in nx) {
      ctl[[i]] <- ct[grep(paste0(i, "."), rownames(ct), fixed = TRUE), , drop = FALSE]
      rownames(ctl[[i]]) <- gsub(paste0(i, ".p."), "", rownames(ctl[[i]]), fixed = TRUE)
    }
  }
  sg <- object[c("call", "family", "df", "nobs", "logLik", "dev.reduction", "iterations", "elapsed")]
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

## For MCMC samples.
summary.gamlss2.mcmc <- function(object, thres = 0.01, ...)
{
  df.res <- object$nobs - object$df
  par <- coef(object, full = FALSE, dropall = FALSE)
  ct <- t(apply(object$samples[, names(par), drop = FALSE], 2, function(x) {
    c("Mean" = mean(x, na.rm =TRUE), quantile(x, probs = c(0.025, 0.975)), mean(abs(x) < thres))
  }))
  colnames(ct)[ncol(ct)] <- paste0("Pr(<|", thres, "|)")
  nx <- unique(sapply(strsplit(names(par), ".", fixed = TRUE),
    function(x) x[1L]))
  ctl <- list()
  for(i in nx) {
    ctl[[i]] <- ct[grep(paste0(i, "."), rownames(ct), fixed = TRUE), , drop = FALSE]
    rownames(ctl[[i]]) <- gsub(paste0(i, ".p."), "", rownames(ctl[[i]]), fixed = TRUE)
  }
  sg <- object[c("call", "family", "df", "nobs", "logLik", "dev.reduction", "iterations", "elapsed")]
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

  ## Collect linear coefficients and specials.
  lin_coef <- specials <- list()
  pn <- x$family$names
  for(i in seq_along(pn)) {
    if(length(x$coefficients[[pn[i]]])) {
      lin_coef[[i]] <- x$coefficients[[pn[i]]]
      rownames(lin_coef[[i]]) <- paste0(pn[i], ".", rownames(lin_coef[[i]]))
    }
    if(length(x$specials[[pn[i]]])) {
      specials[[i]] <- x$specials[[pn[i]]]
      rownames(specials[[i]]) <- paste0(pn[i], ".", rownames(specials[[i]]))
    }
  }

  cat("*--------\n")
  if(length(lin_coef)) {
    lin_coef <- do.call("rbind", lin_coef)
    cat("Coefficients:\n")
    printCoefmat(lin_coef, digits = digits,
      signif.stars = signif.stars, na.print = "NA",
      signif.legend = TRUE, ...)
  }

  if(length(specials)) {
    specials <- do.call("rbind", specials)
    cat("---\nSmooth terms:\n")
    printCoefmat(t(specials))
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
    paste0("Null Dev. Red. = ", round(x$dev.reduction * 100, digits = 2), "%")
  )

  rt <- x$elapsed
  rt <- if(rt > 60) {
    paste(formatC(format(round(rt / 60, 2), nsmall = 2), width = 5), "min", sep = "")
  } else paste(formatC(format(round(rt, 2), nsmall = 2), width = 5), "sec", sep = "")
  info3 <- c(
    paste("AIC =", round(-2 * x$logLik + 2*x$df, digits = 4)),
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

## R2.
.R2 <- function(object, type = c("Cox Snell", "Cragg Uhler", "both", "simple"), newdata = NULL, ...)
{
  type <- match.arg(type)

  y <- if (!is.null(newdata)) {
    model.response(model.frame(object, data = newdata, keepresponse = TRUE))
  } else {
    model.response(model.frame(object, keepresponse = TRUE))
  }

  ## FIXME: bd?

  if(type == "simple") {
    if(!is.null(family(object)$type)) {
      if(family(object)$type != "continuous")
        stop("R-squared only for continuous responses!")
    }

    par <- fitted(object, newdata = newdata, type = "parameter")

    if(is.null(family(object)$mean)) {
      fit <- family(object)$quantile(0.5, par)
    } else {
      fit <- family(object)$mean(par)
    }

    nobs <- length(fit)
    rsq <- 1 - var(y - fit) * (nobs - 1)/(var(y - mean(y)) * (nobs - object$df))
  } else {
    ll0 <- object$null.deviance / -2
    ll1 <- object$deviance / -2
    rsq1 <- 1 - exp((2/object$nobs) * (ll0 - ll1))
    rsq2 <- rsq1 / (1 - exp((2/object$nobs) * ll0))
    if(type == "Cox Snell")
      rsq <- rsq1
    if(type == "Cragg Uhler")
      rsq <- rsq2
    if(type == "both")
      rsq <- list("CoxSnell" = rsq1, "CraggUhler" = rsq2)
  }

  return(rsq)
}

## New r-squared function.
Rsq <- function(object, ..., type = c("Cox Snell", "Cragg Uhler", "both", "simple"), newdata = NULL)
{
  objs <- list(object, ...)
  rsq <- NULL
  drop <- NULL
  for(j in 1:length(objs)) {
    if(inherits(objs[[j]], c("gamlss", "gamlss2"))) {
      if(inherits(objs[[j]], "gamlss")) {
        type2 <- unique(gsub("simple", "Cox Snell", type))
        rsq <- c(rsq, Rsq_gamlss(objs[[j]], type = type2))
      } else {
        rsq <- c(rsq, Rsq_gamlss2(objs[[j]], type = type, newdata = newdata))
      }
    } else {
       drop <- c(drop, j)
    }
  }
  if(length(drop)) {
    rsq <- rsq[-drop]
  }
  if(length(rsq) & !is.list(rsq)) {
    if(length(rsq) > 1L) {
      Call <- match.call()
      names(rsq) <- as.character(Call[-1L])
    }
  }
  if(!is.list(rsq)) {
    if(is.vector(rsq))
      rsq <- sort(rsq)
  }
  return(rsq)
}

Rsq_gamlss <- function(object, ...)
{
  utils::getFromNamespace("Rsq", "gamlss")(object, ...)
}

Rsq_gamlss2 <- function(object, ...)
{
  .R2(object, ...)
}

## New GAIC function.
GAIC <- function(object, ..., k = 2, corrected = FALSE)
{
  objs <- list(object, ...)
  gaic <- edf <- NULL
  drop <- NULL
  for(j in 1:length(objs)) {
    if(inherits(objs[[j]], c("gamlss", "gamlss2"))) {
      if(inherits(objs[[j]], "gamlss")) {
        df <- objs[[j]]$df.fit
        N <- objs[[j]]$N
        deviance <- objs[[j]]$G.deviance
      } else {
        df <- objs[[j]]$df
        N <- objs[[j]]$nobs
        deviance <- objs[[j]]$deviance
      }
      Cor <- if((k == 2) && corrected) {
        (2 * df * (df + 1))/(N - df - 1)
      } else {
        0
      }
      gaic <- c(gaic, deviance + df * k + Cor)
      edf <- c(edf, df)
    } else {
      drop <- c(drop, j)
      gaic <- c(gaic, NA)
      edf <- c(edf, NA)
    }
  }

  if(length(drop)) {
    gaic <- gaic[-drop]
    edf <- edf[-drop]
  }
  if(length(gaic) > 1L) {
    Call <- match.call()
    rn <- as.character(Call[-1L])
    if("..1" %in% rn) {
      rn <- as.character(sys.call(-1))[-1L]
    }
    if(length(drop))
      rn <- rn[-drop]
    i <- order(gaic, decreasing = FALSE)
    gaic <- data.frame("AIC" = gaic[i], "df" = edf[i])
    rn <- rn[i]
    if(any(j <- duplicated(rn)))
      rn[j] <- paste0(rn[j], ".", 1:sum(j))
    rownames(gaic) <- rn
  }
  return(gaic)
}

## Generic AIC/BIC method.
AIC.gamlss2 <- function(object, ..., k = 2, corrected = FALSE)
{
  GAIC(object, ..., k = k, corrected = corrected)
}

BIC.gamlss2 <- function(object, ...)
{
  GAIC(object, ..., k = log(object$nobs), corrected = FALSE)
}


get_df2 <- function(object) {
  df <- 0
  if(length(object$fitted.linear))
    df <- df + sum(sapply(object$fitted.linear, function(x) length(x$coefficients)))
  if(length(object$fitted.specials)) {
    for(j in seq_along(object$fitted.specials)) {
      dfj <- sapply(object$fitted.specials[[j]], function(x) x$edf)
      df <- df + sum(unlist(dfj))
    }
  }
  return(df)
}

stepwise <- function(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
{ 
  nx <- family$names

  for(i in nx) {
    if(!any(grepl("(Intercept)", xterms[[i]]))) {
      stop(paste0("intercept is missing for parameter ", i, "!"))
    }
  }

  ## Set optimizer.
  

  ## Penalty for AIC.
  K <- if(is.null(control$K)) log(nrow(x)) else control$K

  ## Extract factor levels. FIXME?
  ## xlev <- attr(xterms, "xlevels")

  xsel <- ssel <- list()
  for(i in nx) {
    xsel[[i]] <- rep(FALSE, length(xterms[[i]]))
    names(xsel[[i]]) <- xterms[[i]]
    ssel[[i]] <- rep(FALSE, length(sterms[[i]]))
    names(ssel[[i]]) <- sterms[[i]]
    xsel[[i]][names(xsel[[i]]) == "(Intercept)"] <- TRUE
  }

  ## (1) Forward linear only.
  ici <- Inf
  ll <- -Inf
  df <- 0
  first <- TRUE
  do <- TRUE
  while(do) {
    ic <- list()
    xterms_sel <- lapply(xsel, function(x) names(x)[x])
    xterms_nosel <- lapply(xsel, function(x) names(x)[!x])
    k <- 0L
    for(i in nx) {
      ic[[i]] <- list()
      if(length(xterms_nosel[[i]])) {
        for(j in xterms_nosel[[i]]) {
          xtermsi <- xterms_sel
          xtermsi[[i]] <- c(xtermsi[[i]], j)
          m <- RS(x, y, specials = NULL, family, offsets,
            weights, start = start, xtermsi, sterms = NULL, control)
          dfj <- get_df2(m)
          ic[[i]][[j]] <- list("AIC" = -2 * m$logLik + K * dfj, "logLik" = m$logLik, "df" = dfj)
          if(first) {
            ici[1L] <- m$null.deviance + K * length(nx)
            ll[1L] <- m$null.deviance / -2
            df[1L] <- length(nx)
            names(ici)[1L] <- names(ll)[1L] <- names(df)[1L] <- "nullmodel"
            first <- FALSE
            start <- m$fitted.values
          }
          k <- k + 1L
        }
      }
    }

    if(k > 0L) {
      i <- lapply(ic, function(x) {
        x <- sapply(x, function(z) z$AIC)
        if(length(x)) {
          return(min(unlist(x)))
        } else {
          return(Inf)
        }
      })
      i <- nx[which.min(unlist(i))]
      j <- names(ic[[i]])[which.min(sapply(ic[[i]], function(z) z$AIC))]

      ici <- c(ici, ic[[i]][[j]]$AIC)
      ll <- c(ll, ic[[i]][[j]]$logLik)
      df <- c(df, ic[[i]][[j]]$df)
      lici <- length(ici)
      names(ici)[lici] <- names(ll)[lici] <- names(df)[lici] <- paste0(i, ".", j)

      if(length(ici) > 1L) {
        tstat <- 2 * (ll[length(ll)] - ll[length(ll) - 1L])
        dfij <- df[length(df)] - df[length(df) - 1L]
        pval <- pchisq(tstat, df = dfij, lower.tail = FALSE)

        if((ici[length(ici)] >= ici[length(ici) - 1L]) | (pval >= 0.05)) {
          ll <- ll[-length(ll)]
          ici <- ici[-length(ici)]
          df <- df[-length(df)]
          do <- FALSE
        }
      }

      if(do) {
        xsel[[i]][j] <- TRUE

        xterms_sel2 <- lapply(xsel, function(x) names(x)[x])

        m <- RS(x, y, specials = NULL, family, offsets,
          weights, start = start, xterms_sel2, sterms = NULL, control)

        start <- m$fitted.values

        cat("  <+> parameter", i, "term", j, "\n")
      }
    } else {
      do <- FALSE
    }
  }

  xterms_sel2 <- lapply(xsel, function(x) names(x)[x])

  m <- RS(x, y, specials = NULL, family, offsets,
    weights, start = start, xterms_sel2, sterms = NULL, control)

  start <- m$fitted.values

  ## (2) Forward nonlinear.
  if(length(sterms)) {
    sterms_sel2 <- NULL
    do <- TRUE
    while(do) {
      ic <- list()
      sterms_sel <- lapply(ssel, function(x) names(x)[x])
      sterms_nosel <- lapply(ssel, function(x) names(x)[!x])
      k <- 0L
      for(i in nx) {
        ic[[i]] <- list()
        if(length(sterms_nosel[[i]])) {
          for(j in sterms_nosel[[i]]) {
            stermsi <- sterms_sel
            stermsi[[i]] <- c(stermsi[[i]], j)
            m <- RS(x, y, specials, family, offsets,
              weights, start = start, xterms_sel2, sterms = stermsi, control)
            dfj <- get_df2(m)
            ic[[i]][[j]] <- list("AIC" = -2 * m$logLik + K * dfj, "logLik" = m$logLik, "df" = dfj)
            k <- k + 1L
          }
        }
      }

      if(k > 0L) {
        i <- lapply(ic, function(x) {
          x <- sapply(x, function(z) z$AIC)
          if(length(x)) {
            return(min(unlist(x)))
          } else {
            return(Inf)
          }
        })
        i <- nx[which.min(unlist(i))]
        j <- names(ic[[i]])[which.min(sapply(ic[[i]], function(z) z$AIC))]

        ici <- c(ici, ic[[i]][[j]]$AIC)
        ll <- c(ll, ic[[i]][[j]]$logLik)
        df <- c(df, ic[[i]][[j]]$df)
        lici <- length(ici)
        names(ici)[lici] <- names(ll)[lici] <- names(df)[lici] <- paste0(i, ".", j)

        if(length(ici) > 1L) {
          tstat <- 2 * (ll[length(ll)] - ll[length(ll) - 1L])
          dfij <- df[length(df)] - df[length(df) - 1L]
          pval <- pchisq(tstat, df = dfij, lower.tail = FALSE)

          if((ici[length(ici)] >= ici[length(ici) - 1L]) | (pval >= 0.05)) {
            ll <- ll[-length(ll)]
            ici <- ici[-length(ici)]
            df <- df[-length(df)]
            do <- FALSE
          }
        }

        if(do) {
          ssel[[i]][j] <- TRUE

          sterms_sel2 <- lapply(ssel, function(x) names(x)[x])

          m <- RS(x, y, specials, family, offsets,
            weights, start = start, xterms_sel2, sterms = sterms_sel2, control)

          start <- m$fitted.values

          cat("  <+> parameter", i, "term", j, "\n")
        }
      } else {
        do <- FALSE
      }
    }
  }

  sterms_sel2 <- lapply(ssel, function(x) names(x)[x])

  m <- RS(x, y, specials, family, offsets,
    weights, start = start, xterms_sel2, sterms = sterms_sel2, control)

  ## (3) Mixed.
  selected <- NULL
  for(i in names(xterms_sel2))
    selected <- c(selected, paste0(i, ".p.", xterms_sel2[[i]]))
  for(i in names(sterms_sel2))
    selected <- c(selected, paste0(i, ".s.", sterms_sel2[[i]]))
 
  allterms <- NULL
  for(i in names(xterms_sel2))
    allterms <- c(allterms, paste0(i, ".p.", xterms[[i]]))
  for(i in names(sterms_sel2))
    allterms <- c(allterms, paste0(i, ".s.", sterms[[i]]))

  notselected <- allterms[!(allterms %in% selected)]

  do <- length(notselected) > 0L
  while(do) {
    anysel <- FALSE

    ## Backward step.
    ic0 <- ici[length(ici)]
    ic <- list()
    for(sel in selected) {
      if(grepl("(Intercept)", sel, fixed = TRUE))
        next

      xterms_sel3 <- xterms_sel2
      if(grepl(".p.", sel, fixed = TRUE)) {
        j <- strsplit(sel, ".p.", fixed = TRUE)[[1]]

        i <- j[1L]
        j <- j[2L]

        xterms_sel3[[i]] <- xterms_sel3[[i]][-grep(j, xterms_sel3[[i]], fixed = TRUE)]
      }

      sterms_sel3 <- sterms_sel2
      if(grepl(".s.", sel, fixed = TRUE)) {
        j <- strsplit(sel, ".s.", fixed = TRUE)[[1]]

        i <- j[1L]
        j <- j[2L]

        sterms_sel3[[i]] <- sterms_sel3[[i]][-grep(j, sterms_sel3[[i]], fixed = TRUE)]
      }

      if(is.na(j))
        next

      m <- RS(x, y, specials, family, offsets,
        weights, start = start, xterms_sel3, sterms = sterms_sel3, control)

      ic[[i]][[j]] <- list("AIC" = -2 * m$logLik + K * dfj, "logLik" = m$logLik, "df" = dfj, "term" = sel)
    }

    i <- lapply(ic, function(x) {
      x <- sapply(x, function(z) z$AIC)
      if(length(x)) {
        return(min(unlist(x)))
      } else {
        return(Inf)
      }
    })

    i <- nx[which.min(unlist(i))]
    j <- names(ic[[i]])[which.min(sapply(ic[[i]], function(z) z$AIC))]

    if(ic[[i]][[j]]$AIC < ic0) {
      anysel <- TRUE
      sel <- ic[[i]][[j]]$term
      if(grepl(".p.", sel, fixed = TRUE)) {
        j <- strsplit(sel, ".p.", fixed = TRUE)[[1]]

        i <- j[1L]
        j <- j[2L]

        xterms_sel2[[i]] <- xterms_sel2[[i]][-grep(j, xterms_sel2[[i]], fixed = TRUE)]
      }

      if(grepl(".s.", sel, fixed = TRUE)) {
        j <- strsplit(sel, ".s.", fixed = TRUE)[[1]]

        i <- j[1L]
        j <- j[2L]

        sterms_sel2[[i]] <- sterms_sel2[[i]][-grep(j, sterms_sel2[[i]], fixed = TRUE)]
      }

      ici <- c(ici, ic[[i]][[j]]$AIC)
      ll <- c(ll, ic[[i]][[j]]$logLik)
      df <- c(df, ic[[i]][[j]]$df)

      lici <- length(ici)
      names(ici)[lici] <- names(ll)[lici] <- names(df)[lici] <- paste0(i, ".", j)

      ic0 <- ici[length(ici)]

      m <- RS(x, y, specials, family, offsets,
        weights, start = start, xterms_sel2, sterms = sterms_sel2, control)

      start <- m$fitted.values

      cat("  <-> parameter", i, "term", j, "\n")
    }

    ## Forward step.
    selected <- NULL
    for(i in names(xterms_sel2))
      selected <- c(selected, paste0(i, ".p.", xterms_sel2[[i]]))
    for(i in names(sterms_sel2))
      selected <- c(selected, paste0(i, ".s.", sterms_sel2[[i]]))
    notselected <- allterms[!(allterms %in% selected)]

    ic <- list()
    for(sel in notselected) {
      xterms_sel3 <- xterms_sel2
      if(grepl(".p.", sel, fixed = TRUE)) {
        j <- strsplit(sel, ".p.", fixed = TRUE)[[1]]

        i <- j[1L]
        j <- j[2L]

        xterms_sel3[[i]] <- c(xterms_sel3[[i]], j)
      }

      sterms_sel3 <- sterms_sel2
      if(grepl(".s.", sel, fixed = TRUE)) {
        j <- strsplit(sel, ".s.", fixed = TRUE)[[1]]

        i <- j[1L]
        j <- j[2L]

        sterms_sel3[[i]] <- c(sterms_sel3[[i]], j)
      }

      if(is.na(j))
        next

      m <- RS(x, y, specials, family, offsets,
        weights, start = start, xterms_sel3, sterms = sterms_sel3, control)

      ic[[i]][[j]] <- list("AIC" = -2 * m$logLik + K * dfj, "logLik" = m$logLik, "df" = dfj, "term" = sel)
    }

    i <- lapply(ic, function(x) {
      x <- sapply(x, function(z) z$AIC)
      if(length(x)) {
        return(min(unlist(x)))
      } else {
        return(Inf)
      }
    })

    i <- nx[which.min(unlist(i))]
    j <- names(ic[[i]])[which.min(sapply(ic[[i]], function(z) z$AIC))]

    if(ic[[i]][[j]]$AIC < ic0) {
      anysel <- TRUE
      sel <- ic[[i]][[j]]$term
      if(grepl(".p.", sel, fixed = TRUE)) {
        j <- strsplit(sel, ".p.", fixed = TRUE)[[1]]

        i <- j[1L]
        j <- j[2L]

        xterms_sel2[[i]] <- xterms_sel2[[i]][-grep(j, xterms_sel2[[i]], fixed = TRUE)]
      }

      if(grepl(".s.", sel, fixed = TRUE)) {
        j <- strsplit(sel, ".s.", fixed = TRUE)[[1]]

        i <- j[1L]
        j <- j[2L]

        sterms_sel2[[i]] <- sterms_sel2[[i]][-grep(j, sterms_sel2[[i]], fixed = TRUE)]
      }

      ici <- c(ici, ic[[i]][[j]]$AIC)
      ll <- c(ll, ic[[i]][[j]]$logLik)
      df <- c(df, ic[[i]][[j]]$df)

      lici <- length(ici)
      names(ici)[lici] <- names(ll)[lici] <- names(df)[lici] <- paste0(i, ".", j)

      ic0 <- ici[length(ici)]

      m <- RS(x, y, specials, family, offsets,
        weights, start = start, xterms_sel2, sterms = sterms_sel2, control)

      start <- m$fitted.values

      cat("  <+> parameter", i, "term", j, "\n")
    }

    do <- anysel
  }

  m <- RS(x, y, specials, family, offsets,
    weights, start = start, xterms_sel2, sterms = sterms_sel2, control)

  m$selection <- list("AIC" = ici, "logLik" = ll, "df" = df)
  m$xterms <- xterms_sel2
  m$sterms <- sterms_sel2
  m$specials <- specials[unique(unlist(sterms_sel2))]

  return(m)
}


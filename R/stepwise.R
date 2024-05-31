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

  ## Models must include intercepts!
  for(i in nx) {
    if(!any(grepl("(Intercept)", xterms[[i]]))) {
      stop(paste0("intercept is missing for parameter ", i, "!"))
    }
  }

  ## Which strategies to use?
  strategy <- control$strategy
  choices <- c("forward.linear", "forward", "backward", "replace")
  if(is.null(strategy)) {
    strategy <- choices
  } else {
    strategy <- match.arg(strategy, choices, several.ok = TRUE)
  }

  ## Inner and outer trace.
  trace <- control$trace
  if(length(trace) < 2L) {
    trace <- c(FALSE, trace)
  }
  control$trace <- trace[1L]

  ## Penalty for AIC.
  K <- if(is.null(control$K)) log(nrow(x)) else control$K

  ## Extract factor levels. FIXME?
  ## xlev <- attr(xterms, "xlevels")

  ## Setup vectors for selected and not selected model terms.
  ## Start with linear terms first.
  selected <- NULL
  notselected <- NULL
  for(i in nx) {
    if(length(xterms[[i]]))
      notselected <- c(notselected, paste0(i, ".p.", xterms[[i]]))
  }
  notselected <- notselected[-grep(".p.(Intercept)", notselected, fixed = TRUE)]

  ## Linear or nonlinear?
  linear <- function(x) {
    grepl(".p.", x, fixed = TRUE)
  }

  ## Helper function for splitting term names.
  splitname <- function(x) {
    split <- if(linear(x)) ".p." else ".s."
    x <- strsplit(x, split, fixed = TRUE)[[1L]]
    return(x)
  }

  ## Helper function to add to xterms/sterms.
  addterm <- function(x, y) {
    x <- splitname(x)
    y[[x[1L]]] <- c(y[[x[1L]]], x[2L])
    return(y)
  }

  ## Helper function to drop from xterms/sterms.
  dropterm <- function(x, y) {
    x <- splitname(x)
    y[[x[1L]]] <- y[[x[1L]]][y[[x[1L]]] != x[2L]]
    return(y)
  }

  ## Helper function to extract GAIC etc.
  modelstats <- function(model) {
    stats <- list("df" = get_df2(model))
    stats$GAIC <- model$deviance + K * stats$df
    stats$logLik <- model$logLik
    return(stats)
  }

  ## Estimate nullmodel first.
  xterms_itcpt <- list()
  for(i in nx)
    xterms_itcpt[[i]] <- "(Intercept)"

  m <- RS(x, y, specials = NULL, family, offsets,
    weights, start = start, xterms_itcpt, sterms = NULL, control)

  ## Update starting values.
  start <- fit <- m$fitted.values

  ## Nullmodel GAIC, logLik, df.
  stats_save <- list()
  stats_save[[1L]] <- modelstats(m)
  stats_save[[1L]]$term <- "nullmodel"

  ## Save initial GAIC.
  gaic <- stats_save[[1L]]$GAIC

  ## Set starting xterms.
  xterms <- xterms_itcpt

  ## Save number of iterations.
  iter <- 1L

  ## (1) Forward linear.
  if(length(notselected) & ("forward.linear" %in% strategy)) {
    if(trace[2L])
      cat("Forward Linear Selection\n")

    ## Start forward linear selection.
    do <- TRUE
    while(do) {
      ## List for term wise stats.
      stats <- list()

      ## Try all terms.
      for(j in notselected) {
        ## Add term.
        xtermsj <- addterm(j, xterms)

        ## Estimate model.
        m <- RS(x, y, specials = NULL, family, offsets,
          weights, start = start, xtermsj, sterms = NULL, control)

        ## Store stats.
        stats[[j]] <- modelstats(m)

        ## If GAIC improves, save fitted values for setting
        ## start later.
        if(stats[[j]]$GAIC < gaic) {
          gaic <- stats[[j]]$GAIC
          fit <- m$fitted.values
        }
      }

      ## Add term with minimum GAIC?
      j <- sapply(stats, function(x) x$GAIC)
      j <- names(j)[which.min(j)]

      ## Check if GAIC improved?
      if(stats[[j]]$GAIC < stats_save[[iter]]$GAIC) {
        ## Update linear terms.
        xterms <- addterm(j, xterms)

        ## Remove from not selected.
        notselected <- notselected[notselected != j]

        ## Add to selected.
        selected <- c(selected, j)

        ## Increase iterator.
        iter <- iter + 1L

        ## Save stats.
        stats_save[[iter]] <- stats[[j]]
        stats_save[[iter]]$term <- j

        ## Set starting values.
        start <- fit

        ## Print info.
        if(trace[2L]) {
          j <- splitname(j)
          cat("  GAIC = ", fmt(stats_save[[iter]]$GAIC, width = 10, digits = 4),
            " <+> parameter ", j[1L], ", term ", j[2L],  "\n", sep = "")
        }
      } else {
        do <- FALSE
      }
    }
  }

  ## (2) Forward all.
  for(i in nx) {
    if(length(sterms[[i]]))
      notselected <- c(notselected, paste0(i, ".s.", sterms[[i]]))
  }

  ## sterms start list().
  sterms <- list()
  for(i in nx)
    sterms[[i]] <- character(0)

  ## Max. iterations.
  n.iter <- control$n.iter
  if(is.null(n.iter))
    n.iter <- Inf
  gaic0 <- gaic
  improved <- TRUE
  k <- 0L
  while(improved & (k < n.iter)) {
    if(length(notselected) & ("forward" %in% strategy)) {
      if(trace[2L])
        cat("Forward Selection\n")

      do <- TRUE
      while(do) {
        ## List for term wise stats.
        stats <- list()

        ## Try all terms.
        for(j in notselected) {
          ## Add term.
          if(linear(j)) {
            xtermsj <- addterm(j, xterms)
            stermsj <- sterms
          } else {
            xtermsj <- xterms
            stermsj <- addterm(j, sterms)
          }

          ## Estimate model.
          m <- RS(x, y, specials = specials, family, offsets,
            weights, start = start, xtermsj, sterms = stermsj, control)

          ## Store stats.
          stats[[j]] <- modelstats(m)

          ## Save fitted values.
          if(stats[[j]]$GAIC < gaic) {
            gaic <- stats[[j]]$GAIC
            fit <- m$fitted.values
          }
        }

        ## Add term with minimum GAIC?
        j <- sapply(stats, function(x) x$GAIC)
        j <- names(j)[which.min(j)]

        ## Check if GAIC improved?
        if(stats[[j]]$GAIC < stats_save[[iter]]$GAIC) {
          ## Update terms.
          if(linear(j)) {
            xterms <- addterm(j, xterms)
          } else {
            sterms <- addterm(j, sterms)
          }

          ## Remove from not selected.
          notselected <- notselected[notselected != j]

          ## Add to selected.
          selected <- c(selected, j)

          ## Increase iterator.
          iter <- iter + 1L

          ## Save stats.
          stats_save[[iter]] <- stats[[j]]
          stats_save[[iter]]$term <- j

          ## Set starting values.
          start <- fit

          ## Print info.
          if(trace[2L]) {
            j <- splitname(j)
            cat("  GAIC = ", fmt(stats_save[[iter]]$GAIC, width = 10, digits = 4),
              " <+> parameter ", j[1L], ", term ", j[2L],  "\n", sep = "")
          }
        } else {
          do <- FALSE
        }
      }
    }

    ## (3) Replacements.
    if(length(notselected) & ("replace" %in% strategy)) {
      if(trace[2L])
        cat("Replacement Step\n")
      do <- TRUE
      while(do) {
        stats <- list()
        for(j in selected) {
          for(r in notselected) {
            ## Exchange terms.
            if(linear(j)) {
              xtermsj <- dropterm(j, xterms)
              stermsj <- sterms
            } else {
              xtermsj <- xterms
              stermsj <- dropterm(j, sterms)
            }

            if(linear(r)) {
              xtermsj <- addterm(r, xterms)
              stermsj <- sterms
            } else {
              xtermsj <- xterms
              stermsj <- addterm(r, sterms)
            }

            ## Estimate model.
            m <- RS(x, y, specials = specials, family, offsets,
              weights, start = start, xtermsj, sterms = stermsj, control)

            ## Store stats.
            stats[[j]] <- modelstats(m)
            stats[[j]]$old <- j
            stats[[j]]$new <- r

            ## Save fitted values.
            if(stats[[j]]$GAIC < gaic) {
              gaic <- stats[[j]]$GAIC
              fit <- m$fitted.values
            }
          }
        }

        ## Exchange term with minimum GAIC?
        j <- sapply(stats, function(x) x$GAIC)
        j <- names(j)[which.min(j)]

        ## Check if GAIC improved?
        if(stats[[j]]$GAIC < stats_save[[iter]]$GAIC) {
          ## Exchange terms.
          if(linear(j)) {
            xterms <- dropterm(j, xterms)
          } else {
            sterms <- dropterm(j, sterms)
          }
          if(linear(stats[[j]]$new)) {
            xterms <- addterm(stats[[j]]$new, xterms)
          } else {
            sterms <- addterm(stats[[j]]$new, sterms)
          }

          ## Remove from not selected.
          notselected <- c(notselected, j)
          notselected <- notselected[notselected != stats[[j]]$new]

          ## Add to selected.
          selected <- c(selected, stats[[j]]$new)
          selected <- selected[selected != j]

          ## Increase iterator.
          iter <- iter + 1L

          ## Save stats.
          stats_save[[iter]] <- stats[[j]]
          stats_save[[iter]]$term <- stats[[j]]$new

          ## Set starting values.
          start <- fit

          ## Print info.
          if(trace[2L]) {
            r <- splitname(stats[[j]]$old)
            j <- splitname(stats[[j]]$new)
            cat("  GAIC = ", fmt(stats_save[[iter]]$GAIC, width = 10, digits = 4),
              " <-> parameter ", r[1L], ", term ", r[2L],  "\n", sep = "")
            cat("                    <+> parameter ", j[1L], ", term ", j[2L],  "\n", sep = "")
          }
        } else {
          do <- FALSE
        }
      }
    }

    ## (4) Backward.
    if(length(selected) & ("backward" %in% strategy)) {
      if(trace[2L])
        cat("Backward Elimination\n")

      do <- TRUE
      while(do) {
        ## List for term wise stats.
        stats <- list()

        ## Try all terms.
        for(j in selected) {
          ## Add term.
          if(linear(j)) {
            xtermsj <- dropterm(j, xterms)
            stermsj <- sterms
          } else {
            xtermsj <- xterms
            stermsj <- dropterm(j, sterms)
          }

          ## Estimate model.
          m <- RS(x, y, specials = specials, family, offsets,
            weights, start = start, xtermsj, sterms = stermsj, control)

          ## Store stats.
          stats[[j]] <- modelstats(m)

          ## Save fitted values.
          if(stats[[j]]$GAIC < gaic) {
            gaic <- stats[[j]]$GAIC
            fit <- m$fitted.values
          }
        }

        ## Remove term with minimum GAIC?
        j <- sapply(stats, function(x) x$GAIC)
        j <- names(j)[which.min(j)]

        ## Check if GAIC improved?
        if(stats[[j]]$GAIC < stats_save[[iter]]$GAIC) {
          ## Update terms.
          if(linear(j)) {
            xterms <- dropterm(j, xterms)
          } else {
            sterms <- dropterm(j, sterms)
          }

          ## Remove from selected.
          selected <- selected[selected != j]

          ## Add to selected.
          notselected <- c(notselected, j)

          ## Increase iterator.
          iter <- iter + 1L

          ## Save stats.
          stats_save[[iter]] <- stats[[j]]
          stats_save[[iter]]$term <- j

          ## Set starting values.
          start <- fit

          ## Print info.
          if(trace[2L]) {
            j <- splitname(j)
            cat("  GAIC = ", fmt(stats_save[[iter]]$GAIC, width = 10, digits = 4),
            " <-> parameter ", j[1L], ", term ", j[2L],  "\n", sep = "")
          }
        } else {
          do <- FALSE
        }
      }
    }

    improved <- gaic < gaic0
    if(improved)
      gaic0 <- gaic

    k <- k + 1L
  }

  ## Final model.
  m <- RS(x, y, specials, family, offsets,
    weights, start = start, xterms, sterms = sterms, control)

  ## Small extractor function.
  ge <- function(j) { sapply(stats_save, function(z) z[[j]]) }

  m$selection <- list("GAIC" = ge("GAIC"), "logLik" = ge("logLik"),
    "df" = ge("df"), "K" = K)
  names(m$selection$GAIC) <- ge("term")
  m$xterms <- xterms
  m$sterms <- sterms
  m$specials <- specials[unique(unlist(sterms))]

  return(m)
}


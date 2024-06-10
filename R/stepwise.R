## Helper functions.
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
modelstats <- function(model, K) {
  stats <- list("df" = get_df2(model))
  stats$GAIC <- model$deviance + K * stats$df
  stats$logLik <- model$logLik
  return(stats)
}

forward_backward_step <- function(selected, notselected, xterms, sterms, strategy, stats_save,
  x, y, specials, family, offsets, weights, start, control, trace, K, maxit = Inf)
{
  ## Which model terms?
  modelterms <- switch(strategy,
    "forward" = notselected,
    "backward" = selected
  )

  ## Get current iterator.
  iter <- length(stats_save)

  ## Multiple core processing?
  cores <- control$cores

  ## Start select loop.
  do <- TRUE
  model <- NULL
  k <- 0L
  while(do && (k < maxit)) {
    if(length(modelterms) > 1L) {
      stats <- list()

      fitfun <- function(j) {
        ## Set model terms.
        xtermsj <- xterms
        stermsj <- sterms

        ## Forward or backward?
        if(strategy == "forward") {
          if(linear(j)) {
            xtermsj <- addterm(j, xterms)
          } else {
            stermsj <- addterm(j, sterms)
          }
        } else {
          if(linear(j)) {
            xtermsj <- dropterm(j, xterms)
          } else {
            stermsj <- dropterm(j, sterms)
          }
        }

        ## Estimate model.
        m <- RS(x, y, specials = specials, family, offsets,
          weights, start = start, xtermsj, sterms = stermsj, control)

        ## Return stats.
        return(modelstats(m, K))
      }

      if(cores < 2L) {
        stats <- lapply(modelterms, fitfun)
      } else {
        stats <- parallel::mclapply(modelterms, fitfun, mc.cores = cores)
      }

      ## Term with minimum GAIC?
      names(stats) <- modelterms
      j <- sapply(stats, function(x) x$GAIC)
      j <- names(j)[which.min(j)]

      ## Check if GAIC improved?
      if(stats[[j]]$GAIC < stats_save[[iter]]$GAIC) {
        ## Update terms.
        if(strategy == "forward") {
          if(linear(j)) {
            xterms <- addterm(j, xterms)
          } else {
            sterms <- addterm(j, sterms)
          }

          ## Remove from not selected.
          notselected <- notselected[notselected != j]

          ## Add to selected.
          selected <- c(selected, j)
        } else {
          if(linear(j)) {
            xterms <- dropterm(j, xterms)
          } else {
            sterms <- dropterm(j, sterms)
          }

          ## Remove from selected.
          selected <- selected[selected != j]

          ## Add to not selected.
          notselected <- c(notselected, j)
        }

        ## Increase iterator.
        iter <- iter + 1L

        ## Save stats.
        stats_save[[iter]] <- stats[[j]]
        stats_save[[iter]]$term <- j

        modelterms <- switch(strategy,
          "forward" = notselected,
          "backward" = selected
        )

        ## Print info.
        if(trace[2L]) {
          j <- splitname(j)
          cat("  GAIC = ", fmt(stats_save[[iter]]$GAIC, width = 10, digits = 4),
            if(strategy == "forward") " <+> parameter " else " <-> parameter ",
            j[1L], ", term ", j[2L],  "\n", sep = "")
        }
      } else {
        do <- FALSE
      }
    } else {
      do <- FALSE
    }

    k <- k + 1L
  }

  model <- RS(x, y, specials = specials, family, offsets,
    weights, start = start, xterms, sterms = sterms, control)

  rval <- list("selected" = selected, "notselected" = notselected,
    "stats_save" = stats_save, "xterms" = xterms, "sterms" = sterms,
    "model" = model)

  return(rval)
}

replace_step <- function(selected, notselected, xterms, sterms, stats_save,
  x, y, specials, family, offsets, weights, start, control, trace, K, maxit = Inf)
{
  ## Get current iterator.
  iter <- length(stats_save)

  ## Multiple core processing?
  cores <- control$cores

  ## Start select loop.
  do <- TRUE
  model <- NULL
  k <- 0L
  while(do & (k < maxit)) {
    if(length(selected) > 1L) {
      stats <- NULL

      ## Set model terms.
      xtermsj <- xterms
      stermsj <- sterms

      for(j in selected) {
        ## Drop model term.
        if(linear(j)) {
          xtermsj <- dropterm(j, xterms)
        } else {
          stermsj <- dropterm(j, sterms)
        }

        fitfun <- function(i) {
          ## Add model term.
          if(linear(i)) {
            xtermsj <- addterm(i, xterms)
          } else {
            stermsj <- addterm(i, sterms)
          }

          ## Estimate model.
          m <- RS(x, y, specials = specials, family, offsets,
            weights, start = start, xtermsj, sterms = stermsj, control)

          ## Store stats.
          ms <- modelstats(m, K)
          ms$drop <- j
          ms$add <- i

          return(ms)
        }

        if(cores < 2L) {
          si <- lapply(notselected, fitfun)
        } else {
          si <- parallel::mclapply(notselected, fitfun, mc.cores = cores)
        }

        stats <- c(stats, si)
      }

      ## Term with minimum GAIC?
      k <- which.min(sapply(stats, function(x) x$GAIC))

      ## Check if GAIC improved?
      if(stats[[k]]$GAIC < stats_save[[iter]]$GAIC) {
        ## Update terms.
        j <- stats[[k]]$drop
        i <- stats[[k]]$add

        if(linear(i)) {
          xterms <- addterm(i, xterms)
        } else {
          sterms <- addterm(i, sterms)
        }

        if(linear(j)) {
          xterms <- dropterm(j, xterms)
        } else {
          sterms <- dropterm(j, sterms)
        }

        ## Remove from selected.
        selected <- selected[selected != j]
        notselected <- c(notselected, j)

        ## Add to selected.
        selected <- c(selected, i)
        notselected <- notselected[notselected != i]

        ## Increase iterator.
        iter <- iter + 1L

        ## Save stats.
        stats_save[[iter]] <- stats[[k]]
        stats_save[[iter]]$term <- i
        stats_save[[iter]]$dropped <- j

        ## Print info.
        if(trace[2L]) {
          r <- splitname(j)
          a <- splitname(i)
          cat("  GAIC = ", fmt(stats_save[[iter]]$GAIC, width = 10, digits = 4),
              " <-> parameter ", r[1L], ", term ", r[2L],  "\n", sep = "")
          cat("                    <+> parameter ", a[1L], ", term ", a[2L],  "\n", sep = "")
        }
      } else {
        do <- FALSE
      }
    } else {
      do <- FALSE
    }
    k <- k + 1L
  }

  model <- RS(x, y, specials = specials, family, offsets,
    weights, start = start, xterms, sterms = sterms, control)

  rval <- list("selected" = selected, "notselected" = notselected,
    "stats_save" = stats_save, "xterms" = xterms, "sterms" = sterms,
    "model" = model)

  return(rval)
}

both_step <- function(selected, notselected, xterms, sterms, stats_save,
  x, y, specials, family, offsets, weights, start, control, trace, K)
{
  trace2 <- trace
  trace2[2L] <- FALSE

  gaic <- stats_save[[length(stats_save)]]$GAIC
  model <- NULL

  do <- TRUE
  while(do) {
    fs <- forward_backward_step(selected, notselected, xterms, sterms,
      strategy = "forward", stats_save, x, y, specials, family, offsets, weights, start,
      control, trace2, K, maxit = 1L)

    bs <- forward_backward_step(selected, notselected, xterms, sterms,
      strategy = "backward", stats_save, x, y, specials, family, offsets, weights, start,
      control, trace2, K, maxit = 1L)

    i1 <- length(fs$stats_save)
    i2 <- length(bs$stats_save)

    if((fs$stats_save[[i1]]$GAIC < gaic) || (bs$stats_save[[i2]]$GAIC < gaic)) {
      if((fs$stats_save[[i1]]$GAIC < bs$stats_save[[i2]]$GAIC) && (fs$stats_save[[i1]]$GAIC < gaic)) {
        ## Print info.
        if(trace[2L]) {
          j <- splitname(fs$stats_save[[i1]]$term)
          cat("  GAIC = ", fmt(fs$stats_save[[i1]]$GAIC, width = 10, digits = 4),
            " <+> parameter ", j[1L], ", term ", j[2L],  "\n", sep = "")
        }
        selected <- fs$selected
        notselected <- fs$notselected
        xterms <- fs$xterms
        sterms <- fs$sterms
        stats_save <- fs$stats_save
        model <- fs$model
        gaic <- fs$stats_save[[i1]]$GAIC
      }

      if((bs$stats_save[[i2]]$GAIC < fs$stats_save[[i1]]$GAIC) && (bs$stats_save[[i2]]$GAIC < gaic)) {
        ## Print info.
        if(trace[2L]) {
          j <- splitname(bs$stats_save[[i2]]$term)
          cat("  GAIC = ", fmt(bs$stats_save[[i2]]$GAIC, width = 10, digits = 4),
            " <-> parameter ", j[1L], ", term ", j[2L],  "\n", sep = "")
        }
        selected <- bs$selected
        notselected <- bs$notselected
        xterms <- bs$xterms
        sterms <- bs$sterms
        stats_save <- bs$stats_save
        model <- bs$model
        gaic <- bs$stats_save[[i2]]$GAIC
      }
    } else {
      do <- FALSE
    }
  }

  rval <- list("selected" = selected, "notselected" = notselected,
    "stats_save" = stats_save, "xterms" = xterms, "sterms" = sterms,
    "model" = model)

  return(rval)
}

all_step <- function(selected, notselected, xterms, sterms, stats_save,
  x, y, specials, family, offsets, weights, start, control, trace, K)
{
  trace2 <- trace
  trace2[2L] <- FALSE

  gaic <- stats_save[[length(stats_save)]]$GAIC
  model <- NULL

  do <- TRUE
  while(do) {
    fs <- forward_backward_step(selected, notselected, xterms, sterms,
      strategy = "forward", stats_save, x, y, specials, family, offsets, weights, start,
      control, trace2, K, maxit = 1L)

    bs <- forward_backward_step(selected, notselected, xterms, sterms,
      strategy = "backward", stats_save, x, y, specials, family, offsets, weights, start,
      control, trace2, K, maxit = 1L)

    rs <- replace_step(selected, notselected, xterms, sterms,
      stats_save, x, y, specials = specials, family, offsets, weights, start,
      control, trace, K, maxit = 1L)

    i1 <- length(fs$stats_save)
    i2 <- length(bs$stats_save)
    i3 <- length(rs$stats_save)

    fs_step <- (fs$stats_save[[i1]]$GAIC < bs$stats_save[[i2]]$GAIC) &&
      (fs$stats_save[[i1]]$GAIC < rs$stats_save[[i3]]$GAIC) &&
      (fs$stats_save[[i1]]$GAIC < gaic)

    bs_step <- (bs$stats_save[[i2]]$GAIC < fs$stats_save[[i1]]$GAIC) &&
      (bs$stats_save[[i2]]$GAIC < rs$stats_save[[i3]]$GAIC) &&
      (bs$stats_save[[i2]]$GAIC < gaic)

    rs_step <- (rs$stats_save[[i3]]$GAIC < fs$stats_save[[i1]]$GAIC) &&
      (rs$stats_save[[i3]]$GAIC < bs$stats_save[[i2]]$GAIC) &&
      (rs$stats_save[[i3]]$GAIC < gaic)

    if(any(c(fs_step, bs_step, rs_step))) {
      if(fs_step) {
        ## Print info.
        if(trace[2L]) {
          j <- splitname(fs$stats_save[[i1]]$term)
          cat("  GAIC = ", fmt(fs$stats_save[[i1]]$GAIC, width = 10, digits = 4),
            " <+> parameter ", j[1L], ", term ", j[2L],  "\n", sep = "")
        }
        selected <- fs$selected
        notselected <- fs$notselected
        xterms <- fs$xterms
        sterms <- fs$sterms
        stats_save <- fs$stats_save
        model <- fs$model
        gaic <- fs$stats_save[[i1]]$GAIC
      }

      if(bs_step) {
        ## Print info.
        if(trace[2L]) {
          j <- splitname(bs$stats_save[[i2]]$term)
          cat("  GAIC = ", fmt(bs$stats_save[[i2]]$GAIC, width = 10, digits = 4),
            " <-> parameter ", j[1L], ", term ", j[2L],  "\n", sep = "")
        }
        selected <- bs$selected
        notselected <- bs$notselected
        xterms <- bs$xterms
        sterms <- bs$sterms
        stats_save <- bs$stats_save
        model <- bs$model
        gaic <- bs$stats_save[[i2]]$GAIC
      }

      if(rs_step) {
        ## Print info.
        if(trace[2L]) {
          r <- splitname(rs$stats_save[[i3]]$dropped)
          a <- splitname(rs$stats_save[[i3]]$term)
          cat("  GAIC = ", fmt(rs$stats_save[[i3]]$GAIC, width = 10, digits = 4),
              " <-> parameter ", r[1L], ", term ", r[2L],  "\n", sep = "")
          cat("                    <+> parameter ", a[1L], ", term ", a[2L],  "\n", sep = "")
        }
        selected <- rs$selected
        notselected <- rs$notselected
        xterms <- rs$xterms
        sterms <- rs$sterms
        stats_save <- rs$stats_save
        model <- rs$model
        gaic <- rs$stats_save[[i3]]$GAIC
      }
    } else {
      do <- FALSE
    }
  }

  rval <- list("selected" = selected, "notselected" = notselected,
    "stats_save" = stats_save, "xterms" = xterms, "sterms" = sterms,
    "model" = model)

  return(rval)
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
  choices <- c("forward.linear", "forward", "backward", "backward.linear",
     "replace", "replace.linear", "both", "both.linear", "all", "all.linear")
  if(is.null(strategy)) {
    strategy <- c("both.linear", "both")
  } else {
    strategy <- match.arg(strategy, choices, several.ok = TRUE)
  }

  ## Should the order of the parameters be retained?
  keeporder <- isTRUE(control$keeporder)

  ## Inner and outer trace.
  trace <- control$trace
  if(length(trace) < 2L) {
    trace <- c(FALSE, trace)
  }
  control$trace <- trace[1L]

  ## Penalty for AIC.
  K <- if(is.null(control$K)) 2 else control$K

  ## Parallel loops?
  if(is.null(control$cores)) {
    control$cores <- 1L
  }

  ## Estimate nullmodel first.
  xterms_itcpt <- list()
  for(i in nx)
    xterms_itcpt[[i]] <- "(Intercept)"

  m <- RS(x, y, specials = NULL, family, offsets,
    weights, start = start, xterms_itcpt, sterms = NULL, control)

  ## Save model stats and number of iterations.
  stats_save <- list()
  stats_save[[1L]] <- modelstats(m, K)

  ## Save initial GAIC.
  gaic <- stats_save[[1L]]$GAIC

  ## Setup vectors for selected and not selected model terms.
  ## Start with linear terms first.
  selected <- NULL
  notselected <- NULL
  for(i in nx) {
    if(length(xterms[[i]]))
      notselected <- c(notselected, paste0(i, ".p.", xterms[[i]]))
  }
  notselected <- notselected[-grep(".p.(Intercept)", notselected, fixed = TRUE)]

  ## Set starting xterms.
  xterms <- xterms_itcpt

  ## Start linear only stepwise algorithm.
  do <- TRUE
  iter <- 1L
  while(do) {
    if("forward.linear" %in% strategy) {
      if(trace[2L] && (iter < 2L))
        cat("Forward Linear Selection\n")
      if(keeporder) {
        for(i in nx) {
          ## Get only linear terms
          isel <- grep(paste0(i, ".p."), selected, fixed = TRUE, value = TRUE)
          inot <- grep(paste0(i, ".p."), notselected, fixed = TRUE, value = TRUE)

          ## Run forward selection.
          si <- forward_backward_step(isel, inot, xterms, sterms = NULL, strategy = "forward",
            stats_save, x, y, specials = NULL, family, offsets, weights, start, control, trace, K)

          ## Update.
          if(length(si$model)) {
            selected <- c(selected[!grepl(paste0(i, ".p."), selected, fixed = TRUE)], si$selected)
            notselected <- c(notselected[!grepl(paste0(i, ".p."), notselected, fixed = TRUE)], si$notselected)
            xterms <- si$xterms
            stats_save <- si$stats_save
            model <- si$model
          }
        }
      } else {
        ## Run forward selection.
        si <- forward_backward_step(selected, notselected, xterms, sterms = NULL, strategy = "forward",
          stats_save, x, y, specials = NULL, family, offsets, weights, start, control, trace, K)

        ## Update.
        if(length(si$model)) {
          selected <- si$selected
          notselected <- si$notselected
          xterms <- si$xterms
          stats_save <- si$stats_save
          model <- si$model
        }
      }
    }

    if("both.linear" %in% strategy) {
      if(trace[2L] && (iter < 2L))
        cat("Bidirectional Linear Selection\n")
      if(keeporder) {
        for(i in nx) {
          ## Get only linear terms
          isel <- grep(paste0(i, ".p."), selected, fixed = TRUE, value = TRUE)
          inot <- grep(paste0(i, ".p."), notselected, fixed = TRUE, value = TRUE)

          ## Run bidirectional selection.
          si <- both_step(isel, inot, xterms, sterms = NULL,
            stats_save, x, y, specials = NULL, family, offsets, weights, start, control, trace, K)

          ## Update.
          if(length(si$model)) {
            selected <- c(selected[!grepl(paste0(i, ".p."), selected, fixed = TRUE)], si$selected)
            notselected <- c(notselected[!grepl(paste0(i, ".p."), notselected, fixed = TRUE)], si$notselected)
            xterms <- si$xterms
            stats_save <- si$stats_save
            model <- si$model
          }
        }
      } else {
        ## Bidirectional selection.
        si <- both_step(selected, notselected, xterms, sterms = NULL,
          stats_save, x, y, specials = NULL, family, offsets, weights, start, control, trace, K)

        ## Update.
        if(length(si$model)) {
          selected <- si$selected
          notselected <- si$notselected
          xterms <- si$xterms
          stats_save <- si$stats_save
          model <- si$model
        }
      }
    }

    if("all.linear" %in% strategy) {
      if(trace[2L] && (iter < 2L))
        cat("All Directions Linear Selection\n")
      if(keeporder) {
        for(i in nx) {
          ## Get only linear terms
          isel <- grep(paste0(i, ".p."), selected, fixed = TRUE, value = TRUE)
          inot <- grep(paste0(i, ".p."), notselected, fixed = TRUE, value = TRUE)

          ## Run bidirectional selection.
          si <- all_step(isel, inot, xterms, sterms = NULL,
            stats_save, x, y, specials = NULL, family, offsets, weights, start, control, trace, K)

          ## Update.
          if(length(si$model)) {
            selected <- c(selected[!grepl(paste0(i, ".p."), selected, fixed = TRUE)], si$selected)
            notselected <- c(notselected[!grepl(paste0(i, ".p."), notselected, fixed = TRUE)], si$notselected)
            xterms <- si$xterms
            stats_save <- si$stats_save
            model <- si$model
          }
        }
      } else {
        ## All directions selection.
        si <- all_step(selected, notselected, xterms, sterms = NULL,
          stats_save, x, y, specials = NULL, family, offsets, weights, start, control, trace, K)

        ## Update.
        if(length(si$model)) {
          selected <- si$selected
          notselected <- si$notselected
          xterms <- si$xterms
          stats_save <- si$stats_save
          model <- si$model
        }
      }
    }

    if("replace.linear" %in% strategy) {
      if(trace[2L] && (iter < 2L))
        cat("Replace Linear Selection\n")
      if(keeporder) {
        for(i in nx) {
          ## Get only linear terms
          isel <- grep(paste0(i, ".p."), selected, fixed = TRUE, value = TRUE)
          inot <- grep(paste0(i, ".p."), notselected, fixed = TRUE, value = TRUE)

          ## Run replace selection.
          si <- replace_step(isel, inot, xterms, sterms = NULL,
            stats_save, x, y, specials = NULL, family, offsets, weights, start, control, trace, K)

          ## Update.
          if(length(si$model)) {
            selected <- c(selected[!grepl(paste0(i, ".p."), selected, fixed = TRUE)], si$selected)
            notselected <- c(notselected[!grepl(paste0(i, ".p."), notselected, fixed = TRUE)], si$notselected)
            xterms <- si$xterms
            stats_save <- si$stats_save
            model <- si$model
          }
        }
      } else {
        ## Run replace selection.
        si <- replace_step(selected, notselected, xterms, sterms = NULL,
          stats_save, x, y, specials = NULL, family, offsets, weights, start, control, trace, K)

        ## Update.
        if(length(si$model)) {
          selected <- si$selected
          notselected <- si$notselected
          xterms <- si$xterms
          stats_save <- si$stats_save
          model <- si$model
        }
      }
    }

    if("backward.linear" %in% strategy) {
      if(trace[2L] && (iter < 2L))
        cat("Backward Linear Selection\n")
      if(keeporder) {
        for(i in rev(nx)) {
          ## Get only linear terms
          isel <- grep(paste0(i, ".p."), selected, fixed = TRUE, value = TRUE)
          inot <- grep(paste0(i, ".p."), notselected, fixed = TRUE, value = TRUE)

          ## Run backward elimination.
          si <- forward_backward_step(isel, inot, xterms, sterms = NULL, strategy = "backward",
            stats_save, x, y, specials = NULL, family, offsets, weights, start, control, trace, K)

          ## Update.
          if(length(si$model)) {
            selected <- c(selected[!grepl(paste0(i, ".p."), selected, fixed = TRUE)], si$selected)
            notselected <- c(notselected[!grepl(paste0(i, ".p."), notselected, fixed = TRUE)], si$notselected)
            xterms <- si$xterms
            stats_save <- si$stats_save
            model <- si$model
          }
        }
      } else {
        ## Run backward elimination.
        si <- forward_backward_step(selected, notselected, xterms, sterms = NULL, strategy = "backward",
          stats_save, x, y, specials = NULL, family, offsets, weights, start, control, trace, K)

        ## Update.
        if(length(si$model)) {
          selected <- si$selected
          notselected <- si$notselected
          xterms <- si$xterms
          stats_save <- si$stats_save
          model <- si$model
        }
      }
    }

    do <- stats_save[[length(stats_save)]]$GAIC < gaic
    if(do) {
      gaic <- stats_save[[length(stats_save)]]$GAIC
      iter <- iter + 1L
    }
  }

  ## Using all terms now.
  for(i in nx) {
    if(length(sterms[[i]]))
      notselected <- c(notselected, paste0(i, ".s.", sterms[[i]]))
  }

  ## sterms start list().
  sterms <- list()
  for(i in nx)
    sterms[[i]] <- character(0)

  do <- TRUE
  iter <- 1L
  while(do) {
    if(trace[2L] && (iter > 1L))
      cat("Continue\n")

    if("forward" %in% strategy) {
      if(trace[2L] && (iter < 2L))
        cat("Forward Selection\n")
      if(keeporder) {
        for(i in nx) {
          ## Get selected terms.
          isel <- c(
            grep(paste0(i, ".p."), selected, fixed = TRUE, value = TRUE),
            grep(paste0(i, ".s."), selected, fixed = TRUE, value = TRUE)
          )

          ## Get not selected terms.
          inot <- c(
            grep(paste0(i, ".p."), notselected, fixed = TRUE, value = TRUE),
            grep(paste0(i, ".s."), notselected, fixed = TRUE, value = TRUE)
          )

          ## Run forward selection.
          si <- forward_backward_step(isel, inot, xterms, sterms, strategy = "forward",
            stats_save, x, y, specials, family, offsets, weights, start, control, trace, K)

          ## Update.
          if(length(si$model)) {
            jd <- !grepl(paste0(i, ".p."), selected, fixed = TRUE) & !grepl(paste0(i, ".s."), selected, fixed = TRUE)
            selected <- c(selected[jd], si$selected)
            jd <- !grepl(paste0(i, ".p."), notselected, fixed = TRUE) & !grepl(paste0(i, ".s."), notselected, fixed = TRUE)
            notselected <- c(notselected[jd], si$notselected)
            selected <- unique(c(selected, si$selected))
            notselected <- unique(c(notselected, si$notselected))
            xterms <- si$xterms
            sterms <- si$sterms
            stats_save <- si$stats_save
            model <- si$model
          }
        }
      } else {
        ## Run forward selection.
        si <- forward_backward_step(selected, notselected, xterms, sterms, strategy = "forward",
          stats_save, x, y, specials, family, offsets, weights, start, control, trace, K)

        ## Update.
        if(length(si$model)) {
          selected <- si$selected
          notselected <- si$notselected
          xterms <- si$xterms
          sterms <- si$sterms
          stats_save <- si$stats_save
          model <- si$model
        }
      }
    }

    if("both" %in% strategy) {
      if(trace[2L] && (iter < 2L))
        cat("Bidirectional Selection\n")
      if(keeporder) {
        for(i in nx) {
          ## Get selected terms.
          isel <- c(
            grep(paste0(i, ".p."), selected, fixed = TRUE, value = TRUE),
            grep(paste0(i, ".s."), selected, fixed = TRUE, value = TRUE)
          )

          ## Get not selected terms.
          inot <- c(
            grep(paste0(i, ".p."), notselected, fixed = TRUE, value = TRUE),
            grep(paste0(i, ".s."), notselected, fixed = TRUE, value = TRUE)
          )

          ## Run bidirectional selection.
          si <- both_step(isel, inot, xterms, sterms,
            stats_save, x, y, specials, family, offsets, weights, start, control, trace, K)

          ## Update.
          if(length(si$model)) {
            jd <- !grepl(paste0(i, ".p."), selected, fixed = TRUE) & !grepl(paste0(i, ".s."), selected, fixed = TRUE)
            selected <- c(selected[jd], si$selected)
            jd <- !grepl(paste0(i, ".p."), notselected, fixed = TRUE) & !grepl(paste0(i, ".s."), notselected, fixed = TRUE)
            notselected <- c(notselected[jd], si$notselected)
            xterms <- si$xterms
            sterms <- si$sterms
            stats_save <- si$stats_save
            model <- si$model
          }
        }
      } else {
        ## Bidirectional selection.
        si <- both_step(selected, notselected, xterms, sterms,
          stats_save, x, y, specials, family, offsets, weights, start, control, trace, K)

        ## Update.
        if(length(si$model)) {
          selected <- si$selected
          notselected <- si$notselected
          xterms <- si$xterms
          sterms <- si$sterms
          stats_save <- si$stats_save
          model <- si$model
        }
      }
    }

    if("all" %in% strategy) {
      if(trace[2L] && (iter < 2L))
        cat("All Directions Selection\n")
      if(keeporder) {
        for(i in nx) {
          ## Get selected terms.
          isel <- c(
            grep(paste0(i, ".p."), selected, fixed = TRUE, value = TRUE),
            grep(paste0(i, ".s."), selected, fixed = TRUE, value = TRUE)
          )

          ## Get not selected terms.
          inot <- c(
            grep(paste0(i, ".p."), notselected, fixed = TRUE, value = TRUE),
            grep(paste0(i, ".s."), notselected, fixed = TRUE, value = TRUE)
          )

          ## Run bidirectional selection.
          si <- all_step(isel, inot, xterms, sterms,
            stats_save, x, y, specials, family, offsets, weights, start, control, trace, K)

          ## Update.
          if(length(si$model)) {
            jd <- !grepl(paste0(i, ".p."), selected, fixed = TRUE) & !grepl(paste0(i, ".s."), selected, fixed = TRUE)
            selected <- c(selected[jd], si$selected)
            jd <- !grepl(paste0(i, ".p."), notselected, fixed = TRUE) & !grepl(paste0(i, ".s."), notselected, fixed = TRUE)
            notselected <- c(notselected[jd], si$notselected)
            xterms <- si$xterms
            sterms <- si$sterms
            stats_save <- si$stats_save
            model <- si$model
          }
        }
      } else {
        ## Bidirectional selection.
        si <- all_step(selected, notselected, xterms, sterms,
          stats_save, x, y, specials, family, offsets, weights, start, control, trace, K)

        ## Update.
        if(length(si$model)) {
          selected <- si$selected
          notselected <- si$notselected
          xterms <- si$xterms
          sterms <- si$sterms
          stats_save <- si$stats_save
          model <- si$model
        }
      }
    }

    if("replace" %in% strategy) {
      if(trace[2L] && (iter < 2L))
        cat("Replace Selection\n")
      if(keeporder) {
        for(i in nx) {
          ## Get selected terms.
          isel <- c(
            grep(paste0(i, ".p."), selected, fixed = TRUE, value = TRUE),
            grep(paste0(i, ".s."), selected, fixed = TRUE, value = TRUE)
          )

          ## Get not selected terms.
          inot <- c(
            grep(paste0(i, ".p."), notselected, fixed = TRUE, value = TRUE),
            grep(paste0(i, ".s."), notselected, fixed = TRUE, value = TRUE)
          )

          ## Run replace selection.
          si <- replace_step(isel, inot, xterms, sterms,
            stats_save, x, y, specials, family, offsets, weights, start, control, trace, K)

          ## Update.
          if(length(si$model)) {
            jd <- !grepl(paste0(i, ".p."), selected, fixed = TRUE) & !grepl(paste0(i, ".s."), selected, fixed = TRUE)
            selected <- c(selected[jd], si$selected)
            jd <- !grepl(paste0(i, ".p."), notselected, fixed = TRUE) & !grepl(paste0(i, ".s."), notselected, fixed = TRUE)
            notselected <- c(notselected[jd], si$notselected)
            xterms <- si$xterms
            sterms <- si$sterms
            stats_save <- si$stats_save
            model <- si$model
          }
        }
      } else {
        ## Run replace selection.
        si <- replace_step(selected, notselected, xterms, sterms,
          stats_save, x, y, specials, family, offsets, weights, start, control, trace, K)

        ## Update.
        if(length(si$model)) {
          selected <- si$selected
          notselected <- si$notselected
          xterms <- si$xterms
          sterms <- si$sterms
          stats_save <- si$stats_save
          model <- si$model
        }
      }
    }

    if("backward" %in% strategy) {
      if(trace[2L] && (iter < 2L))
        cat("Backward Selection\n")
      if(keeporder) {
        for(i in rev(nx)) {
          ## Get selected terms.
          isel <- c(
            grep(paste0(i, ".p."), selected, fixed = TRUE, value = TRUE),
            grep(paste0(i, ".s."), selected, fixed = TRUE, value = TRUE)
          )

          ## Get not selected terms.
          inot <- c(
            grep(paste0(i, ".p."), notselected, fixed = TRUE, value = TRUE),
            grep(paste0(i, ".s."), notselected, fixed = TRUE, value = TRUE)
          )

          ## Run backward elimination.
          si <- forward_backward_step(isel, inot, xterms, sterms, strategy = "backward",
            stats_save, x, y, specials, family, offsets, weights, start, control, trace, K)

          ## Update.
          if(length(si$model)) {
            jd <- !grepl(paste0(i, ".p."), selected, fixed = TRUE) & !grepl(paste0(i, ".s."), selected, fixed = TRUE)
            selected <- c(selected[jd], si$selected)
            jd <- !grepl(paste0(i, ".p."), notselected, fixed = TRUE) & !grepl(paste0(i, ".s."), notselected, fixed = TRUE)
            notselected <- c(notselected[jd], si$notselected)
            xterms <- si$xterms
            sterms <- si$sterms
            stats_save <- si$stats_save
            model <- si$model
          }
        }
      } else {
        ## Run backward elimination.
        si <- forward_backward_step(selected, notselected, xterms, sterms, strategy = "backward",
          stats_save, x, y, specials, family, offsets, weights, start, control, trace, K)

        ## Update.
        if(length(si$model)) {
          selected <- si$selected
          notselected <- si$notselected
          xterms <- si$xterms
          sterms <- si$sterms
          stats_save <- si$stats_save
          model <- si$model
        }
      }
    }

    do <- stats_save[[length(stats_save)]]$GAIC < gaic
    if(do) {
      gaic <- stats_save[[length(stats_save)]]$GAIC
      iter <- iter + 1L
    }
  }

  ## Small extractor function.
  ge <- function(j) { sapply(stats_save, function(z) z[[j]]) }

  model$selection <- list("GAIC" = ge("GAIC"), "logLik" = ge("logLik"),
    "df" = ge("df"), "K" = K, "formula" = xs2formula(xterms, sterms))
  names(model$selection$GAIC) <- gsub(".s.", ".", gsub(".p.", ".", ge("term"), fixed = TRUE), fixed = TRUE)
  model$xterms <- xterms
  model$sterms <- sterms
  model$specials <- specials[unique(unlist(sterms))]

  if(trace[2L]) {
    cat("Final Model\n")
    for(j in names(model$selection$formula)) {
      cat(paste0("$", j, "\n.. "))
      print(model$selection$formula[[j]])
    }
  }

  return(model)
}

xs2formula <- function(x, s)
{
  f <- list()
  for(j in unique(names(x), names(s))) {
    f[[j]] <- paste(c(unique(x[[j]]), unique(s[[j]])), collapse = "+")
    f[[j]] <- gsub("(Intercept)", "1", f[[j]], fixed = TRUE)
    f[[j]] <- as.formula(paste("~", f[[j]]), env = .GlobalEnv)
  }
  return(f)
}

stepGAMLSS <- function(formula, ..., K = 2,
  strategy = c("both.linear", "both"), keeporder = FALSE,
  cores = 1L)
{
  gamlss2(formula, ..., optimizer = stepwise, K = K, strategy = strategy,
    keeporder = keeporder, cores = cores)
}


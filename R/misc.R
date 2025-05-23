varimp <- function(object, newdata = NULL, scale = TRUE, nrep = 20, ...)
{
  vn <- unique(all.vars(formula(object$fake_formula, lhs = 0)))
  if(is.null(newdata))
    newdata <- model.frame(object)
  ll0 <- logLik(object, newdata = newdata)
  ll <- rep(NA, length(vn))
  names(ll) <- vn
  for(j in seq_along(vn)) {
    xj <- newdata[[vn[j]]]
    vi <- rep(NA, nrep)
    for(i in 1:nrep) {
      newdata[[vn[j]]] <- sample(xj)
      vi[i] <- logLik(object, newdata = newdata) - ll0
    }
    ll[vn[j]] <- mean(vi, na.rm = TRUE)
    newdata[[vn[j]]] <- xj
  }
  ll <- ll / ll0
  if(scale) {
    ll <- (ll - min(ll)) / diff(range(ll))
  }
  return(sort(ll))
}

available_families <- function(type = c("continuous", "discrete"), families = NULL, ...)
{
  stopifnot(requireNamespace("gamlss.dist"))
  if(!("package:gamlss.dist" %in% search()))
    attachNamespace("gamlss.dist")
  funs <- if(is.null(families)) {
    ls("package:gamlss.dist")
  } else {
    as.character(families)
  }
  type <- match.arg(type)
  d <- list()
  tf <- tempfile()
  tf2 <- tempfile()
  warn <- getOption("warn")
  options("warn" = -1)
  for(j in seq_along(funs)) {
    fj0 <- get(funs[j])    
    png(tf)
    capture.output(fj <- try(fj0(), silent = TRUE), file = tf2)
    dev.off()
    if(!inherits(fj, "try-error")) {
      if(inherits(fj, "gamlss.family")) {
        if(tolower(fj$type) == tolower(type)) {
          args <- list(...)
          i <- names(args) %in% names(formals(fj0))
          if(any(i)) {
            args  <- args[i]
            d[[funs[j]]] <- do.call("fj0", args)
          } else {
            d[[funs[j]]] <- fj0()
          }
        }
      }
    }
  }
  unlink(tf)
  unlink(tf2)
  options("warn" = warn)
  return(d)
}

find_family <- function(y, families = NULL, k = 2, verbose = TRUE, ...) {
  if(is.null(families)) {
    families <- if(is.numeric(y)) {
      available_families(type = "continuous")
    } else {
      available_families(type = "discrete")
    }
  }
  if(is.character(families)) {
    families <- available_families(families = families)
  }

  if(is.null(names(families))) {
    nf <- sapply(families, function(x) {
      if(is.function(x))
        x <- x()
      x$family[1L]
    })
    names(families) <- nf
  }

  ic <- rep(NA, length(families))
  names(ic) <- names(families)

  warn <- options("warn")$warn
  options("warn" = -1)
  on.exit(options("warn" = warn))

  for(j in names(families)) {
    if(verbose) {
      cat("..", j, "family\n")
    }

    warning_occurred <- FALSE

    b <- try({
      withCallingHandlers(
        gamlss2(y ~ 1, family = families[[j]], trace = FALSE, ...),
        warning = function(w) {
          warning_occurred <<- TRUE
          invokeRestart("muffleWarning")
        }
      )
    }, silent = TRUE)

    if(!inherits(b, "try-error") && !warning_occurred) {
      ic[j] <- GAIC(b, k = k)
      cat(".. .. IC =", round(ic[j], 4), "\n")
    } else {
      cat(".. .. error\n")
    }
  }

  return(sort(ic, decreasing = TRUE))
}

fit_family <- function(y, family = NO, plot = TRUE, ...)
{
  if(is.character(family))
    family <- available_families(families = family[1L])[[1L]]

  y <- na.omit(y)

  b <- gamlss2(y ~ 1, family = family, ...)

  if(plot) {
    par <- predict(b, type = "parameter", drop = FALSE)

    main <- list(...)$main
    xlab <- list(...)$xlab
    if(is.null(xlab))
      xlab <- "Response"
    ylim <- list(...)$ylim
    xlim <- list(...)$xlim

    if(is.null(b$family$type))
      b$family$type <- "continuous"

    if((b$family$type == "continuous") || (b$family$type == "mixed")) {
      if(is.null(main)) {
        k <- list(...)$k
        if(is.null(k))
          k <- 2
        main <- paste("Histogram and estimated",  b$family$family,
          "density\nGAIC =", round(GAIC(b, k = k), 4))
      }

      dy <- b$family$pdf(y, par)

      h <- hist(y, breaks = "Scott", plot = FALSE)

      if(is.null(ylim))
        ylim <- range(c(h$density, dy))
      if(is.null(xlim))
        xlim <- range(y, na.rm = TRUE)

      hist(y, breaks = "Scott", freq = FALSE,
        ylim = ylim, xlim = xlim,
        xlab = xlab, main = main)
      i <- order(y)
      lines(dy[i] ~ y[i], col = 4, lwd = 2)
    }

    if(b$family$type == "discrete") {
      if(is.null(main)) {
        k <- list(...)$k
        if(is.null(k))
          k <- 2
        main <- paste("Proportions and estimated",  b$family$family,
          "probabilties\nGAIC =", round(GAIC(b, k = k), 4))
      }

      dy <- b$family$pdf(y, par)

      Y <- unique(cbind(y, dy))
      Y <- Y[order(Y[, 1]), ]

      tab <- prop.table(table(y))

      if(is.null(ylim))
        ylim <- range(c(0, tab, dy * 1.1))

      ylab <- list(...)$ylab
      if(is.null(ylab))
        ylab <- "Probability"

      bp <- barplot(tab, ylim = ylim, xlim = xlim,
        xlab = xlab, main = main, ylab = ylab)

      bp <- as.numeric(bp)
      #xr <- min(diff(bp))
      #bp <- bp - 0.5 * xr

      lines(Y[, 2] ~ bp, col = 4, lwd = 2, type = "h")
      points(bp, Y[, 2], col = 4, pch = 16)
      points(bp, Y[, 2], col = rgb(0.1, 0.1, 0.1, alpha = 0.6))
    }

    legend <- list(...)$legend
    if(is.null(legend))
      legend <- TRUE

    if(legend) {
      np <- names(par)
      par <- par[1L, ]

      pos <- list(...)$pos
      if(is.null(pos))
        pos <- "topright"

      legend(pos,
        paste0(np, " = ", round(par, 4)),
        title = "Parameters:", title.font = 2,
        lwd = 1, col = NA, bty = "n")
    }

    return(invisible(b))
  } else {
    return(b)
  }  
}

find_gamlss2 <- function(formula, families = NULL, k = 2,
  select = FALSE, verbose = TRUE, ...)
{
  if(is.null(families))
    families <- available_families(...)

  if(is.character(families)) {
    families <- available_families(families = families, ...)
  }

  if(is.null(names(families))) {
    nf <- sapply(families, function(x) {
      if(is.function(x))
        x <- x()
      x$family[1L]
    })
    names(families) <- nf
  }

  engine <- if(select) select_gamlss2 else gamlss2

  ic <- rep(NA, length(families))
  names(ic) <- names(families)

  m <- NULL

  warn <- options("warn")$warn
  options("warn" = -1)
  on.exit(options("warn" = warn))

  for(j in names(families)) {
    if(verbose) {
      cat("..", j, "family\n")
    }

    warning_occurred <- FALSE

    b <- try({
      withCallingHandlers(
        engine(formula, family = families[[j]], trace = FALSE, ...),
        warning = function(w) {
          warning_occurred <<- TRUE
          invokeRestart("muffleWarning")
        }
      )
    }, silent = TRUE)

    if(!inherits(b, "try-error") && !warning_occurred) {
      ic[j] <- GAIC(b, k = k)
      cat(".. .. IC =", round(ic[j], 4), "\n")
      if(ic[j] == min(ic, na.rm = TRUE)) {
        m <- b
      }
    } else {
      cat(".. .. error\n")
    }
  }

  if(!is.null(m))
    m$ic <- sort(ic, decreasing = TRUE)

  return(m)
}


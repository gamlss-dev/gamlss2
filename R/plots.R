## A plotting method.
plot.gamlss2 <- function(x, parameter = NULL,
  which = "effects", terms = NULL,
  scale = TRUE, spar = TRUE, ...)
{
  ## What should be plotted?
  which.match <- c("effects", "hist-resid", "qq-resid", "wp-resid", "scatter-resid", "selection")
  if(!is.character(which)) {
    if(any(which > 5L))
      which <- which[which <= 5L]
    which <- which.match[which]
  } else which <- which.match[grep(tolower(which), which.match, fixed = TRUE)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")

  ## If effects, which parameters to plot?
  if(is.null(parameter)) {
    parameter <- list(...)$what
    if(is.null(parameter))
    parameter <- list(...)$model
    if(is.null(parameter))
      parameter <- x$family$names
  }
  if(!is.character(parameter))
    parameter <- x$family$names[parameter]
  parameter <- x$family$names[pmatch(parameter, x$family$names)]

  parameter <- parameter[!is.na(parameter)]
  if(length(parameter) < 1L)
    stop("argument parameter is specified wrong!")

  ## Check for any effect plots.
  if(length(which) < 2L) {
    if(which == "effects" & length(x$results$effects) < 1L)
      which  <- c("hist-resid", "qq-resid", "wp-resid", "scatter-resid")
  }


  ask <- list(...)$ask
  if(is.null(ask))
    ask <- FALSE
  spare <- spar

  ## Effect plots.
  if("effects" %in% which) {
    if(is.null(x$results))
      x$results <- results(x)

    en <- grep2(parameter, names(x$results$effects), fixed = TRUE, value = TRUE)

    if(length(en) < 2L)
      spare <- FALSE

    if(spare) {
      owd <- par(no.readonly = TRUE)
      on.exit(par(owd))
    }

    if(!is.null(terms)) {
      if(is.character(terms)) {
        en <- grep2(terms, en, fixed = TRUE, value = TRUE)
      } else {
        if(max(terms) > length(en))
          stop("argument terms is specified wrong!")
        en <- en[terms]
      }
    }

    if(length(x$results$effects)) {
      if(spare) {
        if(isTRUE(ask)) {
          par(ask = TRUE)
        } else {
          nplts <- length(en)
          par(mfrow = if(nplts <= 4) c(1, nplts) else n2mfrow(nplts))
        }
      }

      ylim <- list(...)$ylim

      if(scale > 0) {
        ylim <- list()
        for(i in parameter) {
          gok <- grep(i, en, fixed = TRUE, value = TRUE)
          for(j in gok) {
            if("lower" %in% colnames(x$results$effects[[j]])) {
              ylim[[i]] <- c(ylim[[i]], range(x$results$effects[[j]][, c("lower", "upper")]))
            } else {
              ylim[[i]] <- c(ylim[[i]], range(x$results$effects[[j]][, "fit"]))
            }
          }
          if(length(gok))
            ylim[[i]] <- range(ylim[[i]])
        }
      }

      for(j in en) {
        p <- strsplit(j, ".", fixed = TRUE)[[1L]][1L]
        if(!is.factor(x$results$effects[[j]][[1L]])) {
          xn <- colnames(x$results$effects[[j]])
          xn <- xn[!(xn %in% c("lower", "upper", "fit"))]
          if(length(xn) < 2) {
            plot_smooth_effect(x$results$effects[[j]], ylim = ylim[[p]], ...)
          } else {
            if(any(sapply(x$results$effects[[j]], is.factor))) {
              plot_smooth_effect(x$results$effects[[j]], ylim = ylim[[p]], ...)
            } else {
              plot_smooth_effect_2d(x$results$effects[[j]], ylim = ylim[[p]], ...)
            }
          }
        } else {
          plot_factor_effect(x$results$effects[[j]], ylim = ylim[[p]], ...)
        }
      }
    }

    ## No further plotting.
    return(invisible(NULL))
  }

  spare <- spar

  ## Residual plot.
  if(any(grepl("resid", which))) {
    if(length(which) < 2L)
      spare <- FALSE

    if(spare) {
      owd <- par(no.readonly = TRUE)
      on.exit(par(owd))
    }

    ## Number of plots.
    if(spare)
      par(mfrow = n2mfrow(length(which)))

    ## Compute residuals.
    resids <- residuals(x, type = "quantile", ...)

    if("hist-resid" %in% which) {
      plot_hist(resids, ...)
    }

    if("qq-resid" %in% which) {
      plot_qq(resids, ...)
    }

    if("wp-resid" %in% which) {
      plot_wp(resids, ...)
    }

    if("scatter-resid" %in% which) {
      p <- predict(x, type = "parameter", drop = FALSE, ...)
      m <- family(x)$quantile(0.5, p)
      n <- if(is.null(dim(p))) length(p) else nrow(p)
      if(length(m) != n) {
        if(is.null(family(x)$mean)) {
          m <- p[[1L]]
        } else {
          m <- family(x)$mean(p)
        }
      }
      plot_sr(m, resids, ...)
    }
  }

  ## Stepwise model selection.
  if(("selection" %in% which) & !is.null(x$selection)) {
    ylim <- list(...)$ylim
    labels <- list(...)$labels
    if(is.null(labels))
      labels <- TRUE
    xlim <- list(...)$xlim
    pch <- list(...)$pch
    if(is.null(pch))
      pch <- 16
    type <- list(...)$type
    if(is.null(type))
      type <- "b"
    col <- list(...)$col
    if(is.null(col))
      col <- 4
    xlab <- list(...)$xlab
    if(is.null(xlab))
      xlab <- "Stepwise Iteration"
    ylab <- list(...)$ylab
    if(is.null(ylab))
      ylab <- paste("-2 * logLik +", round(x$selection$K, 4), "* df")
    main <- list(...)$main
    if(is.null(main))
      main <- "GAIC Path"

    n <- length(x$selection$GAIC) - 1L

    plot(0:n, x$selection$GAIC, type = type, xlab = xlab,
      ylab = ylab, main = main, col = col, pch = pch,
      ylim = ylim, xlim = xlim)

    if(type == "b") {
      points(0:n, x$selection$GAIC, type = type)
      points(0:n, x$selection$GAIC)
    }

    if(is.logical(labels)) {
      if(labels) {
        nl <- names(x$selection$GAIC)
        nl[1L] <- "nullmodel"
        pos <- rep(3L, length(nl))
        pos[1L] <- 4L
        text(0:n, x$selection$GAIC, nl, pos = pos)
      }
    }
  }
}

## Plot univariate smooth effects.
plot_smooth_effect <- function(x, col = NULL, ncol = 20L,
  xlab = NULL, ylab = NULL, main = NULL,
  xlim = NULL, ylim = NULL, ...)
{
  if(is.null(col)) {
    col <- gray.colors(ncol, start = 0.3, end = 1)
  } else {
    if(is.function(col))
      col <- col(ncol)
  }
  ncol <- length(col)
  if(!("lower" %in% names(x)))
    x$lower <- x$upper <- x$fit
  polyl <- apply(x[, c("lower", "fit", "upper")], 1,
    function(x)  {
      seq(x[1], x[2], length = ncol)
    }
  )
  polyu <- apply(x[, c("lower", "fit", "upper")], 1,
    function(x)  {
      rev(seq(x[2], x[3], length = ncol))
    }
  )
  px <- x[, 1L]
  fx <- x[, "fit"]
  colr <- rev(col)
  if(is.null(xlim))
    xlim <- range(x[, 1L])
  if(is.null(ylim))
    ylim <- range(x[, c("fit", "lower", "upper")])
  if(is.null(xlab))
    xlab <- colnames(x)[1L]
  if(is.null(ylab))
    ylab <- attr(x, "label")
  add <- list(...)$add
  if(is.null(add))
    add <- FALSE
  if(!add) {
    plot(1, 1, xlim = xlim, ylim = ylim,
      xlab = xlab, ylab = ylab, type = "n", main = main)
  }
  for(i in 1:nrow(polyl)) {
    polygon(cbind(c(px, rev(px)), c(polyl[i, ], rev(fx))),
      col = colr[i], border = NA)
    polygon(cbind(c(px, rev(px)), c(fx, rev(polyu[i, ]))),
      col = colr[i], border = NA)
  }
  addline <- list(...)$addline
  if(is.null(addline))
    addline <- TRUE
  if(addline)
    lines(x[, "fit"] ~ x[, 1L], lwd = 1)
}

## Plot bivariate smooth effects.
plot_smooth_effect_2d <- function(x, col = NULL, ncol = 20L,
  xlab = NULL, ylab = NULL, zlab = NULL, main = NULL,
  xlim = NULL, ylim = NULL,
  persp = FALSE, contour = TRUE, symmetric = TRUE,
  theta = 40, phi = 40, expand = 0.9, ticktype = "simple", ...)
{
  n <- sqrt(nrow(x))
  x1 <- seq(min(x[, 1L]), max(x[, 1L]), length = n)
  x2 <- seq(min(x[, 2L]), max(x[, 2L]), length = n)
  m <- matrix(x$fit, n, n)

  if(is.null(col))
    col <- colorspace::diverge_hcl

  if(is.null(xlab))
    xlab <- names(x)[1L]
  if(is.null(ylab))
    ylab <- names(x)[2L]
  if(is.null(ylim)) {
    ylim <- range(x[, "fit"])
  }
  ylim0 <- ylim
  if(symmetric)
    ylim <- c(-max(abs(ylim)), max(abs(ylim)))

  image <- list(...)$image
  if(!is.null(image)) {
    if(image)
      persp <- FALSE
  }

  if(persp) {
    pal <- make_pal(col, ncol = ncol, data = m,
      range = ylim, symmetric = symmetric)
    col <- pal$map(m)
  } else {
    pal <- make_pal(col, ncol = ncol, data = x[, "fit"],
      range = ylim, symmetric = symmetric)
    col <- pal$colors
  }

  if(persp) {
    if(is.null(zlab))
      zlab <- attr(x, "label")
    pmat <- persp(x1, x2, m, xlab = xlab, ylab = ylab, col = col,
      theta = theta, phi = phi, zlab = zlab, expand = expand,
      ticktype = ticktype, zlim = ylim0, border = NA)
    cl <- contourLines(x1, x2, m)
    if(length(cl)) {
      for(i in 1:length(cl)) {
        lines(trans3d(cl[[i]]$x, cl[[i]]$y, cl[[i]]$level, pmat), col = "black")
        tco <- trans3d(mean(cl[[i]]$x), mean(cl[[i]]$y), cl[[i]]$level, pmat)
        text(tco[1L], tco[2L], cl[[i]]$level, col = "black",
          pos = 3, cex = 0.6, offset = 0.5)
      }
    }
  } else {
    if(is.null(main)) {
      main <- attr(x, "label")
    }
    image(x1, x2, m, main = main,
      xlab = xlab, ylab = ylab, zlim = ylim0, col = col,
      xlim = range(x1), ylim = range(x2), breaks = pal$breaks)
    contour(x1, x2, m, add = TRUE)
  }
}

## Factor effects.
plot_factor_effect <- function(x, col = NULL, ncol = -1L, width = 0.6,
  xlab = NULL, ylab = NULL, main = NULL,
  xlim = NULL, ylim = NULL, ...)
{
  if(ncol < 0L) {
    ncol <- if(nrow(x) > 10L) 1L else 50L
  }
  if(is.null(col)) {
    col <- gray.colors(ncol, start = 0.3, end = 1)
  } else {
    if(is.function(col))
      col <- col(ncol)
  }
  ncol <- length(col)
  px <- t(x[, c("lower", "fit", "upper")])
  colnames(px) <- levels(x[[1L]])
  pos <- 1:ncol(px)
  if(is.null(xlim))
    xlim <- c(0.5, ncol(px) + 0.5)
  if(is.null(ylim))
    ylim <- range(px)
  if(is.null(xlab))
    xlab <- colnames(x)[1L]
  if(is.null(ylab))
    ylab <- attr(x, "label")
  plot(pos, xlim = xlim, ylim = ylim,
    xlab = xlab, ylab = ylab,
    axes = FALSE, type = "n")
  axis(1, at = pos, labels = colnames(px))
  width <- width/2
  coll <- rev(col)
  colu <- col
  if(ncol > 1) {
    for(j in 1:ncol(px)) {
      xr <- range(px[, j])
      xrl <- seq(xr[1], px["fit", j], length = ncol)
      xru <- rev(seq(px["fit", j], xr[2], length = ncol))
      for(i in 1:ncol) {
        rect(pos[j]-width, xrl[i], pos[j]+width, px["fit", j],
          col = coll[i], border = NA)
        rect(pos[j]-width, px["fit", j], pos[j]+width, xru[i],
          col = coll[i], border = NA)
      }
    }
  }
  for(j in 1:ncol(px))
    lines(c(pos[j]-width, pos[j]+width), rep(px["fit", j], 2L), lwd = 1)
  box()
  axis(2)
}

## Histogram and density plot.
plot_hist <- function(x, ...)
{
  x <- na.omit(x)
  h <- hist(x, breaks = "Scott", plot = FALSE)
  d <- density(x)
  ylim <- list(...)$ylim
  if(is.null(ylim))
    ylim <- range(c(h$density, d$y))
  main <- list(...)$main
  if(is.null(main))
    main <- "Histogram and Density"
  xlab <- list(...)$xlab
  if(is.null(xlab))
    xlab <- "Quantile Residuals"
  hist(x, breaks = "Scott", freq = FALSE, ylim = ylim,
    xlab = xlab, main = main, ...)
  lines(d, lwd = 2, col = 4)
  rug(x, col = rgb(0.1, 0.1, 0.1, alpha = 0.3))
}

## Q-Q plot.
plot_qq <- function(x, ...)
{
  z <- qnorm(ppoints(length(x)))
  pch <- list(...)$pch
  if(is.null(pch))
    pch <- 19
  col <- list(...)$col
  if(is.null(col))
    col <- rgb(0.1, 0.1, 0.1, alpha = 0.3)
  qqnorm(x, col = col, pch = pch)
  lines(z, z, lwd = 2, col = 4)
}

## Wormplot.
plot_wp <- function(x, ...)
{
  x <- na.omit(x)
  d <- qqnorm(x, plot = FALSE)
  probs <- c(0.25, 0.75)
  y3 <- quantile(x, probs, type = 7, na.rm = TRUE)
  x3 <- qnorm(probs)
  slope <- diff(y3)/diff(x3)
  int <- y3[1L] - slope * x3[1L]
  d$y <- d$y - (int + slope * d$x)

  xlim <- list(...)$xlim
  if(is.null(xlim)) {
    xlim <- range(d$x)
    xlim <- xlim + c(-0.1, 0.1) * diff(xlim)
  }

  ylim <- list(...)$ylim
  if(is.null(ylim)) {
    ylim <- range(d$y)
    ylim <- ylim + c(-0.3, 0.3) * diff(ylim)
  }

  main <- list(...)$main
  if(is.null(main))
    main <- "Worm Plot"
  xlab <- list(...)$xlab
  if(is.null(xlab))
    xlab <- "Theoretical Quantiles"
  ylab <- list(...)$ylab
  if(is.null(ylab))
    ylab <- "Deviation"
  pch <- list(...)$pch
  if(is.null(pch))
    pch <- 19
  col <- list(...)$col
  if(is.null(col))
    col <- rgb(0.1, 0.1, 0.1, alpha = 0.3)

  plot(d$x, d$y, xlim = xlim, ylim = ylim,
    xlab = xlab, ylab = ylab, main = main,
    col = col, pch = pch)
  grid(lty = "solid")

  dz <- 0.25
  z <- seq(xlim[1L], xlim[2L], dz)
  p <- pnorm(z)
  se <- (1/dnorm(z)) * (sqrt(p * (1 - p)/length(d$y)))
  low <- qnorm((1 - 0.95)/2) * se
  high <- -low
  lines(z, low, lty = 2)
  lines(z, high, lty = 2)

  fit <- lm(d$y ~ d$x + I(d$x^2) + I(d$x^3))
  i <- order(d$x)
  lines(d$x[i], fitted(fit)[i], col = 4, lwd = 2)
}

## Fitted vs. residuals.
plot_sr <- function(f, x, ...) {
  if(any(is.na(x))) {
    i <- which(is.na(x))
    x <- x[-i]
    f <- f[-i]
  }
  main <- list(...)$main
  if(is.null(main))
    main <- "Against Fitted Values"
  xlab <- list(...)$xlab
  if(is.null(xlab))
    xlab <- "Fitted Values (Median)"
  ylab <- list(...)$ylab
  if(is.null(ylab))
    ylab <- "Quantile Residuals"
  pch <- list(...)$pch
  if(is.null(pch))
    pch <- 19
  col <- list(...)$col
  if(is.null(col))
    col <- rgb(0.1, 0.1, 0.1, alpha = 0.3)
  plot(f, x, xlab = xlab, ylab = ylab, main = main,
    col = col, pch = pch)
  abline(h = 0, col = "lightgray")

  m <- lm(x ~ f + I(f^2) + I(f^3))
  i <- order(f)
  lines(f[i], fitted(m)[i], col = 4, lwd = 2)
}

## Create color palette.
make_pal <- function(col, ncol = NULL, data = NULL, range = NULL, 
  breaks = NULL, swap = FALSE, symmetric = TRUE) 
{
  if(is.null(symmetric))
    symmetric <- TRUE
  if(is.null(col))
    col <- colorspace::diverge_hcl
  if(is.null(ncol) && is.null(breaks))
    ncol <- 99L
  if(is.null(ncol) && !is.null(breaks))
    ncol <- length(breaks) - 1L
  if(!is.null(ncol) && symmetric) {
    if(ncol %% 2 == 0)
      ncol <- ncol + 1L
  }
  if(is.function(col))
    col <- col(ncol)    
  else 
    ncol <- length(col)
  if(swap) 
    col <- rev(col)
  if(all(is.null(data), is.null(range), is.null(breaks))) 
    stop("at least one needs to be specified")
  if(is.null(breaks)) {
    if(is.null(range)) {
      range <- range(data, na.rm = TRUE)
      if(symmetric) { 
        mar <- max(abs(range))
        range <- c(0 - mar, mar)
      }
    }
    if(diff(range) == 0)
      breaks <- seq(min(range) - 1, min(range) + 1, length.out = ncol + 1L)
    else
      breaks <- seq(range[1L], range[2L], length.out = ncol + 1L)
  } else stopifnot(length(breaks) == ncol + 1L)
  if(is.matrix(data)) {
    obs2col <- function(x) {
      hgt <- (x[-1L, -1L] + x[-1L, -(ncol(x) - 1L)] + 
        x[-(nrow(x) -1L), -1L] + x[-(nrow(x) -1L), -(ncol(x) - 1L)])/4
      c(col[1L], col, col[ncol])[cut(hgt, c(-Inf, breaks, Inf), include.lowest = TRUE)]
    }
  } else {
    obs2col <- function(x) c(col[1L], col, col[ncol])[cut(x, c(-Inf, breaks, Inf))]
  }

  return(list(colors = col, breaks = breaks, map = obs2col))
}

## Plot a list of objects.
plot.gamlss2.list <- function(x, parameter = NULL, which = "effects", terms = NULL, spar = TRUE, legend = TRUE, ...)
{
  ## What should be plotted?
  which.match <- c("effects", "hist-resid", "qq-resid", "wp-resid", "scatter-resid")
  if(!is.character(which)) {
    if(any(which > 5L))
      which <- which[which <= 5L]
    which <- which.match[which]
  } else which <- which.match[grep(tolower(which), which.match, fixed = TRUE)]
  if(length(which) > length(which.match) || !any(which %in% which.match))
    stop("argument which is specified wrong!")

  if("effects" %in% which) {
    ok <- sapply(x, function(z) { !is.null(z$results$effects) })
    x <- x[ok]
    if(length(x)) {
      effects <- unique(unlist(sapply(x, function(z) { names(z$results$effects) })))
      e2 <- strsplit(effects, ",", fixed = TRUE)
      e2 <- sapply(e2, function(x) { paste0(x[-grep(")", x, fixed = TRUE)], collapse = ",") })
      e2 <- unique(e2)
      eff <- list()
      for(j in e2) {
        jn <- paste0(j, ")")
        for(i in seq_along(x)) {
          jj <- grep(j, names(x[[i]]$results$effects), fixed = TRUE)
          if(length(jj)) {
            xn <- colnames(x[[i]]$results$effects[[jj]])
            xn <- xn[!(xn %in% c("fit", "lower", "upper"))]
            if(length(xn) < 2L) {
              if(length(eff[[jn]]) < 1L) {
                eff[[jn]] <- x[[i]]$results$effects[[jj]][, c(xn, "fit")]
                colnames(eff[[jn]]) <- c(xn, names(x)[i])
              } else {
                fj <- data.frame(x[[i]]$results$effects[[jj]][, "fit"])
                colnames(fj) <- names(x)[i]
                eff[[jn]] <- cbind(eff[[jn]], "fit" = fj)
              }
            }
          }
        }
      }
    }

    if(is.null(parameter)) {
      parameter <- list(...)$what
      if(is.null(parameter))
      parameter <- list(...)$model
      if(is.null(parameter))
        parameter <- x$family$names
    }

    if(!is.null(parameter)) {
      if(is.character(parameter)) {
        ok <- grep2(parameter, names(eff), fixed = TRUE, value = TRUE)
      } else {
        ok <- parameter
      }
      eff <- eff[ok]   
    }
    if(length(eff) > 0L) {
      if(!is.null(terms)) {
        if(is.character(terms)) {
          ok <- grep2(terms, names(eff), fixed = TRUE, value = TRUE)
        } else {
          ok <- terms
        }
        eff <- eff[ok]
      }
    }
    if(length(eff) > 0L) {
      if(spar) {
        owd <- par(no.readonly = TRUE)
        on.exit(par(owd))
        par(mfrow = n2mfrow(length(eff)))
      }
      if(legend) {
        pos <- list(...)$pos
        if(is.null(pos))
          pos <- "topright"
        pos <- rep(pos, length.out = length(eff))
      }
      col <- list(...)$col
      if(is.null(col))
        col <- 1L
      lwd <- list(...)$lwd
      if(is.null(lwd))
        lwd <- 2L
      for(j in seq_along(eff)) {
        nm <- ncol(eff[[j]]) - 1L
        o <- order(eff[[j]][, 1L])
        eff[[j]] <- eff[[j]][o, ]
        lty <- list(...)$lty
        if(is.null(lty))
          lty <- 1:nm
        matplot(eff[[j]][, 1L], eff[[j]][, -1L],
          type = "l", col = col, lwd = lwd, lty = lty,
          xlab = colnames(eff[[j]])[1L], ylab = names(eff)[j])
        if(legend) {
          legend(pos[j], colnames(eff[[j]])[-1L],
            lwd = lwd, lty = lty, bty = "n", col = col)
        }
      }
    }
  }

  return(invisible(NULL))
}


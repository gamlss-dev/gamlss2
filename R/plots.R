## A plotting method.
plot.gamlss2 <- function(x, parameter = NULL,
  which = "effects", terms = NULL,
  scale = TRUE, spar = TRUE, ...)
{
  if(spar) {
    owd <- par(no.readonly = TRUE)
    on.exit(par(owd))
  }

  ## Which parameters to plot?
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

  ## Effect plots.
  if("effects" %in% which) {
    if(is.null(x$results))
      x$results <- results(x)

    en <- grep2(parameter, names(x$results$effects), fixed = TRUE, value = TRUE)
    if(!is.null(terms))
      en <- grep2(terms, en, fixed = TRUE, value = TRUE)

    if(length(x$results$effects)) {
      if(spar)
        par(mfrow = n2mfrow(length(en)))

      ylim <- list(...)$ylim

      if(scale > 0) {
        ylim <- list()
        for(i in parameter) {
          for(j in grep(i, en, fixed = TRUE, value = TRUE))
            ylim[[i]] <- c(ylim[[i]], range(x$results$effects[[j]][, c("lower", "upper")]))
          ylim[[i]] <- range(ylim[[i]])
        }
      }

      for(j in en) {
        p <- strsplit(j, ".", fixed = TRUE)[[1L]][1L]
        if(!is.factor(x$results$effects[[j]][[1L]])) {
          plot_smooth_effect(x$results$effects[[j]], ylim = ylim[[p]], ...)
        } else {
          plot_factor_effect(x$results$effects[[j]], ylim = ylim[[p]], ...)
        }
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
  plot(1, 1, xlim = xlim, ylim = ylim,
    xlab = xlab, ylab = ylab, type = "n", main = main)
  for(i in 1:nrow(polyl)) {
    polygon(cbind(c(px, rev(px)), c(polyl[i, ], rev(fx))),
      col = colr[i], border = NA)
    polygon(cbind(c(px, rev(px)), c(fx, rev(polyu[i, ]))),
      col = colr[i], border = NA)
  }
  lines(x[, "fit"] ~ x[, 1L], lwd = 1)
}

## Factor effects.
plot_factor_effect <- function(x, col = NULL, ncol = 20L, width = 0.6,
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
    axes = FALSE)
  axis(1, at = pos, labels = colnames(px))
  width <- width/2
  coll <- rev(col)
  colu <- col
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
  for(j in 1:ncol(px))
    lines(c(pos[j]-width, pos[j]+width), rep(px["fit", j], 2L), lwd = 1)
  box()
  axis(2)
}


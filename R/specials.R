## Function simply evaluates the
## special terms in the model formula
## and assigns appropriate model fitting
## functions for the backfitting steps.
special_terms <- function(x, data, ...)
{
  sterms <- list()

  if(length(x)) {
    for(j in unlist(x)) {
      sj <- eval(parse(text = j), envir = data)
      if(any(grepl(".smooth.spec", class(sj)))) {
        stopifnot(requireNamespace("mgcv"))
        knots <- list(...)$knots
        sj <- mgcv::smoothCon(sj, data = data, knots = knots,
          absorb.cons = TRUE, scale.penalty = TRUE)
        for(i in 1:length(sj))
          sj[[i]]$orig.label <- j
        sjn <- sapply(sj, function(x) x$label)
        names(sj) <- sjn
        sterms <- c(sterms, sj)
      } else {
        sterms[[j]] <- sj
      }
    }
  }

  return(sterms)
}

## Special term fit function, works with gamlss model terms, too.
special.wfit <- function(x, z, w, y, eta, j, family, control, ...)
{
  if(inherits(x, "smooth")) {
    call <- attr(x, "call")
    call[[2]] <- quote(x)
    fe <- eval(call)
    if(!is.null(fe$y))
      fe$fitted.values <- fe$y
    fit <- list(
      "fitted.values" = as.numeric(fe$fitted.values),
      "coefficients" = fe$coefSmo,
      "lambdas" = fe$lambda,
      "edf" = fe$nl.df,
      "df" = length(z) - fe$nl.df,
      "model" = fe$model
    )
  } else {
    if(inherits(x, "special")) {
      fit <- special_fit(x, z, w, y, eta, j, family, control, ...)
    } else {
      ff <- if(is.null(x$special.wfit)) {
        smooth.construct.wfit
      } else {
        x$special.wfit
      }
      fit <- ff(x, z, w, y, eta, j, family, control)
    }
  }
  return(fit)
}

## Fitting function for mgcv smooth terms.
smooth.construct.wfit <- function(x, z, w, y, eta, j, family, control)
{
  ## Number of observations.
  n <- nrow(x$X)

  ## Set up smoothing parameters.
  lambdas <- x$lambdas
  if(is.null(lambdas))
    lambdas <- 10
  lambdas <- rep(lambdas, x$dim)

  ## Pre compute matrices.
  XW <- x$X * w
  XWX <- crossprod(XW, x$X)
  XWz <- crossprod(XW, z)
  S <- diag(1e-05, ncol(x$X))
  
  ## Function to search for smoothing parameters using GCV.
  fl <- function(l, rf = FALSE) {
    for(j in 1:length(x$S))
      S <- S + l[j] * x$S[[j]]

    P <- try(chol2inv(chol(XWX + S)), silent = TRUE)
    if(inherits(P, "try-error"))
      P <- solve(XWX + S)

    b <- drop(P %*% XWz)

    fit <- drop(x$X %*% b)
    edf <- sum(diag(XWX %*% P))

    if(rf) {
      return(list("coefficients" = b, "fitted.values" = fit, "edf" = edf,
        "lambdas" = l, "vcov" = P, "df" = nrow(x$X) - edf))
    } else {
      if(is.null(control$criterion))
        control$criterion <- "gcv"

      rss <- sum(w * (z - fit)^2)

      rval <- switch(tolower(control$criterion),
        "gcv" = rss * n / (n - edf)^2,
        "aic" = rss + 2 * edf,
        "bic" = rss + log(n) * edf
      )

      return(rval)
    }
  }

  opt <- nlminb(lambdas, objective = fl)

  return(fl(opt$par, rf = TRUE))
}

## A method for fitting special terms.
special_fit <- function(x, ...)
{
  UseMethod("special_fit")
}

## A method for predicting special terms.
special_predict <- function(x, ...)
{
  UseMethod("special_predict")
}

## Default method.
special_predict.default <- function(x, data, ...)
{
  if(is.null(x$model))
    return(predict(x, newdata = data))
  else
    return(predict(x$model, newdata = data))
}


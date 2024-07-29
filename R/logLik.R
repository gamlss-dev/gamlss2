.logLik.gamlss2 <- function(object, newdata = NULL)
{
  if(is.null(newdata)) {
    ll <- object$logLik
    nobs <- object$nobs
  } else {
    par <- predict(object, type = "parameter", newdata = newdata)
    y <- if(!is.null(newdata)) {
      model.response(model.frame(object, data = newdata, keepresponse = TRUE))
    } else {
      model.response(model.frame(object, keepresponse = TRUE))
    }
    ll <- family(object)$loglik(y, par)
    nobs <- length(y)
  }
  attr(ll, "nobs") <- nobs
  attr(ll, "df") <- object$df
  class(ll) <- "logLik"
  return(ll)
}

logLik.gamlss2 <- function(object, ..., newdata = NULL)
{
  objs <- list(object, ...)

  ll <- list()
  for(j in 1:length(objs)) {
    ll[[j]] <- .logLik.gamlss2(objs[[j]], newdata = newdata)
  }

  if(length(ll) > 1L) {
    Call <- match.call()
    rn <- as.character(Call[-1L])
    if("..1" %in% rn) {
      rn <- as.character(sys.call(-1))[-1L]
    }
    i <- order(unlist(ll), decreasing = TRUE)
    ll <- ll[i]
    rn <- rn[i]
    if(any(j <- duplicated(rn)))
      rn[j] <- paste0(rn[j], ".", 1:sum(j))
    names(ll) <- rn
  } else {
    ll <- ll[[1L]]
  }

  return(ll)
}

get_df <- function(object)
{
  df <- 0
  if(length(object$xterms)) {
    if(is.null(attr(object$fitted.linear, "edf"))) {
      df <- df + sum(sapply(object$fitted.linear, function(x) length(x$coefficients)))
    } else {
      df <- df + sum(attr(object$fitted.linear, "edf"))
    }
  }
  if(length(object$sterms)) {
    for(j in seq_along(object$fitted.specials)) {
      dfj <- sapply(object$fitted.specials[[j]], function(x) x$edf)
      df <- df + sum(unlist(dfj))
    }
  }
  return(df)
}

deviance.gamlss2 <- function(object, ..., newdata = NULL)
{
  objs <- list(object, ...)

  devs <- NULL
  for(j in 1:length(objs)) {
    devs <- c(devs, -2 * as.numeric(logLik(objs[[j]], newdata = newdata)))
  }

  if(length(devs) > 1L) {
    Call <- match.call()
    rn <- as.character(Call[-1L])
    if("..1" %in% rn) {
      rn <- as.character(sys.call(-1))[-1L]
    }
    i <- order(devs, decreasing = FALSE)
    devs <- devs[i]
    rn <- rn[i]
    if(any(j <- duplicated(rn)))
      rn[j] <- paste0(rn[j], ".", 1:sum(j))
    names(devs) <- rn
  }

  return(devs)
}

response_name <- function(formula) {
  if(is.list(formula)) {
    if(!is.null(formula$formula)) {
      formula <- formula$formula
    } else {
      formula <- do.call("as.Formula", formula)
    }
  }
  formula <- as.Formula(formula)
  formula <- formula(formula, rhs = 0, collapse = TRUE)
  rn <- all.vars(formula)
  rn <- rn[rn != "."]
  return(rn)
}


logLik.gamlss2 <- function(object, newdata = NULL, ...)
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

get_df <- function(object)
{
  df <- 0
  if(length(object$xterms))
    df <- df + sum(sapply(object$fitted.linear, function(x) length(x$coefficients)))
  if(length(object$sterms)) {
    for(j in seq_along(object$fitted.specials)) {
      dfj <- sapply(object$fitted.specials[[j]], function(x) x$edf)
      df <- df + sum(unlist(dfj))
    }
  }
  return(df)
}

deviance.gamlss2 <- function(object, ...)
{
  -2 * as.numeric(logLik(object))
}


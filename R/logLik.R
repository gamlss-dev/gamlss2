logLik.gamlss2 <- function(object, newdata = NULL, ...)
{
  par <- predict(object, type = "parameter", newdata = newdata)
  y <- if(!is.null(newdata)) {
    model.response(model.frame(object, data = newdata, keepresponse = TRUE))
  } else {
    model.response(model.frame(object, keepresponse = TRUE))
  }
  ll <- family(object)$loglik(y, par)
  df <- get_df(object)
  attr(ll, "nobs") <- length(y)
  attr(ll, "df") <- df
  class(ll) <- "logLik"
  return(ll)
}

get_df <- function(object)
{
  df <- 0
  if(length(object$xterms))
    df <- df + sum(sapply(object$fitted.linear, function(x) length(x$coefficients)))
  if(length(object$sterms)) {
    for(j in seq_along(object$fitted.specials))
      df <- df + sum(sapply(object$fitted.specials[[j]], function(x) x$edf))
  }
  return(df)
}


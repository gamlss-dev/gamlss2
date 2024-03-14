## S3 method for extracting fitted/predicted distributions3 objects
## associated methods are in gamlss.dist (as well as distributions3, topmodels, etc.)
prodist.gamlss2 <- function(object, ...) {
  ## extract fitted parameters
  d <- predict(object, ..., type = "parameter", drop = FALSE)

  ## set class to general GAMLSS distribution (distributions3 object)
  class(d) <- c("GAMLSS", "distribution")

  ## include family information
  ## (FIXME: For now just use the character label, but moving forward maybe
  ## the functions in object$family should be used)
  attr(d, "family") <- object$family$family

  ## return distributions3 object
  return(d)
}

## newresponse() method for topmodels.
newresponse.gamlss2 <- function(object, newdata, ...)
{
  y <- if(missing(newdata) || is.null(newdata)) {
    model.response(model.frame(object, keepresponse = TRUE))
  } else {
    model.response(model.frame(object, keepresponse = TRUE, data = newdata))
  }
  return(y)
}


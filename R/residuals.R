## Residuals.
residuals.gamlss2 <- function(object,
  type = c("quantile", "response", "parameter"), newdata = NULL, ...)
{
  family <- family(object)

  ## If there is a residuals function in the family
  ## object, use this function.
  if(!is.null(family$residuals)) {
    res <- family$residuals(object, type = type, ...)
    if(length(class(res)) < 2) {
      if(inherits(res, "numeric"))
        class(res) <- c("gamlss2.residuals", class(res))
    }
  } else {
    type <- match.arg(type)

    ## Get the response.
    object$y <- if(is.null(newdata)) {
      model.response(model.frame(object, keepresponse = TRUE))
    } else {
      model.response(model.frame(object, keepresponse = TRUE, data = newdata))
    }

    if(is.null(object$y))
      stop("response variable is missing, cannot compute residuals!")

    nobs <- nrow(object$y)
    y <- if(is.data.frame(object$y)) {
      if(ncol(object$y) < 2) {
        object$y[[1]]
      } else object$y
    } else {
      object$y
    }

    ## Predict parameters.
    par <- predict(object, drop = FALSE, newdata = newdata, type = "parameter", ...)

    ## Randomized quantile residuals.
    if(type == "quantile") {
      if(!is.null(family$rqres)) {
        res <- family$rqres(y, par)
      } else {
        if(is.null(family$p)) {
          type <- "response"
          warning(paste("no $p() function in family '", family$family,
            "', cannot compute quantile residuals, computing response resdiuals instead!", sep = ""))
        } else {
          discrete <- FALSE
          if(!is.null(family$type)) {
            if(tolower(family$type) == "discrete")
              discrete <- TRUE
          }
          if(family$family == "binomial")
            discrete <- TRUE

          ## Discrete case.
          if(discrete) {
            ymin <- min(y, na.rm = TRUE)
            a <- family$p(ifelse(y == ymin, y, y - 1), par)
            a <- ifelse(y == ymin, 0, a)
            b <- family$p(y, par)
            u <- runif(length(y), a, b)
            u <- ifelse(u > 0.999999, u - 1e-16, u)
            u <- ifelse(u < 1e-06, u + 1e-16, u)
            res <- qnorm(u)
          } else {
            if(family$type == "mixed") {
              mass.p <- list(...)$mass.p
              prob.mp <- list(...)$prob.mp

              if(is.null(dim(mass.p)))
                mass.p <- matrix(mass.p, ncol = 1L)
              if(is.null(dim(prob.mp)))
                prob.mp <- matrix(prob.mp, ncol = 1L)

              if(ncol(mass.p) != ncol(prob.mp))
                stop()

              res <- family$p(y, par)
              for(j in mass.p) {
                i <- which(y == j)
                if(j == 1) {
                  pj <- 1
                }
              }
            } else {
              ## Continuous case.
              prob <- family$p(y, par)
              thres <- 0.999999999999999
              prob[prob > thres] <- thres
              prob[prob < (1 - thres)] <- 1 - thres
              res <- qnorm(prob)
            }
          }
        }
      }

      if(any(isnf <- !is.finite(res))) {
        warning("non finite quantiles from probabilities, set to NA!")
        res[isnf] <- NA
      }

      attr(res, "type") <- "quantile"
    }

    ## Response residuals.
    if(type == "response") {
      mu <- if(is.null(family$mu)) {
        function(par, ...) { par[[1]] }
      } else family$mu
      res <- y - mu(par)
      attr(res, "type") <- "response"
    }

    ## FIXME: residuals for one parameter.
    if(type == "parameter") {
      stop("Not implemented yet!")
    }
   
    class(res) <- c("gamlss2.residuals", class(res))
  }

  if(any(j <- !is.finite(res)))
    res[j] <- NA

  return(res)
}


## Function simply evaluates the
## sepcial terms in the model formula
## and assigns appropriate model fitting
## functions for the backfitting steps.
special_terms <- function(x, data, ...)
{
  sterms <- list()

  if(length(x)) {
    for(j in unlist(x)) {
      sterms[[j]] <- eval(parse(text = j), envir = data)
      if(any(grepl(".smooth.spec", class(sterms[[j]])))) {
        stopifnot(requireNamespace("mgcv"))
        knots <- list(...)$knots
        sterms[[j]] <- mgcv::smoothCon(sterms[[j]], data = data, knots = knots)[[1L]]
      }
    }
  }

  return(sterms)
}

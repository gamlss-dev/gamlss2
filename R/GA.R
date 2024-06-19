genetic <- function(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
{
  ## Extract parameter names.
  nx <- family$names

  ## Models must include intercepts!
  for(i in nx) {
    if(!any(grepl("(Intercept)", xterms[[i]]))) {
      stop(paste0("intercept is missing for parameter ", i, "!"))
    }
  }

  ## Inner and outer trace.
  trace <- control$trace
  if(length(trace) < 2L) {
    trace <- c(FALSE, trace)
  }
  control$trace <- trace[1L]

  ## Penalty for AIC.
  K <- if(is.null(control$K)) 2 else control$K

  ## All possible model terms.
  modelterms <- NULL
  for(i in nx) {
    if(length(xterms[[i]]))
      modelterms <- c(modelterms, paste0(i, ".p.", xterms[[i]]))
    if(length(sterms[[i]]))
      modelterms <- c(modelterms, paste0(i, ".s.", sterms[[i]]))
  }
  modelterms <- modelterms[-grep(".p.(Intercept)", modelterms, fixed = TRUE)]

  ## Helper function to create xterms and sterms.
  make_xs_terms <- function(x) {
    xterms <- sterms <- list()
    for(i in nx) {
      xterms[[i]] <- "(Intercept)"
      sterms[[i]] <- character(0)
    }
    for(j in x) {
      js <- splitname(j)
      if(linear(j)) {
        xterms[[js[1L]]] <- c(xterms[[js[1L]]], js[2L])
      } else {
        sterms[[js[1L]]] <- c(sterms[[js[1L]]], js[2L])
      }
    }

    return(list("xterms" = xterms, "sterms" = sterms))
  }

  ## The GAMLSS fitness function.
  fitness <- function(string) {
    i <- which(string == 1)
    xs <- make_xs_terms(modelterms[i])
    m <- RS(x, y, specials[unique(unlist(xs$sterms))], family, offsets,
      weights, start = start, xterms = xs$xterms, sterms = xs$sterms, control)
    return(-1 * modelstats(m, K)$GAIC)
  }

  ## ga() specific controls.
  gac <- control$ga
  if(!is.list(gac))
    gac <- as.list(gac)
  if(is.null(gac$maxiter))
    gac$maxiter <- 100L
  if(is.null(gac$parallel))
    gac$parallel <- FALSE

  ## Start the GA.
  b <- GA::ga("binary", fitness = fitness, nBits = length(modelterms),
    maxiter = gac$maxiter, parallel = gac$parallel,
    monitor = trace[2L])

  ## Optimum solution.
  i <- which(b@solution == 1)
  xs <- make_xs_terms(modelterms[i])

  ## Reestimate final model.
  model <- RS(x, y, specials[unique(unlist(xs$sterms))], family, offsets,
    weights, start = start, xterms = xs$xterms, sterms = xs$sterms, control)

  ## Collect all info.
  model$selection <- list("GAIC" = modelstats(model, K)$GAIC, "logLik" = model$logLik,
    "df" = model$df, "K" = K, "formula" = xs2formula(xs$xterms, xs$sterms))
  model$xterms <- xs$xterms
  model$sterms <- xs$sterms
  model$specials <- specials[unique(unlist(xs$sterms))]
  model$genetic <- b

  if(trace[2L]) {
    cat("Final Model\n")
    for(j in names(model$selection$formula)) {
      cat(paste0("$", j, "\n.. "))
      print(model$selection$formula[[j]])
    }
  }

  return(model)
}

gaGAMLSS <- function(formula, ..., K = 2, maxiter = 100L, cores = 1L)
{
  gamlss2(formula, ..., optimizer = genetic, K = K,
    ga = list("parallel" = if(cores < 2L) FALSE else cores,
      "maxiter" = maxiter))
}


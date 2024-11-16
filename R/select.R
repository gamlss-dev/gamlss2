## Model term selection based on null space penalties.
select_gamlss2 <- function(formula, ..., criterion = "BIC", thres = c(0.9, 0.2))
{
  criterion <- rep(criterion, length.out = 2L)

  m <- match.call()
  m[[1L]] <- as.name("gamlss2")
  m["criterion"] <- criterion[1L]
  m["criterion_refit"] <- criterion[2L]
  m["select"] <- TRUE
  m["thres"] <- parse(text = paste("c(", thres[1L], ",", thres[2L], ")"))
  m["optimizer"] <- parse(text = "sRS")

  model <- eval(m, parent.frame())

  return(model)
}

## Internal selection optimizer function.
sRS <- function(x, y, specials, family, offsets, weights,
  start, xterms, sterms, control)
{
  trace <- control$trace
  control$trace <- FALSE

  if(trace)
    cat(".. selection step\n")

  m <- RS(x, y, specials, family, offsets, weights, start, xterms, sterms, control)

  refit <- control$refit
  if(is.null(refit))
    refit <- TRUE

  if(length(m$fitted.specials) & refit) {
    drop <- fr <- list()

    thres <- control$thres
    if(length(thres) < 2)
      thres <- c(thres, 0.2)

    ## edf threshold.
    for(j in names(m$fitted.specials)) {
      if(length(m$fitted.specials[[j]])) {
        for(i in names(m$fitted.specials[[j]])) {
          fr[[j]] <- c(fr[[j]], range(m$fitted.specials[[j]][[i]]$fitted.values))
          if(m$fitted.specials[[j]][[i]]$edf <= thres[1L])
            if(inherits(specials[[i]], "mgcv.smooth"))
              drop[[j]] <- unique(c(drop[[j]], i))
        }
      }
    }

    ## range threshold.
    fr <- lapply(fr, function(x) diff(range(x)))
    for(j in names(m$fitted.specials)) {
      if(length(m$fitted.specials[[j]])) {
        for(i in names(m$fitted.specials[[j]])) {
          fri <- diff(range(m$fitted.specials[[j]][[i]]$fitted.values))
          if(fri/fr[[j]] <= thres[2L]) {
            if(inherits(specials[[i]], "mgcv.smooth"))
              drop[[j]] <- unique(c(drop[[j]], i))
          }
        }
      }
    }

    if(length(drop)) {
      for(j in names(drop))
        sterms[[j]] <- sterms[[j]][!(sterms[[j]] %in% drop[[j]])]

      specials <- specials[names(specials) %in% unlist(sterms)]

      for(j in names(specials)) {
        if(inherits(specials[[j]], "mgcv.smooth")) {
          specials[[j]]$S <- specials[[j]]$S[1:(length(specials[[j]]$S) - 1L)]
        }
      }

      if(trace)
        cat(".. refitting step\n")

      control$trace <- trace
      control$criterion <- control$criterion_refit

      m <- RS(x, y, specials, family, offsets, weights, start, xterms, sterms, control)

      m$selection <- list("formula" = xs2formula(xterms, sterms), "select" = TRUE)
    } else {
      if(control$criterion_refit != control$criterion) {
        control$trace <- trace
        control$criterion <- control$criterion_refit
        m <- RS(x, y, specials, family, offsets, weights, start, xterms, sterms, control)
      }
    }
  }

  class(m) <- c(class(m), "select")

  return(m)
}

new_formula.select <- function(object, ...) {
  f <- object$selection$formula
  rn <- response_name(object)
  f[[1L]] <- eval(parse(text = paste("update(f[[1L]], ", rn, " ~ .)")))
  return(f)
}


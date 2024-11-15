## Model term selection based on null space penalties.
select_gamlss2 <- function(formula, ..., criterion = "BIC", thres = 0.2)
{
  criterion <- rep(criterion, length.out = 2L)

  m <- match.call()
  m[[1L]] <- as.name("gamlss2")
  m["criterion"] <- criterion[1L]
  m["criterion_refit"] <- criterion[2L]
  m["select"] <- TRUE
  m["thres"] <- thres
  m["optimizer"] <- expression(.select_gamlss2)

  model <- eval(m, parent.frame())

  return(model)
}

## Internal selection optimizer function.
.select_gamlss2 <- function(x, y, specials, family, offsets, weights,
  start, xterms, sterms, control)
{
  trace <- control$trace
  control$trace <- FALSE

  if(trace)
    cat(".. selection step\n")

  m <- RS(x, y, specials, family, offsets, weights, start, xterms, sterms, control)

  if(length(m$fitted.specials)) {
    drop <- list()
    for(j in names(m$fitted.specials)) {
      if(length(m$fitted.specials[[j]])) {
        for(i in names(m$fitted.specials[[j]])) {
          if(m$fitted.specials[[j]][[i]]$edf <= control$thres)
            if(inherits(specials[[i]], "mgcv.smooth"))
              drop[[j]] <- c(drop[[j]], i)
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

  return(m)
}


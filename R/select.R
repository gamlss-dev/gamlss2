## Model term selection based on null space penalties.
select_gamlss2 <- function(formula, ..., criterion = "BIC", thres = 0.1)
{
  m <- match.call()
  m[[1L]] <- as.name("gamlss2")
  m["criterion"] <- criterion
  m["select"] <- TRUE
  m["thres"] <- thres
  m["optimizer"] <- expression(.select_gamlss2)

  model <- eval(m, parent.frame())

  return(model)
}

.select_gamlss2 <- function(x, y, specials, family, offsets, weights,
  start, xterms, sterms, control)
{
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

      m <- RS(x, y, specials, family, offsets, weights, start, xterms, sterms, control)

      m$selection <- list("formula" = xs2formula(xterms, sterms), "select" = TRUE)
    }
  }

  return(m)
}


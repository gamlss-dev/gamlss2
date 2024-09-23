## Second make.link function.
make.link2 <- function(link)
{
  if(is.null(link))
    link <- "identity"
  if(is.function(link)) {
    rval <- link()
    if(!all(c("linkfun", "linkinv", "mu.eta", "valideta", "name") %in% names(rval)))
      stop("link is spefified wrong!")
  } else {
    link0 <- link
    if(link0 == "tanhalf"){
      rval <- list(
        "linkfun" = function (mu) {
          tan(mu/2)},
        "linkinv" = function(eta) {
          2 * atan(eta)},
        "mu.eta" = function(eta) {
          2 / (eta^2 + 1)},
        "mu.eta2" = function(eta) {
          (-4 * eta ) / (eta^2 + 1)^2},
        "valideta" = function(eta) TRUE,
        "name" = "tanhalf"
      )
    } else {
      mu.eta2 <- function(x) {
        if(link0 == "identity") {
          x$mu.eta2 <- function(eta) rep.int(0, length(eta))
          return(x)
        }
        if(link0 == "log") {
          x$mu.eta2 <- function(eta) exp(eta)
          return(x)
        }
        if(link0 == "logit") {
          x$mu.eta2 <- function(eta) {
            eta <- exp(eta)
            return(-eta * (eta - 1) / (eta + 1)^3)
          }
          return(x)
        }
        if(link0 == "probit") {
          x$mu.eta2 <- function(eta) {
            -eta * dnorm(eta, mean = 0, sd = 1)
          }
          return(x)
        }
        if(link0 == "inverse") {
          x$mu.eta2 <- function(eta) {
            2 / (eta^3)
          }
          return(x)
        }
        if(link0 == "1/mu^2") {
          x$mu.eta2 <- function(eta) {
            0.75 / eta^(2.5)
          }
          return(x)
        }
        if(link0 == "sqrt") {
          x$mu.eta2 <- function(eta) { rep(2, length = length(eta)) }
          return(x)
        }
        x$mu.eta2 <- function(eta) rep.int(0, length(eta))
        ## warning(paste('higher derivatives of link "', link, '" not available!', sep = ''))
        return(x)
      }

      if(link %in% c("logit", "probit", "cauchit", "cloglog", "identity",
                     "log", "sqrt", "1/mu^2", "inverse")) {
        rval <- make.link(link)
      } else {
        rval <- switch(link,
          "rhogit" = list(
            "linkfun" = function(mu) { mu / sqrt(1 - mu^2) },
            "linkinv" = function(eta) {
                rval <- eta / sqrt(1 + eta^2)
                rval <- (abs(rval) - .Machine$double.eps) * sign(rval)
                rval
            },
            "mu.eta" = function(eta) { 1 / (1 + eta^2)^1.5 }
          ),
          "cloglog2" = list(
            "linkfun" = function(mu) { log(-log(mu)) },
            "linkinv" = function(eta) {
              pmax(pmin(1 - expm1(-exp(eta)), .Machine$double.eps), .Machine$double.eps)
            },
            "mu.eta" = function(eta) {
              eta <- pmin(eta, 700)
              pmax(-exp(eta) * exp(-exp(eta)), .Machine$double.eps)
            }
          ),
          "sigmoid" = list(
            "linkfun" = function(mu) {
              i <- mu <= -1
              if(any(i))
                mu[i] <- mu[i] <- -0.9999
              i <- mu >= 1
              if(any(i))
                mu[i] <- mu[i] <- 0.9999 
              -log(2/(mu + 1) - 1)
            },
            "linkinv" = function(eta) {
              tanh(eta/2)
            },
            "mu.eta" = function(eta) {
              0.5 / cosh(eta * 0.5)^2
            },
            "mu.eta2" = function(eta) {
              eta2 <- eta * 0.5
              -(0.5 * (2 * (sinh(eta2) * 0.5 * cosh(eta2)))/(cosh(eta2)^2)^2)
            }
          )
        )
      }

      rval <- mu.eta2(rval)
    }
    rval$name <- link
  }

  rval
}

## Parsing links helper function.
parse_links <- function(links, default.links, ...)
{
  dots <- list(...)
  nl <- names(default.links)
  if(length(dots))
    links <- as.character(dots)
  if(is.null(names(links)))
    names(links) <- rep(nl, length.out = length(links))
  links <- as.list(links)
  for(j in nl) {
    if(is.null(links[[j]]))
      links[[j]] <- default.links[j]
  }
  links <- links[nl]
  links <- as.character(links)
  names(links) <- nl
  links
}

## Function takes a gamlss family and sets it up
## a bit different in order to support more than
## 4 parameter models.
tF <- function(x, ...)
{
  if(is.function(x))
    x <- x()
  if(!inherits(x, "gamlss.family"))
    stop('only "gamlss.family" objects can be transformed!')

  args <- list(...)
  pr <- args$range
  check_range <- function(par) {
    for(j in names(par)) {
      if(!is.null(pr[[j]])) {
        par[[j]][par[[j]] < min(pr[[j]])] <- min(pr[[j]])
        par[[j]][par[[j]] > max(pr[[j]])] <- max(pr[[j]])
      }
    }
    par
  }
  nx <- names(x$parameters)
  score <- hess <- initialize <- list()

  make_call <- function(fun) {
    fn <- deparse(substitute(fun), backtick = TRUE, width.cutoff = 500)
    nf <- names(formals(fun))
    if(length(nf) < 1) {
      call <- paste(fn, "()", sep = "")
    } else {
      call <- paste(fn, "(", if("y" %in% nf) "y," else "", sep = "")
      np <- nx[nx %in% nf]
      if(!length(np)) {
        call <- paste(fn, "(y", sep = "")
      } else {
        call <- paste(call, paste(np, '=', 'par$', np, sep = '', collapse = ','), sep = "")
      }
      if("bd" %in% nf) {
        call <- paste(call, ",bd", sep = "")
      }
    }
    call <- parse(text = paste(call, ")", sep = ""))
    return(call)
  }

  if("mu" %in% nx) {
    mu.link <- make.link2(x$mu.link)
    mu.cs <- make_call(x$dldm)
    mu.hs <- make_call(x$d2ldm2)
    score$mu  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(mu.cs) * mu.link$mu.eta(mu.link$linkfun(par$mu))
      if(!is.null(dim(res))) {
        if(length(dim(res)) > 1)
          res <- res[, 1]
      }
      res
    }
    hess$mu <- function(y, par, ...) {
      par <- check_range(par)
      ## score <- eval(mu.cs)
      hess <- -1 * eval(mu.hs)
      eta <- mu.link$linkfun(par$mu)
      res <- drop(hess * mu.link$mu.eta(eta)^2)
      if(!is.null(dim(res))) {
        if(length(dim(res)) > 1)
          res <- res[, 1]
      }
      res
    }
    if(!is.null(x$mu.initial)) {
      initialize$mu <- function(y, ...) {
        if(!is.null(attr(y, "contrasts"))) {
          if(!is.null(dim(y)))
            y <- y[, ncol(y)]
        }
        res <- eval(x$mu.initial)
        if(!is.null(dim(res))) {
          if(length(dim(res)) > 1)
            res <- res[, 1]
        }
        res
      }
    }
  }

  if("sigma" %in% nx) {
    sigma.link <- make.link2(x$sigma.link)
    sigma.cs <- make_call(x$dldd)
    sigma.hs <- make_call(x$d2ldd2)
    score$sigma  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(sigma.cs) * sigma.link$mu.eta(sigma.link$linkfun(par$sigma))
      if(!is.null(dim(res))) {
        if(length(dim(res)) > 1)
          res <- res[, 1]
      }
      res
    }
    hess$sigma <- function(y, par, ...) {
      par <- check_range(par)
      ## score <- eval(sigma.cs)
      hess <- -1 * eval(sigma.hs)
      eta <- sigma.link$linkfun(par$sigma)
      ## res <- drop(score * sigma.link$mu.eta2(eta) + hess * sigma.link$mu.eta(eta)^2)
      res <- drop(hess * sigma.link$mu.eta(eta)^2)
      if(!is.null(dim(res))) {
        if(length(dim(res)) > 1)
          res <- res[, 1]
      }
      res
    }
    if(!is.null(x$sigma.initial)) {
      initialize$sigma <- function(y, ...) {
        res <- eval(x$sigma.initial)
        if(!is.null(dim(res))) {
          if(length(dim(res)) > 1)
            res <- res[, 1]
        }
        res
      }
    }
  }

  if("nu" %in% nx) {
    nu.link <- make.link2(x$nu.link)
    nu.cs <- make_call(x$dldv)
    nu.hs <- make_call(x$d2ldv2)
    score$nu  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(nu.cs) * nu.link$mu.eta(nu.link$linkfun(par$nu))
      if(!is.null(dim(res))) {
        if(length(dim(res)) > 1)
          res <- res[, 1]
      }
      res
    }
    hess$nu <- function(y, par, ...) {
      par <- check_range(par)
      ## score <- eval(nu.cs)
      hess <- -1 * eval(nu.hs)
      eta <- nu.link$linkfun(par$nu)
      ## res <- drop(score * nu.link$mu.eta2(eta) + hess * nu.link$mu.eta(eta)^2)
      res <- drop(hess * nu.link$mu.eta(eta)^2)
      if(!is.null(dim(res))) {
        if(length(dim(res)) > 1)
          res <- res[, 1]
      }
      res
    }
    if(!is.null(x$nu.initial)) {
      initialize$nu <- function(y, ...) {
        res <- eval(x$nu.initial)
        if(!is.null(dim(res))) {
          if(length(dim(res)) > 1)
            res <- res[, 1]
        }
        res
      }
    }
  }

  if("tau" %in% nx) {
    tau.link <- make.link2(x$tau.link)
    tau.cs <- make_call(x$dldt)
    tau.hs <- make_call(x$d2ldt2)
    score$tau  <- function(y, par, ...) {
      par <- check_range(par)
      res <- eval(tau.cs) * tau.link$mu.eta(tau.link$linkfun(par$tau))
      if(!is.null(dim(res))) {
        if(length(dim(res)) > 1)
          res <- res[, 1]
      }
      res
    }
    hess$tau <- function(y, par, ...) {
      par <- check_range(par)
      ## score <- eval(tau.cs)
      hess <- -1 * eval(tau.hs)
      eta <- tau.link$linkfun(par$tau)
      ## res <- drop(score * tau.link$mu.eta2(eta) + hess * tau.link$mu.eta(eta)^2)
      res <- drop(hess * tau.link$mu.eta(eta)^2)
      if(!is.null(dim(res))) {
        if(length(dim(res)) > 1)
          res <- res[, 1]
      }
      res
    }
    if(!is.null(x$tau.initial)) {
      initialize$tau <- function(y, ...) {
        res <- eval(x$tau.initial)
        if(!is.null(dim(res))) {
          if(length(dim(res)) > 1)
            res <- res[, 1]
        }
        res
      }
    }
  }

  ## For CG algorithm.
  if(all(c("mu", "sigma") %in% nx)) {
    mu.sigma.hs <- make_call(x$d2ldmdd)
    hess$mu.sigma <- function(y, par, ...) {
      par <- check_range(par)
      eta.mu <- mu.link$linkfun(par$mu)
      eta.sigma <- sigma.link$linkfun(par$sigma)
      hess <- -1 * eval(mu.sigma.hs)
      hess * mu.link$mu.eta(eta.mu) * sigma.link$mu.eta(eta.sigma)
    }
    hess$sigma.mu <- hess$mu.sigma
  }
  if("nu" %in% nx) {
    mu.nu.hs <- make_call(x$d2ldmdv)
    hess$mu.nu <- function(y, par, ...) {
      par <- check_range(par)
      eta.mu <- mu.link$linkfun(par$mu)
      eta.nu <- nu.link$linkfun(par$nu)
      hess <- -1 * eval(mu.nu.hs)
      hess * mu.link$mu.eta(eta.mu) * nu.link$mu.eta(eta.nu)
    }
    hess$nu.mu <- hess$mu.nu

    sigma.nu.hs <- make_call(x$d2ldddv)
    hess$sigma.nu <- function(y, par, ...) {
      par <- check_range(par)
      eta.sigma <- sigma.link$linkfun(par$sigma)
      eta.nu <- nu.link$linkfun(par$nu)
      hess <- -1 * eval(sigma.nu.hs)
      hess * sigma.link$mu.eta(eta.sigma) * nu.link$mu.eta(eta.nu)
    }
    hess$nu.sigma <- hess$sigma.nu
  }

  if("tau" %in% nx) {
    mu.tau.hs <- make_call(x$d2ldmdt)
    hess$mu.tau <- function(y, par, ...) {
      par <- check_range(par)
      eta.mu <- mu.link$linkfun(par$mu)
      eta.tau <- tau.link$linkfun(par$tau)
      hess <- -1 * eval(mu.tau.hs)
      hess * mu.link$mu.eta(eta.mu) * tau.link$mu.eta(eta.tau)
    }
    hess$tau.mu <- hess$mu.tau

    sigma.tau.hs <- make_call(x$d2ldddt)
    hess$sigma.tau <- function(y, par, ...) {
      par <- check_range(par)
      eta.sigma <- sigma.link$linkfun(par$sigma)
      eta.tau <- tau.link$linkfun(par$tau)
      hess <- -1 * eval(sigma.tau.hs)
      hess * sigma.link$mu.eta(eta.sigma) * tau.link$mu.eta(eta.tau)
    }
    hess$tau.sigma <- hess$sigma.tau

    nu.tau.hs <- make_call(x$d2ldvdt)
    hess$nu.tau <- function(y, par, ...) {
      par <- check_range(par)
      eta.nu <- nu.link$linkfun(par$nu)
      eta.tau <- tau.link$linkfun(par$tau)
      hess <- -1 * eval(nu.tau.hs)
      hess * nu.link$mu.eta(eta.nu) * tau.link$mu.eta(eta.tau)
    }
    hess$tau.nu <- hess$nu.tau
  }

  dfun <- get(paste("d", x$family[1], sep = ""))
  pfun <- try(get(paste("p", x$family[1], sep = "")), silent = TRUE)
  qfun <- try(get(paste("q", x$family[1], sep = "")), silent = TRUE)
  rfun <- try(get(paste("r", x$family[1], sep = "")), silent = TRUE)

  nf <- names(formals(dfun))
  bdc <- "bd" %in% nf

  dc <- parse(text = paste('dfun(y,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',log=log,...',
    if(bdc) paste0(",bd") else NULL, ")", sep = ""))
  pc <- parse(text = paste('pfun(q,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',log=log,...',
    if(bdc) paste0(",bd") else NULL, ")", sep = ""))
  qc <- parse(text = paste('qfun(p,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',log=log,...',
    if(bdc) paste0(",bd") else NULL, ")", sep = ""))
  rc <- parse(text = paste('rfun(n,', paste(paste(nx, 'par$', sep = "="),
    nx, sep = '', collapse = ','), ',...',
    if(bdc) paste0(",bd") else NULL, ")", sep = ""))

  rval <- list(
    "family" = x$family[1],
    "names" = nx,
    "links" = unlist(x[paste(nx, "link", sep = ".")]),
    "score" = score,
    "hess" = hess,
    "d" = function(y, par, log = FALSE, ...) {
       par <- check_range(par)
       d <- eval(dc)
       if(log)
         d[!is.finite(d)] <- -100
       return(d)
    },
    "p" = if(!inherits(pfun, "try-error")) function(q, par, log = FALSE, ...) {
      par <- check_range(par)
      p <- eval(pc)
      if(length(p) < length(par[[1L]])) {
        q <- rep(q, length.out = length(par[[1L]]))
        p <- eval(pc)
      }
      return(p)
    } else NULL,
    "q" = if(!inherits(qfun, "try-error")) function(p, par, log = FALSE, ...) {
      par <- check_range(par)
      q <- eval(qc)
      if(length(q) < length(par[[1L]])) {
        p <- rep(p, length.out = length(par[[1L]]))
        q <- eval(qc)
      }
      return(q)
    } else NULL,
    "r" = if(!inherits(rfun, "try-error")) function(n, par, ...) {
      par <- check_range(par)
      return(eval(rc))
    } else NULL
  )
  names(rval$links) <- nx
  rval$valid.response <- x$y.valid
  rval$initialize <- initialize
  rval$type <- tolower(x$type)

  if(!is.null(x$mean)) {
    meanc <- make_call(x$mean)
    rval$mean  <- function(par, ...) {
      par <- check_range(par)
      res <- eval(meanc)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
  }

  if(!is.null(x$variance)) {
    varc <- make_call(x$variance)
    rval$variance  <- function(par, ...) {
      par <- check_range(par)
      res <- eval(varc)
      if(!is.null(dim(res)))
        res <- res[, 1]
      res
    }
  }

  linkinv <- list()
  for(j in rval$names) {
    link <- make.link2(rval$links[[j]])
    linkinv[[j]] <- link$linkinv
  }

  rval$map2par <- function(eta) {
    for(j in names(eta)) {
      eta[[j]] <- linkinv[[j]](eta[[j]])
      eta[[j]][is.na(eta[[j]])] <- 0
      if(any(jj <- eta[[j]] == Inf))
        eta[[j]][jj] <- 10
      if(any(jj <- eta[[j]] == -Inf))
        eta[[j]][jj] <- -10
    }
    return(eta)
  }

  rval$loglik <- function(y, par) {
    par <- check_range(par)
    log <- TRUE
    d <- try(eval(dc), silent = TRUE)
    if(inherits(d, "try-error")) {
      warning("problems evaluating the log-density of the model, set log-likelihood to -Inf")
      return(-Inf)
    }
    d[!is.finite(d)] <- -100
    return(sum(d, na.rm = TRUE))
  }

  if(!is.null(x$rqres)) {
    rqres <- utils::getFromNamespace("rqres", "gamlss")
    nenv <- new.env()
    assign("rqres", utils::getFromNamespace("rqres", "gamlss"), envir = nenv)

    rval$rqres <- function(y, par, ...) {
      assign("y", y, envir = nenv)
      for(i in nx)
        assign(i, par[[i]], envir = nenv)
      eval(x$rqres, envir = nenv)
    }
  }

  class(rval) <- "gamlss2.family"
  rval
}

## Complete a family object, e.g.,
## if derivatives are not supplied they
## will be approximated numerically.
complete_family <- function(family)
{
  if(is.character(family)) {
    family <- get(family)
  }
    
  if(is.function(family)) {
    family <- family()
  }

  if(inherits(family, "gamlss.family")) {
    return(tF(family))
  }

  if(is.null(family$family)) {
    family$family <- "No family name supplied!"
  }

  if(is.null(names(family$links)))
    names(family$links) <- family$names

  linkinv <- linkfun <- mu.eta <- list()
  for(j in family$names) {
    link <- make.link2(family$links[[j]])
    linkinv[[j]] <- link$linkinv
    linkfun[[j]] <- link$linkfun
    mu.eta[[j]] <- link$mu.eta
  }

  if(is.null(family$map2par)) {
    family$map2par <- function(eta) {
      for(j in names(eta)) {
        eta[[j]] <- linkinv[[j]](eta[[j]])
        eta[[j]][is.na(eta[[j]])] <- 0
        if(any(jj <- eta[[j]] == Inf))
          eta[[j]][jj] <- 10
        if(any(jj <- eta[[j]] == -Inf))
          eta[[j]][jj] <- -10
      }
      return(eta)
    }
  }

  if(is.null(family$mu.eta)) {
    family$mu.eta <- mu.eta
  }

  if(is.null(family$mu)) {
    family$mu <- function(par) { par[[1]] }
  }

  if(is.null(family$mean)) {
    family$mean <- function(par) { par[[1]] }
  }

  if(is.null(family$loglik)) {
    if(!is.null(family$d)) {
      family$loglik <- function(y, par, ...) {
        logdens <- try(family$d(y, par, log = TRUE), silent = TRUE)
        if(inherits(logdens, "try-error")) {
          warning("problems evaluating the log-density of the model, set log-likelihood to -Inf")
          return(-Inf)
        }
        if(any(i <- !is.finite(logdens))) {
          logdens[i] <- -100
        }
        return(sum(logdens, na.rm = TRUE))
      }
    } else {
      stop("the family object does not have a $d() function!")
    }
  }

  err01 <- .Machine$double.eps^(1/3)
  err02 <- err01 * 2
  err11 <- .Machine$double.eps^(1/4)
  err12 <- err11 * 2

  if(is.null(family$score) & !is.null(family$d))
    family$score <- list()
  for(i in family$names) {
    if(is.null(family$score[[i]]) & !is.null(family$d)) {
      fun <- c(
        "function(y, par, ...) {",
        paste("  eta <- linkfun[['", i, "']](par[['", i, "']]);", sep = ""),
        paste("  par[['", i, "']] <- linkinv[['", i, "']](eta + err01);", sep = ""),
        "  d1 <- family$d(y, par, log = TRUE);",
        paste("  par[['", i, "']] <- linkinv[['", i, "']](eta - err01);", sep = ""),
        "  d2 <- family$d(y, par, log = TRUE);",
        "  return((d1 - d2) / err02)",
        "}"
      )
      family$score[[i]] <- eval(parse(text = paste(fun, collapse = "")))
      attr(family$score[[i]], "dnum") <- TRUE
    }
  }

  if(is.null(family$hess) & !is.null(family$d))
    family$hess <- list()
  for(i in family$names) {
    if(is.null(family$hess[[i]]) & !is.null(family$d)) {
      fun <- if(!is.null(attr(family$score[[i]], "dnum"))) {
        c(
          "function(y, par, ...) {",
          paste("  eta <- linkfun[['", i, "']](par[['", i, "']]);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta + err11);", sep = ""),
          paste("  d1 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta - err11);", sep = ""),
          paste("  d2 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          "  return(-1 * (d1 - d2) / err12)",
          "}"
        )
      } else {
        c(
          "function(y, par, ...) {",
          paste("  eta <- linkfun[['", i, "']](par[['", i, "']]);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta + err01);", sep = ""),
          paste("  d1 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          paste("  par[['", i, "']] <- linkinv[['", i, "']](eta - err01);", sep = ""),
          paste("  d2 <- family$score[['", i, "']](y, par, ...);", sep = ""),
          "  return(-1 * (d1 - d2) / err02)",
          "}"
        )
      }

#      fun <- c(
#        "function(y, par, ...) {",
#        paste0("sc <- family$score[['", i, "']](y, par, ...); return(sc^2)"),
#        "}"
#      )

      family$hess[[i]] <- eval(parse(text = paste(fun, collapse = "")))
    }
  }
  for(i in seq_along(family$names)) {
    for(j in seq_along(family$names)) {
      if(i < j) {
        hij <- paste0(family$names[i], ".", family$names[j])
        if(is.null(family$hess[[hij]])) {
          ni <- family$names[i]
          nj <- family$names[j]

          fun <- c(
            "function(y, par, ...) {",
            paste("  eta <- linkfun[['", ni, "']](par[['", ni, "']]);", sep = ""),
            paste("  par[['", ni, "']] <- linkinv[['", ni, "']](eta + err01);", sep = ""),
            paste("  d1 <- family$score[['", nj, "']](y, par, ...);", sep = ""),
            paste("  par[['", ni, "']] <- linkinv[['", ni, "']](eta - err01);", sep = ""),
            paste("  d2 <- family$score[['", nj, "']](y, par, ...);", sep = ""),
            "  return(-1 * (d1 - d2) / err02)",
            "}"
          )

#          fun <- c(
#            "function(y, par, ...) {",
#            paste0("family$score[['", ni, "']](y, par, ...)*family$score[['", nj, "']](y, par, ...)"),
#            "}"
#          )

          family$hess[[hij]] <- eval(parse(text = paste(fun, collapse = "")))
          hji <- paste0(family$names[j], ".", family$names[i])
          family$hess[[hji]] <- family$hess[[hij]]
        }
      }
    }
  }

  return(family)
}

## A simple print method.
print.gamlss2.family <- function(x, full = TRUE, ...)
{
  cat("Family:", x$family, if(!is.null(x$full.name)) paste0("(", x$full.name, ")") else NULL,  "\n")
  links <- paste(names(x$links), x$links, sep = " = ")
  links <- paste(links, collapse = ", ")
  cat(if(length(links) > 1) "Link functions:" else "Link function:", links, sep = " ")
  cat("\n")
  if(full) {
    nfun <- names(x[c("transform", "optimizer", "sampler", "results", "predict")])
    if(!all(is.na(nfun))) {
      nfun <- nfun[!is.na(nfun)]
      cat("---\nFamily specific functions:\n")
      for(j in nfun)
        cat(" ..$ ", j, "\n", sep = "")
    }
    nfun <- names(x[c("score", "hess")])
    if(!all(is.na(nfun))) {
      nfun <- nfun[!is.na(nfun)]
      cat("---\nDerivative functions:\n")
      for(j in nfun) {
        cat(" ..$ ", j, "\n", sep = "")
        for(i in names(x[[j]]))
          cat(" .. ..$ ", i, "\n", sep = "")
      }
    }
  }
}

## Some example families.
Gaussian <- function(...)
{
  links <- c(mu = "identity", sigma = "log")

  rval <- list(
    "family" = "gaussian",
    "names" = c("mu", "sigma"),
    "links" = parse_links(links, c(mu = "identity", sigma = "log"), ...),
    "score" = list(
      "mu" = function(y, par, ...) { drop((y - par$mu) / (par$sigma^2)) },
      "sigma" = function(y, par, ...) { drop(-1 + (y - par$mu)^2 / (par$sigma^2)) }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) { drop(1 / (par$sigma^2)) },
      "sigma" = function(y, par, ...) { rep(2, length(y)) }
    ),
    "loglik" = function(y, par, ...) {
      sum(dnorm(y, par$mu, par$sigma, log = TRUE))
    },
    "mu" = function(par, ...) {
      par$mu
    },
    "d" = function(y, par, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "p" = function(y, par, ...) {
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    },
    "r" = function(n, par) {
      rnorm(n, mean = par$mu, sd = par$sigma)
    },
    "q" = function(p, par) {
      qnorm(p, mean = par$mu, sd = par$sigma)
    },
    "crps" = function(y, par, ...) {
      sum(scoringRules::crps_norm(y, mean = par$mu, sd = par$sigma), na.rm = TRUE)
    },
    "initialize" = list(
      "mu"    = function(y, ...) { (y + mean(y)) / 2 },
      "sigma" = function(y, ...) { rep(sd(y), length(y)) }
    ),
    "mean"      = function(par) par$mu,
    "variance"  = function(par) par$sigma^2,
    "valid.response" = function(x) {
      if(is.factor(x) | is.character(x))
        stop("the response should be numeric!")
      return(TRUE)
    }
  )

  class(rval) <- "gamlss2.family"
  rval
}

Weibull <- function(...)
{
  rval <- list(
    "family" = "Weibull",
    "names" = c("mu", "sigma"),
    "links" = c(mu = "identity", sigma = "log"),
    "d" = function(y, par, log = FALSE, ...) {
      delta <- y[, "status"]
      y <- log(y[, "time"])
      yms <- (y - par$mu) / par$sigma
      fy <- delta * (yms - par$sigma - exp(yms))
      Sy <- (1 - delta) * -exp(yms)
      d <- fy + Sy
      if(!log)
        d <- exp(d)
      return(d)
    },
    "p" = function(y, par, ...) {
      delta <- y[, "status"]
      y <- log(y[, "time"])
      p1 <- 1 - exp(-exp((y - par$mu) / par$sigma))
      p2 <- runif(length(y), p1, 1)
      prob <- ifelse(delta > 0, p1, p2)
      return(prob)
    },
    "q" = function(p, par, ...) {
      lambda <- exp(-par$mu/par$sigma)
      alpha <- 1 / par$sigma
      q <- lambda * (-log(1 - p))^(1 / alpha)
      return(q)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        delta <- y[, "status"]
        y <- log(y[, "time"])
        eyms <- exp((y - par$mu)/par$sigma)
        s1 <- 1/par$sigma
        eymss1 <- eyms * s1
        a <- -(delta * (s1 - eymss1))
        b <- (1 - delta) * (eymss1)
        return(a + b)
      },
      "sigma" = function(y, par, ...) {
        delta <- y[, "status"]
        y <- log(y[, "time"])
        yms <- (y - par$mu)/par$sigma
        eyms <- exp(yms)
        eyms2 <-  eyms * yms
        a <- -(delta * (yms + par$sigma - eyms2))
        b <- (1 - delta) * eyms2
        return(a + b)
      }
    ),
    "hess" <- list(
      "mu" = function(y, par, ...) {
        delta <- y[, "status"]
        y <- log(y[, "time"])
        eyms <- exp((y - par$mu)/par$sigma) * 1/par$sigma^2
        a <- -(delta * eyms)
        b <- -((1 - delta) * eyms)
        return(-(a + b))
      },
      "sigma" = function(y, par, ...) {
        delta <- y[, "status"]
        y <- log(y[, "time"])
        yms <- (y - par$mu)/par$sigma
        eyms <- exp(yms)
        a <- -(delta * (-yms + par$sigma - (eyms * (-yms - eyms * yms^2))))
        b <- (1 - delta) * (eyms * -yms - eyms * yms^2)
        return(-(a + b))
      }
    ),
    "valid.response" = function(x) {
      if(!inherits(x, "Surv"))
        stop("the response should be a survival object!")
      return(TRUE)
    }
  )

  class(rval) <- "gamlss2.family"

  rval
}

## From VGAM.
is.Numeric <- function (x, length.arg = Inf, integer.valued = FALSE, positive = FALSE) {
  if (all(is.numeric(x)) && all(is.finite(x)) && (if (is.finite(length.arg)) length(x) == 
    length.arg else TRUE) && (if (integer.valued) all(x == round(x)) else TRUE) && 
    (if (positive) all(x > 0) else TRUE)) TRUE else FALSE
}

## Yeo-Johnson transform family. From VGAM.
YJt <- function(y, lambda = 1, derivative = 0,
  epsilon = sqrt(.Machine$double.eps), inverse = FALSE) {
   if(!is.Numeric(derivative, length.arg = 1, integer.valued = TRUE) || 
        derivative < 0) 
        stop("argument 'derivative' must be a non-negative integer")
    ans <- y
    if(!is.Numeric(epsilon, length.arg = 1, positive = TRUE)) 
        stop("argument 'epsilon' must be a single positive number")
    L <- max(length(lambda), length(y))
    if(length(y) != L) 
        y <- rep_len(y, L)
    if(length(lambda) != L) 
        lambda <- rep_len(lambda, L)
    if(inverse) {
        if(derivative != 0) 
            stop("argument 'derivative' must 0 when inverse = TRUE")
        if(any(index <- y >= 0 & abs(lambda) > epsilon)) 
            ans[index] <- (y[index] * lambda[index] + 1)^(1/lambda[index]) - 
                1
        if(any(index <- y >= 0 & abs(lambda) <= epsilon)) 
            ans[index] <- expm1(y[index])
        if(any(index <- y < 0 & abs(lambda - 2) > epsilon)) 
            ans[index] <- 1 - (-(2 - lambda[index]) * y[index] + 
                1)^(1/(2 - lambda[index]))
        if(any(index <- y < 0 & abs(lambda - 2) <= epsilon)) 
            ans[index] <- -expm1(-y[index])
        return(ans)
    }
    if(derivative == 0) {
        if(any(index <- y >= 0 & abs(lambda) > epsilon)) 
            ans[index] <- ((y[index] + 1)^(lambda[index]) - 1)/lambda[index]
        if(any(index <- y >= 0 & abs(lambda) <= epsilon)) 
            ans[index] <- log1p(y[index])
        if(any(index <- y < 0 & abs(lambda - 2) > epsilon)) 
            ans[index] <- -((-y[index] + 1)^(2 - lambda[index]) - 
                1)/(2 - lambda[index])
        if(any(index <- y < 0 & abs(lambda - 2) <= epsilon)) 
            ans[index] <- -log1p(-y[index])
    } else {
        psi <- Recall(y = y, lambda = lambda, derivative = derivative - 
            1, epsilon = epsilon, inverse = inverse)
        if(any(index <- y >= 0 & abs(lambda) > epsilon)) 
            ans[index] <- ((y[index] + 1)^(lambda[index]) * (log1p(y[index]))^(derivative) - 
                derivative * psi[index])/lambda[index]
        if(any(index <- y >= 0 & abs(lambda) <= epsilon)) 
            ans[index] <- (log1p(y[index]))^(derivative + 1)/(derivative + 
                1)
        if(any(index <- y < 0 & abs(lambda - 2) > epsilon)) 
            ans[index] <- -((-y[index] + 1)^(2 - lambda[index]) * 
                (-log1p(-y[index]))^(derivative) - derivative * 
                psi[index])/(2 - lambda[index])
        if(any(index <- y < 0 & abs(lambda - 2) <= epsilon)) 
            ans[index] <- (-log1p(-y[index]))^(derivative + 1)/(derivative + 
                1)
    }
    ans
}

YJ <- function(...) {
  fam <- list(
    "family" = "Yeo-Johnson",
    "names" = c("mu", "sigma", "lambda"),
    "links" = c(mu = "identity", sigma = "log", lambda = "identity"),
    "d" = function(y, par, log = FALSE, ...) {
      psi <- YJt(y, par$lambda)
      d <- -0.918938533204675 - log(par$sigma) - 0.5 * ((psi - par$mu)/par$sigma)^2 +
        (par$lambda - 1) * sign(y) * log1p(abs(y))
      if(!log)
        d <- exp(d)
      return(d)
    },
    "p" = function(y, par) {
      psi <- YJt(y, par$lambda)
      pnorm(psi, mean = par$mu, sd = par$sigma)
    },
    "q" = function(p, par) {
      q <- qnorm(p, mean = par$mu, sd = par$sigma)
      YJt(q, par$lambda, inverse = TRUE)
    },
    "score" = list(
      "mu" = function(y, par, ...) {
        psi <- YJt(y, par$lambda)
        (psi - par$mu)/(par$sigma^2)
      },
      "sigma" = function(y, par, ...) {
        psi <- YJt(y, par$lambda)
        -1/par$sigma + (psi - par$mu)^2/(par$sigma^3)
      },
      "lambda" = function(y, par, ...) {
        psi <- YJt(y, par$lambda)
        -(psi - par$mu) / (par$sigma^2) * YJt(y, par$lambda, derivative = 1) + sign(y) * log1p(abs(y))
      }
    ),
    "hess" = list(
      "mu" = function(y, par, ...) {
        1 / par$sigma^2
      },
      "sigma" = function(y, par, ...) {
        2 / (par$sigma^2)
      }
    ),
    "mean" = function(par) {
      YJt(par$mu, par$lambda, inverse = TRUE)
    }
  )
  class(fam) <- "gamlss2.family"
  return(fam)
}

## For binomial families.
.bi.list <- c("BI", "Binomial", "BB", "Beta Binomial", "ZIBI", "ZIBB", 
  "ZABI", "ZABB", "DBI", "BItr", "BBtr",  "ZIBItr", "ZIBBtr", 
  "ZABItr", "ZABBtr", "DBItr")

get_y_bd <- function(Y) {
  if(is.null(Y))
    return(list(y = 1, bd = 1))
  if(NCOL(Y) == 1) {
    y <- if(is.factor(Y))  Y != levels(Y)[1] else Y
    bd <- if(is.null(dim(Y))) rep(1, length(Y)) else rep(1, nrow(Y))
    if(any(y < 0 | y > 1))
      stop("y values must be 0 <= y <= 1")
  } else if(NCOL(Y) == 2) {
    if(any(abs(Y - round(Y)) > 0.001)) {
      warning("non-integer counts in a binomial GAMLSS!")
    }
    bd <- Y[,1] + Y[,2]
    y <-  Y[,1]
    if (any(y < 0 | y > bd)) stop("y values must be 0 <= y <= N") # MS Monday, October 17, 2005 
  } else {
    stop(paste("For the binomial family, Y must be", 
      "a vector of 0 and 1's or a 2 column", "matrix where col 1 is no. successes", 
      "and col 2 is no. failures"))
  }
  return(data.frame(y = y, bd = bd))
}

## softplus link object.
softplus <- function(a = 1) {
  link <- list(
    linkfun = function(mu) {
      eta <- mu + log(1 - exp(-abs(a * mu)))/a
      eta[mu < log(2)/a] <- log(expm1(a * mu[mu < log(2)/a]))/a
      return(eta)
    },
    linkinv = function(eta) pmax(0, eta) + log1p(exp(-abs(a * eta)))/a,
    mu.eta = function(eta) 1/(1 + exp(-a * eta)),
    dmu.eta = function(eta) a * exp(-a * eta)/(1 + exp(-a * eta)),
    valideta = function(eta) TRUE,
    name = sprintf("softplus(%s)", format(a, digits = 3))
  )
  class(link) <- "link-glm"
  return(link)
}

ologit4 <- function(...) {
  fam <- list(
    "family" = "Ordered Logit",
    "names" = c("mu", "r1", "r2", "r3"),
    "links" = c(mu = "identity", r1 = "identity", r2 = "identity", r3 = "identity"),
    "d" = function(y, par, log = FALSE, ...) {
      e1 <- exp(par$mu - par$r1) / (1 + exp(par$mu - par$r1))
      e2 <- exp(par$mu - par$r2) / (1 + exp(par$mu - par$r2))
      e3 <- exp(par$mu - par$r3) / (1 + exp(par$mu - par$r3))

      p1 <- 1 - e1
      p2 <- e1 - e2
      p3 <- e2 - e3
      p4 <- e3

      d <- rep(NA, length(y))

      d[y == 1L] <- p1[y == 1L]
      d[y == 2L] <- p2[y == 2L]
      d[y == 3L] <- p3[y == 3L]
      d[y == 4L] <- p4[y == 4L]

      d[d < 1e-08 | is.na(d)] <- 1e-08

      if(log) {
        d <- log(d)
        d[is.na(d)] <- -1e+10
      }

      return(d)
    }
  )
  class(fam) <- c("gamlss2.family", "family.bamlss")
  return(fam)
}

##if(FALSE) {
#library("bamlss")
#library("gamlss2")

#set.seed(123)

#n <- 2000
#x <- runif(n, -3, 3)
#y <- sin(x) + rnorm(n, sd = 0.3)
#yf <- cut(y, breaks = seq(min(y), max(y), length = 5), include.lowest = TRUE)
#yi <- as.integer(yf)

#b <- gamlss2(yi ~ s(x), family = ologit4, step = 0.1, maxit = c(400, 400))
#plot(b)

#m <- bamlss(yi ~ s(x), family = ologit4)
#plot(m)
##}

## Shifted log-link.
log1 <- function(shift = 1) {
  linkfun <- function(mu) log(mu - shift)
  linkinv <- function(eta) exp(eta) + shift
  mu.eta <- function(eta) exp(eta)
  valideta <- function(eta) TRUE
  
  structure(
    list(
      linkfun = linkfun,
      linkinv = linkinv,
      mu.eta = mu.eta,
      valideta = valideta,
      name = paste0("exp(x) +", shift)
    ),
    class = "link-glm"
  )
}

## Kumaraswamy distriubtion.
Kumaraswamy <- function(a.link = log1, b.link = log1, ...) {
  lfa <- make.link2(a.link)
  lfb <- make.link2(b.link)

  fam <- list(
    "family" = "Kumaraswamy",
    "names" = c("a", "b"),
    "links" = c("a" = a.link, "b" = b.link),
    "d" = function(y, par, log = FALSE, ...) {
      d <- log(par$a) + log(par$b) + (par$a - 1) * log(y) + (par$b - 1) * log(1 - y^(par$a))
      if(!log)
        d <- exp(d)
      return(d)
    },
    "score" = list(
      "a" = function(y, par, ...) {
        ly <- log(y)
        ya <- y^par$a
        (1/par$a + ly - (par$b - 1) * (ya * ly/(1 - ya))) * lfa$mu.eta(lfa$linkfun(par$a))
      },
      "b" = function(y, par, ...) {
        (1/par$b + log(1 - y^par$a)) * lfb$mu.eta(lfb$linkfun(par$b))
      }
    ),
    "hess" = list(
      "a" = function(y, par, ...) {
        ya <- y^par$a
        ly <- log(y)
        y1a <- 1 - ya
        ly2 <- ly^2

        (1/par$a^2 + (par$b - 1) * (ya * ly2/(y1a) + ya^2 * ly2 /(y1a)^2)) * lfa$mu.eta(lfa$linkfun(par$a))^2
      },
      "b" = function(y, par, ...) {
        1/par$b^2 * lfb$mu.eta(lfb$linkfun(par$b))^2
      }
    ),
    "p" = function(y, par) {
      1 - (1 - y^par$a)^par$b
    },
    "q" = function(p, par) {
      (1 - (1 - p)^par$b)^(1 / par$a)
    },
    "r" = function(n, par) {
      par <- as.data.frame(par)
      rn <- apply(par, 1, function(p2) {
        p <- runif(n)
        (1 - (1 - p)^p2["b"])^(1 / p2["a"])
      })
      if(!is.null(dim(rn)))
        rn <- t(rn)
      return(rn)
    },
    "mean" = function(par) {
      par$b * gamma(1 + 1/par$a) * gamma(par$b) / gamma(1 + 1/par$a + par$b)
    },
    "mode" = function(par) {
      ((par$a - 1) / (par$a * par$b - 1))^(1/par$a)
    },
    "valid.response" = function(x) {
      if(any(x < 0) | any(x > 1))
        stop("the response should be in (0,1)!")
      return(TRUE)
    }
  )

  class(fam) <- "gamlss2.family"
  return(fam)
}


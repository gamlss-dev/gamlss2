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
    if(inherits(link, "link-glm")) {
      class(link) <- c("link-gamlss2", "link-glm")
      return(link)
    }
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

  if(is.null(rval$linkinv) | is.null(rval$linkfun))
    rval <- gamlss.dist::make.link.gamlss(as.character(rval$name))

  class(rval) <- c("link-gamlss2", "link-glm")

  rval
}

## Helper function.
"c.link-gamlss2" <- function(...) {
  return(list(...))
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
.gamlss2_family_cache <- new.env(parent = emptyenv())

tF <- function(x, ...)
{
  if(is.function(x)) x <- x()
  if(!inherits(x, "gamlss.family")) return(x)

  dots <- list(...)
  pr <- dots$range
  nx <- names(x$parameters)[unlist(x$parameters)]
  score <- hessian <- initialize <- list()

  ## Find make.link2() in the package namespace or the search path.
  make_link2 <- if(exists("make.link2", mode = "function", inherits = TRUE)) {
    get("make.link2", mode = "function", inherits = TRUE)
  } else {
    utils::getFromNamespace("make.link2", "gamlss2")
  }

  ## Code-generation helpers.
  qstr <- function(z) encodeString(z, quote = '"')
  pref <- function(z) paste0("par[[", qstr(z), "]]" )

  mkfun <- function(src, bindings = list()) {
    e <- list2env(bindings, parent = baseenv())
    f <- eval(parse(text = src), envir = e)
    ## Compile generated functions once during family setup.
    compiler::cmpfun(f)
  }

  ## Precompute parameter bounds used by the generated functions.
  rnames <- nx[vapply(nx, function(j) !is.null(pr[[j]]), logical(1L))]
  if(length(rnames)) {
    rlo <- vapply(rnames, function(j) min(pr[[j]]), numeric(1L))
    rhi <- vapply(rnames, function(j) max(pr[[j]]), numeric(1L))
  } else {
    rlo <- rhi <- numeric()
  }

  numtxt <- function(z) {
    if(is.na(z)) return("NA_real_")
    if(is.infinite(z)) return(if(z > 0) "Inf" else "-Inf")
    sprintf("%.17g", z)
  }

  range_src <- if(length(rnames)) {
    paste(vapply(seq_along(rnames), function(k) {
      j <- rnames[[k]]
      lo <- numtxt(rlo[[k]])
      hi <- numtxt(rhi[[k]])
      pj <- pref(j)
      paste0(
        "z__ <- ", pj, "\n",
        "ii__ <- z__ < ", lo, "\n",
        "chg__ <- any(ii__)\n",
        "if(chg__) z__[ii__] <- ", lo, "\n",
        "ii__ <- z__ > ", hi, "\n",
        "if(any(ii__)) { z__[ii__] <- ", hi, "; chg__ <- TRUE }\n",
        "if(chg__) ", pj, " <- z__"
      )
    }, character(1L)), collapse = "\n")
  } else ""

  deriv_call <- function(fun, binding = ".fun") {
    nf <- names(formals(fun))
    if(!length(nf)) return(paste0(binding, "()"))

    aa <- character()
    if("y" %in% nf) aa <- c(aa, "y")

    np <- nx[nx %in% nf]
    if(length(np)) {
      aa <- c(aa, paste0(np, "=", vapply(np, pref, character(1L))))
    } else if(!("y" %in% nf)) {
      ## Keep the fallback used by make_call().
      aa <- c(aa, "y")
    }

    if("bd" %in% nf)
      aa <- c(aa, 'bd=attr(y, "bd", exact=TRUE)')

    paste0(binding, "(", paste(aa, collapse = ","), ")")
  }

  ## Derivative of a parameter with respect to its linear predictor.
  link_objects <- setNames(lapply(nx, function(j) {
    make_link2(x[[paste0(j, ".link")]])
  }), nx)

  link_name <- function(j) {
    z <- link_objects[[j]]$name
    if(is.null(z) || length(z) != 1L) "" else as.character(z)
  }

  mult_expr <- function(j, tag = "") {
    p <- pref(j)
    switch(link_name(j),
      identity = "1",
      log      = p,
      inverse  = paste0("-(", p, " * ", p, ")"),
      sqrt     = paste0("2 * sqrt(", p, ")"),
      logit    = paste0(p, " * (1 - ", p, ")"),
      paste0(".", tag, "mueta(.", tag, "linkfun(", p, "))")
    )
  }

  hot_bindings <- function(fun, j, tag = "") {
    out <- list(.fun = fun)
    ## Bind link functions for the generated expression.
    out[[paste0(".", tag, "linkfun")]] <- link_objects[[j]]$linkfun
    out[[paste0(".", tag, "mueta")]] <- link_objects[[j]]$mu.eta
    out
  }

  make_score <- function(fun, j) {
    dc <- deriv_call(fun)
    mm <- mult_expr(j)
    rr <- if(identical(mm, "1")) dc else paste0("(", dc, ") * (", mm, ")")
    src <- paste0(
      "function(par, y, ...) {\n", range_src, "\n",
      "res <- ", rr, "\n",
      "if(length(dim(res)) > 1L) res <- res[, 1L]\n",
      "res\n}"
    )
    mkfun(src, hot_bindings(fun, j))
  }

  make_hessian <- function(fun, j) {
    dc <- deriv_call(fun)
    mm <- mult_expr(j)
    rr <- if(identical(mm, "1")) {
      paste0("-", dc)
    } else {
      paste0("-", dc, " * (", mm, ") * (", mm, ")")
    }
    src <- paste0(
      "function(par, y, ...) {\n", range_src, "\n",
      "res <- drop(", rr, ")\n",
      "if(length(dim(res)) > 1L) res <- res[, 1L]\n",
      "res\n}"
    )
    mkfun(src, hot_bindings(fun, j))
  }

  make_cross_hessian <- function(fun, j1, j2) {
    dc <- deriv_call(fun)
    m1 <- mult_expr(j1, "a_")
    m2 <- mult_expr(j2, "b_")
    src <- paste0(
      "function(par, y, ...) {\n", range_src, "\n",
      "-", dc, " * (", m1, ") * (", m2, ")\n}"
    )
    b <- list(.fun = fun)
    b$.a_linkfun <- link_objects[[j1]]$linkfun
    b$.a_mueta   <- link_objects[[j1]]$mu.eta
    b$.b_linkfun <- link_objects[[j2]]$linkfun
    b$.b_mueta   <- link_objects[[j2]]$mu.eta
    mkfun(src, b)
  }

  ## Score and Hessian functions.
  if("mu" %in% nx) {
    score$mu   <- make_score(x$dldm, "mu")
    hessian$mu <- make_hessian(x$d2ldm2, "mu")

    if(!is.null(x$mu.initial)) {
      mu.initial <- x$mu.initial
      initialize$mu <- function(y, ...) {
        if(!is.null(attr(y, "contrasts")) && !is.null(dim(y)))
          y <- y[, ncol(y)]
        bd <- attr(y, "bd", exact = TRUE)
        res <- eval(mu.initial)
        if(length(dim(res)) > 1L) res <- res[, 1L]
        res
      }
    }
  }

  if("sigma" %in% nx) {
    score$sigma   <- make_score(x$dldd, "sigma")
    hessian$sigma <- make_hessian(x$d2ldd2, "sigma")

    if(!is.null(x$sigma.initial)) {
      sigma.initial <- x$sigma.initial
      initialize$sigma <- function(y, ...) {
        res <- eval(sigma.initial)
        if(length(dim(res)) > 1L) res <- res[, 1L]
        res
      }
    }
  }

  if("nu" %in% nx) {
    score$nu   <- make_score(x$dldv, "nu")
    hessian$nu <- make_hessian(x$d2ldv2, "nu")

    if(!is.null(x$nu.initial)) {
      nu.initial <- x$nu.initial
      initialize$nu <- function(y, ...) {
        res <- eval(nu.initial)
        if(length(dim(res)) > 1L) res <- res[, 1L]
        res
      }
    }
  }

  if("tau" %in% nx) {
    score$tau   <- make_score(x$dldt, "tau")
    hessian$tau <- make_hessian(x$d2ldt2, "tau")

    if(!is.null(x$tau.initial)) {
      tau.initial <- x$tau.initial
      initialize$tau <- function(y, ...) {
        res <- eval(tau.initial)
        if(length(dim(res)) > 1L) res <- res[, 1L]
        res
      }
    }
  }

  ## Cross derivatives for CG.
  if(all(c("mu", "sigma") %in% nx)) {
    hessian[["mu:sigma"]] <- make_cross_hessian(x$d2ldmdd, "mu", "sigma")
    hessian[["sigma:mu"]] <- hessian[["mu:sigma"]]
  }
  if("nu" %in% nx) {
    hessian[["mu:nu"]] <- make_cross_hessian(x$d2ldmdv, "mu", "nu")
    hessian[["nu:mu"]] <- hessian[["mu:nu"]]

    hessian[["sigma:nu"]] <- make_cross_hessian(x$d2ldddv, "sigma", "nu")
    hessian[["nu:sigma"]] <- hessian[["sigma:nu"]]
  }
  if("tau" %in% nx) {
    hessian[["mu:tau"]] <- make_cross_hessian(x$d2ldmdt, "mu", "tau")
    hessian[["tau:mu"]] <- hessian[["mu:tau"]]

    hessian[["sigma:tau"]] <- make_cross_hessian(x$d2ldddt, "sigma", "tau")
    hessian[["tau:sigma"]] <- hessian[["sigma:tau"]]

    hessian[["nu:tau"]] <- make_cross_hessian(x$d2ldvdt, "nu", "tau")
    hessian[["tau:nu"]] <- hessian[["nu:tau"]]
  }

  ## Working response and weights for RS.
  diag_score_fun <- list()
  diag_hess_fun <- list()
  if("mu" %in% nx) {
    diag_score_fun$mu <- x$dldm
    diag_hess_fun$mu <- x$d2ldm2
  }
  if("sigma" %in% nx) {
    diag_score_fun$sigma <- x$dldd
    diag_hess_fun$sigma <- x$d2ldd2
  }
  if("nu" %in% nx) {
    diag_score_fun$nu <- x$dldv
    diag_hess_fun$nu <- x$d2ldv2
  }
  if("tau" %in% nx) {
    diag_score_fun$tau <- x$dldt
    diag_hess_fun$tau <- x$d2ldt2
  }

  make_update <- function() {
    zb <- list()
    branches <- vapply(seq_along(nx), function(k) {
      j <- nx[[k]]
      sf <- diag_score_fun[[j]]
      hf <- diag_hess_fun[[j]]
      if(is.null(sf) || is.null(hf))
        stop("cannot build fused update() for parameter '", j, "'")

      sb <- paste0(".zscore", k)
      hb <- paste0(".zhess", k)
      sc <- deriv_call(sf, binding = sb)
      hc <- deriv_call(hf, binding = hb)

      tag <- paste0("zw", k, "_")
      mm <- mult_expr(j, tag = tag)
      if(identical(mm, "1")) {
        msrc <- ""
        score_expr <- sc
        hess_expr <- paste0("-", hc)
      } else {
        msrc <- paste0("m__ <- ", mm, "\n")
        score_expr <- paste0("(", sc, ") * m__")
        hess_expr <- paste0("-", hc, " * m__ * m__")
      }

      zb[[sb]] <<- sf
      zb[[hb]] <<- hf
      zb[[paste0(".", tag, "linkfun")]] <<- link_objects[[j]]$linkfun
      zb[[paste0(".", tag, "mueta")]] <<- link_objects[[j]]$mu.eta

      ## Match deriv_checks() used by RS().
      paste0(
        qstr(j), " = {\n",
        msrc,
        "score__ <- ", score_expr, "\n",
        "if(length(dim(score__)) > 1L) score__ <- score__[, 1L]\n",
        "hess__ <- drop(", hess_expr, ")\n",
        "if(length(dim(hess__)) > 1L) hess__ <- hess__[, 1L]\n",
        "score__[is.na(score__)] <- 1.490116e-08\n",
        "score__[score__ > 1e10] <- 1e10\n",
        "score__[score__ < -1e10] <- -1e10\n",
        "hess__[is.na(hess__)] <- 1.490116e-08\n",
        "hess__[hess__ > 1e10] <- 1e10\n",
        "ii__ <- (hess__ == 0) | !is.finite(hess__)\n",
        "hess__[ii__] <- 1.490116e-08\n",
        "ii__ <- hess__ < 0\n",
        "hess__[ii__] <- -hess__[ii__]\n",
        "hess__[hess__ < 1e-10] <- 1e-10\n",
        "list(eta = eta + score__ / hess__, weights = hess__)\n",
        "}"
      )
    }, character(1L))

    src <- paste0(
      "function(par, y, eta, which) {\n",
      ## Clip parameters in the local copy when ranges are supplied.
      range_src, "\n",
      "switch(which,\n", paste(branches, collapse = ",\n"),
      ",\nstop(\"unknown parameter in update(): \", which, call. = FALSE))\n",
      "}"
    )
    mkfun(src, zb)
  }

  update_fun <- make_update()

  ## Distribution functions.
  fam <- x$family[1L]
  dfun <- get(paste0("d", fam), mode = "function", inherits = TRUE)
  pfun <- get0(paste0("p", fam), mode = "function", inherits = TRUE)
  qfun <- get0(paste0("q", fam), mode = "function", inherits = TRUE)
  rfun <- get0(paste0("r", fam), mode = "function", inherits = TRUE)

  dist_call <- function(fun, first, with_log = FALSE, with_dots = TRUE,
                        binding = ".fun") {
    aa <- c(first, paste0(nx, "=", vapply(nx, pref, character(1L))))
    nf <- names(formals(fun))
    if(with_log) {
      log_argument <- if("log" %in% nf) "log" else
        if("log.p" %in% nf) "log.p" else "log"
      aa <- c(aa, paste0(log_argument, "=log"))
    }
    if(with_dots) aa <- c(aa, "...")

    ## Add bd only for functions that use a response argument.
    if("bd" %in% nf && identical(first, "y"))
      aa <- c(aa, 'bd=attr(y, "bd", exact=TRUE)')

    paste0(binding, "(", paste(aa, collapse = ","), ")")
  }

  dcall <- dist_call(dfun, "y", with_log = TRUE, with_dots = TRUE)
  pdf_fun <- mkfun(
    paste0("function(par, y, log=FALSE, ...) ", dcall),
    list(.fun = dfun)
  )

  cdf_fun <- NULL
  if(!is.null(pfun)) {
    pcall <- dist_call(pfun, "y", with_log = TRUE, with_dots = TRUE)
    cdf_fun <- mkfun(paste0(
      "function(par, y, log=FALSE, ...) {\n",
      "p__ <- ", pcall, "\n",
      "n__ <- length(par[[1L]])\n",
      "if(length(p__) < n__) { y <- rep(y, length.out=n__); p__ <- ", pcall, " }\n",
      "p__\n}"
    ), list(.fun = pfun))
  }

  quantile_fun <- NULL
  support_fun <- x[["support"]]
  if(!is.null(qfun)) {
    qcall <- dist_call(qfun, "p", with_log = TRUE, with_dots = TRUE)
    quantile_fun <- mkfun(paste0(
      "function(par, p, log=FALSE, ...) {\n",
      "ii__ <- p <= 1e-10; if(any(ii__)) p[ii__] <- 1e-10\n",
      "ii__ <- p >= (1 - 1e-10); if(any(ii__)) p[ii__] <- 1 - 1e-10\n",
      "q__ <- ", qcall, "\n",
      "n__ <- length(par[[1L]])\n",
      "if(length(q__) < n__) { p <- rep(p, length.out=n__); q__ <- ", qcall, " }\n",
      "q__\n}"
    ), list(.fun = qfun))

    ## The prediction quantile wrapper clips probabilities away from zero and
    ## one. Use the original distribution function for support endpoints so
    ## that unbounded supports remain infinite.
    if(is.null(support_fun)) {
      qcall_support <- dist_call(qfun, "p", with_log = FALSE,
        with_dots = TRUE)
      support_fun <- mkfun(paste0(
        "function(par, ...) {\n",
        "n__ <- length(par[[1L]])\n",
        "p <- rep(0, length.out=n__); lo__ <- ", qcall_support, "\n",
        "p <- rep(1, length.out=n__); hi__ <- ", qcall_support, "\n",
        "cbind(min=lo__, max=hi__)\n}"
      ), list(.fun = qfun))
      attr(support_fun, "inferred") <- "quantile"
    }
  }

  random_fun <- NULL
  if(!is.null(rfun)) {
    rcall <- dist_call(rfun, "n", with_log = FALSE, with_dots = TRUE)
    random_fun <- mkfun(
      paste0("function(par, n, ...) ", rcall),
      list(.fun = rfun)
    )
  }

  rval <- list(
    family = fam,
    names = nx,
    links = unlist(x[paste(nx, "link", sep = ".")]),
    score = score,
    hessian = hessian,
    update = update_fun,
    pdf = pdf_fun,
    cdf = cdf_fun,
    quantile = quantile_fun,
    random = random_fun,
    support = support_fun
  )
  names(rval$links) <- nx
  rval$valid.response <- x$y.valid
  rval$initialize <- initialize
  rval$type <- tolower(x$type)

  ## Moments.
  moment_call <- function(fun) {
    nf <- names(formals(fun))
    if(!length(nf)) return(".fun()")
    np <- nx[nx %in% nf]
    if(!length(np)) return(".fun()")
    paste0(".fun(", paste0(np, "=", vapply(np, pref, character(1L)),
      collapse = ","), ")")
  }

  if(!is.null(x$mean)) {
    mc <- moment_call(x$mean)
    rval$mean <- mkfun(paste0(
      "function(par, ...) { res <- ", mc,
      "; if(length(dim(res)) > 1L) res <- res[,1L]; res }"
    ), list(.fun = x$mean))
  }

  if(!is.null(x$variance)) {
    vc <- moment_call(x$variance)
    rval$variance <- mkfun(paste0(
      "function(par, ...) { res <- ", vc,
      "; if(length(dim(res)) > 1L) res <- res[,1L]; res }"
    ), list(.fun = x$variance))
  }

  ## Map linear predictors to parameters.
  map_blocks <- vapply(seq_along(nx), function(k) {
    j <- nx[[k]]
    ej <- paste0("eta[[", qstr(j), "]]" )
    nm <- link_name(j)
    inv <- switch(nm,
      identity = ej,
      log = paste0("exp(", ej, ")"),
      paste0(".li", k, "(", ej, ")")
    )
    paste0(
      "if(!is.null(", ej, ")) {\n",
      "z__ <- ", inv, "\n",
      "bad__ <- !is.finite(z__)\n",
      "if(any(bad__)) {\n",
      "  ii__ <- is.na(z__); if(any(ii__)) z__[ii__] <- 0\n",
      "  ii__ <- z__ == Inf; if(any(ii__)) z__[ii__] <- 10\n",
      "  ii__ <- z__ == -Inf; if(any(ii__)) z__[ii__] <- -10\n",
      "}\n",
      ej, " <- z__\n",
      "}"
    )
  }, character(1L))

  li_bindings <- setNames(lapply(seq_along(nx), function(k) {
    link_objects[[nx[[k]]]]$linkinv
  }), paste0(".li", seq_along(nx)))

  rval$map2par <- mkfun(paste0(
    "function(eta) {\n", paste(map_blocks, collapse = "\n"), "\neta\n}"
  ), li_bindings)

  ## Log-likelihood.
  dll <- dist_call(dfun, "y", with_log = TRUE, with_dots = FALSE, binding = ".dfun")
  ## Use log = TRUE directly in the generated call.
  dll <- sub("log=log", "log=TRUE", dll, fixed = TRUE)
  rval$log_likelihood <- mkfun(paste0(
    "function(par, y) {\n",
    "d__ <- ", dll, "\n",
    "ii__ <- !is.finite(d__); if(any(ii__)) d__[ii__] <- -100\n",
    "sum(d__, na.rm=TRUE)\n}"
  ), list(.dfun = dfun))

  ## Randomized quantile residuals.
  if(!is.null(x$rqres)) {
    rqres_expr <- x$rqres
    nenv <- new.env(parent = baseenv())
    assign("rqres", utils::getFromNamespace("rqres", "gamlss"), envir = nenv)

    rval$rqres <- function(par, y, ...) {
      assign("y", y, envir = nenv)
      for(i in nx) assign(i, par[[i]], envir = nenv)
      eval(rqres_expr, envir = nenv)
    }
  }

  ## Record the original generated callbacks for RS memoization. Arbitrary
  ## user callbacks are not opted in, and replacing any recorded function
  ## invalidates reuse of that function. These wrappers bind the standard
  ## gamlss.dist density and construct their links from character names.
  if(identical(environment(dfun), asNamespace("gamlss.dist")) &&
      all(vapply(x[paste(nx, "link", sep = ".")], is.character, logical(1L))))
    attr(rval, "rs.cache") <- rval[c("map2par", "pdf", "log_likelihood")]

  rval <- complete_family_support(rval)
  rval <- complete_family_cdf(rval)
  rval <- complete_family_quantile(rval)
  rval <- complete_family_moments(rval)
  rval <- complete_family_random(rval)
  class(rval) <- "gamlss2.family"
  rval
}

## Factory for numerically approximated score functions.
make_numeric_score <- function(parameter, pdf, linkfun, linkinv,
  step = .Machine$double.eps^(1/3))
{
  force(parameter)
  force(pdf)
  force(linkfun)
  force(linkinv)
  force(step)

  function(par, y, ...) {
    eta <- linkfun(par[[parameter]])

    par[[parameter]] <- linkinv(eta + step)
    upper <- pdf(par = par, y = y, log = TRUE)

    par[[parameter]] <- linkinv(eta - step)
    lower <- pdf(par = par, y = y, log = TRUE)

    (upper - lower) / (2 * step)
  }
}

## Factory for numerically approximated negative Hessian functions.
make_numeric_hessian <- function(parameter, score, linkfun, linkinv,
  step = .Machine$double.eps^(1/3))
{
  force(parameter)
  force(score)
  force(linkfun)
  force(linkinv)
  force(step)

  function(par, y, ...) {
    eta <- linkfun(par[[parameter]])

    par[[parameter]] <- linkinv(eta + step)
    upper <- score(par = par, y = y, ...)

    par[[parameter]] <- linkinv(eta - step)
    lower <- score(par = par, y = y, ...)

    -(upper - lower) / (2 * step)
  }
}

## Factory for joint numerical score and negative Hessian updates.
make_numeric_update <- function(pdf, linkinv,
  step = .Machine$double.eps^(1/4))
{
  force(pdf)
  force(linkinv)
  force(step)

  function(par, y, eta, which) {
    par_plus <- par_minus <- par

    par_plus[[which]] <- linkinv[[which]](eta + step)
    upper <- pdf(par = par_plus, y = y, log = TRUE)

    center <- pdf(par = par, y = y, log = TRUE)

    par_minus[[which]] <- linkinv[[which]](eta - step)
    lower <- pdf(par = par_minus, y = y, log = TRUE)

    score <- deriv_checks(
      (upper - lower) / (2 * step),
      is.weight = FALSE
    )
    hessian <- deriv_checks(
      -(upper - 2 * center + lower) / step^2,
      is.weight = TRUE
    )

    list(
      eta = eta + score / hessian,
      weights = hessian
    )
  }
}

## Number of elementwise distributions represented by a parameter object.
family_parameter_rows <- function(par)
{
  if(is.data.frame(par) || is.matrix(par))
    return(nrow(par))
  if(is.list(par)) {
    n <- lengths(par)
    if(!length(n)) return(0L)
    nr <- max(n)
    if(any(!n %in% c(1L, nr)))
      stop("family parameters have incompatible lengths.")
    return(nr)
  }
  1L
}

## Normalize the family-level support contract to a two-column matrix.
normalize_family_support <- function(value, par)
{
  nr <- family_parameter_rows(par)
  if(is.list(value) && !is.data.frame(value)) {
    nms <- names(value)
    lower <- if("min" %in% nms) value[["min"]] else
      if("lower" %in% nms) value[["lower"]] else NULL
    upper <- if("max" %in% nms) value[["max"]] else
      if("upper" %in% nms) value[["upper"]] else NULL
    if(is.null(lower) || is.null(upper))
      stop("family support must contain 'min' and 'max' endpoints.")
    value <- cbind(min = lower, max = upper)
  } else if(is.null(dim(value))) {
    if(!is.numeric(value) || length(value) != 2L)
      stop("family support must be two numeric endpoints or a two-column matrix.")
    value <- matrix(rep(value, each = nr), nrow = nr, ncol = 2L)
  } else {
    value <- as.matrix(value)
  }

  if(!is.numeric(value) || ncol(value) != 2L)
    stop("family support must be a numeric matrix with two columns.")
  if(nrow(value) == 1L && nr != 1L)
    value <- value[rep.int(1L, nr), , drop = FALSE]
  if(nrow(value) != nr)
    stop("family support must provide one row per parameter combination.")
  colnames(value) <- c("min", "max")
  if(anyNA(value))
    stop("family support endpoints must not be missing.")
  if(any(value[, "min"] > value[, "max"]))
    stop("family support lower endpoint exceeds its upper endpoint.")

  rn <- if(is.data.frame(par) || is.matrix(par)) rownames(par) else NULL
  if(length(rn) == nr) rownames(value) <- rn
  value
}

## Infer the interval containing a continuous density's positive mass. Using
## log-densities avoids mistaking ordinary floating-point underflow in the
## tails for a finite endpoint. This remains a heuristic for arbitrary user
## densities; an explicit support specification is always authoritative.
make_density_support <- function(pdf, max_steps = 20L, refine = 64L)
{
  force(pdf)
  force(max_steps)
  force(refine)

  support <- function(par, ...) {
    nr <- family_parameter_rows(par)
    if(!nr) stop("family parameters must not be empty.")
    parameters <- if(is.matrix(par)) as.list(as.data.frame(par)) else
      as.list(par)
    parameters <- lapply(parameters, rep, length.out = nr)
    dots <- list(...)
    bounds <- matrix(NA_real_, nrow = nr, ncol = 2L,
      dimnames = list(NULL, c("min", "max")))

    for(i in seq_len(nr)) {
      pari <- lapply(parameters, function(z) z[i])
      evaluate <- function(y) {
        ans <- tryCatch(
          suppressWarnings(if(!length(dots)) {
            pdf(par = pari, y = y, log = TRUE)
          } else {
            do.call(pdf, c(list(par = pari, y = y, log = TRUE), dots))
          }),
          error = function(e) NA_real_
        )
        if(length(ans) != 1L) return(NA_real_)
        as.numeric(ans)
      }
      inside <- function(y) {
        z <- evaluate(y)
        !is.na(z) && z > -Inf
      }

      ## Include conventional response-scale anchors, parameter values, their
      ## midpoints, and local offsets around the first parameter (usually a
      ## location or mean). These find common fixed and parameter-dependent
      ## boundaries before the outward search is needed.
      pv <- suppressWarnings(as.numeric(unlist(pari, use.names = FALSE)))
      pv <- pv[is.finite(pv)]
      anchors <- unique(c(0, pv))
      midpoints <- numeric()
      if(length(anchors) > 1L) {
        pairs <- utils::combn(anchors, 2L)
        midpoints <- pairs[1L, ] / 2 + pairs[2L, ] / 2
      }
      local <- numeric()
      if(length(pv)) {
        scales <- unique(c(0.5, 1, abs(pv[-1L]),
          abs(pv[-1L] - pv[1L])))
        scales <- scales[is.finite(scales) & scales > 0]
        local <- c(pv[1L] - scales, pv[1L] + scales)
      }
      candidates <- unique(c(-1, -0.5, 0, 0.5, 1, pv,
        midpoints, local))
      candidates <- sort(candidates[is.finite(candidates)])
      log_density <- vapply(candidates, evaluate, numeric(1L))
      in_support <- !is.na(log_density) & log_density > -Inf
      if(!any(in_support))
        stop("could not infer support from the family density at row ", i,
          "; please supply 'family$support'.", call. = FALSE)

      ## Start near the mode among the probed points. Positive infinity is a
      ## valid log-density at an integrable boundary singularity.
      x0 <- candidates[which.max(replace(log_density, !in_support, -Inf))]
      distances <- abs(candidates - x0)
      distances <- distances[is.finite(distances) & distances > 0]
      step <- if(length(distances)) min(distances) else max(1, abs(x0))
      step <- max(step, 16 * .Machine$double.eps * max(1, abs(x0)))

      refine_boundary <- function(xin, xout) {
        for(k in seq_len(refine)) {
          midpoint <- xin / 2 + xout / 2
          if(!is.finite(midpoint) || midpoint == xin || midpoint == xout)
            break
          if(inside(midpoint)) xin <- midpoint else xout <- midpoint
        }
        endpoint <- xin / 2 + xout / 2

        ## Recover exact, meaningful endpoints such as zero or a parameter
        ## value when bisection has converged to one of the probe anchors.
        nearest <- candidates[which.min(abs(candidates - endpoint))]
        tolerance <- 128 * .Machine$double.eps *
          max(1, abs(endpoint), abs(nearest))
        if(abs(nearest - endpoint) <= tolerance) endpoint <- nearest
        endpoint
      }

      find_boundary <- function(direction) {
        side <- if(direction < 0) {
          which(candidates < x0)
        } else {
          which(candidates > x0)
        }
        if(length(side)) {
          side <- side[order(candidates[side], decreasing = direction < 0)]
          terminal <- vapply(seq_along(side), function(k) {
            !in_support[side[k]] && all(!in_support[side[k:length(side)]])
          }, logical(1L))
          if(any(terminal))
            return(refine_boundary(x0,
              candidates[side[which(terminal)[1L]]]))
        }

        ## No nearby zero-density region was found. Expand geometrically; a
        ## density that stays positive over this wide range is treated as
        ## having an infinite endpoint.
        last_inside <- x0
        for(k in 0:(max_steps - 1L)) {
          probe <- x0 + direction * step * 2^k
          if(!is.finite(probe)) return(direction * Inf)
          if(inside(probe)) {
            last_inside <- probe
          } else {
            farther <- x0 + direction * step * 2^(k + 1L)
            if(!is.finite(farther) || !inside(farther))
              return(refine_boundary(last_inside, probe))
            last_inside <- farther
          }
        }
        direction * Inf
      }

      bounds[i, ] <- c(find_boundary(-1), find_boundary(1))

      ## Most distribution families have parameter-independent support. After
      ## inferring the first row, verify its boundary pattern for every row in
      ## a few vectorized log-density evaluations. This avoids repeating the
      ## search for large prediction vectors. If the density is not vectorized
      ## or any pattern differs, continue with row-wise inference.
      if(i == 1L && nr > 1L) {
        scale <- max(1, abs(x0), abs(bounds[i, is.finite(bounds[i, ])]))
        delta <- 128 * .Machine$double.eps * scale
        probes <- expected <- numeric()
        if(is.finite(bounds[i, "min"])) {
          probes <- c(probes, bounds[i, "min"] - delta,
            bounds[i, "min"] + delta)
          expected <- c(expected, 0, 1)
        } else {
          probes <- c(probes, x0 - step * 2^(max_steps - 1L))
          expected <- c(expected, 1)
        }
        if(is.finite(bounds[i, "max"])) {
          probes <- c(probes, bounds[i, "max"] - delta,
            bounds[i, "max"] + delta)
          expected <- c(expected, 1, 0)
        } else {
          probes <- c(probes, x0 + step * 2^(max_steps - 1L))
          expected <- c(expected, 1)
        }

        same_support <- TRUE
        for(k in seq_along(probes)) {
          ans <- tryCatch(
            suppressWarnings(if(!length(dots)) {
              pdf(par = parameters, y = rep.int(probes[k], nr), log = TRUE)
            } else {
              do.call(pdf, c(list(par = parameters,
                y = rep.int(probes[k], nr), log = TRUE), dots))
            }),
            error = function(e) NULL
          )
          if(is.null(ans) || length(ans) != nr) {
            same_support <- FALSE
            break
          }
          observed <- !is.na(ans) & as.numeric(ans) > -Inf
          if(any(observed != as.logical(expected[k]))) {
            same_support <- FALSE
            break
          }
        }
        if(same_support) {
          bounds[-1L, ] <- bounds[rep.int(1L, nr - 1L), , drop = FALSE]
          break
        }
      }
    }
    bounds
  }
  attr(support, "inferred") <- "density"
  support
}

## Preserve explicit support specifications and otherwise infer endpoints from
## an existing quantile function. Without a quantile, use the canonical
## nonnegative-integer support for a count family or probe the log-density for
## a continuous family.
complete_family_support <- function(family)
{
  support <- family[["support"]]
  if(is.function(support) &&
      isTRUE(attr(support, "gamlss2.normalized", exact = TRUE)))
    return(family)
  inferred <- NULL
  if(is.null(support)) {
    quantile <- family[["quantile"]]
    if(is.function(quantile)) {
      support <- function(par, ...) {
        cbind(
          min = quantile(par = par, p = 0),
          max = quantile(par = par, p = 1)
        )
      }
      inferred <- "quantile"
    } else {
      type <- family[["type"]]
      if(is.null(type)) type <- "continuous"
      type <- tolower(type[1L])
      if(identical(type, "discrete")) {
        endpoints <- c(0, Inf)
        support <- function(par, ...) endpoints
        inferred <- "count"
      } else {
        if(!identical(type, "continuous") ||
            !is.function(family[["pdf"]]))
          return(family)
        support <- make_density_support(family[["pdf"]])
        inferred <- "density"
      }
    }
  } else if(!is.function(support)) {
    if(!is.numeric(support) || length(support) != 2L || anyNA(support) ||
        support[1L] > support[2L])
      stop("'family$support' must be a function or two ordered numeric endpoints.")
    endpoints <- unname(support)
    support <- function(par, ...) endpoints
  }

  support0 <- support
  support <- function(par, ...) {
    normalize_family_support(support0(par, ...), par)
  }
  if(is.null(inferred)) inferred <- attr(support0, "inferred", exact = TRUE)
  if(!is.null(inferred)) attr(support, "inferred") <- inferred
  attr(support, "gamlss2.normalized") <- TRUE
  family$support <- support
  family
}

## Numerically integrate a continuous density over its declared support.
make_numeric_cdf <- function(pdf, support)
{
  cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE,
    rel.tol = 1e-8, abs.tol = 0, subdivisions = 100L, ...)
  {
    if(!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail))
      stop("'lower.tail' must be TRUE or FALSE.")
    if(!is.logical(log.p) || length(log.p) != 1L || is.na(log.p))
      stop("'log.p' must be TRUE or FALSE.")
    if(!is.numeric(y)) stop("'y' must be numeric.")
    if(!length(y)) return(numeric())
    if(!is.numeric(rel.tol) || length(rel.tol) != 1L ||
        !is.finite(rel.tol) || rel.tol <= 0)
      stop("'rel.tol' must be a positive finite number.")
    if(!is.numeric(abs.tol) || length(abs.tol) != 1L ||
        !is.finite(abs.tol) || abs.tol < 0)
      stop("'abs.tol' must be a nonnegative finite number.")
    if(!is.numeric(subdivisions) || length(subdivisions) != 1L ||
        !is.finite(subdivisions) || subdivisions < 1 ||
        subdivisions != floor(subdivisions))
      stop("'subdivisions' must be a positive integer.")
    subdivisions <- as.integer(subdivisions)

    nr <- family_parameter_rows(par)
    ny <- length(y)
    if(!nr) stop("family parameters must not be empty.")
    n <- max(nr, ny)
    if(!nr %in% c(1L, n) || !ny %in% c(1L, n))
      stop("lengths of family parameters and 'y' are incompatible.")

    parameters <- if(is.matrix(par)) as.list(as.data.frame(par)) else
      as.list(par)
    parameters <- lapply(parameters, rep, length.out = n)
    yy <- rep(y, length.out = n)
    dots <- list(...)
    bounds <- if(!length(dots)) {
      support(parameters)
    } else {
      do.call(support, c(list(par = parameters), dots))
    }
    value <- rep(NA_real_, n)

    for(i in seq_len(n)) {
      if(is.na(yy[i])) next
      lo <- bounds[i, "min"]
      hi <- bounds[i, "max"]
      if(lower.tail && yy[i] <= lo || !lower.tail && yy[i] >= hi) {
        value[i] <- 0
        next
      }
      if(lower.tail && yy[i] >= hi || !lower.tail && yy[i] <= lo) {
        value[i] <- 1
        next
      }

      pari <- lapply(parameters, function(z) z[i])
      evaluate <- function(z) {
        if(!length(dots)) {
          pdf(par = pari, y = z, log = FALSE)
        } else {
          do.call(pdf, c(list(par = pari, y = z, log = FALSE), dots))
        }
      }

      ## stats::integrate() supplies vector-valued batches of abscissae. Most
      ## family densities are vectorized; retain a scalar fallback for custom
      ## densities that are not.
      vectorized <- TRUE
      integrand <- function(z) {
        density <- NULL
        if(vectorized) {
          density <- try(evaluate(z), silent = TRUE)
          if(inherits(density, "try-error") || length(density) != length(z)) {
            vectorized <<- FALSE
            density <- NULL
          }
        }
        if(is.null(density)) {
          density <- vapply(z, function(zz) {
            ans <- evaluate(zz)
            if(length(ans) != 1L)
              stop("the family density must return one value per response.")
            as.numeric(ans)
          }, numeric(1L))
        }
        density <- as.numeric(density)
        if(anyNA(density) || any(!is.finite(density)))
          stop("the family density returned non-finite values during integration.")
        if(any(density < 0))
          stop("the family density returned negative values during integration.")
        density
      }

      limits <- if(lower.tail) c(lo, yy[i]) else c(yy[i], hi)
      result <- tryCatch(
        stats::integrate(integrand, lower = limits[1L], upper = limits[2L],
          subdivisions = subdivisions, rel.tol = rel.tol,
          abs.tol = abs.tol, stop.on.error = FALSE),
        error = function(e) e
      )
      if(inherits(result, "error"))
        stop("numerical CDF integration failed at row ", i, ": ",
          conditionMessage(result), call. = FALSE)
      if(!identical(result$message, "OK"))
        stop("numerical CDF integration failed at row ", i, ": ",
          result$message, call. = FALSE)

      probability <- result$value
      tolerance <- max(10 * rel.tol, 10 * abs.tol,
        100 * .Machine$double.eps)
      if(!is.finite(probability) || probability < -tolerance ||
          probability > 1 + tolerance)
        stop("numerical CDF integration produced a value outside [0, 1] at row ",
          i, ".", call. = FALSE)
      value[i] <- min(max(probability, 0), 1)
    }

    if(log.p) value <- log(value)
    value
  }
  attr(cdf, "dnum") <- TRUE
  cdf
}

## Numerically sum a count density over its declared integer support.
make_numeric_count_cdf <- function(pdf, support)
{
  cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE,
    rel.tol = 1e-8, abs.tol = 0, max.terms = 1e6L, ...)
  {
    if(!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail))
      stop("'lower.tail' must be TRUE or FALSE.")
    if(!is.logical(log.p) || length(log.p) != 1L || is.na(log.p))
      stop("'log.p' must be TRUE or FALSE.")
    if(!is.numeric(y)) stop("'y' must be numeric.")
    if(!length(y)) return(numeric())
    if(!is.numeric(rel.tol) || length(rel.tol) != 1L ||
        !is.finite(rel.tol) || rel.tol <= 0)
      stop("'rel.tol' must be a positive finite number.")
    if(!is.numeric(abs.tol) || length(abs.tol) != 1L ||
        !is.finite(abs.tol) || abs.tol < 0)
      stop("'abs.tol' must be a nonnegative finite number.")
    if(!is.numeric(max.terms) || length(max.terms) != 1L ||
        !is.finite(max.terms) || max.terms < 1 ||
        max.terms != floor(max.terms))
      stop("'max.terms' must be a positive integer.")
    max.terms <- as.double(max.terms)

    nr <- family_parameter_rows(par)
    ny <- length(y)
    if(!nr) stop("family parameters must not be empty.")
    n <- max(nr, ny)
    if(!nr %in% c(1L, n) || !ny %in% c(1L, n))
      stop("lengths of family parameters and 'y' are incompatible.")

    parameters <- if(is.matrix(par)) as.list(as.data.frame(par)) else
      as.list(par)
    parameters <- lapply(parameters, rep, length.out = n)
    yy <- rep(y, length.out = n)
    dots <- list(...)
    bounds <- if(!length(dots)) {
      support(parameters)
    } else {
      do.call(support, c(list(par = parameters), dots))
    }
    value <- rep(NA_real_, n)
    batch.size <- 256L
    probability.tolerance <- max(10 * rel.tol, 10 * abs.tol,
      100 * .Machine$double.eps)

    for(i in seq_len(n)) {
      if(is.na(yy[i])) next

      lo <- ceiling(bounds[i, "min"])
      hi <- floor(bounds[i, "max"])
      if(!is.finite(lo))
        stop("count-family support must have a finite lower endpoint at row ",
          i, ".", call. = FALSE)
      if(lo > hi)
        stop("count-family support contains no integers at row ", i, ".",
          call. = FALSE)

      if(yy[i] < lo) {
        value[i] <- if(lower.tail) 0 else 1
        next
      }
      if(is.finite(hi) && yy[i] >= hi) {
        value[i] <- if(lower.tail) 1 else 0
        next
      }
      if(yy[i] == Inf) {
        value[i] <- if(lower.tail) 1 else 0
        next
      }

      pari <- lapply(parameters, function(z) z[i])
      evaluate <- function(z) {
        ans <- try(if(!length(dots)) {
          pdf(par = pari, y = z, log = FALSE)
        } else {
          do.call(pdf, c(list(par = pari, y = z, log = FALSE), dots))
        }, silent = TRUE)
        if(inherits(ans, "try-error") || length(ans) != length(z)) {
          ans <- vapply(z, function(zz) {
            zz <- if(!length(dots)) {
              pdf(par = pari, y = zz, log = FALSE)
            } else {
              do.call(pdf, c(list(par = pari, y = zz, log = FALSE), dots))
            }
            if(length(zz) != 1L)
              stop("the family density must return one value per response.")
            as.numeric(zz)
          }, numeric(1L))
        }
        ans <- as.numeric(ans)
        if(anyNA(ans) || any(!is.finite(ans)))
          stop("the family density returned non-finite values during summation.")
        if(any(ans < 0))
          stop("the family density returned negative values during summation.")
        ans
      }

      ## Sum in modest batches so very large thresholds do not require a huge
      ## temporary integer vector. On an infinite upper support, convergence
      ## is assessed only after two decreasing batches. This direct upper-tail
      ## sum avoids cancellation when the requested tail is small.
      sum_mass <- function(from, to, infinite = FALSE) {
        total <- 0
        used <- 0
        previous <- Inf
        current <- from
        repeat {
          remaining <- if(infinite) Inf else to - current + 1
          if(remaining <= 0) break
          take <- min(batch.size, remaining, max.terms - used)
          if(take < 1)
            stop("numerical CDF summation exceeded 'max.terms' at row ", i,
              ".", call. = FALSE)
          take <- as.integer(take)
          points <- current + seq_len(take) - 1L
          chunk <- sum(evaluate(points))
          total <- total + chunk
          used <- used + take

          if(!is.finite(total) || total > 1 + probability.tolerance)
            stop("numerical CDF summation produced a value outside [0, 1] at row ",
              i, ".", call. = FALSE)

          tolerance <- max(abs.tol, rel.tol * total)
          if(1 - total <= tolerance) break
          if(infinite && used >= 2L * batch.size && chunk <= previous &&
              chunk <= tolerance)
            break
          if(!infinite && take >= remaining) break

          next.current <- points[take] + 1
          if(!is.finite(next.current) || next.current == current)
            stop("integer support is too large to enumerate at row ", i, ".",
              call. = FALSE)
          current <- next.current
          previous <- chunk
        }
        min(max(total, 0), 1)
      }

      cutoff <- floor(yy[i])
      if(lower.tail) {
        probability <- sum_mass(lo, min(cutoff, hi))
      } else if(is.finite(hi)) {
        probability <- sum_mass(cutoff + 1, hi)
      } else {
        lower.probability <- sum_mass(lo, cutoff)
        probability <- if(lower.probability <= 0.5) {
          1 - lower.probability
        } else {
          sum_mass(cutoff + 1, Inf, infinite = TRUE)
        }
      }
      value[i] <- probability
    }

    if(log.p) value <- log(value)
    value
  }
  attr(cdf, "dnum") <- TRUE
  cdf
}

## Add a deterministic numerical CDF only when a family has a density and
## support. Continuous densities are integrated; discrete densities are
## summed over the integer support.
complete_family_cdf <- function(family)
{
  if(is.function(family[["cdf"]])) return(family)
  type <- family[["type"]]
  if(is.null(type)) type <- "continuous"
  type <- tolower(type[1L])
  if(is.na(type) || !type %in% c("continuous", "discrete") ||
      !is.function(family[["pdf"]]) ||
      !is.function(family[["support"]]))
    return(family)
  family$cdf <- if(type == "continuous") {
    make_numeric_cdf(family$pdf, family$support)
  } else {
    make_numeric_count_cdf(family$pdf, family$support)
  }
  family
}

## Numerically invert a CDF over the declared family support.
make_numeric_quantile <- function(cdf, support, type, pdf = NULL)
{
  type <- tolower(type[1L])
  force(cdf)
  force(support)
  force(type)
  force(pdf)

  quantile <- function(par, p, lower.tail = TRUE, log.p = FALSE,
    tol = sqrt(.Machine$double.eps), maxiter = 1000L,
    max.terms = 1e6L, ...)
  {
    if(!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail))
      stop("'lower.tail' must be TRUE or FALSE.")
    if(!is.logical(log.p) || length(log.p) != 1L || is.na(log.p))
      stop("'log.p' must be TRUE or FALSE.")
    if(!is.numeric(p)) stop("'p' must be numeric.")
    if(!length(p)) return(numeric())
    if(!is.numeric(tol) || length(tol) != 1L ||
        !is.finite(tol) || tol <= 0)
      stop("'tol' must be a positive finite number.")
    if(!is.numeric(maxiter) || length(maxiter) != 1L ||
        !is.finite(maxiter) || maxiter < 1 ||
        maxiter != floor(maxiter))
      stop("'maxiter' must be a positive integer.")
    maxiter <- as.integer(maxiter)
    if(!is.numeric(max.terms) || length(max.terms) != 1L ||
        !is.finite(max.terms) || max.terms < 1 ||
        max.terms != floor(max.terms))
      stop("'max.terms' must be a positive integer.")
    max.terms <- as.double(max.terms)

    probability <- p
    if(log.p) {
      if(any(probability > 0, na.rm = TRUE))
        stop("log probabilities must not be greater than zero.")
      probability <- if(lower.tail) {
        exp(probability)
      } else {
        -expm1(probability)
      }
    } else {
      if(any(probability < 0 | probability > 1, na.rm = TRUE))
        stop("'p' must contain probabilities in [0, 1].")
      if(!lower.tail) probability <- 1 - probability
    }

    nr <- family_parameter_rows(par)
    np <- length(probability)
    if(!nr) stop("family parameters must not be empty.")
    n <- max(nr, np)
    if(!nr %in% c(1L, n) || !np %in% c(1L, n))
      stop("lengths of family parameters and 'p' are incompatible.")

    parameters <- if(is.matrix(par)) as.list(as.data.frame(par)) else
      as.list(par)
    parameters <- lapply(parameters, rep, length.out = n)
    probability <- rep(probability, length.out = n)
    dots <- list(...)
    bounds <- if(!length(dots)) {
      support(parameters)
    } else {
      do.call(support, c(list(par = parameters), dots))
    }
    value <- rep(NA_real_, n)

    lo <- bounds[, "min"]
    hi <- bounds[, "max"]
    if(identical(type, "discrete")) {
      lo <- ceiling(lo)
      hi <- floor(hi)
    }
    relevant <- which(!is.na(probability))
    empty <- relevant[lo[relevant] > hi[relevant]]
    if(length(empty))
      stop("family support is empty at row ", empty[1L], ".", call. = FALSE)

    at.zero <- relevant[probability[relevant] <= 0]
    at.one <- relevant[probability[relevant] >= 1]
    value[at.zero] <- lo[at.zero]
    value[at.one] <- hi[at.one]
    interior <- relevant[
      probability[relevant] > 0 & probability[relevant] < 1
    ]
    if(!length(interior)) return(value)

    ## When the CDF itself was generated from a count PMF, invert that PMF in
    ## one forward pass. Repeated CDF bisection would rescan the same mass many
    ## times and is substantially slower.
    if(identical(type, "discrete") &&
        isTRUE(attr(cdf, "dnum")) && is.function(pdf)) {
      bad <- interior[!is.finite(lo[interior])]
      if(length(bad))
        stop("count-family support must have a finite lower endpoint at row ",
          bad[1L], ".", call. = FALSE)

      pdf.dots <- dots
      pdf.dots[c("par", "y", "log", "lower.tail", "log.p")] <- NULL
      batch.size <- 64L
      probability.tolerance <- 100 * .Machine$double.eps

      for(i in interior) {
        pari <- lapply(parameters, function(z) z[i])
        evaluate <- function(y) {
          args <- c(list(par = pari, y = y, log = FALSE), pdf.dots)
          mass <- try(do.call(pdf, args), silent = TRUE)
          if(inherits(mass, "try-error") || length(mass) != length(y)) {
            mass <- vapply(y, function(yy) {
              args$y <- yy
              ans <- tryCatch(
                do.call(pdf, args),
                error = function(e) e
              )
              if(inherits(ans, "error"))
                stop("numerical quantile inversion failed at row ", i, ": ",
                  conditionMessage(ans), call. = FALSE)
              if(!is.numeric(ans) || length(ans) != 1L)
                stop("the family density must return one value per response.")
              as.numeric(ans)
            }, numeric(1L))
          } else {
            mass <- as.numeric(mass)
          }
          if(anyNA(mass) || any(!is.finite(mass)))
            stop("the family density returned non-finite values during ",
              "quantile inversion at row ", i, ".", call. = FALSE)
          if(any(mass < 0))
            stop("the family density returned negative values during ",
              "quantile inversion at row ", i, ".", call. = FALSE)
          mass
        }

        current <- lo[i]
        cumulative <- 0
        used <- 0
        repeat {
          remaining <- if(is.finite(hi[i])) hi[i] - current + 1 else Inf
          if(remaining <= 0)
            stop("the family density and support do not reach probability at row ",
              i, ".", call. = FALSE)
          take <- min(batch.size, remaining, max.terms - used)
          if(take < 1)
            stop("numerical quantile inversion exceeded 'max.terms' at row ",
              i, ".", call. = FALSE)
          take <- as.integer(take)
          points <- current + seq_len(take) - 1L
          cumulative.mass <- cumulative + cumsum(evaluate(points))
          hit <- which(cumulative.mass >= probability[i])[1L]
          if(!is.na(hit)) {
            value[i] <- points[hit]
            break
          }

          cumulative <- cumulative.mass[take]
          if(!is.finite(cumulative) ||
              cumulative > 1 + probability.tolerance)
            stop("numerical quantile inversion produced cumulative mass ",
              "outside [0, 1] at row ", i, ".", call. = FALSE)
          used <- used + take
          if(take >= remaining)
            stop("the family density and support do not reach probability at row ",
              i, ".", call. = FALSE)

          next.current <- points[take] + 1
          if(!is.finite(next.current) || next.current == current)
            stop("integer support is too large to enumerate at row ", i, ".",
              call. = FALSE)
          current <- next.current
        }
      }
      return(value)
    }

    ## Start bracketing near the first finite parameter, which is generally a
    ## location or mean. Construct this without a row-wise list conversion.
    center <- rep(NA_real_, n)
    for(parameter in parameters) {
      candidate <- suppressWarnings(as.numeric(parameter))
      replace <- !is.finite(center) & is.finite(candidate)
      center[replace] <- candidate[replace]
    }
    center[!is.finite(center)] <- 0
    finite.lower <- is.finite(lo)
    finite.upper <- is.finite(hi)
    center[finite.lower] <- pmax(center[finite.lower], lo[finite.lower])
    center[finite.upper] <- pmin(center[finite.upper], hi[finite.upper])

    ## Evaluate all active rows in one CDF call. If an arbitrary user CDF is
    ## scalar-only, remember that fact and retain the compatible row fallback.
    cdf.formals <- names(formals(cdf))
    cdf.dots <- dots
    cdf.dots[c("par", "y", "lower.tail", "log.p", "log")] <- NULL
    vectorized <- TRUE

    cdf_call <- function(index, y) {
      args <- list(
        par = lapply(parameters, function(z) z[index]),
        y = y
      )
      if("lower.tail" %in% cdf.formals) args$lower.tail <- TRUE
      if("log.p" %in% cdf.formals) {
        args$log.p <- FALSE
      } else if("log" %in% cdf.formals) {
        args$log <- FALSE
      }
      do.call(cdf, c(args, cdf.dots))
    }

    cdf_values <- function(index, y) {
      if(!length(index)) return(numeric())
      ans <- NULL
      if(vectorized) {
        ans <- try(cdf_call(index, y), silent = TRUE)
        if(inherits(ans, "try-error") || length(ans) != length(index)) {
          vectorized <<- FALSE
          ans <- NULL
        }
      }
      if(is.null(ans)) {
        ans <- vapply(seq_along(index), function(k) {
          row <- index[k]
          out <- tryCatch(
            cdf_call(row, y[k]),
            error = function(e) e
          )
          if(inherits(out, "error"))
            stop("numerical quantile inversion failed at row ", row, ": ",
              conditionMessage(out), call. = FALSE)
          if(!is.numeric(out) || length(out) != 1L)
            stop("the family CDF must return one value per response.")
          as.numeric(out)
        }, numeric(1L))
      } else {
        ans <- as.numeric(ans)
      }

      bad <- which(is.na(ans) | !is.finite(ans))
      if(length(bad))
        stop("the family CDF returned a non-finite value at row ",
          index[bad[1L]], ".", call. = FALSE)
      cdf.tolerance <- 100 * .Machine$double.eps
      bad <- which(ans < -cdf.tolerance | ans > 1 + cdf.tolerance)
      if(length(bad))
        stop("the family CDF returned a value outside [0, 1] at row ",
          index[bad[1L]], ".", call. = FALSE)
      pmin(pmax(ans, 0), 1)
    }

    if(identical(type, "discrete")) {
      bad <- interior[!is.finite(lo[interior])]
      if(length(bad))
        stop("count-family support must have a finite lower endpoint at row ",
          bad[1L], ".", call. = FALSE)

      center <- floor(center)
      center <- pmax(center, lo)
      center[finite.upper] <- pmin(center[finite.upper], hi[finite.upper])
      distribution <- cdf_values(interior, center[interior])
      left <- right <- rep(NA_real_, n)
      step <- rep(1, n)

      ## Search downward together for rows whose starting CDF is already above
      ## the target.
      pending <- interior[distribution >= probability[interior]]
      right[pending] <- center[pending]
      for(iteration in seq_len(maxiter)) {
        if(!length(pending)) break
        candidate <- pmax(lo[pending], center[pending] - step[pending])
        candidate.distribution <- cdf_values(pending, candidate)
        below <- candidate.distribution < probability[pending]

        if(any(below)) {
          found <- pending[below]
          left[found] <- candidate[below] + 1
        }
        stay <- pending[!below]
        if(length(stay)) {
          stay.candidate <- candidate[!below]
          right[stay] <- stay.candidate
          at.lower <- stay.candidate <= lo[stay]
          if(any(at.lower)) {
            found <- stay[at.lower]
            left[found] <- right[found]
          }
          stay <- stay[!at.lower]
        }
        pending <- stay
        step[pending] <- step[pending] * 2
      }
      if(length(pending))
        stop("could not bracket the requested quantile at row ",
          pending[1L], ".", call. = FALSE)

      ## Search upward together for rows whose starting CDF is below target.
      pending <- interior[distribution < probability[interior]]
      left[pending] <- center[pending] + 1
      step[pending] <- 1
      for(iteration in seq_len(maxiter)) {
        if(!length(pending)) break
        candidate <- center[pending] + step[pending]
        candidate[!is.finite(candidate)] <- .Machine$double.xmax
        bounded <- is.finite(hi[pending])
        candidate[bounded] <- pmin(
          candidate[bounded], hi[pending][bounded]
        )
        candidate.distribution <- cdf_values(pending, candidate)
        above <- candidate.distribution >= probability[pending]

        if(any(above)) {
          found <- pending[above]
          right[found] <- candidate[above]
        }
        stay <- pending[!above]
        if(length(stay)) {
          stay.candidate <- candidate[!above]
          next.left <- stay.candidate + 1
          stuck <- stay.candidate >= hi[stay] |
            !is.finite(next.left) | next.left <= stay.candidate
          if(any(stuck))
            stop("could not bracket the requested quantile at row ",
              stay[which(stuck)[1L]], ".", call. = FALSE)
          left[stay] <- next.left
        }
        pending <- stay
        step[pending] <- step[pending] * 2
      }
      if(length(pending))
        stop("could not bracket the requested quantile at row ",
          pending[1L], ".", call. = FALSE)

      ## Batched integer binary search for the generalized inverse.
      pending <- interior[left[interior] < right[interior]]
      for(iteration in seq_len(maxiter)) {
        if(!length(pending)) break
        midpoint <- floor(left[pending] / 2 + right[pending] / 2)
        adjacent <- midpoint <= left[pending]

        if(any(adjacent)) {
          rows <- pending[adjacent]
          distribution <- cdf_values(rows, left[rows])
          use.left <- distribution >= probability[rows]
          right[rows[use.left]] <- left[rows[use.left]]
          left[rows[!use.left]] <- right[rows[!use.left]]
        }

        rows <- pending[!adjacent]
        if(length(rows)) {
          distribution <- cdf_values(rows, midpoint[!adjacent])
          move.right <- distribution >= probability[rows]
          right[rows[move.right]] <- midpoint[!adjacent][move.right]
          left[rows[!move.right]] <- midpoint[!adjacent][!move.right] + 1
        }
        pending <- pending[left[pending] < right[pending]]
      }
      if(length(pending))
        stop("numerical quantile inversion did not converge at row ",
          pending[1L], ".", call. = FALSE)
      value[interior] <- left[interior]
    } else {
      lower <- lo
      upper <- hi
      unbounded <- interior[
        !is.finite(lo[interior]) | !is.finite(hi[interior])
      ]
      center.distribution <- rep(NA_real_, n)
      center.distribution[unbounded] <- cdf_values(
        unbounded, center[unbounded]
      )
      initial.step <- pmax(1, abs(center))

      ## Batch the geometric search for all unbounded lower endpoints.
      pending <- unbounded[
        !is.finite(lo[unbounded]) &
          center.distribution[unbounded] >= probability[unbounded]
      ]
      upper[pending] <- center[pending]
      step <- initial.step
      previous <- center
      for(iteration in seq_len(maxiter)) {
        if(!length(pending)) break
        candidate <- center[pending] - step[pending]
        candidate[!is.finite(candidate)] <- -.Machine$double.xmax
        stuck <- candidate == previous[pending]
        if(any(stuck))
          stop("could not bracket the requested quantile at row ",
            pending[which(stuck)[1L]], ".", call. = FALSE)
        candidate.distribution <- cdf_values(pending, candidate)
        below <- candidate.distribution < probability[pending]
        if(any(below)) {
          found <- pending[below]
          lower[found] <- candidate[below]
        }
        stay <- pending[!below]
        if(length(stay)) {
          upper[stay] <- candidate[!below]
          previous[stay] <- candidate[!below]
          step[stay] <- step[stay] * 2
        }
        pending <- stay
      }
      if(length(pending))
        stop("could not bracket the requested quantile at row ",
          pending[1L], ".", call. = FALSE)
      rows <- unbounded[
        !is.finite(lo[unbounded]) &
          center.distribution[unbounded] < probability[unbounded]
      ]
      lower[rows] <- center[rows]

      ## Batch the corresponding search for unbounded upper endpoints.
      pending <- unbounded[
        !is.finite(hi[unbounded]) &
          center.distribution[unbounded] < probability[unbounded]
      ]
      lower[pending] <- center[pending]
      step <- initial.step
      previous <- center
      for(iteration in seq_len(maxiter)) {
        if(!length(pending)) break
        candidate <- center[pending] + step[pending]
        candidate[!is.finite(candidate)] <- .Machine$double.xmax
        stuck <- candidate == previous[pending]
        if(any(stuck))
          stop("could not bracket the requested quantile at row ",
            pending[which(stuck)[1L]], ".", call. = FALSE)
        candidate.distribution <- cdf_values(pending, candidate)
        above <- candidate.distribution >= probability[pending]
        if(any(above)) {
          found <- pending[above]
          upper[found] <- candidate[above]
        }
        stay <- pending[!above]
        if(length(stay)) {
          lower[stay] <- candidate[!above]
          previous[stay] <- candidate[!above]
          step[stay] <- step[stay] * 2
        }
        pending <- stay
      }
      if(length(pending))
        stop("could not bracket the requested quantile at row ",
          pending[1L], ".", call. = FALSE)
      rows <- unbounded[
        !is.finite(hi[unbounded]) &
          center.distribution[unbounded] >= probability[unbounded]
      ]
      upper[rows] <- center[rows]

      lower.distribution <- cdf_values(interior, lower[interior])
      at.lower <- interior[
        lower.distribution >= probability[interior]
      ]
      value[at.lower] <- lower[at.lower]
      bracket <- setdiff(interior, at.lower)
      if(length(bracket)) {
        upper.distribution <- cdf_values(bracket, upper[bracket])
        bad <- bracket[upper.distribution < probability[bracket]]
        if(length(bad))
          stop("the family CDF and support do not bracket probability at row ",
            bad[1L], ".", call. = FALSE)

        bracket.width <- upper - lower
        bracket.scale <- ifelse(
          is.finite(bracket.width),
          pmax(1, bracket.width),
          pmax(1, abs(lower), abs(upper))
        )
        pending <- bracket
        for(iteration in seq_len(maxiter)) {
          if(!length(pending)) break
          width <- upper[pending] - lower[pending]
          rounding.tolerance <- 4 * .Machine$double.eps *
            pmax(1, abs(lower[pending]), abs(upper[pending]))
          converged <- width <= pmax(
            tol * bracket.scale[pending], rounding.tolerance
          )
          pending <- pending[!converged]
          if(!length(pending)) break

          midpoint <- lower[pending] / 2 + upper[pending] / 2
          stuck <- midpoint == lower[pending] |
            midpoint == upper[pending]
          pending <- pending[!stuck]
          midpoint <- midpoint[!stuck]
          if(!length(pending)) break

          distribution <- cdf_values(pending, midpoint)
          move.upper <- distribution >= probability[pending]
          upper[pending[move.upper]] <- midpoint[move.upper]
          lower[pending[!move.upper]] <- midpoint[!move.upper]
        }

        if(length(pending)) {
          width <- upper[pending] - lower[pending]
          rounding.tolerance <- 4 * .Machine$double.eps *
            pmax(1, abs(lower[pending]), abs(upper[pending]))
          midpoint <- lower[pending] / 2 + upper[pending] / 2
          converged <- width <= pmax(
            tol * bracket.scale[pending], rounding.tolerance
          ) | midpoint == lower[pending] | midpoint == upper[pending]
          pending <- pending[!converged]
        }
        if(length(pending))
          stop("numerical quantile inversion did not converge at row ",
            pending[1L], ".", call. = FALSE)
        value[bracket] <- lower[bracket] / 2 + upper[bracket] / 2
      }
    }

    value
  }
  attr(quantile, "qnum") <- TRUE
  quantile
}
## Add a deterministic numerical quantile when the CDF and support are
## available. Existing analytical quantile functions are retained.
complete_family_quantile <- function(family)
{
  if(is.function(family[["quantile"]])) return(family)
  type <- family[["type"]]
  if(is.null(type)) type <- "continuous"
  type <- tolower(type[1L])
  if(is.na(type) || !type %in% c("continuous", "discrete") ||
      !is.function(family[["cdf"]]) ||
      !is.function(family[["support"]]))
    return(family)
  family$quantile <- make_numeric_quantile(
    family$cdf, family$support, type, family$pdf
  )
  family
}

## Numerically integrate the mean or centered second moment of a continuous
## density. The scalar-density fallback mirrors make_numeric_cdf().
make_numeric_continuous_moment <- function(pdf, support, center = NULL,
  quantile = NULL)
{
  force(pdf)
  force(support)
  force(center)
  force(quantile)

  moment <- function(par, rel.tol = 1e-8, abs.tol = 0,
    subdivisions = 100L, ...)
  {
    if(!is.numeric(rel.tol) || length(rel.tol) != 1L ||
        !is.finite(rel.tol) || rel.tol <= 0)
      stop("'rel.tol' must be a positive finite number.")
    if(!is.numeric(abs.tol) || length(abs.tol) != 1L ||
        !is.finite(abs.tol) || abs.tol < 0)
      stop("'abs.tol' must be a nonnegative finite number.")
    if(!is.numeric(subdivisions) || length(subdivisions) != 1L ||
        !is.finite(subdivisions) || subdivisions < 1 ||
        subdivisions != floor(subdivisions))
      stop("'subdivisions' must be a positive integer.")
    subdivisions <- as.integer(subdivisions)

    nr <- family_parameter_rows(par)
    if(!nr) stop("family parameters must not be empty.")
    parameters <- if(is.matrix(par)) as.list(as.data.frame(par)) else
      as.list(par)
    parameters <- lapply(parameters, rep, length.out = nr)
    dots <- list(...)
    bounds <- if(!length(dots)) {
      support(parameters)
    } else {
      do.call(support, c(list(par = parameters), dots))
    }

    centers <- NULL
    if(is.function(center)) {
      center.dots <- dots
      center.dots[c("par", "rel.tol", "abs.tol", "subdivisions")] <- NULL
      if(isTRUE(attr(center, "mnum", exact = TRUE))) {
        centers <- do.call(center, c(list(par = parameters,
          rel.tol = rel.tol, abs.tol = abs.tol,
          subdivisions = subdivisions), center.dots))
      } else {
        centers <- do.call(center, c(list(par = parameters), center.dots))
      }
      centers <- as.numeric(centers)
      if(length(centers) != nr || anyNA(centers) || any(!is.finite(centers)))
        stop("the family mean must return one finite value per parameter row.")
    }

    value <- numeric(nr)
    for(i in seq_len(nr)) {
      lo <- bounds[i, "min"]
      hi <- bounds[i, "max"]
      if(lo >= hi) {
        value[i] <- if(is.null(centers)) lo else 0
        next
      }

      pari <- lapply(parameters, function(z) z[i])
      evaluate <- function(z) {
        if(!length(dots)) {
          pdf(par = pari, y = z, log = FALSE)
        } else {
          do.call(pdf, c(list(par = pari, y = z, log = FALSE), dots))
        }
      }
      vectorized <- TRUE
      density <- function(z) {
        ans <- NULL
        if(vectorized) {
          ans <- try(evaluate(z), silent = TRUE)
          if(inherits(ans, "try-error") || length(ans) != length(z)) {
            vectorized <<- FALSE
            ans <- NULL
          }
        }
        if(is.null(ans)) {
          ans <- vapply(z, function(zz) {
            res <- evaluate(zz)
            if(length(res) != 1L)
              stop("the family density must return one value per response.")
            as.numeric(res)
          }, numeric(1L))
        }
        ans <- as.numeric(ans)
        if(anyNA(ans) || any(!is.finite(ans)))
          stop("the family density returned non-finite values during integration.")
        if(any(ans < 0))
          stop("the family density returned negative values during integration.")
        ans
      }

      ## Splitting an infinite interval near its probability mass prevents
      ## integrate() from missing a density located far from zero. Prefer an
      ## analytical median, then a known mean, and finally ordinary support
      ## and parameter anchors.
      anchor <- NA_real_
      if(is.function(quantile) &&
          !isTRUE(attr(quantile, "qnum", exact = TRUE))) {
        qdots <- dots
        qdots[c("par", "p", "lower.tail", "log.p")] <- NULL
        median <- try(do.call(quantile,
          c(list(par = pari, p = 0.5), qdots)), silent = TRUE)
        if(!inherits(median, "try-error") && length(median) == 1L &&
            is.finite(median) && median >= lo && median <= hi)
          anchor <- as.numeric(median)
      }
      if(!is.finite(anchor) && !is.null(centers) &&
          centers[i] >= lo && centers[i] <= hi)
        anchor <- centers[i]
      if(!is.finite(anchor)) {
        pv <- suppressWarnings(as.numeric(unlist(pari, use.names = FALSE)))
        midpoint <- if(is.finite(lo) && is.finite(hi)) lo / 2 + hi / 2 else
          numeric()
        candidates <- unique(c(midpoint, 0, pv))
        candidates <- candidates[
          is.finite(candidates) & candidates >= lo & candidates <= hi
        ]
        if(length(candidates)) {
          log_density <- vapply(candidates, function(z) {
            ans <- try(suppressWarnings(if(!length(dots)) {
              pdf(par = pari, y = z, log = TRUE)
            } else {
              do.call(pdf, c(list(par = pari, y = z, log = TRUE), dots))
            }), silent = TRUE)
            if(inherits(ans, "try-error") || length(ans) != 1L || is.na(ans))
              return(-Inf)
            as.numeric(ans)
          }, numeric(1L))
          if(any(log_density > -Inf))
            anchor <- candidates[which.max(log_density)]
        }
      }

      transform <- if(is.null(centers)) {
        function(z) z
      } else {
        centeri <- centers[i]
        function(z) (z - centeri)^2
      }
      integrand <- function(z) {
        ans <- transform(z) * density(z)
        if(anyNA(ans) || any(!is.finite(ans)))
          stop("the requested numerical moment is not finite.")
        ans
      }
      integrate_one <- function(fun) {
        label <- if(is.null(centers)) "mean" else "variance"
        limits <- list(c(lo, hi))
        if(is.finite(anchor) && anchor > lo && anchor < hi)
          limits <- list(c(lo, anchor), c(anchor, hi))
        values <- vapply(limits, function(limitsi) {
          result <- tryCatch(
            stats::integrate(fun, lower = limitsi[1L], upper = limitsi[2L],
              subdivisions = subdivisions, rel.tol = rel.tol,
              abs.tol = abs.tol, stop.on.error = FALSE),
            error = function(e) e
          )
          if(inherits(result, "error"))
            stop("numerical ", label, " integration failed at row ", i,
              ": ", conditionMessage(result), call. = FALSE)
          if(!identical(result$message, "OK"))
            stop("numerical ", label, " integration failed at row ", i,
              ": ", result$message, call. = FALSE)
          if(!is.finite(result$value))
            stop("numerical ", label, " is not finite at row ", i, ".",
              call. = FALSE)
          result$value
        }, numeric(1L))
        sum(values)
      }

      value[i] <- integrate_one(integrand)

      ## A signed first-moment integral can appear finite through cancellation
      ## even when the expectation does not exist (for example, a Cauchy
      ## density). Verify absolute integrability whenever support crosses zero.
      if(is.null(centers) && lo < 0 && hi > 0)
        integrate_one(function(z) abs(z) * density(z))
    }

    if(!is.null(centers)) {
      tolerance <- max(10 * abs.tol, 100 * .Machine$double.eps)
      if(any(value < -tolerance))
        stop("numerical variance is negative.", call. = FALSE)
      value <- pmax(value, 0)
    }
    value
  }
  attr(moment, if(is.null(center)) "mnum" else "vnum") <- TRUE
  moment
}

## Numerically sum count moments in batches. Variance uses a weighted online
## update, avoiding cancellation in E[X^2] - E[X]^2.
make_numeric_count_moment <- function(pdf, support, second = FALSE)
{
  force(pdf)
  force(support)
  force(second)

  moment <- function(par, rel.tol = 1e-8, abs.tol = 0,
    max.terms = 1e6L, ...)
  {
    if(!is.numeric(rel.tol) || length(rel.tol) != 1L ||
        !is.finite(rel.tol) || rel.tol <= 0)
      stop("'rel.tol' must be a positive finite number.")
    if(!is.numeric(abs.tol) || length(abs.tol) != 1L ||
        !is.finite(abs.tol) || abs.tol < 0)
      stop("'abs.tol' must be a nonnegative finite number.")
    if(!is.numeric(max.terms) || length(max.terms) != 1L ||
        !is.finite(max.terms) || max.terms < 1 ||
        max.terms != floor(max.terms))
      stop("'max.terms' must be a positive integer.")
    max.terms <- as.double(max.terms)

    nr <- family_parameter_rows(par)
    if(!nr) stop("family parameters must not be empty.")
    parameters <- if(is.matrix(par)) as.list(as.data.frame(par)) else
      as.list(par)
    parameters <- lapply(parameters, rep, length.out = nr)
    dots <- list(...)
    bounds <- if(!length(dots)) {
      support(parameters)
    } else {
      do.call(support, c(list(par = parameters), dots))
    }
    value <- numeric(nr)
    batch.size <- 256L
    probability.tolerance <- max(10 * rel.tol, 10 * abs.tol,
      100 * .Machine$double.eps)

    for(i in seq_len(nr)) {
      lo <- ceiling(bounds[i, "min"])
      hi <- floor(bounds[i, "max"])
      if(!is.finite(lo))
        stop("count-family support must have a finite lower endpoint at row ",
          i, ".", call. = FALSE)
      if(lo > hi)
        stop("count-family support contains no integers at row ", i, ".",
          call. = FALSE)

      pari <- lapply(parameters, function(z) z[i])
      evaluate <- function(z) {
        ans <- try(if(!length(dots)) {
          pdf(par = pari, y = z, log = FALSE)
        } else {
          do.call(pdf, c(list(par = pari, y = z, log = FALSE), dots))
        }, silent = TRUE)
        if(inherits(ans, "try-error") || length(ans) != length(z)) {
          ans <- vapply(z, function(zz) {
            res <- if(!length(dots)) {
              pdf(par = pari, y = zz, log = FALSE)
            } else {
              do.call(pdf, c(list(par = pari, y = zz, log = FALSE), dots))
            }
            if(length(res) != 1L)
              stop("the family density must return one value per response.")
            as.numeric(res)
          }, numeric(1L))
        }
        ans <- as.numeric(ans)
        if(anyNA(ans) || any(!is.finite(ans)))
          stop("the family density returned non-finite values during summation.")
        if(any(ans < 0))
          stop("the family density returned negative values during summation.")
        ans
      }

      mass <- mean.value <- m2 <- absolute.first <- raw.second <- 0
      previous <- c(mass = Inf, first = Inf, second = Inf)
      used <- 0
      current <- lo
      converged <- FALSE
      repeat {
        remaining <- if(is.finite(hi)) hi - current + 1 else Inf
        if(remaining <= 0) {
          converged <- TRUE
          break
        }
        take <- min(batch.size, remaining, max.terms - used)
        if(take < 1) break
        take <- as.integer(take)
        points <- current + seq_len(take) - 1L
        probabilities <- evaluate(points)
        chunk.mass <- sum(probabilities)
        chunk.first <- sum(abs(points) * probabilities)
        chunk.second <- if(second) sum(points^2 * probabilities) else 0
        if(any(!is.finite(c(chunk.mass, chunk.first, chunk.second))))
          stop("the requested numerical moment is not finite at row ", i,
            ".", call. = FALSE)

        if(chunk.mass > 0) {
          chunk.mean <- sum(points * probabilities) / chunk.mass
          chunk.m2 <- if(second)
            sum((points - chunk.mean)^2 * probabilities) else 0
          combined.mass <- mass + chunk.mass
          delta <- chunk.mean - mean.value
          if(second)
            m2 <- m2 + chunk.m2 + delta^2 * mass * chunk.mass / combined.mass
          mean.value <- mean.value + delta * chunk.mass / combined.mass
          mass <- combined.mass
        }
        absolute.first <- absolute.first + chunk.first
        raw.second <- raw.second + chunk.second
        used <- used + take

        if(!is.finite(mean.value) || !is.finite(absolute.first) ||
            (second && (!is.finite(m2) || !is.finite(raw.second))))
          stop("the requested numerical moment is not finite at row ", i,
            ".", call. = FALSE)
        if(!is.finite(mass) || mass > 1 + probability.tolerance)
          stop("numerical moment summation produced probability mass outside ",
            "[0, 1] at row ", i, ".", call. = FALSE)

        if(is.finite(hi) && take >= remaining) {
          converged <- TRUE
          break
        }
        mass.tol <- max(abs.tol, rel.tol)
        first.tol <- max(abs.tol,
          rel.tol * max(1, absolute.first))
        second.tol <- max(abs.tol,
          rel.tol * max(1, raw.second))
        if(used >= 2L * batch.size && chunk.mass <= previous["mass"] &&
            chunk.first <= previous["first"] &&
            (!second || chunk.second <= previous["second"]) &&
            abs(1 - mass) <= mass.tol && chunk.first <= first.tol &&
            (!second || chunk.second <= second.tol)) {
          converged <- TRUE
          break
        }

        next.current <- points[take] + 1
        if(!is.finite(next.current) || next.current == current)
          stop("integer support is too large to enumerate at row ", i, ".",
            call. = FALSE)
        current <- next.current
        previous <- c(mass = chunk.mass, first = chunk.first,
          second = chunk.second)
      }

      label <- if(second) "variance" else "mean"
      if(!converged)
        stop("numerical ", label, " summation exceeded 'max.terms' at row ",
          i, ".", call. = FALSE)
      if(mass <= 0 || abs(1 - mass) > probability.tolerance)
        stop("numerical ", label, " summation did not recover unit ",
          "probability mass at row ", i, ".", call. = FALSE)
      value[i] <- if(second) m2 / mass else mean.value
    }
    value
  }
  attr(moment, if(second) "vnum" else "mnum") <- TRUE
  moment
}

## Add numerical mean and variance functions without replacing analytical
## implementations.
complete_family_moments <- function(family)
{
  type <- family[["type"]]
  if(is.null(type)) type <- "continuous"
  type <- tolower(type[1L])
  if(is.na(type) || !type %in% c("continuous", "discrete") ||
      !is.function(family[["pdf"]]) ||
      !is.function(family[["support"]]))
    return(family)

  if(!is.function(family[["mean"]])) {
    family$mean <- if(type == "continuous") {
      make_numeric_continuous_moment(
        family$pdf, family$support, quantile = family$quantile
      )
    } else {
      make_numeric_count_moment(family$pdf, family$support)
    }
  }
  if(!is.function(family[["variance"]])) {
    family$variance <- if(type == "continuous") {
      make_numeric_continuous_moment(
        family$pdf, family$support, center = family$mean,
        quantile = family$quantile
      )
    } else {
      make_numeric_count_moment(family$pdf, family$support, second = TRUE)
    }
  }
  family
}

## Generate random values by inverse transform when a quantile function is
## available. One vectorized quantile call handles every requested draw.
make_numeric_random <- function(quantile)
{
  force(quantile)
  random <- function(par, n, ...)
  {
    if(!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 ||
        n != floor(n) || n > .Machine$integer.max)
      stop("'n' must be a nonnegative integer.", call. = FALSE)
    n <- as.integer(n)
    if(n == 0L) return(numeric())

    parameters <- if(is.matrix(par)) as.list(as.data.frame(par)) else
      as.list(par)
    nr <- family_parameter_rows(parameters)
    if(!nr) stop("family parameters must not be empty.")
    parameters <- lapply(parameters, function(z) {
      rep(rep(z, length.out = nr), times = n)
    })
    probability <- stats::runif(nr * n)
    dots <- list(...)
    dots[c("par", "p", "lower.tail", "log.p")] <- NULL
    value <- as.numeric(do.call(quantile,
      c(list(par = parameters, p = probability), dots)))
    if(length(value) != nr * n)
      stop("the family quantile must return one value per probability.")

    value <- matrix(value, nrow = nr, ncol = n)
    if(nr == 1L) return(as.vector(value))
    if(n == 1L) return(value[, 1L])
    colnames(value) <- paste0("r_", seq_len(n))
    value
  }
  attr(random, "rnum") <- TRUE
  random
}

complete_family_random <- function(family)
{
  if(!is.function(family[["random"]]) &&
      is.function(family[["quantile"]]))
    family$random <- make_numeric_random(family$quantile)
  family
}

## Complete a family object, e.g.,
## if derivatives are not supplied they
## will be approximated numerically.
complete_family <- function(family, .links = NULL)
{
  if(!is.null(attr(family, "family"))) {
    family <- attr(family, "family")
  }

  if(is.character(family)) {
    family <- get(family)
  }
    
  family_constructor <- NULL
  if(is.function(family)) {
    family_constructor <- family
    family <- family()
  }

  if(inherits(family, "gamlss.family")) {
    ## The generated family functions are costly to byte-compile. Cache only
    ## zero-argument constructors bound under their family name in a locked
    ## package namespace; user closures and modified family objects must retain
    ## their existing per-call semantics.
    cache_key <- NULL
    if(!is.null(family_constructor)) {
      constructor_env <- environment(family_constructor)
      family_name <- family$family[1L]
      if(isNamespace(constructor_env) &&
          exists(family_name, envir = constructor_env, inherits = FALSE) &&
          bindingIsLocked(family_name, constructor_env) &&
          identical(get(family_name, envir = constructor_env, inherits = FALSE),
            family_constructor)) {
        cache_key <- paste0(environmentName(constructor_env), "::", family_name)
      }
    }

    family_signature <- if(is.null(cache_key)) NULL else
      try(serialize(family, connection = NULL, version = 3L), silent = TRUE)
    if(inherits(family_signature, "try-error")) {
      cache_key <- NULL
      family_signature <- NULL
    }
    if(!is.null(cache_key) &&
        exists(cache_key, envir = .gamlss2_family_cache, inherits = FALSE)) {
      cached <- get(cache_key, envir = .gamlss2_family_cache, inherits = FALSE)
      if(identical(cached$signature, family_signature))
        return(cached$family)
    }

    family <- tF(family)
    if(!is.null(cache_key))
      assign(cache_key, list(signature = family_signature, family = family),
        envir = .gamlss2_family_cache)
    return(family)
  }

  if(is.null(family$family)) {
    family$family <- "No family name supplied!"
  }

  if(inherits(family, "distribution")) {
    fn <- class(family)[1L]
    ff <- get(fn)
    np <- names(formals(ff))
    if(is.null(.links)) {
      stop(paste0("no links for parameters (",
        paste0(np, collapse = ", "),
        ") supplied!"))
    }
    names(.links) <- np
    family <- family(family, links = .links)
  }

  if(!is.null(family[["d"]])) {
    family[["pdf"]] <- family[["d"]]
  }
  if(!is.null(family[["p"]])) {
    family[["cdf"]] <- family[["p"]]
  }
  if(!is.null(family[["q"]])) {
    family[["quantile"]] <- family[["q"]]
  }
  if(!is.null(family[["r"]])) {
    family[["random"]] <- family[["r"]]
  }
  family[c("d", "p", "q", "r")] <- NULL

  family <- complete_family_support(family)

  if(is.null(family$pdf))
    stop("the family needs a $pdf() function!")

  family <- complete_family_cdf(family)
  family <- complete_family_quantile(family)
  family <- complete_family_moments(family)
  family <- complete_family_random(family)

  use_numeric_update <-
    is.null(family[["update"]]) &&
    is.null(family[["score"]]) &&
    is.null(family[["hessian"]]) &&
    is.null(family[["hess"]])

  if(is.null(family$log_likelihood)) {
    family$log_likelihood <- function(par, y, ...) {
      sum(family$pdf(par = par, y = y, log = TRUE, ...), na.rm = TRUE)
    }
    family["logLik"] <- NULL
  }

  if(!is.list(family$links))
    family$links <- as.list(family$links)
  if(is.null(names(family$links)))
    names(family$links) <- family$names

  linkinv <- linkfun <- mu.eta <- list()
  for(j in family$names) {
    link <- make.link2(family$links[[j]])
    linkinv[[j]] <- link$linkinv
    linkfun[[j]] <- link$linkfun
    mu.eta[[j]] <- link$mu.eta
  }

  if(use_numeric_update) {
    family$update <- make_numeric_update(
      pdf = family$pdf,
      linkinv = linkinv
    )
  }

  if(is.null(family$map2par)) {
    family$map2par <- function(eta) {
      for(j in names(eta)) {
        z <- linkinv[[j]](eta[[j]])
        if(any(!is.finite(z))) {
          if(any(jj <- is.na(z)))
            z[jj] <- 0
          if(any(jj <- z == Inf))
            z[jj] <- 10
          if(any(jj <- z == -Inf))
            z[jj] <- -10
        }
        eta[[j]] <- z
      }
      return(eta)
    }
  }

  if(is.null(family$mu.eta)) {
    family$mu.eta <- mu.eta
  }

  if(is.null(family$mean)) {
    family$mean <- function(par) { par[[1]] }
  }

  if(is.null(family$log_likelihood)) {
    if(!is.null(family$pdf)) {
      family$log_likelihood <- function(par, y, ...) {
        logdens <- try(family$pdf(par = par, y = y, log = TRUE), silent = TRUE)
        if(inherits(logdens, "try-error")) {
          warning("problems evaluating the log-density of the model, set log-likelihood to -Inf")
          return(-Inf)
        }
#        if(any(is.na(logdens))) {
#          warning("NA log-density values!")
#        }
        if(any(i <- !is.finite(logdens))) {
          ## warning("non finite log-density values, set to -100!")
          logdens[i] <- -100
        }
        return(sum(logdens, na.rm = TRUE))
      }
    } else {
      stop("the family object does not have a $pdf() function!")
    }
  }

  err01 <- .Machine$double.eps^(1/3)
  err11 <- .Machine$double.eps^(1/4)

  if(is.null(family$score) && !is.null(family$pdf))
    family$score <- list()
  for(parameter in family$names) {
    if(is.null(family$score[[parameter]]) && !is.null(family$pdf)) {
      family$score[[parameter]] <- make_numeric_score(
        parameter = parameter,
        pdf = family$pdf,
        linkfun = linkfun[[parameter]],
        linkinv = linkinv[[parameter]],
        step = err01
      )
      attr(family$score[[parameter]], "dnum") <- TRUE
    }
  }

  if(is.null(family[["hessian"]]) && !is.null(family[["hess"]])) {
    family[["hessian"]] <- family[["hess"]]
    family[["hess"]] <- NULL
  }

  if(is.null(family$hessian) && !is.null(family$pdf))
    family$hessian <- list()
  for(parameter in family$names) {
    if(is.null(family$hessian[[parameter]]) && !is.null(family$pdf)) {
      score <- family$score[[parameter]]
      step <- if(isTRUE(attr(score, "dnum"))) err11 else err01
      family$hessian[[parameter]] <- make_numeric_hessian(
        parameter = parameter,
        score = score,
        linkfun = linkfun[[parameter]],
        linkinv = linkinv[[parameter]],
        step = step
      )
    }
  }
  for(i in seq_along(family$names)) {
    for(j in seq_along(family$names)) {
      if(i < j) {
        hij <- paste0(family$names[i], ":", family$names[j])
        if(is.null(family$hessian[[hij]])) {
          ni <- family$names[i]
          nj <- family$names[j]

          family$hessian[[hij]] <- make_numeric_hessian(
            parameter = ni,
            score = family$score[[nj]],
            linkfun = linkfun[[ni]],
            linkinv = linkinv[[ni]],
            step = err01
          )
        }
        hji <- paste0(family$names[j], ":", family$names[i])
        if(is.null(family$hessian[[hji]]))
          family$hessian[[hji]] <- family$hessian[[hij]]
      }
    }
  }

  if(is.null(family$type)) {
    family$type <- "continuous"
  }

  class(family) <- "gamlss2.family"

  return(family)
}

family.gamlss2.family <- function(object, ...) {
  complete_family(object, ...)
}

## A simple print method.
print.gamlss2.family <- function(x, full = TRUE, ...)
{
  cat("Family:", x$family, if(!is.null(x$full.name)) paste0("(", x$full.name, ")") else NULL,  "\n")
  if(!is.character(x$links)) {
    if(inherits(x$links, c("link-gamlss2", "link-glm")))
      links <- x$links$name
    else
      links <- sapply(x$links, function(x) {
        if(is.character(x)) return(x) else return(x$name)
      })
  } else {
    links <- x$links
  }
  links <- paste(links, collapse = ", ")
  if(links != "") {
    cat(if(length(x$links) > 1) "Link functions:" else "Link function:", links, sep = " ")
    cat("\n")
  }
  if(full) {
    nfun <- names(x[c("transform", "optimizer", "sampler", "results", "predict")])
    if(!all(is.na(nfun))) {
      nfun <- nfun[!is.na(nfun)]
      cat("---\nFamily specific functions:\n")
      for(j in nfun)
        cat(" ..$ ", j, "\n", sep = "")
    }
    nfun <- names(x[c("score", "hessian")])
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
  return(invisible(x))
}

## Some example families.
Gaussian <- function(...)
{
  links <- c(mu = "identity", sigma = "log")

  rval <- list(
    "family" = "Gaussian",
    "names" = c("mu", "sigma"),
    "links" = parse_links(links, c(mu = "identity", sigma = "log"), ...),
    "score" = list(
      "mu" = function(par, y, ...) { drop((y - par$mu) / (par$sigma^2)) },
      "sigma" = function(par, y, ...) { drop(-1 + (y - par$mu)^2 / (par$sigma^2)) }
    ),
    "hessian" = list(
      "mu" = function(par, y, ...) { drop(1 / (par$sigma^2)) },
      "sigma" = function(par, y, ...) { rep(2, length(y)) }
    ),
    "log_likelihood" = function(par, y, ...) {
      sum(dnorm(y, par$mu, par$sigma, log = TRUE))
    },
    "mu" = function(par, ...) {
      par$mu
    },
    "pdf" = function(par, y, log = FALSE) {
      dnorm(y, mean = par$mu, sd = par$sigma, log = log)
    },
    "cdf" = function(par, y, ...) {
      pnorm(y, mean = par$mu, sd = par$sigma, ...)
    },
    "random" = function(n, par) {
      rnorm(n, mean = par$mu, sd = par$sigma)
    },
    "quantile" = function(par, p) {
      qnorm(p, mean = par$mu, sd = par$sigma)
    },
    "crps" = function(par, y, ...) {
      sum(scoringRules::crps_norm(y, mean = par$mu, sd = par$sigma), na.rm = TRUE)
    },
    "initialize" = list(
      "mu"    = function(y, ...) { (y + mean(y)) / 2 },
      "sigma" = function(y, ...) { rep(sd(y), length(y)) }
    ),
    "mean"      = function(par) par$mu,
    "variance"  = function(par) par$sigma^2,
    "skewness" = function(par) { rep(0, length(par$mu)) },
    "kurtosis" = function(par) { rep(3, length(par$mu)) },
    "valid.response" = function(x) {
      if(is.factor(x) | is.character(x))
        stop("the response should be numeric!")
      return(TRUE)
    }
  )

#  rval$update <- function(par, y, eta, which) {
#    score <- deriv_checks(rval$score[[which]](par = par, y = y, id = which), is.weight = FALSE)
#    hessian <- deriv_checks(rval$hessian[[which]](par = par, y = y, id = which), is.weight = TRUE)
#    z <- eta + 1 / hessian * score
#    return(list("eta" = z, "weights" = hessian))
#  }

  rval$update <- update_Gaussian
  rval$type <- "continuous"

  class(rval) <- "gamlss2.family"
  rval
}

update_Gaussian <- function(par, y, eta, which) {
  .Call("update_Gaussian", par, y, eta, which, PACKAGE = "gamlss2")
}

Weibull <- function(...)
{
  rval <- list(
    "family" = "Weibull",
    "names" = c("mu", "sigma"),
    "links" = c(mu = "identity", sigma = "log"),
    "pdf" = function(par, y, log = FALSE, ...) {
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
    "cdf" = function(par, y, ...) {
      delta <- y[, "status"]
      y <- log(y[, "time"])
      p1 <- 1 - exp(-exp((y - par$mu) / par$sigma))
      p2 <- runif(length(y), p1, 1)
      prob <- ifelse(delta > 0, p1, p2)
      return(prob)
    },
    "quantile" = function(par, p, ...) {
      lambda <- exp(-par$mu/par$sigma)
      alpha <- 1 / par$sigma
      q <- lambda * (-log(1 - p))^(1 / alpha)
      return(q)
    },
    "score" = list(
      "mu" = function(par, y, ...) {
        delta <- y[, "status"]
        y <- log(y[, "time"])
        eyms <- exp((y - par$mu)/par$sigma)
        s1 <- 1/par$sigma
        eymss1 <- eyms * s1
        a <- -(delta * (s1 - eymss1))
        b <- (1 - delta) * (eymss1)
        return(a + b)
      },
      "sigma" = function(par, y, ...) {
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
    "hessian" <- list(
      "mu" = function(par, y, ...) {
        delta <- y[, "status"]
        y <- log(y[, "time"])
        eyms <- exp((y - par$mu)/par$sigma) * 1/par$sigma^2
        a <- -(delta * eyms)
        b <- -((1 - delta) * eyms)
        return(-(a + b))
      },
      "sigma" = function(par, y, ...) {
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
    "pdf" = function(par, y, log = FALSE, ...) {
      psi <- YJt(y, par$lambda)
      d <- -0.918938533204675 - log(par$sigma) - 0.5 * ((psi - par$mu)/par$sigma)^2 +
        (par$lambda - 1) * sign(y) * log1p(abs(y))
      if(!log)
        d <- exp(d)
      return(d)
    },
    "cdf" = function(par, y) {
      psi <- YJt(y, par$lambda)
      pnorm(psi, mean = par$mu, sd = par$sigma)
    },
    "quantile" = function(par, p) {
      q <- qnorm(p, mean = par$mu, sd = par$sigma)
      YJt(q, par$lambda, inverse = TRUE)
    },
    "score" = list(
      "mu" = function(par, y, ...) {
        psi <- YJt(y, par$lambda)
        (psi - par$mu)/(par$sigma^2)
      },
      "sigma" = function(par, y, ...) {
        psi <- YJt(y, par$lambda)
        -1/par$sigma + (psi - par$mu)^2/(par$sigma^3)
      },
      "lambda" = function(par, y, ...) {
        psi <- YJt(y, par$lambda)
        -(psi - par$mu) / (par$sigma^2) * YJt(y, par$lambda, derivative = 1) + sign(y) * log1p(abs(y))
      }
    ),
    "hessian" = list(
      "mu" = function(par, y, ...) {
        1 / par$sigma^2
      },
      "sigma" = function(par, y, ...) {
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
    "pdf" = function(par, y, log = FALSE, ...) {
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

OL <- function(k) {
  if(length(k) != 1L || !is.numeric(k) || is.na(k) || !is.finite(k) ||
      k < 2 || k != floor(k) || k > .Machine$integer.max) {
    stop("argument k must be a single integer greater than or equal to 2.")
  }
  k <- as.integer(k)

  ## Parameter names: location and delta-encoded cutpoints.
  delta_indices <- if(k > 2L) seq.int(2L, k - 1L) else integer(0L)
  threshold_names <- c(
    "theta1", if(length(delta_indices)) paste0("delta", delta_indices)
  )
  par_names <- c("location", threshold_names)
  categories <- seq_len(k)
  category_levels <- as.character(categories)

  ## Identity links for now.
  links <- rep("identity", length(par_names))
  names(links) <- par_names

  as_categories <- function(y, argument = "y") {
    if(is.factor(y)) {
      if(!identical(levels(y), category_levels)) {
        stop(
          argument, " must have factor levels ",
          paste(category_levels, collapse = ", "), "."
        )
      }
      y <- as.integer(y)
    } else {
      if(!is.numeric(y))
        stop(argument, " must be numeric or a factor.")
      if(anyNA(y) || any(!is.finite(y)))
        stop(argument, " must not contain missing or non-finite values.")
      if(any(y != floor(y)))
        stop(argument, " must contain integer category values.")
      y <- as.integer(y)
    }

    if(any(y < 1L | y > k))
      stop(argument, " must contain category values in 1, ..., ", k, ".")

    y
  }

  compute_components <- function(par, n = NULL) {
    missing_parameters <- setdiff(par_names, names(par))
    if(length(missing_parameters)) {
      stop(
        "missing distribution parameter(s): ",
        paste(missing_parameters, collapse = ", "), "."
      )
    }

    parameter_lengths <- vapply(par[par_names], length, integer(1L))
    if(any(parameter_lengths < 1L))
      stop("distribution parameters must not be empty.")

    if(is.null(n))
      n <- max(parameter_lengths)
    else
      n <- max(c(n, parameter_lengths))

    if(any(!parameter_lengths %in% c(1L, n)))
      stop("distribution parameters must have length 1 or a common length.")

    ## Build increasing cutpoints.
    cuts <- matrix(NA_real_, nrow = n, ncol = k - 1L)
    cuts[, 1L] <- rep_len(par$theta1, n)
    if(k > 2L) {
      for(j in delta_indices) {
        cuts[, j] <- cuts[, j - 1L] +
          exp(rep_len(par[[paste0("delta", j)]], n))
      }
    }

    ## Cumulative probabilities: c_j = P(Y > j).
    cum_probs <- do.call(
      cbind,
      lapply(seq_len(k - 1L), function(j) {
        plogis(rep_len(par$location, n) - cuts[, j])
      })
    )

    ## Category probabilities.
    probs <- matrix(NA_real_, nrow = n, ncol = k)
    probs[, 1L] <- 1 - cum_probs[, 1L]
    if(k > 2L) {
      for(j in delta_indices) {
        probs[, j] <- cum_probs[, j - 1L] - cum_probs[, j]
      }
    }
    probs[, k] <- cum_probs[, k - 1L]

    list(
      cuts = cuts,
      cum_probs = cum_probs,
      probs = probs
    )
  }

  initial_cutpoints <- function(y) {
    y <- as_categories(y, "response")
    if(!length(y))
      stop("response must not be empty.")

    counts <- tabulate(y, nbins = k)
    n <- sum(counts)
    cumulative_probs <- cumsum(counts)[seq_len(k - 1L)] / n

    ## Keep empirical cumulative probabilities away from 0 and 1 so
    ## initial cutpoints remain finite when boundary categories are empty.
    eps <- 0.5 / (n + 1)
    cumulative_probs <- pmin(pmax(cumulative_probs, eps), 1 - eps)
    qlogis(cumulative_probs)
  }

  fam <- list(
    family = paste0("Ordered Logit (", k, " categories)"),
    names = par_names,
    links = links,

    pdf = function(par, y, log = FALSE, ...) {
      comps <- compute_components(par, n = length(y))
      n <- nrow(comps$probs)
      if(length(y) == 1L && n > 1L)
        y <- rep(y, n)
      if(length(y) != n)
        stop("y must have length 1 or match the distribution parameters.")
      y_int <- as_categories(y)

      p <- comps$probs[cbind(seq_len(n), y_int)]
      p[p < 1e-8 | is.na(p)] <- 1e-8

      if(log) {
        p <- log(p)
        p[is.na(p)] <- -1e10
      }

      p
    },

    initialize = {
      init_list <- list()

      ## Start with zero latent location so empirical cumulative logits
      ## directly initialize the cutpoints.
      init_list$location <- function(y, ...) {
        y <- as_categories(y, "response")
        rep(0, length(y))
      }

      ## Initialize theta1.
      init_list$theta1 <- function(y, ...) {
        cuts <- initial_cutpoints(y)
        rep(cuts[1L], length(y))
      }

      ## Initialize deltas as log-spacings between adjacent cutpoints.
      for(j in delta_indices) {
        init_list[[paste0("delta", j)]] <- local({
          jj <- j
          function(y, ...) {
            cuts <- initial_cutpoints(y)
            val <- log(max(cuts[jj] - cuts[jj - 1L], 1e-4))
            rep(val, length(y))
          }
        })
      }

      init_list
    }
  )

  ## Probabilities on response scale.
  fam$probabilities <- function(par, ...) {
    probs <- compute_components(par)$probs
    colnames(probs) <- paste0("Pr(Y=", categories, ")")
    probs
  }

  ## Moments of the numeric category labels on the response scale.
  fam$mean <- function(par, ...) {
    probs <- fam$probabilities(par)
    drop(probs %*% categories)
  }

  fam$variance <- function(par, ...) {
    probs <- fam$probabilities(par)
    ey <- drop(probs %*% categories)
    drop(probs %*% categories^2) - ey^2
  }

  fam$transition <- function(par, ...) {
    tp <- compute_components(par)$cum_probs
    colnames(tp) <- paste0("Pr(Y>", seq_len(k - 1L), ")")
    tp
  }

  fam$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
    probs <- compute_components(par, n = length(y))$probs
    n <- nrow(probs)

    if(length(y) == 1L && n > 1L)
      y <- rep(y, n)
    if(length(y) != n)
      stop("y must have length 1 or match the distribution parameters.")
    y_int <- as_categories(y)

    cprobs <- t(apply(probs, 1L, cumsum))
    ans <- cprobs[cbind(seq_len(n), y_int)]

    if(!lower.tail)
      ans <- 1 - ans

    if(log.p)
      ans <- log(ans)

    ans
  }

  fam$quantile <- function(par, p, ...) {
    if(!is.numeric(p))
      stop("p must be numeric.")

    probs <- compute_components(par, n = length(p))$probs
    n <- nrow(probs)

    if(length(p) == 1L && n > 1L)
      p <- rep.int(p, n)

    if(length(p) != n)
      stop("p must have length 1 or match the distribution parameters.")
    if(anyNA(p) || any(!is.finite(p)))
      stop("p must not contain missing or non-finite values.")
    if(any(p < 0 | p > 1))
      stop("p must be in [0, 1].")

    cprobs <- t(apply(probs, 1L, cumsum))

    q <- integer(n)
    for(i in seq_len(n)) {
      idx <- which(cprobs[i, ] >= p[i])[1L]
      if(is.na(idx))
        idx <- k
      q[i] <- idx
    }

    q
  }

  fam$log_likelihood <- function(par, y, ...) {
    sum(fam$pdf(par, y, log = TRUE))
  }

  fam$valid.response <- function(x) {
    as_categories(x, "response")
    TRUE
  }

  fam$residuals <- function(object, ...) {
    rqres_ologit(object, ...)
  }

  fam$type <- "discrete"

  class(fam) <- c("gamlss2.family", "family.bamlss")
  fam
}

rqres_ologit <- function(object, ...) {
  fam <- family(object)
  if(!grepl("^Ordered Logit", fam$family))
    stop("OL family required.")

  mf <- model.frame(object)
  y <- stats::model.response(mf)

  par <- predict(object)
  probs <- fam$probabilities(par)

  K <- ncol(probs)
  y_int <- if(is.factor(y)) {
    match(as.character(y), as.character(seq_len(K)))
  } else {
    as.integer(y)
  }
  n <- length(y_int)

  if(any(y_int < 1L | y_int > K | is.na(y_int))) {
    bad_vals <- sort(unique(y[is.na(y_int) | y_int < 1L | y_int > K]))
    stop(
      "response has values outside 1..", K, " implied by OL(k).\n",
      "Offending values: ", paste(bad_vals, collapse = ", ")
    )
  }

  cprobs <- t(apply(probs, 1L, cumsum))
  F_upper <- cprobs[cbind(seq_len(n), y_int)]

  F_lower <- numeric(n)
  idx_gt1 <- which(y_int > 1L)
  if(length(idx_gt1) > 0L) {
    F_lower[idx_gt1] <- cprobs[cbind(idx_gt1, y_int[idx_gt1] - 1L)]
  }

  u <- stats::runif(n, min = F_lower, max = F_upper)
  u[u <= 0] <- 1e-12
  u[u >= 1] <- 1 - 1e-12
  stats::qnorm(u)
}

ologit <- function(k) {
  OL(k = k)
}

#if(FALSE) {
#library("gamlss2")

### From MASS.
#library("MASS")

#options(contrasts = c("contr.treatment", "contr.poly"))

#m <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
#summary(m)

### Response needs to be integer.
#housing$Satint <- as.integer(housing$Sat)

### Estimate model.
#b <- gamlss2(Satint ~ Infl + Type + Cont, data = housing, weights = Freq, family = ologitK(k = 3))

### Compare.
#coef(m)
#coef(b)

### Predict probabilities.
#pm <- predict(m, type = "p")
#pb <- predict(b)
#pb <- family(b)$probabilities(pb)

#print(head(pm))
#print(head(pb))
#}

## Shifted log-link.
shiftlog <- function(shift = 1) {
  if(length(shift) != 1L || !is.numeric(shift) ||
      is.na(shift) || !is.finite(shift)) {
    stop("argument shift must be a single finite number.")
  }

  linkfun <- function(mu) {
    if(!is.numeric(mu) || anyNA(mu) || any(!is.finite(mu)) ||
        any(mu <= shift)) {
      stop("values must be finite and greater than shift.")
    }
    log(mu - shift)
  }
  linkinv <- function(eta) exp(eta) + shift
  mu.eta <- function(eta) exp(eta)
  valideta <- function(eta) all(is.finite(eta))
  validmu <- function(mu) {
    is.numeric(mu) && all(is.finite(mu)) && all(mu > shift)
  }

  shift_label <- if(shift == 0) {
    ""
  } else {
    paste0(if(shift > 0) " + " else " - ", format(abs(shift)))
  }

  structure(
    list(
      linkfun = linkfun,
      linkinv = linkinv,
      mu.eta = mu.eta,
      mu.eta2 = mu.eta,
      valideta = valideta,
      validmu = validmu,
      name = paste0("exp(x)", shift_label)
    ),
    class = "link-glm"
  )
}

## Kumaraswamy distribution.
Kumaraswamy <- KS <- function(a.link = shiftlog, b.link = shiftlog, ...) {
  lfa <- make.link2(a.link)
  lfb <- make.link2(b.link)

  prepare_parameters <- function(par, n = NULL) {
    if(!is.list(par))
      stop("par must be a list-like object containing a and b.")

    a_value <- par$a
    b_value <- par$b
    missing_parameters <- c(
      if(is.null(a_value)) "a",
      if(is.null(b_value)) "b"
    )
    if(length(missing_parameters)) {
      stop(
        "missing distribution parameter(s): ",
        paste(missing_parameters, collapse = ", "), "."
      )
    }
    if(!is.numeric(a_value) || !is.numeric(b_value))
      stop("parameters a and b must be numeric.")

    parameter_lengths <- c(a = length(a_value), b = length(b_value))
    if(any(parameter_lengths < 1L))
      stop("parameters a and b must not be empty.")

    if(is.null(n))
      n <- max(parameter_lengths)
    else
      n <- max(c(n, parameter_lengths))

    if(any(!parameter_lengths %in% c(1L, n)))
      stop("parameters a and b must have length 1 or a common length.")

    a <- rep_len(as.numeric(a_value), n)
    b <- rep_len(as.numeric(b_value), n)
    if(anyNA(a) || anyNA(b) || any(!is.finite(a)) || any(!is.finite(b)) ||
        any(a <= 0) || any(b <= 0)) {
      stop("parameters a and b must be finite and strictly positive.")
    }

    list(a = a, b = b, n = n)
  }

  recycle_argument <- function(x, n, argument) {
    if(!is.numeric(x))
      stop(argument, " must be numeric.")
    if(length(x) == 1L && n > 1L)
      x <- rep(x, n)
    if(length(x) != n)
      stop(argument, " must have length 1 or match the distribution parameters.")
    as.numeric(x)
  }

  prepare_scores <- function(par, y) {
    par <- prepare_parameters(par, n = length(y))
    y <- recycle_argument(y, par$n, "y")
    if(anyNA(y) || any(!is.finite(y)) || any(y <= 0) || any(y >= 1))
      stop("y must contain finite values strictly between 0 and 1.")
    list(par = par, y = y)
  }

  fam <- list(
    "family" = "Kumaraswamy",
    "names" = c("a", "b"),
    "links" = list("a" = a.link, "b" = b.link),

    "pdf" = function(par, y, log = FALSE, ...) {
      par <- prepare_parameters(par, n = length(y))
      y <- recycle_argument(y, par$n, "y")
      logd <- rep(-Inf, par$n)
      logd[is.na(y)] <- NA_real_

      inside <- !is.na(y) & is.finite(y) & y > 0 & y < 1
      if(any(inside)) {
        ly <- log(y[inside])
        log_one_minus_ya <- log(-expm1(par$a[inside] * ly))
        logd[inside] <- log(par$a[inside]) + log(par$b[inside]) +
          (par$a[inside] - 1) * ly +
          (par$b[inside] - 1) * log_one_minus_ya
      }

      at_zero <- !is.na(y) & y == 0
      if(any(at_zero)) {
        logd[at_zero & par$a < 1] <- Inf
        logd[at_zero & par$a == 1] <- log(par$b[at_zero & par$a == 1])
      }

      at_one <- !is.na(y) & y == 1
      if(any(at_one)) {
        logd[at_one & par$b < 1] <- Inf
        logd[at_one & par$b == 1] <- log(par$a[at_one & par$b == 1])
      }

      if(log) logd else exp(logd)
    },

    "score" = list(
      "a" = function(par, y, ...) {
        z <- prepare_scores(par, y)
        par <- z$par
        ly <- log(z$y)
        ya <- exp(par$a * ly)
        one_minus_ya <- -expm1(par$a * ly)
        score <- 1/par$a + ly -
          (par$b - 1) * ya * ly/one_minus_ya
        score * lfa$mu.eta(lfa$linkfun(par$a))
      },
      "b" = function(par, y, ...) {
        z <- prepare_scores(par, y)
        par <- z$par
        log_one_minus_ya <- log(-expm1(par$a * log(z$y)))
        score <- 1/par$b + log_one_minus_ya
        score * lfb$mu.eta(lfb$linkfun(par$b))
      }
    ),

    "hessian" = list(
      "a" = function(par, y, ...) {
        z <- prepare_scores(par, y)
        par <- z$par
        ly <- log(z$y)
        ya <- exp(par$a * ly)
        one_minus_ya <- -expm1(par$a * ly)
        hessian <- 1/par$a^2 +
          (par$b - 1) * ya * ly^2/one_minus_ya^2
        hessian * lfa$mu.eta(lfa$linkfun(par$a))^2
      },
      "b" = function(par, y, ...) {
        z <- prepare_scores(par, y)
        par <- z$par
        1/par$b^2 * lfb$mu.eta(lfb$linkfun(par$b))^2
      }
    ),

    "cdf" = function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
      par <- prepare_parameters(par, n = length(y))
      y <- recycle_argument(y, par$n, "y")
      p <- rep(NA_real_, par$n)

      not_missing <- !is.na(y)
      p[not_missing & y <= 0] <- if(lower.tail) 0 else 1
      p[not_missing & y >= 1] <- if(lower.tail) 1 else 0

      inside <- not_missing & y > 0 & y < 1
      if(any(inside)) {
        log_survival <- par$b[inside] *
          log(-expm1(par$a[inside] * log(y[inside])))
        p[inside] <- if(lower.tail) -expm1(log_survival) else exp(log_survival)
      }

      if(log.p)
        p <- log(p)
      p
    },

    "quantile" = function(par, p, lower.tail = TRUE, log.p = FALSE, ...) {
      if(!is.numeric(p))
        stop("p must be numeric.")

      par <- prepare_parameters(par, n = length(p))
      p <- recycle_argument(p, par$n, "p")
      if(anyNA(p))
        stop("p must not contain missing values.")

      if(log.p) {
        if(any(p > 0))
          stop("log probabilities must not be greater than 0.")
        p <- exp(p)
      } else if(any(!is.finite(p))) {
        stop("p must contain finite probabilities.")
      }

      if(!lower.tail)
        p <- 1 - p
      if(any(p < 0 | p > 1))
        stop("p must be in [0, 1].")

      inner <- -expm1(log1p(-p)/par$b)
      inner^(1/par$a)
    },

    "random" = function(par, n) {
      if(length(n) != 1L || !is.numeric(n) || is.na(n) ||
          !is.finite(n) || n < 0 || n != floor(n) ||
          n > .Machine$integer.max) {
        stop("n must be a single non-negative integer.")
      }
      n <- as.integer(n)
      if(n == 0L)
        return(numeric(0L))

      par <- prepare_parameters(par)
      draws <- matrix(runif(par$n * n), nrow = par$n, ncol = n)
      parameters <- list(a = par$a, b = par$b)
      for(j in seq_len(n))
        draws[, j] <- fam$quantile(parameters, draws[, j])

      if(par$n == 1L)
        return(as.vector(draws))
      if(n == 1L)
        return(draws[, 1L])

      colnames(draws) <- paste0("r_", seq_len(n))
      draws
    },

    "mean" = function(par, ...) {
      par <- prepare_parameters(par)
      exp(log(par$b) + lbeta(1 + 1/par$a, par$b))
    },

    "variance" = function(par, ...) {
      par <- prepare_parameters(par)
      mean <- exp(log(par$b) + lbeta(1 + 1/par$a, par$b))
      second_moment <- exp(log(par$b) + lbeta(1 + 2/par$a, par$b))
      pmax(second_moment - mean^2, 0)
    },

    "mode" = function(par, ...) {
      par <- prepare_parameters(par)
      mode <- rep(NA_real_, par$n)
      interior <- par$a > 1 & par$b > 1
      mode[interior] <- (
        (par$a[interior] - 1) /
          (par$a[interior] * par$b[interior] - 1)
      )^(1/par$a[interior])
      mode
    },

    "valid.response" = function(x) {
      if(!is.numeric(x))
        stop("the response must be numeric.")
      if(!length(x) || anyNA(x) || any(!is.finite(x)) ||
          any(x <= 0) || any(x >= 1)) {
        stop("the response must contain finite values strictly between 0 and 1.")
      }
      TRUE
    },

    "type" = "continuous"
  )

  class(fam) <- "gamlss2.family"
  fam
}

## Compatibility wrapper for the former internal duplicate.
LKS <- function(a.link = shiftlog, b.link = shiftlog, ...) {
  Kumaraswamy(a.link = a.link, b.link = b.link, ...)
}

discretize <- function(family = NO) {
  if(is.function(family))
    family <- family()

  if(inherits(family, "gamlss.family"))
    family <- tF(family)

  fam <- list(
    "family" = paste("discretized", family$family),
    "names" = family$names,
    "links" = family$links,
    "valid.response" = function(x) {
      if(is.factor(x))
        return(FALSE)
      if(!(ok <- all(x >= 0)))
        stop("response values smaller than 0 not allowed!", call. = FALSE)
      ok
    }
  )

  fam$pdf <- function(par, y, log = FALSE, ...) {
    n <- length(y)

    par <- lapply(par, function(x) rep(x, length.out = n))
    par <- as.data.frame(par)

    F0 <- family$cdf(par = par, y = rep(0, n), ...)
    S0 <- 1 - F0

    d <- family$cdf(par = par, y = y + 1, ...) - family$cdf(par = par, y = y, ...)
    d <- d / S0

    if(log)
      d <- log(d)

    d
  }

  fam$cdf <- function(par, y, log = FALSE, ...) {
    par <- as.data.frame(par)

    np <- nrow(par)
    ny <- length(y)
    n <- max(ny, np)

    y <- rep(y, length.out = n)

    par <- lapply(par, function(x) rep(x, length.out = n))
    par <- as.data.frame(par)

    yy <- floor(y)
    yy[yy < 0] <- -1

    F0 <- family$cdf(par = par, y = rep(0, n), ...)
    S0 <- 1 - F0

    p <- numeric(n)

    ii <- yy >= 0
    p[ii] <- (
      family$cdf(par = par[ii, , drop = FALSE], y = yy[ii] + 1, ...) -
        F0[ii]
    ) / S0[ii]

    p <- pmin(pmax(p, 0), 1)

    if(log)
      p <- log(p)

    p
  }

  fam$quantile <- function(par, p, ...) {
    par <- as.data.frame(par)

    np <- nrow(par)
    n <- max(length(p), np)

    p <- rep(p, length.out = n)

    par <- lapply(par, function(x) rep(x, length.out = n))
    par <- as.data.frame(par)

    if(any(is.na(p)))
      stop("p must not contain NA.", call. = FALSE)

    if(any(p < 0 | p > 1))
      stop("p must be in [0, 1].", call. = FALSE)

    F0 <- family$cdf(par = par, y = rep(0, n), ...)
    S0 <- 1 - F0

    pp <- F0 + p * S0

    pp[pp <= 0] <- 0
    pp[pp >= 1] <- 1

    qc <- family$quantile(par, pp, ...)

    q <- ceiling(qc) - 1

    q[p <= 0] <- 0
    q[q < 0] <- 0

    q
  }

  fam$type <- "discrete"
  class(fam) <- "gamlss2.family"

  fam
}

MN <- function(k)
{
  stopifnot(k >= 2)

  pn <- paste0("pi", 2:k)
  links <- rep("log", k - 1)
  names(links) <- pn

  rval <- list(
    family = "Multinomial Logit",
    names  = pn,
    links  = links,

    valid.response = function(x) {
      if(!is.factor(x))
        stop("the response must be a factor!")
      if(nlevels(x) != k)
        stop("number of levels of the response and argument k differ in MN()!")
      TRUE
    },

    pdf = function(par, y, log = FALSE) {
      y_int <- as.integer(y)
      w <- do.call("cbind", par)
      denom <- 1 + rowSums(w)

      logp <- numeric(length(y_int))
      is_ref <- (y_int == 1L)
      logp[is_ref] <- -log(denom[is_ref])

      if(any(!is_ref)) {
        jj <- y_int[!is_ref] - 1L
        wsub <- w[!is_ref, , drop = FALSE]
        logp[!is_ref] <- log(wsub[cbind(seq_len(nrow(wsub)), jj)]) - log(denom[!is_ref])
      }

      if(!log) exp(logp) else logp
    },

    log_likelihood = function(par, y, ...) sum(rval$pdf(par, y, log = TRUE), na.rm = TRUE),

    type = "discrete"
  )

  rval$score <- setNames(vector("list", k - 1), pn)
  for(j in seq_len(k - 1)) {
    id <- pn[j]
    rval$score[[id]] <- local({
      jj <- j
      idd <- id
      function(par, y, ...) {
        y_int <- as.integer(y)
        w <- do.call("cbind", par)
        denom <- 1 + rowSums(w)
        p_j <- par[[idd]] / denom
        as.numeric(y_int == (jj + 1L)) - p_j
      }
    })
  }

  rval$hessian <- list()
  for(j in seq_len(k - 1)) {
    idj <- pn[j]

    rval$hessian[[idj]] <- local({
      idd <- idj
      function(par, y, ...) {
        w <- do.call("cbind", par)
        denom <- 1 + rowSums(w)
        p_j <- par[[idd]] / denom
        p_j * (1 - p_j)
      }
    })

    for(m in seq_len(k - 1)) if(m != j) {
      idm <- pn[m]
      nm  <- paste0(idj, ".", idm)

      rval$hessian[[nm]] <- local({
        idd_j <- idj
        idd_m <- idm
        function(par, y, ...) {
          w <- do.call("cbind", par)
          denom <- 1 + rowSums(w)
          p_j <- par[[idd_j]] / denom
          p_m <- par[[idd_m]] / denom
          -1 * p_j * p_m
        }
      })
    }
  }

  rval$probabilities <- function(par, numeric = TRUE, ...) {
    w <- do.call("cbind", par)
    denom <- 1 + rowSums(w)
    p <- cbind(1/denom, w/denom)
    colnames(p) <- c("pi1", names(par))
    as.data.frame(p)
  }

  rval$cdf <- function(par, y, lower.tail = TRUE, log.p = FALSE, ...) {
    probs <- rval$probabilities(par)
    P <- as.matrix(probs)
    n <- nrow(P)
    K <- ncol(P)

    if(length(y) == 1L) y <- rep.int(y, n)
    y_int <- if(is.factor(y)) as.integer(y) else as.integer(y)

    if(anyNA(y_int))
      stop("missing values in y are not allowed in cdf().", call. = FALSE)
    if(any(y_int < 1L | y_int > K)) {
      bad <- sort(unique(y_int[y_int < 1L | y_int > K]))
      stop("y has values outside 1..", K, ". Offending values: ",
        paste(bad, collapse = ", "), call. = FALSE)
    }

    cP <- P
    cP[] <- t(apply(P, 1L, cumsum))

    ans <- cP[cbind(seq_len(n), y_int)]
    if(!lower.tail) ans <- 1 - ans
    if(log.p) ans <- log(ans)
    ans
  }

  rval$quantile <- function(par, p, ...) {
    probs <- rval$probabilities(par)
    P <- as.matrix(probs)
    n <- nrow(P)
    K <- ncol(P)

    if(length(p) == 1L) p <- rep.int(p, n)
    if(length(p) != n)
      stop("length(p) must be 1 or equal to the number of observations.", call. = FALSE)
    if(anyNA(p)) stop("p must not contain NA.", call. = FALSE)
    if(any(p < 0 | p > 1)) stop("p must be in [0, 1].", call. = FALSE)

    cP <- P
    cP[] <- t(apply(P, 1L, cumsum))

    q <- integer(n)
    for(i in seq_len(n)) {
      q[i] <- which(cP[i, ] >= p[i])[1L]
      if(is.na(q[i])) q[i] <- K
    }
    q
  }

  rval$residuals <- function(object, ...) {
    rqres_mn(object, ...)
  }

  class(rval) <- "gamlss2.family"
  rval
}

rqres_mn <- function(object, ...) {
  fam <- family(object)
  if(!identical(fam$family, "Multinomial Logit"))
    stop("MN() family required.", call. = FALSE)

  mf <- model.frame(object)
  y <- stats::model.response(mf)

  if(!is.factor(y))
    stop("response must be a factor for MN().", call. = FALSE)

  y_int <- as.integer(y)
  par <- predict(object)
  probs <- fam$probabilities(par)
  P <- as.matrix(probs)

  n <- length(y_int)
  K <- ncol(P)

  if(length(y_int) != nrow(P))
    stop("length(response) and number of predicted rows differ.", call. = FALSE)
  if(any(y_int < 1L | y_int > K))
    stop("response contains invalid category indices.", call. = FALSE)

  cP <- P
  cP[] <- t(apply(P, 1L, cumsum))

  lower <- numeric(n)
  upper <- numeric(n)
  idx_ref <- (y_int == 1L)
  lower[idx_ref] <- 0
  upper[idx_ref] <- cP[cbind(which(idx_ref), 1L)]

  if(any(!idx_ref)) {
    ii <- which(!idx_ref)
    yy <- y_int[ii]
    lower[ii] <- cP[cbind(ii, yy - 1L)]
    upper[ii] <- cP[cbind(ii, yy)]
  }

  u <- stats::runif(n, min = lower, max = upper)
  stats::qnorm(u)
}

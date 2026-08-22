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
    if(with_log) aa <- c(aa, "log=log")
    if(with_dots) aa <- c(aa, "...")

    ## Add bd only for functions that use a response argument.
    nf <- names(formals(fun))
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
    random = random_fun
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
    
  if(is.function(family)) {
    family <- family()
  }

  if(inherits(family, "gamlss.family")) {
    return(tF(family))
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

  if(is.null(family$pdf))
    stop("the family needs a $pdf() function!")

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
          hji <- paste0(family$names[j], ".", family$names[i])
          family$hessian[[hji]] <- family$hessian[[hij]]
        }
      }
    }
  }

  if(is.null(family$type)) {
    family$type <- "continuous"
  }

  class(family) <- "gamlss2.family"

  return(family)
}

## A simple print method.
print.gamlss2.family <- function(x, full = TRUE, ...)
{
  cat("Family:", x$family, if(!is.null(x$full.name)) paste0("(", x$full.name, ")") else NULL,  "\n")
  if(!is.character(x$links)) {
    if(inherits(x$links, c("link-gamlss2", "link-glm")))
      links <- x$links$name
    else
      links <- sapply(x$links, function(x) x$name)
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

## From stats.
safe_pchisq <- function (q, df, ...)
{
  df[df <= 0] <- NA
  pchisq(q = q, df = df, ...)
}

anova.gamlss2 <- function(object, ..., test = "Chisq") {
  ## Allow logical FALSE or "none" for no test.
  if(is.logical(test) && identical(test, FALSE)) {
    test <- "none"
  }
  test <- match.arg(test, choices = c("Chisq", "none"))

  ## Collect models.
  models <- list(object, ...)

  if(length(models) < 2L)
    return(drop1(models[[1L]]))

  model.names <- as.character(match.call())[-1L]

  if(length(models) < 2) {
    stop("Need at least two models for comparison.")
  }

  ## Extract log-likelihood, degrees of freedom, AIC.
  ll <- sapply(models, logLik)
  df <- sapply(models, function(m) attr(logLik(m), "df"))
  aic <- sapply(models, AIC)

  res <- data.frame(
    Df = df,
    AIC = aic,
    logLik = ll
  )

  if(test == "Chisq") {
    res$LRT <- c(NA, diff(2 * ll))
    res$Df <- c(NA, diff(df))
    res[["Pr(>Chi)"]] <- c(NA, safe_pchisq(res$LRT[-1], df = res$Df[-1], lower.tail = FALSE))
  }

  rownames(res) <- model.names
  class(res) <- c("anova", "data.frame")
  return(res)
}

drop1.gamlss2 <- function(object, scope = NULL, test = c("Chisq", "none"),
  parameter = NULL, ...)
{
  ## Allow logical FALSE or "none" for no test.
  if(is.logical(test) && identical(test, FALSE)) {
    test <- "none"
  }
  test <- match.arg(test, choices = c("Chisq", "none"))

  control <- object$control
  control$trace <- FALSE
  args <- list(...)
  control[names(args)] <- args

  xterms <- xterms0 <- object$xterms
  sterms <- sterms0 <- object$sterms

  ##terms <- object$terms
  xl <- object$xlevels

  if(!is.null(scope)) {
    ## Process variables and special term information.
    Xterms <- list()
    ff <- as.Formula(fake_formula(scope, nospecials = TRUE))
    for(i in 1:(length(ff)[2])) {
      Xterms[[i]] <- attr(terms(formula(ff, rhs = i, lhs = 0)), "term.labels")
      if(length(xl)) {
        for(l in names(xl)) {
          if(l %in% Xterms[[i]]) {
            xil <- as.list(Xterms[[i]])
            names(xil) <- Xterms[[i]]
            xil[[l]] <- paste0(l, xl[[l]])
            Xterms[[i]] <- as.character(unlist(xil))
          }
        }
      }
    }
    Sterms <- fake_formula(scope, onlyspecials = TRUE)
    if(length(Xterms)) {
      names(Xterms) <- names(xterms)[1:length(Xterms)]
    }
    if(length(Sterms)) {
      names(Sterms) <- names(sterms)[1:length(Sterms)]
    }
    xterms <- Xterms
    sterms <- Sterms
  }

  offsets <- model.offset(model.frame(object))
  weights <- model.weights(model.frame(object))

  ll0 <- logLik(object)
  dev0 <- deviance(object)
  df0 <- attr(ll0, "df")

  if(is.null(parameter)) {
    parameter <- list(...)$what
    if(is.null(parameter))
    parameter <- list(...)$model
    if(is.null(parameter))
      parameter <- object$family$names
  }
  if(!is.character(parameter))
    parameter <- object$family$names[parameter]
  parameter <- object$family$names[pmatch(parameter, object$family$names)]

  parameter <- parameter[!is.na(parameter)]
  if(length(parameter) < 1L || all(is.na(parameter)))
    stop("Argument parameter is specified wrong!")

  res <- list()
  if(length(xterms)) {
    for(i in parameter) {
      xi <- xterms[[i]][xterms[[i]] != "(Intercept)"]

      if(length(xi)) {
        j <- xi
        names(j) <- xi
        if(length(xl)) {
          xli <- lapply(1:length(xl), function(i) paste0(names(xl)[i], xl[[i]]))
          names(xli) <- names(xl)
          for(l in names(xli)) {
            for(ll in xli[[l]]) {
              names(j)[j == ll] <- l
            }
          }
        }
        rn <- NULL
        for(l in unique(names(j))) {
          xtl <- xterms0
          xtl[[i]] <- xtl[[i]][!(xtl[[i]] %in% j[names(j) == l])]

          xlab <- l
          if(!is.null(object$xlev)) {
            if(length(object$xlev)) {
              xlev <- sapply(names(object$xlev[[i]]), function(ii) {
                l %in% paste0(ii, object$xlev[[i]][[ii]]) }
              )
              if(any(xlev)) {
                xx <- names(xlev)[xlev]
                xlab <- xx
                xx <- paste0(xx, object$xlev[[i]][[xx]])
                xtl[[i]] <- xtl[[i]][!(xtl[[i]] %in% xx)]
              }
            }
          }

          m <- RS(x = object$x, y = object$y, specials = object$specials,
            family = object$family,
            offsets = offsets, weights = weights, ## start = coef(object),
            xterms = xtl, sterms = object$sterms, control = control)

          ll1 <- m$logLik
          df1 <- get_df2(m)
          dfs <- df0 - df1

          dev <- 2 * (as.numeric(ll0) - as.numeric(ll1))
          pval <- safe_pchisq(dev, dfs, lower.tail = FALSE)

          aod <- data.frame(
            "Df" = dfs,
            "logLik" = ll1,
            "AIC" = -2 * ll1 + 2 * df1,
            "LRT" = dev,
            "Pr(>Chi)" = pval,
            row.names = xlab, check.names = FALSE
          )

          rn <- c(rn, xlab)

          res[[i]] <- rbind(res[[i]], aod)
        }
        res[[i]] <- res[[i]][!duplicated(rn), , drop = FALSE]
      }
    }
  }

  if(length(sterms)) {
    for(i in parameter) {
      si <- sterms[[i]]
      if(length(si)) {
        j <- sterms0
        for(l in si) {
          j[[i]] <- sterms[[i]][sterms[[i]] != l]

          m <- RS(x = object$x, y = object$y, specials = object$specials,
            family = object$family,
            offsets = offsets, weights = weights, ## start = coef(object),
            xterms = object$xterms, sterms = j, control = control)

          ll1 <- m$logLik
          df1 <- get_df2(m)
          dfs <- df0 - df1

          dev <- 2 * (as.numeric(ll0) - as.numeric(ll1))
          pval <- safe_pchisq(dev, dfs, lower.tail = FALSE)

          aod <- data.frame(
            "Df" = dfs,
            "logLik" = ll1,
            "AIC" = -2 * ll1 + 2 * df1,
            "LRT" = dev,
            "Pr(>Chi)" = pval,
            row.names = l, check.names = FALSE
          )

          res[[i]] <- rbind(res[[i]], aod)
        }
      }
    }
  }

  for(i in 1:length(res))
    class(res[[i]]) <- c("anova", "data.frame")

  class(res) <- c("drop1.gamlss2", "list")
  attr(res, "formula") <- formula(object)
  attr(res, "logLik") <- ll0

  return(res)
}

print.drop1.gamlss2 <- function(x, digits = max(getOption("digits") - 2L, 3L), ...)
{
  ll <- attr(x, "logLik")
  aic <- -2 * ll + 2 * attr(ll, "df")

  cat("Model:\n")
  print(attr(x, "formula"))
  cat("df =", round(attr(ll, "df"), digits),
      "logLik =", round(ll, digits),
      "AIC =", round(aic, digits), "\n")
  nx <- names(x)
  for(i in seq_along(nx)) {
    cat("*--------\nParameter:", nx[i], "\n---\n")
    print(x[[i]], digits = digits, signif.legend = i == length(nx), ...)
  }
  return(invisible(NULL))
}

#if(FALSE) {
#  data("rent", package = "gamlss.data")

#  m1 <- gamlss2(R ~ Fl, data = rent)
#  m2 <- gamlss2(R ~ Fl + A, data = rent)
#  m3 <- gamlss2(R ~ s(Fl) + A + loc, data = rent)
#  m4 <- gamlss2(R ~ s(Fl) + A + loc | s(Fl) + s(A) + loc, data = rent)

#  anova(m1, m2, m3, m4)
#  anova(m4)

#  drop1(m4)
#  drop1(m4, scope =~ A + loc | s(Fl))
#  drop1(m4, scope =~ loc)
#  drop1(m4, scope =~ 1 | loc)
#}


## Special conditional inference ctree constructor.
tree <- function(formula, ...)
{
  stopifnot(requireNamespace("rpart"))
  st <- list()
  ctr <- list(...)
  st$control <- do.call(rpart::rpart.control, ctr)
  st$formula <- formula
  st$term <- all.vars(formula)
  st$label <- paste0("tree(", paste0(gsub(" ", "", as.character(formula)), collapse = ""), ")")
  st$data <- model.frame(formula)
  class(st) <- c("special", "tree")
  return(st)
}

## ctree fitting function for the backfitting algorithm.
special_fit.tree <- function(x, z, w, y, eta, j, family, control, ...)
{
  f <- update(x$formula, response_z ~ .)
  x$data$response_z <- z
  x$data$w <- w
  rval <- list(
    "model" = rpart::rpart(formula = f, data = x$data, weights = w,
      control = x$control)
  )
  rval$fitted.values <- predict(rval$model)
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift
  frame <- rval$model$frame
  leaves <- frame$var == "<leaf>"
  size <- sum(leaves)
  edf <- 2 * size - 1
  rval$edf <- edf
  class(rval) <- "tree.fitted"
  return(rval)
}

## A ct predict method.
special_predict.tree.fitted <- function(x, data, ...)
{
  p <- predict(x$model, newdata = data, type = "vector")
  p <- p - x$shift
  return(p)
}

## Special conditional inference forest constructor.
cf <- function(formula, ...)
{
  stopifnot(requireNamespace("partykit"))
  st <- list()
  ctr <- list(...)
  ntree <- ctr$ntree
  if(is.null(ntree))
    ntree <- 100L
  ctr$ntree <- NULL
  st$control <- do.call(partykit::ctree_control, ctr)
  st$formula <- formula
  st$term <- all.vars(formula)
  st$label <- paste0("cf(", paste0(gsub(" ", "", as.character(formula)), collapse = ""), ")")
  st$data <- model.frame(formula)
  st$ntree <- ntree
  class(st) <- c("special", "cf")
  return(st)
}

## The fitting function for the backfitting algorithm.
special_fit.cf <- function(x, z, w, y, eta, j, family, control, ...)
{
  f <- update(x$formula, response_z ~ .)
  x$data$response_z <- z
  rval <- list(
    "model" = partykit::cforest(formula = f, data = x$data, weights = w,
      control = x$control, ntree = x$ntree)
  )
  rval$fitted.values <- predict(rval$model)
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift
  rval$edf <- x$ntree
  class(rval) <- "cf.fitted"
  return(rval)
}

## A predict method.
special_predict.cf.fitted <- function(x, data, se.fit = FALSE, ...)
{
  if(se.fit) {
    ct <- lapply(seq_along(x$model$nodes), function(i) partykit::as.constparty(
      partykit::party(x$model$nodes[[i]], data = x$model$data, terms = x$model$terms,
      fitted = data.frame(
        `(response)` = x$model$fitted[["(response)"]],
        `(weights)` = x$model$weights[[i]],
        check.names = FALSE))
    ))
    p <- sapply(ct, predict, newdata = data)
    p <- apply(p, 1, quantile, prob = c(0.05, 0.5, 0.95))
    p <- t(p)
    colnames(p) <- c("lower", "fit", "upper")
  } else {
    p <- predict(x$model, newdata = data, type = "response")
  }
  p <- p - x$shift
  return(p)
}

## ps2() & pb2() wrapper function using s().
ps <- pb <- function(x, k = 20, ...)
{
  sx <- s(x, bs = "ps", k = k, ...)
  sx$term <- deparse(substitute(x))
  sx$label <- paste0("pb2(", sx$term, ")")
  sx$control <- list("criterion" = "ml")
  sx$localML <- TRUE
  return(sx)
}

## The constructor function is used in the formula
## when calling gamlss2().
n <- function(formula, ...)
{
  stopifnot(requireNamespace("nnet"))

  ## List for setting up the special model term.
  st <- list()

  ## List of control arguments.
  ctr <- list(...)
  if(is.null(ctr$size))
    ctr$size <- 50
  if(is.null(ctr$maxit))
    ctr$maxit <- 1000
  if(is.null(ctr$decay))
    ctr$decay <- 0.1
  if(is.null(ctr$trace))
    ctr$trace <- FALSE
  if(is.null(ctr$MaxNWts))
    ctr$MaxNWts <- 10000
  if(is.null(ctr$scale))
    ctr$scale <- TRUE

  ## Put all information together.
  st$control <- ctr
  st$formula <- formula
  st$term <- all.vars(formula)
  st$label <- paste0("n(", paste0(gsub(" ", "", as.character(formula)), collapse = ""), ")")
  st$data <- model.frame(formula)

  ## Scale per default!
  if(ctr$scale) {
    sx <- list()
    for(j in colnames(st$data)) {
      if(!is.factor(st$data[[j]])) {
        sx[[j]] <- range(st$data[[j]])
        st$data[[j]] <- (st$data[[j]] - sx[[j]][1]) / diff(sx[[j]])
      }
    }
    st$scalex <- sx
  }

  ## Assign the "special" class and the new class "n".
  class(st) <- c("special", "n")

  return(st)
}

## Set up the special "n" model term fitting function
special_fit.n <- function(x, z, w, control, ...)
{
  ## Model formula needs to be updated.
  .fnns <- update(x$formula, response_z ~ .)

  ## Assign current working response.
  x$data$response_z <- z
  x$data$weights_w <- w

  ## Possible weights from last iteration.
  Wts <- list(...)$transfer$Wts

  ## Estimate model.
  nnc <- parse(text = paste0('nnet::nnet(formula = .fnns, data = x$data, weights = weights_w,',
      'size = x$control$size, maxit = x$control$maxit, decay = x$control$decay,',
      'trace = x$control$trace, MaxNWts = x$control$MaxNWts, linout = TRUE',
      if(!is.null(Wts)) ', Wts = Wts)' else ')'))

  rval <- list("model" = eval(nnc))

  ## Get the fitted.values.
  rval$fitted.values <- predict(rval$model)

  ## Transferring the weights for the next backfitting iteration.
  ## Note, "transfer" can be used to transfer anything from one
  ## iteration to the next.
  rval$transfer <- list("Wts" = rval$model$wts)

  ## Center fitted values.
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift

  ## Degrees of freedom.
  rval$edf <- length(coef(rval$model))

  ## Possible scaling.
  rval$scalex <- x$scalex

  ## Assign class for predict method.
  class(rval) <- "n.fitted"

  return(rval)
}

## Finally, the predict method.
special_predict.n.fitted <- function(x, data, se.fit = FALSE, ...)
{
  if(!is.null(x$scalex)) {
    for(j in names(x$scalex)) {
      data[[j]] <- (data[[j]] - x$scalex[[j]][1]) / diff(x$scalex[[j]])
    }
  }
  p <- predict(x$model, newdata = data, type = "raw")
  p <- p - x$shift
  if(se.fit)
    p <- data.frame("fit" = p)
  return(p)
}

## Linear model terms wrapper.
lin <- function(x, ..., ridge = FALSE)
{
  x <- lab <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
  f <- try(as.formula(x), silent = TRUE)
  is_f <- TRUE
  if(inherits(f, "try-error")) {
    v <- as.list(substitute(list(...)))[-1]
    if(length(v)) {
      lab <- paste0(x, ",", as.character(unlist(v)))
      v <- c(x, as.character(unlist(v)))
    } else {
      v <- x
    }
    f <- as.formula(paste("~", paste(v, collapse = "+")))
    is_f <- FALSE
  }
  if(ridge) {
    lab <- paste0(lab,",ridge=", ridge)
  }
  lab <- gsub(" ", "", lab)
  v <- all.vars(f)
  sx <- list("formula" = f, "term" = v,
    "label" = paste0("lin(", lab, ")"),
    "by" = "NA", "dim" = length(v))
  if(!ridge) {
    sx$sp <- 1e-10
  }
  class(sx) <- "lin.smooth.spec"
  return(sx)
}

smooth.construct.lin.smooth.spec <- function(object, data, knots)
{
  object$X <- model.matrix(object$formula, data = if (is.list(data)) 
    data[all.vars(reformulate(names(data))) %in% all.vars(object$formula)]
    else data)
  if(any(grepl("(Intercept)", colnames(object$X), fixed = TRUE))) {
    object$X <- object$X[, -1L, drop = FALSE]
  }
  object$bs.dim <- ncol(object$X)
  object$S <- list(diag(object$bs.dim))
  object$rank <- object$bs.dim
  object$null.space.dim <- 0
  object$C <- matrix(0, 0, ncol(object$X))
  object$side.constrain <- FALSE
  class(object) <- "lin.effect"
  return(object)
}

Predict.matrix.lin.effect <- function(object, data) 
{
  if(is.list(data)) 
    data <- data[all.vars(reformulate(names(data))) %in% all.vars(object$formula)]
  X <- model.matrix(object$formula, model.frame(object$formula, data, na.action = na.pass))
  if(any(grepl("(Intercept)", colnames(X), fixed = TRUE))) {
    X <- X[, -1L, drop = FALSE]
  }
  X[!is.finite(X)] <- 0
  return(X)
}

## Random effects.
re <- function(fixed = ~ 1, random = NULL, ...)
{
  stopifnot(requireNamespace("nlme"))

  call <- match.call(expand.dots = FALSE)

  ## List for setting up the special model term.
  st <- list()

  ## List of control arguments.
  ctr <- list(...)
  if(is.null(ctr$method))
    ctr$method <- "ML"

  ## Put all information together.
  st$control <- ctr$control
  st$fixed <- fixed
  ## Update fixed formula for working response z.
  st$fixed <- update(st$fixed, response_z ~ .)
  st$random <- random
  st$term <- c(all.vars(fixed), all.vars(random))
  st$label <- gsub(" ", "", deparse(call))
  st$label <- gsub("re(", "re2(", st$label, fixed = TRUE)
  st$method <- ctr$method
  ctr$method <- NULL
  st$correlation <- ctr$correlation
  ctr$correlation <- NULL

  ## Assign data.
  rm_v <- function(formula) {
    env <- environment(formula)
    formula <- formula(as.Formula(formula), drop = TRUE, collapse = TRUE)
    environment(formula) <- env
    return(formula)
  }

  df <- model.frame(rm_v(fixed))
  if(inherits(random, "formula"))
    dr <- model.frame(rm_v(random))
  else
    dr <- data.frame()

  if(nrow(df) < 1 && nrow(dr) > 0)
    st$data <- dr
  if(nrow(df) > 0 && nrow(dr) < 1)
    st$data <- df
  if(nrow(df) > 0 && nrow(dr) > 0)
    st$data <- cbind(df, dr)

  ## Assign the "special" class and the new class "n".
  class(st) <- c("special", "re")

  return(st)
}

## Set up the special "re" model term fitting function
special_fit.re <- function(x, z, w, control, ...)
{
  ## Assign current working response.
  x$data$response_z <- z
  x$data$weights_w <- w

  ## Estimate model.
  rem <- parse(text = paste0(
    'nlme::lme(fixed = x$fixed, data = x$data, random = x$random, weights = varFixed(~weights_w),',
    'method="', x$method, '",control=x$control,correlation=x$correlation,keep.data=FALSE)'))

  rval <- list("model" = eval(rem))

  ## Get the fitted.values.
  rval$fitted.values <- fitted(rval$model)

  ## Degrees of freedom.
  N <- sum(w != 0)
  rval$edf <- sum(w * (z - rval$fitted.values)^2)/(rval$model$sigma^2)

  ## Assign class for predict method.
  class(rval) <- "re.fitted"

  return(rval)
}

## The re() predict method.
special_predict.re.fitted <- function(x, data, se.fit = FALSE, ...)
{
  p <- predict(x$model, newdata = data, level = 0)
  if(se.fit)
    p <- data.frame("fit" = p)
  return(p)
}

## Loess smoother.
lo <- function(formula, ...) 
{
  ## Ensure it's a formula.
  if(!inherits(formula, "formula")) {
    formula <- as.character(substitute(formula))
    formula <- as.formula(paste("~", formula))
    environment(formula) <- sys.frame(-1)
  }

  ## List for setting up the special model term. 
  st <- list()

  st$control <- list(...)
  st$term <- all.vars(formula) 
  st$label <- paste0("lo(", paste0(gsub(" ", "",
    as.character(formula)), collapse = ""), ")") 
  st$data <- model.frame(formula)

  ## New model formula used for fitting.
  st$formula <- update(formula, response_z ~ .)

  ## Assign the "special" class and the new class "n".
  class(st) <- c("special", "lo")

  return(st) 
}

special_fit.lo <- function(x, z, w, control, ...)
{
  ## Assign current working response.
  x$data$response_z <- z
  x$data$weights_w <- w

  ## Set up loess call.
  call <- "loess(formula = x$formula, data = x$data, weights = weights_w"

  ## Add optional control parameters.
  if(!is.null(x$control)) {
    for(j in names(x$control))
      call <- paste0(call, ", ", j, "= x$control$", j)
  }

  call <- paste0(call, ")")

  ## Estimate model.
  rval <- list("model" = eval(parse(text = call)))

  ## Get the fitted.values.
  rval$fitted.values <- fitted(rval$model) 

  ## Center fitted values. 
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift 

  ## Degrees of freedom.
  rval$edf <-  rval$model$trace.hat

  ## Assign class for predict method. 
  class(rval) <- "lo.fitted" 

  return(rval) 
}

## Loess predict method.
special_predict.lo.fitted <- function(x, data, se.fit = FALSE, ...) 
{
  p <- as.numeric(predict(x$model, newdata = data))
  p <- p - x$shift
  if(se.fit)
    p <- data.frame("fit" = p)
  return(p)
}

## Lasso.
lasso <- function(formula, ...)
{
  stopifnot(requireNamespace("glmnet"))

  ## Ensure it's a formula.
  if(inherits(formula, "matrix"))
    stop("only formulas are allowed!")
  if(!inherits(formula, "formula")) {
    formula <- as.character(substitute(formula))
    formula <- as.formula(paste("~", formula))
    environment(formula) <- sys.frame(-1)
  }

  ## List for setting up the special model term. 
  st <- list()
  st$control <- list(...)
  if(is.null(st$control$criterion))
    st$control$criterion <- "bic"
  st$term <- all.vars(formula) 
  st$label <- paste0("la(", paste0(gsub(" ", "",
    as.character(formula)), collapse = ""), ")") 
  st$X <- model.matrix(formula)
  if(length(j <- grep("(Intercept)", colnames(st$X), fixed = TRUE))) {
    st$X <- st$X[, -j, drop = FALSE]
  }
  st$formula <- formula

  ## Assign the "special" class and the new class "n".
  class(st) <- c("special", "glmnet")

  return(st) 
}

special_fit.glmnet <- function(x, z, w, control, ...)
{
  ## Set up glmnet call.
  call <- "glmnet::glmnet(x = x$X, y = z, weights = w"

  ## Add optional control parameters.
  if(!is.null(x$control)) {
    nc <- names(x$control)
    nc <- nc[nc != "criterion"]
    for(j in nc)
      call <- paste0(call, ", ", j, "= x$control$", j)
  }

  call <- paste0(call, ")")

  ## Estimate model.
  rval <- list("model" = eval(parse(text = call)))

  ## Get optimum lambda using IC.
  p <- predict(rval$model, newx = x$X)
  cm <- coef(rval$model)
  rss <- apply(p, 2, function(f) {
    sum(w * (z - f)^2)
  })
  edf <- apply(cm , 2, function(b) { sum(abs(b) >  1e-10) })
  n <- length(z)
  ic <- switch(tolower(x$control$criterion),
    "gcv" = rss * n / (n - edf)^2,
    "aic" = rss + 2 * edf,
    "gaic" = rss + K * edf,
    "aicc" = rss + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1),
    "bic" = rss + log(n) * edf
  )
  i <- which.min(ic)

  ## Save optimum lambda.
  rval$lambda <- rval$model$lambda[i]

  ## Get the fitted.values.
  rval$fitted.values <- p[, i]

  ## Save coefficients.
  rval$coefficients <- cm[, i]

  ## Center fitted values. 
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift

  ## Degrees of freedom.
  rval$edf <-  edf[i]

  ## Formula, needed for prediction.
  rval$formula <- x$formula

  ## IC and full edfs for plotting.
  rval$ic <- list("criterion" = x$control$criterion, "value" = ic, "edf" = edf)

  ## Assign class for predict method. 
  class(rval) <- "glmnet.fitted" 

  return(rval) 
}

## glmnet predict method.
special_predict.glmnet.fitted <- function(x, data, se.fit = FALSE, ...) 
{
  X <- model.matrix(x$formula, data = data)
  if(length(j <- grep("(Intercept)", colnames(X), fixed = TRUE))) {
    X <- X[, -j, drop = FALSE]
  }
  p <- as.numeric(predict(x$model, newx = X, s = x$lambda))
  p <- p - x$shift
  if(se.fit)
    p <- data.frame("fit" = p)
  return(p)
}

## Matrix block standardization.
blockstand <- function(x, n)
{
  cn <- colnames(x)
  decomp <- qr(x)
  if(decomp$rank < ncol(x))
    stop("block standardization cannot be computed, matrix is not of full rank!")
  scale <- qr.R(decomp) * 1 / sqrt(n)
  x <- qr.Q(decomp) * sqrt(n)
  attr(x, "blockscale") <- scale
  colnames(x) <- cn
  x
}

## Special lasso from Groll et al.
la <- function(x, type = 1, const = 1e-05, ...)
{
  call <- match.call()

  xn <- substitute(x)

  formula <- try(as.formula(xn), silent = TRUE)

  if(!inherits(formula, "try-error")) {
    xn <- all.vars(formula)
    environment(formula) <- environment(x)
  } else {
    xn <- as.character(xn)
    formula <- NULL
  }

  ## List for setting up the special model term. 
  st <- list()
  st$control <- list(...)
  if(is.null(st$control$criterion))
    st$control$criterion <- "bic"
  st$term <- xn 
  st$label <- gsub(" ", "", paste0("la(", as.character(deparse(call[[2]])), ")"))

  if(!is.null(formula)) {
    st$X <- model.matrix(formula,
      contrasts.arg = st$control$contrasts.arg,
      xlev = st$control$xlev)
  } else {
    if(is.factor(x)) {
      st$X <- model.matrix(~ x, contrasts.arg = st$control$contrasts.arg,
        xlev = st$control$xlev)
      colnames(st$X) <- gsub("x", "", colnames(st$X))
      st$is_factor <- TRUE
    } else {
      st$X <- x
    }
  }

  if(length(j <- grep("(Intercept)", colnames(st$X), fixed = TRUE))) {
    st$X <- st$X[, -j, drop = FALSE]
  }

  ## !FIXME
  if(isTRUE(st$is_factor) & FALSE) {
    st$blockscale <- attr(blockstand(st$X, n = nrow(st$X)), "blockscale")
    st$X <- st$X %*% st$blockscale
  }

  st$formula <- formula
  st$control$const <- const

  lt <- c("normal", "group", "nominal", "ordinal")
  if(!is.character(type))
    type <- lt[type[1L]]
  st$lasso_type <- lt[match(type, lt)]

  ## Assign the "special" class and the new class "n".
  class(st) <- c("special", "lasso")

  return(st) 
}

special_fit.lasso <- function(x, z, w, control, transfer, ...)
{
  if(is.null(control$criterion))
    control$criterion <- x$control$criterion

  ridge <- isTRUE(control$add_ridge)

  k <- ncol(x$X)
  XW <- x$X * w
  XWX <- crossprod(XW, x$X)
  XWz <- crossprod(XW, z)
  n <- length(z)
  logn <- log(n)
  bml <- drop(solve(XWX + diag(1e-08, k), XWz))
  b0 <- transfer$coefficients

  if(is.null(b0))
    b0 <- bml

  ## Penalty function.
  if(x$lasso_type == "normal") {
    pen <- function(b) {
      A <- 1 / sqrt(b^2 + x$control$const)
      A <- A * 1 / abs(bml)
      A <- if(length(A) < 2L) matrix(A, 1, 1) else diag(A)
      A
    }
  }

  if(x$lasso_type == "group") {
    pen <- function(b) {
      df <- ncol(x$X)
      A <- 1 / rep(sqrt(sum(b^2)), df) * 1 / rep(sqrt(sum(bml^2)), df)
      A <- if(length(A) < 2L) matrix(A, 1, 1) else diag(A)
      A
    }
  }

  if(x$lasso_type == "nominal") {
    pen <- function(b) {
      k <- ncol(x$X)
      Af <- matrix(0, ncol = choose(k, 2), nrow = k)
      combis <- combn(k, 2)
      for(ff in 1:ncol(combis)){
        Af[combis[1, ff], ff] <- 1
        Af[combis[2, ff], ff] <- -1
      }
      Af <- cbind(diag(k), Af)
      w <- rep(0, length = ncol(Af))
      df <- colSums(abs(x$X) > 0)
      nref <- n - sum(df)
      for(k in 1:ncol(Af)) {
        ok <- which(Af[, k] != 0)
        w[k] <- if(length(ok) < 2) {
          2 / (k + 1) * sqrt((df[ok[1]] + nref) / n)
        } else {
          2 / (k + 1) * sqrt((df[ok[1]] + df[ok[2]]) / n)
        }
        w[k] <- w[k] * 1 / abs(t(Af[, k]) %*% bml)
      }
      A <- 0
      for(k in 1:ncol(Af)) {
        tAf <- t(Af[, k])
        d <- drop(tAf %*% b)
        A <- A + w[k] / sqrt(d^2 + x$control$const) * Af[, k] %*% tAf
      }
      A
    }
  }

  if(x$lasso_type == "ordinal") {
    pen <- function(b) {
      k <- ncol(x$X)
      Af <- diff(diag(k + 1))
      Af[1, 1] <- 1
      Af <- Af[, -ncol(Af), drop = FALSE]
      w <- rep(0, length = ncol(Af))
      df <- colSums(abs(x$X) > 0)
      nref <- n - sum(df)
      for(k in 1:ncol(Af)) {
        ok <- which(Af[, k] != 0)
        w[k] <- if(length(ok) < 2) {
          sqrt((df[ok[1]] + nref) / n)
        } else {
          sqrt((df[ok[1]] + df[ok[2]]) / n)
        }
        w[k] <- w[k] * 1 / abs(t(Af[, k]) %*% bml)
      }
      A <- 0
      for(k in 1:ncol(Af)) {
        tAf <- t(Af[, k])
        d <- drop(tAf %*% b)
        A <- A + w[k] / sqrt(d^2 + x$control$const) * Af[, k] %*% tAf
      }
      A
    }
  }

  S <- pen(b0)

  if(ridge)
    S <- list(S, diag(ncol(S)))

  fl <- function(l, rf = FALSE, coef = FALSE) {
    if(ridge) {
      P <- try(chol2inv(chol(XWX + l[1]*S[[1]] + l[2]*S[[2]])), silent = TRUE)
    } else {
      P <- try(chol2inv(chol(XWX + l*S)), silent = TRUE)
    }

    if(inherits(P, "try-error"))
      P <- solve(XWX + S)

    b <- drop(P %*% XWz)

    fit <- drop(x$X %*% b)

    edf <- sum(diag(XWX %*% P))

    if(rf) {
      names(b) <- colnames(x$X)
      return(list("coefficients" = b, "fitted.values" = fit, "edf" = edf,
        "lambda" = l, "vcov" = P, "df" = n - edf))
    } else {
      if(isTRUE(control$logLik)) {
        eta[[j]] <- eta[[j]] + fit
        rss <- family$loglik(y, family$map2par(eta))
      } else {
        rss <- sum(w * (z - fit)^2)
      }

      rval <- switch(tolower(control$criterion),
        "gcv" = rss * n / (n - edf)^2,
        "aic" = rss + 2 * edf,
        "gaic" = rss + K * edf,
        "aicc" = rss + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1),
        "bic" = rss + logn * edf
      )

      if(coef) {
        rval <- list("ic" = rval, "coefficients" = b)
        names(rval$coefficients) <- colnames(x$X)
      }

      return(rval)
    }
  }

  ## Set up smoothing parameters.
  lambda <- if(is.null(transfer$lambda)) 10 else transfer$lambda
  if(ridge)
    lambda <- rep(lambda, length.out = 2L)

  opt <- nlminb(lambda, objective = fl, lower = lambda / 100, upper = lambda * 100)

  rval <- fl(opt$par, rf = TRUE)

  ## Tranfer arguments.
  rval$transfer <- list("lambda" = rval$lambda, "coefficients" = rval$coefficients)

  ## Arguments needed for prediction and path plots.
  keep <- c("formula", "term", "blockscale", "X")
  rval[keep] <- x[keep]
  rval$z <- z
  rval$w <- w
  rval$XWX <- XWX
  rval$XWz <- XWz
  rval$S <- S
  rval$criterion <- control$criterion
  rval$label <- x$label

  ## Assign class for predict method. 
  class(rval) <- "lasso.fitted"

  return(rval)
}

## Lasso predict method.
special_predict.lasso.fitted <- function(x, data, se.fit = FALSE, ...) 
{
  if(!is.null(x$formula)) {
    X <- model.matrix(x$formula, data = data,
      contrasts.arg = x$contrasts.arg, xlev = x$xlev)
  } else {
    if(is.factor(data[[x$term]])) {
      X <- model.matrix(~ data[[x$term]], contrasts.arg = x$control$contrasts.arg, xlev = x$control$xlev)
      colnames(X) <- gsub("data[[x$term]]", "", colnames(X), fixed = TRUE)
    } else {
      X <- data[[x$term]]
    }
  }

  if(length(j <- grep("(Intercept)", colnames(X), fixed = TRUE))) {
    X <- X[, -j, drop = FALSE]
  }

  ## !FIXME
  if(!is.null(x$blockscale)) {
    X <- X %*% x$blockscale
  }

  fit <- drop(X %*% x$coefficients)

  if(se.fit) {
    se <- rowSums((X %*% x$vcov) * X)
    se <- 2 * sqrt(se)
    fit <- cbind("fit" = fit, "lower" = fit - se, "upper" = fit + se)
  }

  return(fit)
}

## Lasso plotting function.
plot_lasso <- function(x, terms = NULL,
  which = c("criterion", "coefficients"),
  zoom = c(3, 4), spar = TRUE, ...)
{
  which <- match.arg(which)

  lwd <- list(...)$lwd
  col <- list(...)$col
  if(is.null(lwd))
    lwd <- 1.5
  if(is.null(col))
    col = 1

  if(inherits(x, "gamlss2")) {
    x <- specials(x, drop = FALSE)
    cx <- sapply(x, class)

    if(!any(j <- cx %in% c("lasso.fitted", "glmnet.fitted")))
      return(invisible(NULL))

    if(!is.null(terms)) {
      if(is.character(terms))
        terms <- grep2(terms, names(x), fixed = TRUE)
      x <- x[terms]
    }
 
    cx <- sapply(x, class)
    lmbd <- NULL
    if(any(jj <- cx == "lasso.fitted"))
      lmbd <- sapply(x[jj], function(x) x$lambda[1L])

    if(!is.null(lmbd)) {
      lambdas <- NULL
      if(is.null(zoom))
        zoom <- c(3, 4)
      zoom <- rep(zoom, length.out = 2L)
      zoom <- rev(zoom)
      grid <- list(...)$grid
      if(is.null(grid))
        grid <- 50
      for(l in lmbd) {
        lambdas <- c(lambdas, c(seq(log(l) - abs(log(l)) * zoom[1L], log(l), length = grid), log(l),
          seq(log(l), log(l) + abs(log(l)) * zoom[2L], length = grid)))
      }
      lambdas <- exp(sort(unique(lambdas)))
    }

    nx <- names(x)

    if(spar) {
      opar <- par(no.readonly = TRUE)
      par(mfrow = n2mfrow(length(nx)))
      if(which == "coefficients")
        par(mar = c(4, 4, 4, 8))
      on.exit(par(opar))
    }

    for(i in 1:length(nx))
      plot_lasso(x[[i]], which = which, lambdas = lambdas, label = nx[i], ...)
  } else {
    if(!is.null(x$model)) {
      lambdas <- log(x$model$lambda)
      i <- which.min(x$ic$value)
      if(which == "criterion") {
        plot(lambdas, x$ic$value, type = "l",
          xlab = expression(log(lambda)), ylab = toupper(x$ic$criterion),
          lwd = lwd, col = col, ...)
      } else {
        cm <- as.matrix(coef(x$model))
        matplot(lambdas, t(cm), type = "l", lty = 1,
          lwd = lwd, col = col,
          xlab = expression(log(lambda)), ylab = "Coefficients")
      }

      abline(v = lambdas[i], lty = 2, col = "lightgray")

      axis(3, at = lambdas[i], labels = bquote(lambda == .(round(exp(lambdas[i]), 4)) ~ ", edf =" ~ .(x$ic$edf[i])))
    } else {
      lambdas <- list(...)$lambdas
      n <- length(x$z)
      logn <- log(n)

      if(is.list(x$S))
        x$S <- x$S[[1L]]

      cm <- ic <- edfs <- NULL

      for(l in lambdas) {
        P <- try(chol2inv(chol(x$XWX + l*x$S)), silent = TRUE)

        if(inherits(P, "try-error"))
          P <- solve(x$XWX + x$S)

        b <- drop(P %*% x$XWz)

        fit <- drop(x$X %*% b)

        edf <- sum(diag(x$XWX %*% P))

        rss <- sum(x$w * (x$z - fit)^2)

        icl <- switch(tolower(x$criterion),
          "gcv" = rss * n / (n - edf)^2,
          "aic" = rss + 2 * edf,
          "gaic" = rss + K * edf,
          "aicc" = rss + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1),
          "bic" = rss + logn * edf
        )

        ic <- c(ic, icl)
        edfs <- c(edfs, edf)
        cm <- rbind(cm, b)
      }

      lab <- list(...)$label

      rind <- rev(1:length(ic))
      xlim <- rev(range(log(lambdas)))

      xlab <- list(...)$xlab
      if(is.null(xlab))
        xlab <- expression(log(lambda))

      ylab <- list(...)$ylab
      if(is.null(ylab))
        ylab <- if(which == "criterion") toupper(x$criterion) else "Coefficients"

      if(which == "criterion") {
        plot(log(lambdas), ic, type = "l", lwd = lwd, col = col,
          xlab = xlab, ylab = ylab,
          main = "", axes = FALSE, xlim = xlim)
      } else {
        matplot(log(lambdas), cm,
          type = "l", lty = 1, lwd = lwd, col = col,
          xlab = xlab, ylab = ylab,
          main = "", axes = FALSE, xlim = xlim)

        names <- list(...)$names
        if(is.null(names))
          names <- TRUE
        if(isFALSE(names))
          names <- NULL
        if(!is.null(names)) {
          if(!is.character(names))
            names <- colnames(x$X)
          names <- names[1:ncol(x$X)]
          at <- cm[1, ]

          labs <- labs0 <- names
          plab <- at
          o <- order(plab, decreasing = TRUE)
          labs <- labs[o]
          plab <- plab[o]
          rplab <- diff(range(plab))
          if(length(plab) > 2L) {
            for(i in 1:(length(plab) - 1)) {
              dp <- abs(plab[i] - plab[i + 1]) / rplab
              if(dp <= 0.02) {
                labs[i + 1] <- paste(c(labs[i], labs[i + 1]), collapse = ",")
                labs[i] <- ""
              }
            }
          }
          labs <- labs[order(o)]

          axis(4, at = at, labels = labs, las = 1, gap.axis = -1)
        }
      }

      box()
      axis(1)
      axis(2)

      i <- which.min(ic[rind])
      lo <- log(lambdas)[rind][i]

      abline(v = log(x$lambda[1L]), lty = 2, col = "lightgray")

      main <- list(...)$main
      if(is.null(main))
        main <- TRUE
      if(isTRUE(main) | is.character(main)) {
        if(is.character(main)) {
          mtext(main, side = 3, line = 1.5, font = 2, cex = 1.2)
        } else {
          mtext(bquote(log(lambda) == .(round(log(x$lambda), 3)) ~ " edf =" ~ .(round(x$edf, 2))), side = 3, line = 0.3, cex = 0.8)
          mtext(lab, side = 3, line = 1.5, font = 2, cex = 1.2)
        }
      }
    }
  }
  return(invisible(NULL))
}


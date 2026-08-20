## Special conditional inference tree constructor.
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

## tree fitting function for the backfitting algorithm.
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

## Now ctree.
ct <- function(formula, ...)
{
  stopifnot(requireNamespace("partykit"))
  st <- list()
  ctr <- list(...)
  st$control <- do.call(partykit::ctree_control, ctr)
  st$formula <- formula
  st$term <- all.vars(formula)
  st$label <- paste0("ct(", paste0(gsub(" ", "", as.character(formula)), collapse = ""), ")")
  st$data <- model.frame(formula)
  class(st) <- c("special", "ct")
  return(st)
}

special_fit.ct <- function(x, z, w, y, eta, j, family, control, ...)
{
  f <- update(x$formula, response_z ~ .)
  x$data$response_z <- z
  x$data$w <- w
  rval <- list(
    "model" = partykit::ctree(formula = f, data = x$data, weights = w,
      control = x$control)
  )
  rval$fitted.values <- predict(rval$model)
  rval$shift <- mean(rval$fitted.values)
  rval$fitted.values <- rval$fitted.values - rval$shift
  rval$edf <- length(partykit::nodeids(rval$model, terminal = TRUE))
  class(rval) <- "ct.fitted"
  return(rval)
}

special_predict.ct.fitted <- function(x, data, ...)
{
  p <- predict(x$model, newdata = data, type = "response")
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
lin <- function(x, ..., ridge = FALSE, scale = FALSE)
{
  x_expr <- substitute(x)
  x_txt  <- deparse1(x_expr, backtick = TRUE)

  f <- try(as.formula(x_txt), silent = TRUE)
  is_f <- !inherits(f, "try-error")

  if(!is_f) {
    vdot <- as.list(substitute(list(...)))[-1]
    vars <- c(x_txt, vapply(vdot, deparse1, character(1), backtick = TRUE))
    f <- as.formula(paste("~", paste(vars, collapse = "+")))
  }

  v <- all.vars(f)

  label_terms <- v
  if(ridge && ":" %in% all.names(f, functions = TRUE)) {
    tl <- attr(terms(f), "term.labels")
    if(length(tl))
      label_terms <- tl
  }

  make_label <- function(v, prefix = "lin", max_terms = 3) {
    n <- length(v)
    if(n <= max_terms) {
      inside <- paste(v, collapse = "+")
    } else {
      inside <- paste0(
        paste0(v[1:(max_terms - 1)], collapse = "+"),
        "+... [", n, "]"
      )
    }
    paste0(prefix, "(", inside, ")")
  }

  sx <- list(
    formula = f,
    term    = v,
    label   = make_label(label_terms, if(ridge) "ridge" else "lin"),
    by      = "NA",
    dim     = length(v),
    scale   = scale
  )

  if(!ridge) {
    sx$sp <- 1e-10
  }

  class(sx) <- "lin.smooth.spec"
  sx
}

ridge <- function(...)
{
  lin(..., ridge = TRUE, scale = TRUE)
}

smooth.construct.lin.smooth.spec <- function(object, data, knots)
{
  object$X <- model.matrix(object$formula, data = if(is.list(data)) 
    data[all.vars(reformulate(names(data))) %in% all.vars(object$formula)]
    else data)
  if(any(grepl("(Intercept)", colnames(object$X), fixed = TRUE))) {
    object$X <- object$X[, -1L, drop = FALSE]
  }
  if(object$scale) {
    sx <- list()
    for(j in seq_len(ncol(object$X))) {
      xj <- object$X[, j]
      sdj <- sd(xj)
      if(is.finite(sdj) && sdj > 0) {
        nm <- colnames(object$X)[j]
        muj <- mean(xj)
        sx[[as.character(j)]] <- list("mean" = muj, "sd" = sdj)
        object$X[, j] <- (xj - muj) / sdj
      }
    }
    object$scalex <- sx
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
  if(!is.null(object$scalex)) {
    for(j in as.integer(names(object$scalex))) {
      jc <- as.character(j)
      X[, j] <- (X[, j] - object$scalex[[jc]]$mean) / object$scalex[[jc]]$sd
    }
  }
  X[!is.finite(X)] <- 0
  return(X)
}

## Random effects.
re <- function(random, correlation = NULL, ...)
{
  stopifnot(requireNamespace("nlme"))

  call <- match.call(expand.dots = FALSE)

  ## List for setting up the special model term.
  st <- list()

  ## List of control arguments.
  ctr <- list(...)
  if(is.null(ctr$method))
    ctr$method <- "ML"

  fixed <- if(is.null(ctr$fixed)) ~1 else ctr$fixed
  ctr$fixed <- NULL

  ## Put all information together.
  st$fixed <- fixed
  ## Update fixed formula for working response z.
  st$fixed <- update(st$fixed, response_z ~ .)
  st$random <- random

  vr <- try(stats::formula(random), silent = TRUE)
  if(!inherits(vr, "try-error")) {
    vr <- all.vars(formula(as.Formula(random), lhs = 0, collapse = TRUE))
  } else {
    vr <- all.vars(call[[2L]])
  }

  st$term <- c(all.vars(stats::formula(fixed)), vr)

  if(!is.null(correlation)) {
    st$correlation <- correlation
    cf <- try(stats::formula(correlation), silent = TRUE)
    if(!inherits(cf, "try-error"))
      st$term <- c(st$term, all.vars(cf))
  }

  st$label <- gsub(" ", "", deparse(call))
  st$label <- gsub("re(", "re2(", st$label, fixed = TRUE)
  st$method <- ctr$method
  ctr$method <- NULL
  if(!is.null(ctr$control))
    st$control <- do.call(nlme::lmeControl, ctr$control)
  else
    st$control <- nlme::lmeControl()

  ## Assign data.
  env <- environment(random)
  if(is.null(env))
    env <- parent.frame()
  ff <- as.formula(paste("~", paste(st$term, collapse = "+")), env = env)
  st$data <- model.frame(ff)

  if(!all(i <- (st$term %in% names(st$data)))) {
    fi <- as.formula(paste("~", paste(st$term[!i], collapse = "+")), env = env)
    mf <- try(model.frame(fi), silent = TRUE)
    if(inherits(mf, "try-error")) {
      if(inherits(random, "formula")) {
        mf <- try(model.frame(fi), silent = TRUE)
      }
    } else {
      if(inherits(random, "formula")) {
        mf2 <- try(model.frame(fi), silent = TRUE)
        if(!inherits(mf2, "try-error"))
          mf <- cbind(mf, mf2)
      }
    }
    mf <- mf[, unique(names(mf)), drop = FALSE]
    st$data <- cbind(st$data, mf)
  }

  ## Assign the "special" class and the new class "n".
  class(st) <- c("special", "re")

  return(st)
}

## Set up the special "re" model term fitting function
edf_random_intercept <- function(fit, w) {
  g <- nlme::getGroups(fit)

  vc <- nlme::VarCorr(fit)
  tau2 <- as.numeric(vc[1, "Variance"])
  sigma2 <- fit$sigma^2

  wg <- tapply(w, g, sum)

  h <- tau2 * wg / sigma2
  sum(h / (1 + h))
}

special_fit.re   <- function(x, z, w, control, ...)
{
  ## Assign current working response.
  x$data$response_z <- z
  w <- pmax(w, .Machine$double.eps)
  x$data$weights_w <- w

  args <- list(
    "fixed" = x$fixed,
    "data" = x$data,
    "random" = x$random,
    "weights" = nlme::varFixed(~ weights_w),
    "method" = x$method,
    "control" = x$control,
    "keep.data" = FALSE
  )
  if(!is.null(x$correlation))
    args$correlation <- x$correlation

  ## Estimate model.
  rval <- list("model" = do.call(nlme::lme, args))

  ## Get the fitted.values.
  lvl <- length(rval$model$modelStruct$reStruct)
  p0 <- predict(rval$model, newdata = x$data, level = 0)
  p1 <- predict(rval$model, newdata = x$data, level = lvl)
  fit <- as.numeric(p1 - p0)
  rval$fitted.values <- fit ## fitted(rval$model)
  rval$coefficients <- nlme::ranef(rval$model)

  ## Degrees of freedom.
  ## For random-intercept-only models use shrinkage-aware edf.
  re_df <- nlme::ranef(rval$model)

  if(ncol(re_df) == 1L) {
    rval$edf <- edf_random_intercept(rval$model, 1/w)
  } else {
    ## Fallback: approximate numerical trace, or conservative upper bound.
    rval$edf <- nrow(re_df) * ncol(re_df)
  }

  ## Assign class for predict method.
  class(rval) <- "re.fitted"

  return(rval)
}

## The re() predict method.
special_predict.re.fitted <- function(x, data, se.fit = FALSE, alpha = 0.05, ...) {
  fit <- x$model
  stopifnot(inherits(fit, "lme"))

  ## Empty newdata guard.
  if(nrow(data) == 0) {
    if(!se.fit) return(numeric(0))
    return(data.frame(fit = numeric(0), lower = numeric(0), upper = numeric(0)))
  }

  ## Innermost random-effects level.
  reStruct <- fit$modelStruct$reStruct
  lvl <- length(reStruct)

  ## Random-effects contribution per row.
  p <- as.numeric(
    predict(fit, newdata = data, level = lvl, ...) -
    predict(fit, newdata = data, level = 0,   ...)
  )
  if(!se.fit) return(p)

  ## Grouping variable (must be present in newdata).
  grp_name <- names(reStruct)[lvl]
  if(is.null(grp_name) || !(grp_name %in% names(data))) {
    stop("newdata must contain the grouping variable for the innermost level: ", grp_name)
  }
  g_new   <- factor(data[[grp_name]])
  lev_new <- levels(g_new)

  ## Groups seen in the fitted model.
  g_fit <- levels(nlme::getGroups(fit))

  ## Build Z for the innermost level (keep intercept if present).
  re_form <- formula(reStruct, asList = TRUE)[[lvl]]
  Z <- model.matrix(re_form, data)  ## keeps (Intercept) if in re_form

  ## Determine RE column order from the fitted model (probe any existing group).
  probe <- intersect(lev_new, g_fit)
  if(!length(probe)) {
    ## all groups in newdata are unseen -> BLUPs undefined; keep NA bounds
    return(data.frame(fit = p, lower = NA_real_, upper = NA_real_))
  }
  Vi_probe <- as.matrix(nlme::getVarCov(fit, individual = probe[1], type = "random.effects"))
  re_cols  <- colnames(Vi_probe)

  ## Ensure Z has all required columns; add sensible defaults if missing.
  miss <- setdiff(re_cols, colnames(Z))
  if(length(miss)) {
    for(nm in miss) {
      if(nm == "(Intercept)") {
        Z <- cbind(Z, `(Intercept)` = 1)
      } else {
        Z <- cbind(Z, tmp = 0); colnames(Z)[ncol(Z)] <- nm
      }
    }
  }
  Z <- Z[, re_cols, drop = FALSE]

  ## Fast path: random intercept only.
  intercept_only <- ncol(Z) == 1 && identical(colnames(Z), "(Intercept)")

  ## Row-wise SEs from group-wise Var(b_g | y).
  se <- rep(NA_real_, nrow(data))
  if(intercept_only) {
    for(g in lev_new) {
      idx <- which(g_new == g)
      if(!length(idx) || !(g %in% g_fit)) next
      Vi <- as.matrix(nlme::getVarCov(fit, individual = g, type = "random.effects"))
      se[idx] <- sqrt(Vi[1, 1])
    }
  } else {
    for(g in lev_new) {
      idx <- which(g_new == g)
      if(!length(idx) || !(g %in% g_fit)) next
      Zi <- Z[idx, , drop = FALSE]
      Vi <- as.matrix(nlme::getVarCov(fit, individual = g, type = "random.effects"))
      if(is.null(dim(Vi))) Vi <- matrix(Vi, nrow = ncol(Z), ncol = ncol(Z))
      M <- Zi %*% Vi %*% t(Zi)
      se[idx] <- sqrt(pmax(0, diag(M)))
    }
  }

  z <- qnorm(1 - alpha/2)
  data.frame(fit = p, lower = p - z * se, upper = p + z * se)
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

## glmnet-based model term.
gnet <- function(formula, ...)
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
  st$label <- paste0("gnet(", paste0(gsub(" ", "",
    as.character(formula)), collapse = ""), ")")
  st$X <- model.matrix(formula)
  if(length(j <- grep("(Intercept)", colnames(st$X), fixed = TRUE))) {
    st$X <- st$X[, -j, drop = FALSE]
  }
  st$formula <- formula

  ## Assign the "special" class and the glmnet implementation class.
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

  ## Penalty.
  K <- control$K
  if(is.null(K))
    K <- 2

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

## Normal Lasso scaling (center + per-column RMS scaling).
normal_scale <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  cn <- colnames(X)
  if(is.null(cn)) stop("normal_scale: X must have colnames")

  mu <- setNames(colMeans(X), cn)
  Xc <- sweep(X, 2, mu[cn], "-")
  cn2 <- setNames(colSums(Xc^2), cn)
  s <- setNames(sqrt(n / pmax(cn2, 1e-12)), cn)

  function(X) {
    X <- as.matrix(X)
    if(is.null(colnames(X))) stop("normal_scale: X must have colnames")
    X <- X[, cn, drop = FALSE]
    Xc <- sweep(X, 2, mu[cn], "-")
    Xs <- sweep(Xc, 2, s[cn], "*")
    colnames(Xs) <- cn
    rownames(Xs) <- rownames(X)
    Xs
  }
}

## Fused scaling (scale only, no centering).
fused_scale <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  cn <- colnames(X)
  if(is.null(cn)) stop("fused_scale: X must have colnames")

  cn2 <- colSums(X^2)
  s <- sqrt(n / pmax(mean(cn2), 1e-12))
  s <- as.numeric(s)

  function(X) {
    X <- as.matrix(X)
    if(is.null(colnames(X))) stop("fused_scale: X must have colnames")
    X <- X[, cn, drop = FALSE]
    Xs <- X * s
    colnames(Xs) <- cn
    rownames(Xs) <- rownames(X)
    Xs
  }
}

## Group scaling (center + QR-based orthonormalization).
group_scale <- function(X) {
  X <- as.matrix(X)
  n <- nrow(X)
  cn <- colnames(X)
  if(is.null(cn)) stop("group_scale: X must have colnames")

  mu <- setNames(colMeans(X), cn)
  Xc <- sweep(X, 2, mu[cn], "-")

  decomp <- qr(Xc)
  if(decomp$rank < ncol(Xc)) stop("group_scale: X not full rank after centering")

  R <- qr.R(decomp)
  Tmat <- solve(R) * sqrt(n)

  function(X) {
    X <- as.matrix(X)
    if(is.null(colnames(X))) stop("group_scale: X must have colnames")
    X <- X[, cn, drop = FALSE]
    Xc <- sweep(X, 2, mu[cn], "-")
    Xs <- Xc %*% Tmat
    colnames(Xs) <- cn
    rownames(Xs) <- rownames(X)
    Xs
  }
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
  if(is.null(st$control$scale))
    st$control$scale <- TRUE
  st$term <- xn 
  st$label <- gsub(" ", "", paste0("la(", as.character(deparse(call[[2]])), ")"))

  if(!is.null(formula)) {
    st$X <- model.matrix(formula,
      contrasts.arg = st$control$contrasts.arg,
      xlev = st$control$xlev)
  } else {
    if(is.factor(x)) {
      st$lev <- levels(x)
      st$is_ordered <- inherits(x, "ordered")
      if(st$is_ordered) {
        x <- ordered(as.character(x), levels = st$lev)
      } else {
        x <- factor(as.character(x), levels = st$lev)
      }
      st$X <- model.matrix(~ x, contrasts.arg = st$control$contrasts.arg,
        xlev = st$control$xlev)
      colnames(st$X) <- gsub("x", "", colnames(st$X))
      st$is_factor <- TRUE
    } else {
      st$X <- x
    }
  }

  if(is.null(dim(st$X))) {
    st$X <- matrix(st$X, ncol = 1L)
  }
  if(is.null(colnames(st$X))) {
    colnames(st$X) <- if(length(st$term)) st$term[1L] else "x"
  }

  cn <- colnames(st$X)
  if(!is.null(cn) && length(j <- grep("(Intercept)", cn, fixed = TRUE))) {
    st$X <- st$X[, -j, drop = FALSE]
  }

  st$colnames <- colnames(st$X)
  st$formula <- formula
  st$control$const <- const

  lt <- c("normal", "group", "ordinal", "nominal")
  if(!is.character(type))
    type <- lt[type[1L]]
  st$lasso_type <- lt[match(type, lt)]

  if(st$control$scale && (st$lasso_type == "group")) {
    st$scale_fun <- group_scale(st$X)
  }

  if(st$control$scale && (st$lasso_type %in% c("ordinal", "nominal"))) {
    st$scale_fun <- fused_scale(st$X)
  }

  if(st$control$scale && is.null(st$scale_fun) && st$lasso_type == "normal") {
    st$scale_fun <- normal_scale(st$X)
  }

  if(!is.null(st$scale_fun)) {
    st$X <- st$scale_fun(st$X)
  }

  if(st$lasso_type %in% c("ordinal", "nominal")) {
    if(is.null(st$control$ridge)) {
      st$control$ridge <- TRUE
    }
  }

  ## Assign the "special" class and the new class "n".
  class(st) <- c("special", "lasso")

  return(st) 
}

special_fit.lasso <- function(x, z, w, control, transfer, ...)
{
  if(is.null(control$criterion))
    control$criterion <- x$control$criterion

  ridge <- isTRUE(control$add_ridge) || isTRUE(x$control$ridge)

  K <- control$K
  if(is.null(K))
    K <- 2

  k <- ncol(x$X)
  XW <- x$X * w
  XWX <- crossprod(XW, x$X)
  XWz <- crossprod(XW, z)
  n <- length(z)
  logn <- log(n)

  bml <- try(solve(XWX + diag(1e-08, k), XWz), silent = TRUE)
  if(inherits(bml, "try-error"))
    bml <- try(solve(XWX + diag(1e-07, k), XWz), silent = TRUE)
  if(inherits(bml, "try-error"))
    bml <- try(solve(XWX + diag(1e-06, k), XWz), silent = TRUE)
  if(inherits(bml, "try-error"))
    bml <- try(solve(XWX + diag(1e-05, k), XWz), silent = TRUE)
  if(inherits(bml, "try-error"))
    bml <- try(solve(XWX + diag(1e-04, k), XWz), silent = TRUE)
  if(inherits(bml, "try-error"))
    stop("cannot compute the ML estimator for the lasso term!")

  bml[abs(bml) < 1e-08] <- 1e-08

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
        if(nref < 0) {
          w[k] <- 1
        } else {
          ok <- which(Af[, k] != 0)
          w[k] <- if(length(ok) < 2) {
            2 / (k + 1) * sqrt((df[ok[1]] + nref) / n)
          } else {
            2 / (k + 1) * sqrt((df[ok[1]] + df[ok[2]]) / n)
          }
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
        if(nref < 0) {
          w[k] <- 1
        } else {
          ok <- which(Af[, k] != 0)
          w[k] <- if(length(ok) < 2) {
            sqrt((df[ok[1]] + nref) / n)
          } else {
            sqrt((df[ok[1]] + df[ok[2]]) / n)
          }
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
    A <- if(ridge) {
      XWX + l[1] * S[[1]] + l[2] * S[[2]]
    } else {
      XWX + l * S
    }

    P <- try(chol2inv(chol(A)), silent = TRUE)

    if(inherits(P, "try-error")) {
      dj <- max(1e-10, 1e-12 * mean(diag(A)))
      P <- try(solve(A + diag(dj, nrow(A))), silent = TRUE)
      if(inherits(P, "try-error"))
        return(Inf)
    }

    b <- drop(P %*% XWz)

    fit <- drop(x$X %*% b)

    edf <- sum(diag(XWX %*% P))

    if(rf) {
      names(b) <- colnames(x$X)
      return(list("coefficients" = b, "fitted.values" = fit, "edf" = edf,
        "lambda" = l, "vcov" = P, "df" = n - edf))
    } else {
      rss <- sum(w * (z - fit)^2)

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
  keep <- c("formula", "term", "blockscale", "X", "colnames",
    "scale_fun", "is_factor", "lev", "is_ordered", "control")
  rval[keep] <- x[keep]
  rval$z <- z
  rval$w <- w
  rval$XWX <- XWX
  rval$XWz <- XWz
  rval$S <- S
  rval$K <- K
  rval$criterion <- control$criterion
  rval$label <- x$label

  names(rval$coefficients) <- colnames(rval$X)

  ## Assign class for predict method. 
  class(rval) <- "lasso.fitted"

  return(rval)
}

## Lasso predict method.
special_predict.lasso.fitted <- function(x, data, se.fit = FALSE, ...) 
{
  ## Build design matrix.
  if(!is.null(x$formula)) {
    X <- model.matrix(x$formula, data = data, na.action = na.pass,
      contrasts.arg = x$control$contrasts.arg,
      xlev = x$control$xlev)
  } else {
    if(isTRUE(x$is_factor)) {
      vals <- as.character(data[[x$term]])
      bad <- which(!is.na(vals) & !(vals %in% x$lev))
      if(length(bad)) {
        print(head(unique(vals[bad]), 20))
        stop("Found values not in training levels")
      }
      lev <- x$lev
      if(is.null(lev))
        lev <- levels(data[[x$term]])
      if(isTRUE(x$is_ordered)) {
        data[[x$term]] <- ordered(as.character(data[[x$term]]), levels = lev)
      } else {
        data[[x$term]] <- factor(as.character(data[[x$term]]), levels = lev)
      }
      X <- model.matrix(~ data[[x$term]], na.action = na.pass)
      colnames(X) <- gsub("data[[x$term]]", "", colnames(X), fixed = TRUE)
    } else {
      X <- data[[x$term]]
    }
  }

  if(is.null(dim(X))) {
    X <- matrix(X, ncol = 1L)
  }
  if(is.null(colnames(X))) {
    colnames(X) <- if(length(x$term)) x$term[1L] else "x"
  }

  if(length(j <- grep("(Intercept)", colnames(X), fixed = TRUE))) {
    X <- X[, -j, drop = FALSE]
  }

  if(!is.null(x$colnames)) {
    X <- X[, x$colnames, drop = FALSE]
  }

  ## Scale.
  if(!is.null(x$scale_fun)) {
    X <- x$scale_fun(X)
  } else if(!is.null(x$blockscale)) {
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
    if(any(jj <- cx == "lasso.fitted")) {
      x <- x[jj]
      lmbd <- sapply(x, function(z) z$lambda[1L])
    } else {
      return(invisible(NULL))
    }

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

      cm <- ic <- edfs <- NULL

      for(l in lambdas) {
        if(is.list(x$S)) {
          P <- try(chol2inv(chol(x$XWX + l*x$S[[1]] + x$lambda[2]*x$S[[2]])),
            silent = TRUE)
        } else {
          P <- try(chol2inv(chol(x$XWX + l*x$S)), silent = TRUE)
        }

        if(inherits(P, "try-error"))
          P <- solve(x$XWX + x$S)

        b <- drop(P %*% x$XWz)

        fit <- drop(x$X %*% b)

        edf <- sum(diag(x$XWX %*% P))

        rss <- sum(x$w * (x$z - fit)^2)

        icl <- switch(tolower(x$criterion),
          "gcv" = rss * n / (n - edf)^2,
          "aic" = rss + 2 * edf,
          "gaic" = rss + x$K * edf,
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
          if(!is.character(names)) {
            names <- colnames(x$X)
            if(is.null(names))
              names <- paste0("x", 1:ncol(x$X))
          }

          if(length(names) < ncol(x$X))
            names <- rep(names, length.out = ncol(x$X))
          names <- names[1:ncol(x$X)]

          if(!all(names == "")) {
            at <- cm[1, ]

            labs <- labs0 <- names
            plab <- at
            o <- order(plab, decreasing = TRUE)
            labs <- labs[o]
            plab <- plab[o]
            rplab <- diff(range(plab))
            if(length(plab) > 2L) {
              for(i in 1:(length(plab) - 1)) {
                dp <- abs(plab[i] - plab[i + 1]) / (rplab + 1e-08)
                if((dp <= 0.02) || (rplab <= 0.02)) {
                  labs[i + 1] <- paste(c(labs[i], labs[i + 1]), collapse = ",")
                  labs[i] <- ""
                }
              }
            }
            labs <- labs[order(o)]

            if(!all(labs == "")) {
              axis(4, at = at, labels = labs, las = 1, gap.axis = -1)
            }
          }
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
      if(isTRUE(list(...)$info)) {
        cex.info <- list(...)$cex.info
        if(is.null(cex.info))
          cex.info <- 1
        mtext(bquote(log(lambda) == .(round(log(x$lambda), 3)) ~ " edf =" ~ .(round(x$edf, 2))), side = 3, line = 0.3, cex = cex.info)
      }
    }
  }
  return(invisible(NULL))
}

## Information criterion used by the free-knot term.
fk_ic <- function(rss, n, edf, criterion = "bic", K = 2)
{
  if(!is.finite(rss) || !is.finite(edf) || n < 1L)
    return(Inf)

  switch(criterion,
    "gcv" = if(n > edf) rss * n / (n - edf)^2 else Inf,
    "aic" = rss + 2 * edf,
    "gaic" = rss + K * edf,
    "aicc" = if(n > edf + 1) {
      rss + 2 * edf + (2 * edf * (edf + 1)) / (n - edf - 1)
    } else Inf,
    "bic" = rss + log(n) * edf,
    stop("unknown free-knot selection criterion: ", criterion, call. = FALSE)
  )
}

## Map log gaps to ordered interior knots. Fixing the final log gap at zero
## makes the parameterization identifiable and prevents knot crossings.
fk_theta2knots <- function(theta, boundary)
{
  if(!length(theta))
    return(numeric())

  loggaps <- c(theta, 0)
  loggaps <- loggaps - max(loggaps)
  gaps <- exp(loggaps)
  gaps <- gaps / sum(gaps)
  boundary[1L] + diff(boundary) * cumsum(gaps)[seq_along(theta)]
}

## Inverse transformation, used for supplied and warm-start knot locations.
fk_knots2theta <- function(knots, boundary)
{
  if(!length(knots))
    return(numeric())

  gaps <- diff(c(boundary[1L], knots, boundary[2L])) / diff(boundary)
  if(any(!is.finite(gaps)) || any(gaps <= 0))
    stop("interior knots must be distinct and strictly inside the data range",
      call. = FALSE)
  log(gaps[-length(gaps)]) - log(gaps[length(gaps)])
}

## Construct a centered B-spline basis. Because the rows of the uncentered
## basis sum to one, the final centered column is redundant and is omitted
## during fitting.
fk_basis <- function(x, interior, boundary, degree, center = NULL)
{
  order <- degree + 1L
  knots <- c(rep(boundary[1L], order), interior, rep(boundary[2L], order))
  finite <- is.finite(x)
  if(all(finite)) {
    X <- splines::splineDesign(knots, x, ord = order, outer.ok = FALSE)
  } else {
    if(is.null(center))
      stop("non-finite covariate values in free-knot basis construction",
        call. = FALSE)
    X <- matrix(NA_real_, length(x), length(knots) - order)
    if(any(finite))
      X[finite, ] <- splines::splineDesign(knots, x[finite],
        ord = order, outer.ok = FALSE)
  }
  if(is.null(center))
    center <- colMeans(X)
  X <- sweep(X, 2L, center, FUN = "-")
  keep <- if(ncol(X) > 1L) seq_len(ncol(X) - 1L) else integer()
  list(X = X[, keep, drop = FALSE], X.full = X, center = center,
    knots = knots, keep = keep)
}

## Weighted least-squares fit conditional on a set of interior knots.
fk_wfit <- function(x, y, w, interior, boundary, degree, final = FALSE,
  knot.edf = length(interior), criterion = "bic", K = 2)
{
  bs <- fk_basis(x, interior, boundary, degree)
  X <- bs$X
  n <- sum(w > 0)

  if(!ncol(X)) {
    coefficients <- numeric(ncol(bs$X.full))
    vcov <- matrix(0, ncol(bs$X.full), ncol(bs$X.full))
    fit <- rep(0, length(y))
    basis.edf <- 0
  } else {
    XWX <- crossprod(X * w, X)
    XWy <- crossprod(X, w * y)
    ridge <- max(1e-12, 1e-10 * mean(diag(XWX)))
    A <- XWX + diag(ridge, ncol(X))
    R <- try(chol(A), silent = TRUE)
    if(inherits(R, "try-error"))
      return(NULL)

    coefficients.reduced <- drop(backsolve(R,
      forwardsolve(t(R), XWy)))
    fit <- drop(X %*% coefficients.reduced)

    ## The inverse and conditional EDF are needed only for the final fit.
    if(final) {
      P <- chol2inv(R)
      basis.edf <- sum(diag(XWX %*% P))
      coefficients <- c(coefficients.reduced, 0)
      vcov <- matrix(0, ncol(bs$X.full), ncol(bs$X.full))
      vcov[bs$keep, bs$keep] <- P
    } else {
      basis.edf <- ncol(X)
      coefficients <- vcov <- NULL
    }
  }

  rss <- sum(w * (y - fit)^2)
  edf <- basis.edf + knot.edf
  rval <- list(rss = rss, edf = edf,
    ic = fk_ic(rss, n, edf, criterion = criterion, K = K))

  if(final) {
    rval <- c(rval, list(
      coefficients = coefficients,
      fitted.values = fit,
      basis.edf = basis.edf,
      knot.edf = knot.edf,
      lambdas = 0,
      vcov = vcov,
      df = n - edf,
      knots = bs$knots,
      interior.knots = interior,
      boundary = boundary,
      basis.center = bs$center,
      degree = degree
    ))
  }
  rval
}

## Fit a fixed number of free interior knots. Supplying 'knots' fixes their
## locations; 'start' supplies warm-start locations that remain free.
free_knots <- function(x, y, w = NULL, k = 0L, degree = 1L,
  criterion = "bic", knots = NULL, K = 2, start = NULL,
  optim.control = list(), ...)
{
  x <- as.numeric(x)
  y <- as.numeric(y)
  if(is.null(w))
    w <- rep(1, length(y))
  w <- as.numeric(w)

  if(length(x) != length(y) || length(w) != length(y))
    stop("'x', 'y', and 'w' must have equal lengths", call. = FALSE)
  if(any(!is.finite(x)) || any(!is.finite(y)) || any(!is.finite(w)) ||
      any(w < 0) || !any(w > 0))
    stop("free-knot inputs must be finite with nonnegative, nonzero weights",
      call. = FALSE)

  degree <- as.integer(degree)[1L]
  if(is.na(degree) || degree < 0L)
    stop("'degree' must be a nonnegative integer", call. = FALSE)
  criterion <- match.arg(tolower(criterion),
    c("gcv", "aic", "gaic", "aicc", "bic"))
  K <- as.numeric(K)[1L]
  if(!is.finite(K) || K < 0)
    stop("'K' must be a nonnegative finite number", call. = FALSE)

  xr <- range(x)
  if(diff(xr) <= 0)
    stop("fk() needs at least two distinct covariate values", call. = FALSE)
  boundary <- xr + c(-1, 1) * diff(xr) * 1e-07

  fixed <- !is.null(knots)
  if(fixed) {
    interior <- sort(as.numeric(knots))
    if(any(!is.finite(interior)) || anyDuplicated(interior) ||
        any(interior <= xr[1L]) || any(interior >= xr[2L]))
      stop("'knots' must be distinct, finite, and strictly inside the data range",
        call. = FALSE)
    k <- length(interior)
  } else {
    k <- as.integer(k)[1L]
    if(is.na(k) || k < 0L)
      stop("'k' must be a nonnegative integer", call. = FALSE)
    interior <- NULL
  }

  max.k <- max(0L, length(unique(x[w > 0])) - degree - 1L)
  if(k > max.k)
    stop("too many interior knots for the distinct covariate values",
      call. = FALSE)

  opt <- list(par = numeric(), value = NA_real_, counts = integer(),
    convergence = 0L, message = NULL)
  if(!fixed && k > 0L) {
    valid.start <- !is.null(start) && length(start) == k &&
      all(is.finite(start)) && all(start > xr[1L]) && all(start < xr[2L]) &&
      !anyDuplicated(start)
    if(!valid.start) {
      start <- as.numeric(quantile(x, probs = seq_len(k) / (k + 1),
        names = FALSE, type = 8))
      if(anyDuplicated(start) || any(start <= xr[1L]) || any(start >= xr[2L]))
        start <- seq(xr[1L], xr[2L],
          length.out = k + 2L)[-c(1L, k + 2L)]
    } else {
      start <- sort(start)
    }

    theta <- fk_knots2theta(start, xr)
    objective <- function(theta) {
      interior <- fk_theta2knots(theta, xr)
      fit <- fk_wfit(x, y, w, interior, boundary, degree,
        final = FALSE, knot.edf = k, criterion = criterion, K = K)
      if(is.null(fit)) Inf else fit$rss
    }
    ocontrol <- utils::modifyList(list(maxit = 100L, factr = 1e7),
      optim.control)
    opt <- optim(theta, fn = objective, method = "L-BFGS-B",
      lower = rep(-20, k), upper = rep(20, k), control = ocontrol)

    ## A B-spline objective is only piecewise smooth when a free knot passes an
    ## observed x value. Retry line-search failures with a derivative-free
    ## method while retaining L-BFGS-B as the fast path.
    if(opt$convergence != 0L) {
      retry.control <- list(
        maxit = max(200L, if(is.null(ocontrol$maxit)) 0L else ocontrol$maxit),
        reltol = 1e-08
      )
      retry <- if(k == 1L) {
        optim(opt$par, fn = objective, method = "Brent",
          lower = -20, upper = 20, control = retry.control)
      } else {
        optim(opt$par, fn = objective, method = "Nelder-Mead",
          control = retry.control)
      }
      if(retry$convergence == 0L || retry$value < opt$value)
        opt <- retry
    }
    interior <- fk_theta2knots(opt$par, xr)
  }

  knot.edf <- if(fixed) 0L else k
  rval <- fk_wfit(x, y, w, interior, boundary, degree, final = TRUE,
    knot.edf = knot.edf, criterion = criterion, K = K)
  if(is.null(rval))
    stop("free-knot weighted least-squares fit is singular", call. = FALSE)

  rval$nk <- k
  rval$criterion <- criterion
  rval$K <- K
  rval$fixed.knots <- fixed
  rval$convergence <- opt$convergence
  rval$message <- opt$message
  rval$counts <- opt$counts
  if(opt$convergence != 0L)
    warning("free-knot optimization did not converge",
      if(!is.null(opt$message)) paste0(": ", opt$message), call. = FALSE)
  rval
}

## Search for the optimum number of free interior knots.
n_free_knots <- function(..., nk = 10L, start = NULL)
{
  nk <- as.integer(nk)[1L]
  if(is.na(nk) || nk < 0L)
    stop("'nk' must be a nonnegative integer", call. = FALSE)

  best <- NULL
  search <- vector("list", nk + 1L)
  errors <- character()
  for(k in 0:nk) {
    start.k <- if(!is.null(start) && length(start) == k) start else NULL
    fit <- suppressWarnings(try(
      free_knots(..., k = k, start = start.k), silent = TRUE))
    if(inherits(fit, "try-error")) {
      errors <- c(errors, paste0("k = ", k, ": ",
        conditionMessage(attr(fit, "condition"))))
      next
    }
    if(fit$convergence != 0L) {
      errors <- c(errors, paste0("k = ", k, ": optimization did not converge"))
      next
    }

    search[[k + 1L]] <- data.frame(nk = k, ic = fit$ic,
      convergence = fit$convergence)
    if(is.null(best) || fit$ic < best$ic)
      best <- fit
  }

  if(is.null(best))
    stop("all free-knot candidate fits failed: ", paste(errors, collapse = "; "),
      call. = FALSE)
  best$search <- do.call(rbind, search)
  best$search.errors <- errors
  best
}

## Free-knots model term constructor.
fk <- function(x, k = 10L, nk = NULL, degree = 1L,
  criterion = "bic", knots = NULL, K = NULL, reselect = FALSE,
  optim.control = list(), ...)
{
  xexpr <- substitute(x)
  expr <- deparse1(xexpr, backtick = TRUE)
  term <- all.vars(xexpr)
  if(length(term) != 1L)
    stop("fk() requires an expression involving one numeric covariate",
      call. = FALSE)
  f <- as.formula(paste("~", expr), env = parent.frame())

  ## special_terms() evaluates constructors in the model-frame data mask.
  ## For transformed terms that mask contains a column named by the complete
  ## expression rather than the original variable.
  penv <- parent.frame()
  xv <- if(exists(expr, envir = penv, inherits = FALSE)) {
    get(expr, envir = penv, inherits = FALSE)
  } else x
  if(!is.numeric(xv) || is.matrix(xv))
    stop("fk() requires a numeric covariate", call. = FALSE)

  degree <- as.integer(degree)[1L]
  k <- as.integer(k)[1L]
  if(!is.null(nk))
    nk <- as.integer(nk)[1L]
  if(is.na(k) || k < 0L || (!is.null(nk) && (is.na(nk) || nk < 0L)))
    stop("'k' and 'nk' must be nonnegative integers", call. = FALSE)
  if(is.na(degree) || degree < 0L)
    stop("'degree' must be a nonnegative integer", call. = FALSE)
  criterion <- match.arg(tolower(criterion),
    c("gcv", "aic", "gaic", "aicc", "bic"))
  if(!is.null(K)) {
    K <- as.numeric(K)[1L]
    if(!is.finite(K) || K < 0)
      stop("'K' must be a nonnegative finite number", call. = FALSE)
  }
  if(!is.list(optim.control))
    stop("'optim.control' must be a list", call. = FALSE)
  optim.control <- utils::modifyList(optim.control, list(...))

  if(!is.null(knots)) {
    knots <- as.numeric(knots)
    if(!is.null(nk) && nk != length(knots))
      stop("'nk' must equal length(knots) when both are supplied",
        call. = FALSE)
    nk <- length(knots)
  }

  sx <- list(
    term = term,
    label = paste0("fk(", expr, ")"),
    expression = expr,
    formula = f,
    x = as.numeric(xv),
    control = list(k = k, nk = nk, degree = degree,
      criterion = criterion, knots = knots, K = K,
      reselect = isTRUE(reselect), optim.control = optim.control)
  )
  class(sx) <- c("special", "fk")
  sx
}

## Free-knots fit function.
special_fit.fk <- function(x, z, w, y, eta, j, family, control,
  transfer = NULL, ...)
{
  ctr <- x$control
  K <- if(is.null(ctr$K)) {
    if(is.null(control$K)) 2 else control$K
  } else ctr$K
  criterion <- if(is.null(ctr$criterion)) {
    if(is.null(control$criterion)) "bic" else control$criterion
  } else ctr$criterion

  start <- NULL
  if(!is.null(transfer$interior.knots))
    start <- transfer$interior.knots

  if(!is.null(ctr$knots)) {
    sx <- free_knots(x$x, z, w, degree = ctr$degree,
      criterion = criterion, knots = ctr$knots, K = K,
      optim.control = ctr$optim.control)
  } else if(!is.null(ctr$nk)) {
    sx <- free_knots(x$x, z, w, k = ctr$nk, degree = ctr$degree,
      criterion = criterion, K = K, start = start,
      optim.control = ctr$optim.control)
  } else if(!isTRUE(ctr$reselect) && !is.null(transfer$nk)) {
    sx <- free_knots(x$x, z, w, k = transfer$nk, degree = ctr$degree,
      criterion = criterion, K = K, start = start,
      optim.control = ctr$optim.control)
  } else {
    sx <- n_free_knots(x$x, z, w, nk = ctr$k, degree = ctr$degree,
      criterion = criterion, K = K, start = start,
      optim.control = ctr$optim.control)
  }

  if(is.null(sx$search) && !is.null(transfer$search))
    sx$search <- transfer$search
  if(is.null(sx$search.errors) && !is.null(transfer$search.errors))
    sx$search.errors <- transfer$search.errors
  sx$shift <- 0
  sx$term <- x$term
  sx$formula <- x$formula
  sx$transfer <- list(
    nk = sx$nk,
    interior.knots = sx$interior.knots,
    search = sx$search,
    search.errors = sx$search.errors
  )
  class(sx) <- "fk.fitted"
  sx
}

## Free-knots predict function. Intervals are conditional on the selected knot
## locations; uncertainty from knot search is not included.
special_predict.fk.fitted <- function(x, data, se.fit = FALSE,
  alpha = 0.05, ...)
{
  expr <- deparse1(x$formula[[2L]], backtick = TRUE)
  xv <- if(expr %in% names(data)) {
    data[[expr]]
  } else {
    f <- x$formula
    model.frame(f, data = data, na.action = na.pass)[[1L]]
  }
  bs <- fk_basis(as.numeric(xv), x$interior.knots, x$boundary,
    x$degree, center = x$basis.center)
  X <- bs$X.full
  p <- drop(X %*% coef(x))

  if(se.fit) {
    se <- sqrt(pmax(0, rowSums((X %*% x$vcov) * X)))
    crit <- qnorm(1 - alpha / 2)
    p <- data.frame(
      fit = p,
      lower = p - crit * se,
      upper = p + crit * se
    )
  }
  p
}

## Monotonic P-spline.
ms <- function(x, mono = c("up", "down"), k = 10,
  criterion = "aicc", lambda = NULL,
  lambda.min = 1e-8, lambda.max = 1e8, lambda.tol = 1e-4,
  ...)
{
  mono <- match.arg(mono)

  xexpr <- substitute(x)
  term <- deparse1(xexpr, backtick = TRUE)

  f <- as.formula(paste("~", term), env = parent.frame())
  xv <- model.frame(f, na.action = na.pass)[[1L]]

  if(!is.numeric(xv))
    stop("ms() currently requires a numeric covariate")

  if(length(unique(xv[is.finite(xv)])) < 4L)
    stop("ms() needs at least four distinct finite covariate values")

  sc <- s(x, bs = "ps", k = k, ...)

  sm <- smoothCon(
    sc,
    data = data.frame(x = xv),
    absorb.cons = FALSE
  )[[1L]]

  sx <- list(
    term = term,
    label = paste0("ms(", term, ")"),
    formula = f,
    X = sm$X,
    S = sm$S[[1L]],
    smooth = sm,
    mono = mono,
    control = list(
      criterion = criterion,
      lambda = lambda,
      lambda.min = lambda.min,
      lambda.max = lambda.max,
      lambda.tol = lambda.tol
    )
  )

  class(sx) <- c("special", "ms")
  sx
}

special_fit.ms <- function(x, z, w, y, eta, j, family, control,
  transfer = NULL, iter = NULL, ...)
{
  n <- length(z)
  X <- x$X
  p <- ncol(X)

  if(p < 2L)
    stop("ms() needs at least two spline coefficients")

  if(length(w) != n || any(!is.finite(w)) || any(w < 0))
    stop("invalid working weights in ms()")

  if(is.null(transfer))
    transfer <- list()

  if(!is.null(x$control)) {
    control[names(x$control)] <- x$control
    if(!is.null(control$method))
      control$criterion <- tolower(control$method)
  }

  if(is.null(control$criterion))
    control$criterion <- "aicc"

  criterion <- tolower(control$criterion)

  if(criterion == "ml")
    criterion <- "aicc"

  K <- if(is.null(control$K)) 2 else control$K

  if(!is.null(x$binning)) {
    rw <- numeric(length(x$binning$nodups))
    rz <- numeric(length(x$binning$nodups))

    calc_Xe(
      x$binning$sorted.index,
      w, z, rw, rz,
      x$binning$order
    )

    XWX <- crossprod(X * rw, X)
    XWz <- drop(crossprod(X, rz))
  } else {
    XW <- X * w
    XWX <- crossprod(XW, X)
    XWz <- drop(crossprod(XW, z))
  }

  zWz <- sum(w * z^2)

  ## P-spline penalty.
  S <- if(is.list(x$S)) x$S[[1L]] else x$S
  S <- 0.5 * (S + t(S))

  ## Monotonicity matrix G:
  ##
  ##     G %*% b >= 0.
  ##
  ## For an increasing B-spline, adjacent coefficients are constrained
  ## to be non-decreasing. For decreasing, reverse the sign.
  D <- matrix(0, p - 1L, p)
  ii <- seq_len(p - 1L)

  D[cbind(ii, ii)] <- -1
  D[cbind(ii, ii + 1L)] <- 1

  mono <- if(is.null(x$mono)) "up" else tolower(x$mono[1L])

  if(mono %in% c("up", "increasing", "increase", "+", "1")) {
    G <- D
  } else if(mono %in% c("down", "decreasing", "decrease", "-", "-1")) {
    G <- -D
  } else {
    stop("unknown 'mono' specification in ms(): ", mono)
  }

  warm.active <- transfer$active
  if(is.null(warm.active))
    warm.active <- integer(0L)

  warm.active <- intersect(
    as.integer(warm.active),
    seq_len(nrow(G))
  )

  solve_lambda <- function(lambda, active0 = integer(0L))
  {
    A <- XWX + lambda * S

    R <- try(chol(A), silent = TRUE)

    if(inherits(R, "try-error")) {
      ridge <- max(
        1e-10,
        max(abs(diag(A)), na.rm = TRUE) * 1e-10
      )
      R <- chol(A + diag(ridge, p))
    }

    H <- chol2inv(R)
    b0 <- drop(H %*% XWz)

    active <- active0

    tol <- sqrt(.Machine$double.eps) *
      (1 + max(abs(b0)))

    eqsolve <- function(active)
    {
      if(!length(active)) {
        return(list(
          b = b0,
          mu = numeric(0L),
          P = H
        ))
      }

      C <- G[active, , drop = FALSE]
      HC <- H %*% t(C)
      M <- C %*% HC

      RM <- chol(M)

      mu <- -drop(
        backsolve(
          RM,
          forwardsolve(
            t(RM),
            drop(C %*% b0)
          )
        )
      )

      b <- drop(b0 + HC %*% mu)

      Q <- backsolve(
        RM,
        forwardsolve(
          t(RM),
          t(HC)
        )
      )

      P <- H - HC %*% Q
      P <- 0.5 * (P + t(P))

      list(
        b = b,
        mu = mu,
        P = P
      )
    }


    ans <- eqsolve(active)

    maxit <- max(20L, 4L * p)

    for(it in seq_len(maxit)) {
      gb <- drop(G %*% ans$b)

      bad <- which(gb < -tol)

      if(length(bad)) {
        add <- bad[which.min(gb[bad])]

        if(!(add %in% active))
          active <- c(active, add)

        ## Restore dual feasibility after adding a constraint.
        repeat {
          ans <- eqsolve(active)

          neg <- which(ans$mu < -tol)
          if(!length(neg))
            break

          rem <- neg[which.min(ans$mu[neg])]
          active <- active[-rem]
        }

        next
      }

      ## Primal feasible; verify dual feasibility.
      if(length(active)) {
        ans <- eqsolve(active)

        neg <- which(ans$mu < -tol)

        if(length(neg)) {
          rem <- neg[which.min(ans$mu[neg])]
          active <- active[-rem]
          ans <- eqsolve(active)
          next
        }
      }

      ## tr(X'WX P), using sum(A * t(B)) with symmetric matrices.
      edf <- sum(XWX * ans$P)

      return(list(
        coefficients = ans$b,
        edf = edf,
        vcov = ans$P,
        active = sort(active)
      ))
    }

    stop("active-set algorithm did not converge in special_fit.ms()")
  }

  ## Fast smoothing-parameter criterion.
  ## RSS = z'Wz - 2 b'X'Wz + b'X'WX b.
  objective <- function(rho)
  {
    l <- exp(rho)

    qf <- solve_lambda(
      l,
      active0 = warm.active
    )

    b <- qf$coefficients
    edf <- qf$edf

    rss <- zWz -
      2 * sum(b * XWz) +
      sum(b * (XWX %*% b))

    rss <- max(rss, 0)

    switch(
      criterion,

      "gcv" =
        rss * n / (n - edf)^2,

      "aic" =
        rss + 2 * edf,

      "gaic" =
        rss + K * edf,

      "aicc" = {
        den <- n - edf - 1
        if(den <= 0)
          Inf
        else
          rss + 2 * edf +
            (2 * edf * (edf + 1)) / den
      },

      "bic" =
        rss + log(n) * edf,

      stop(
        "unknown smoothing criterion for ms(): ",
        criterion
      )
    )
  }

  ## Select lambda.
  fixed.lambda <- control$lambda

  if(!is.null(fixed.lambda)) {
    l <- as.numeric(fixed.lambda[1L])

    if(!is.finite(l) || l <= 0)
      stop("'lambda' in ms() must be a positive finite number")
  } else {
    lambda.min <- if(is.null(control$lambda.min))
      1e-8 else as.numeric(control$lambda.min[1L])

    lambda.max <- if(is.null(control$lambda.max))
      1e8 else as.numeric(control$lambda.max[1L])

    lambda.tol <- if(is.null(control$lambda.tol))
      1e-4 else as.numeric(control$lambda.tol[1L])

    if(!is.finite(lambda.min) || !is.finite(lambda.max) ||
       lambda.min <= 0 || lambda.max <= lambda.min)
      stop("invalid lambda search interval in ms()")

    lower <- log(lambda.min)
    upper <- log(lambda.max)

    ## Warm-start the lambda search from the previous backfitting fit.
    lold <- transfer$lambdas

    if(length(lold) && is.finite(lold[1L]) && lold[1L] > 0) {
      r0 <- log(lold[1L])

      lower <- max(lower, r0 - log(10))
      upper <- min(upper, r0 + log(10))

      ## Defensive guard for a degenerate clipped interval.
      if(lower >= upper) {
        lower <- log(lambda.min)
        upper <- log(lambda.max)
      }
    }

    opt <- optimize(
      objective,
      interval = c(lower, upper),
      tol = lambda.tol
    )

    l <- exp(opt$minimum)
  }

  ## Final fit for selected lambda.
  qf <- solve_lambda(
    l,
    active0 = warm.active
  )

  b <- qf$coefficients
  P <- qf$vcov
  edf <- qf$edf

  fit <- drop(X %*% b)

  if(!is.null(x$binning))
    fit <- fit[x$binning$match.index]


  ## Center term for identifiability in additive backfitting.
  shift <- mean(fit)
  fit <- fit - shift

  ## Return the same essential contract as smooth.construct_wfit().
  rval <- list(
    "coefficients" = b,
    "fitted.values" = fit,
    "edf" = edf,
    "lambdas" = l,
    "vcov" = P,
    "df" = n - edf,

    ## ms-specific bookkeeping.
    "shift" = shift,
    "term" = x$term,
    "formula" = x$formula,
    "smooth" = x$smooth,
    "mono" = x$mono,
    "active" = qf$active,

    ## Warm starts for the next gamlss2 backfitting iteration.
    "transfer" = list(
      "lambdas" = l,
      "coefficients" = b,
      "active" = qf$active
    )
  )

  class(rval) <- "ms.fitted"
  rval
}

## Prediction method.
special_predict.ms.fitted <- function(x, data, se.fit = FALSE,
  alpha = 0.05, ...)
{
  if(nrow(data) == 0L) {
    if(!se.fit)
      return(numeric(0L))

    return(data.frame(
      fit = numeric(0L),
      lower = numeric(0L),
      upper = numeric(0L)
    ))
  }

  mf <- model.frame(
    x$formula,
    data = data,
    na.action = na.pass
  )

  xv <- mf[[1L]]

  Xp <- PredictMat(
    x$smooth,
    data.frame(x = xv)
  )

  fit <- drop(Xp %*% x$coefficients) - x$shift

  if(!se.fit)
    return(fit)

  XP <- Xp %*% x$vcov
  se <- sqrt(pmax(0, rowSums(XP * Xp)))

  za <- qnorm(1 - alpha / 2)

  data.frame(
    fit = fit,
    lower = fit - za * se,
    upper = fit + za * se
  )
}


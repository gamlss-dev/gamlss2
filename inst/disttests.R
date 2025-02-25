library("gamlss")
library("gamlss2")
## devtools::install_github("gamlss-dev/gamlss.dist")
library("gamlss.dist")

## Function to evaluate one family.
disttest <- function(family = "NO", nobs = 1000)
{
  fam <- get(family)

  npar <- names(fam()$parameters)[unlist(fam()$parameters)]
  k <- length(npar)

  links <- unlist(fam()[paste0(npar, ".link")])
  linv <- sapply(npar, function(j) make.link2(links[paste0(j, ".link")])$linkinv)

  rfun <- get(paste0("r", family))
  par <- formals(rfun)
  par[c("n", "...")] <- NULL
  par <- lapply(par[npar], function(x) {
    if(is.call(x))
      x <- eval(x)
    return(x)
  })
  y <- rfun(nobs)

  start <- proc.time()
  b0 <- try(gamlss(y ~ 1, family = fam), silent = TRUE)
  e0 <- (proc.time() - start)["elapsed"]
  ok0 <- !inherits(b0, "try-error")

  start <- proc.time()
  b1 <- try(gamlss2(y ~ 1, family = fam), silent = TRUE)
  e1 <- (proc.time() - start)["elapsed"]
  ok1 <- !inherits(b1, "try-error")

  c0 <- if(ok0) unlist(coefAll(b0)) else rep(NA, k)
  c1 <- if(ok1) coef(b1) else rep(NA, k)

  for(i in seq_along(linv)) {
    if(ok0)
      c0[i] <- linv[[i]](c0[i])
    if(ok1)
      c1[i] <- linv[[i]](c1[i])
  }

  d0 <- if(ok0) deviance(b0) else NA
  d1 <- if(ok1) deviance(b1) else NA

  res <- data.frame(
    "family" = family,
    "parameter" = rep(npar, 2),
    "fitfun" = c(rep("gamlss", k), rep("gamlss2", k)),
    "truth" = rep(unlist(par), 2),
    "estimate" = c(c0, c1),
    "deviance" = c(rep(d0, k), rep(d1, k)),
    "elapsed" = c(rep(e0, k), rep(e1, k)),
    "failed" = c(rep(!ok0, k), rep(!ok1, k))
  )

  rownames(res) <- NULL

  return(res)
}

## Run evaluation for all families.
fams <- c(
  available_families(type = "continuous"),
  available_families(type = "discrete")
)

set.seed(123)

e <- lapply(names(fams), function(j) {
  res <- try(disttest(j), silent = TRUE)
  if(inherits(res, "try-error")) {
    res <- rep(NA, 8)
    names(res) <- c("family", "parameter", "fitfun", "truth", "estimate",
      "deviance", "elapsed", "failed")
    res <- as.data.frame(t(res))
    res$family <- j
    return(res)
  } else {
    return(res)
  }
})

ec <- do.call("rbind", e)

ok <- subset(ec, !is.na(ec$deviance))

plot(deviance ~ as.factor(fitfun), data = ok)

ok$error <- ok$truth - ok$estimate

a <- subset(ok, error < -2)

print(a)


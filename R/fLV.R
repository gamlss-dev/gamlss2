LV <- function(lagN, lagNother,
  r.link      = "identity",
  aself.link  = "identity",
  aother.link = "identity",
  sigma.link  = "log")
{
  stopifnot(
    length(lagN) == length(lagNother),
    all(lagN > 0)
  )

  lagN <- as.numeric(lagN)
  lagNother <- as.numeric(lagNother)
  llN <- log(lagN)

  link_objects <- list(
    r      = make.link2(r.link),
    aself  = make.link2(aself.link),
    aother = make.link2(aother.link),
    sigma  = make.link2(sigma.link)
  )

  link_slope <- function(parameter, name) {
    link <- link_objects[[name]]
    link$mu.eta(link$linkfun(parameter))
  }

  lv_mu <- function(par) {
    llN + par$r * (1 - par$aself  * lagN - par$aother * lagNother)
  }

  lv_residual <- function(y, par) {
    y - lv_mu(par)
  }

  fam <- list(
    family = "Lotka-Volterra Model",

    names = c(
      "r",
      "aself",
      "aother",
      "sigma"
    ),

    links = c(
      r      = r.link,
      aself  = aself.link,
      aother = aother.link,
      sigma  = sigma.link
    ),

    pdf = function(y, par, log = FALSE, ...) {
      dnorm(y, mean = lv_mu(par), sd = par$sigma, log = log)
    },
    cdf = function(y, par, ...) {
      pnorm(y, mean = lv_mu(par), sd = par$sigma, ...)
    },
    quantile = function(par, p) {
      qnorm(p, mean = lv_mu(par), sd = par$sigma)
    },

    score = list(
      r = function(par, y, ...) {
        gradient <- 1 - par$aself * lagN - par$aother * lagNother
        lv_residual(y, par) * gradient * link_slope(par$r, "r") /
          par$sigma^2
      },
      aself = function(par, y, ...) {
        gradient <- -par$r * lagN
        lv_residual(y, par) * gradient * link_slope(par$aself, "aself") /
          par$sigma^2
      },
      aother = function(par, y, ...) {
        gradient <- -par$r * lagNother
        lv_residual(y, par) * gradient * link_slope(par$aother, "aother") /
          par$sigma^2
      },
      sigma = function(par, y, ...) {
        residual <- lv_residual(y, par)
        (-1 / par$sigma + residual^2 / par$sigma^3) *
          link_slope(par$sigma, "sigma")
      }
    ),

    initialize = list(
      r = function(y, ...)
        rep(1, length(y)),
      aself = function(y, ...)
        rep(0.01, length(y)),
      aother = function(y, ...)
        rep(0.005, length(y)),
      sigma = function(y, ...)
        rep(sd(y), length(y))
    ),

    valid.response = function(y) {
      if(!is.numeric(y))
        stop("response must be numeric")
      TRUE
    },

    type = "continuous"
  )

  class(fam) <- "gamlss2.family"

  fam
}

## Testing.
if(FALSE) {
set.seed(1328)

## True species-specific parameters on the natural scale.
## For each species, 1 / aself is its single-species carrying capacity.
truth <- data.frame(
  what   = factor(c("S1", "S2")),
  r      = c(1, 0.8),
  aself  = c(1/100, 1/80),
  aother = c(1/100, 0),
  sigma  = c(0.2, 0.1)
)

## Simulate.
n_series <- 100L
time_steps <- 20L

N1 <- N2 <- matrix(
  NA_real_, nrow = time_steps, ncol = n_series
)
N1[1L, ] <- runif(n_series, 5, 80)
N2[1L, ] <- runif(n_series, 5, 70)

for(t in seq_len(time_steps - 1L)) {
  N1[t + 1L, ] <- gamlss.dist::rLOGNO2(
    n_series,
    mu = exp(
      log(N1[t, ]) +
      truth$r[1L] * (
        1 - truth$aself[1L] * N1[t, ] -
        truth$aother[1L] * N2[t, ]
      )
    ),
    sigma = truth$sigma[1L]
  )
  N2[t + 1L, ] <- gamlss.dist::rLOGNO2(
    n_series,
    mu = exp(
      log(N2[t, ]) +
      truth$r[2L] * (
        1 - truth$aself[2L] * N2[t, ] -
        truth$aother[2L] * N1[t, ]
      )
    ),
    sigma = truth$sigma[2L]
  )
}

## Plot one of the simulated communities.
plot(
  seq_len(time_steps), N1[, 1L],
  type = "l", col = 4, ylim = range(N1[, 1L], N2[, 1L]),
  ylab = "Population size", xlab = "Time",
  main = "Lotka-Volterra competition: one simulated community"
)
lines(seq_len(time_steps), N2[, 1L], col = 2)
legend(
  "topright", legend = c("Species 1", "Species 2"),
  col = c(4, 2), lty = 1
)

lag_rows <- seq_len(time_steps - 1L)
response_rows <- lag_rows + 1L
n_transitions <- length(lag_rows) * n_series

d <- data.frame(
  N = c(
    N1[response_rows, , drop = FALSE],
    N2[response_rows, , drop = FALSE]
  ),
  lagN = c(
    N1[lag_rows, , drop = FALSE],
    N2[lag_rows, , drop = FALSE]
  ),
  lagNother = c(
    N2[lag_rows, , drop = FALSE],
    N1[lag_rows, , drop = FALSE]
  ),
  what = factor(
    rep(c("S1", "S2"), each = n_transitions),
    levels = levels(truth$what)
  ),
  series = rep(
    rep(seq_len(n_series), each = length(lag_rows)),
    2L
  ),
  time = rep(
    rep(response_rows, times = n_series),
    2L
  )
)
d$logN <- log(d$N)

## Generate the family.
fLV <- LV(lagN = d$lagN, lagNother = d$lagNother)

## Formula.
f <- logN ~
    what | ## r
    what | ## aself
    what | ## aother
    what   ## sigma

## Estimate model.
m <- gamlss2(f, data = d, family = fLV)
summary(m)
plot(m)

## Estimated coefficients.
truth_coef <- c(
  "r.p.(Intercept)"      = truth$r[1L],
  "r.p.whatS2"           = diff(truth$r),
  "aself.p.(Intercept)"  = truth$aself[1L],
  "aself.p.whatS2"       = diff(truth$aself),
  "aother.p.(Intercept)" = truth$aother[1L],
  "aother.p.whatS2"      = diff(truth$aother),
  "sigma.p.(Intercept)"  = log(truth$sigma[1L]),
  "sigma.p.whatS2"       = diff(log(truth$sigma))
)
estimated_coef <- unclass(coef(m, dropall = FALSE))
coef_recovery <- data.frame(
  true = truth_coef,
  estimate = estimated_coef[names(truth_coef)]
)
coef_recovery$error <- coef_recovery$estimate - coef_recovery$true
print(round(coef_recovery, 4))

## Compare parameters on their natural scales.
estimated <- predict(m, newdata = truth["what"])
parameter_recovery <- data.frame(
  what = rep(truth$what, times = ncol(truth) - 1L),
  parameter = rep(names(truth)[-1L], each = nrow(truth)),
  true = unlist(truth[-1L], use.names = FALSE),
  estimate = unlist(estimated, use.names = FALSE)
)
print(parameter_recovery, digits = 4)

## Full MCMC.
set.seed(6020)
m2 <- mcmc(m, n.iter = 1200, burnin = 200, thin = 1)
summary(m2)
predict(m2, newdata = truth["what"])
}


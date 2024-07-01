

# gamlss2: GAMLSS Infrastructure for Flexible Distributional Regression

## Installation

The development version of *gamlss2* can be installed via

``` r
install.packages("gamlss2",
  repos = c("https://gamlss-dev.R-universe.dev", "https://cloud.R-project.org"))
```

## Overview

The primary purpose of this package is to facilitate the creation of
advanced infrastructures designed to enhance the GAMLSS modeling
framework. Notably, the *gamlss2* package represents a significant
overhaul of its predecessor, *gamlss*, with a key emphasis on improving
estimation speed and incorporating more flexible infrastructures. These
enhancements enable the seamless integration of various algorithms into
GAMLSS, including gradient boosting, Bayesian estimation, regression
trees, and forests, fostering a more versatile and powerful modeling
environment.

Moreover, the package expands its compatibility by supporting all model
terms from the base R *mgcv* package. Additionally, the *gamlss2*
package introduces the capability to accommodate more than four
parameter families. Essentially, this means that users can now specify
any type of model using these new infrastructures, making the package
highly flexible and accommodating to a wide range of modeling
requirements.

-   The main model function is `gamlss2()`.
-   The default optimizer functions is `RS()`. Optimizer functions can
    be exchanged.
-   Most important methods: `summary()`, `plot()`, `predict()`.
-   Easy development of new family objects, see `?family.gamlss2`.
-   User-specific “special” terms are possible, see `?special_terms`.

For examples, please visit the manual pages.

``` r
help(package = "gamlss2")
```

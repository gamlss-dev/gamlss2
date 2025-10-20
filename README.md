
<!-- README.md is generated from README.qmd via: quarto render README.qmd --to gfm -->

<img src="https://gamlss-dev.github.io/gamlss2/gamlss2.png" align="right" alt="gamlss2 logo" width="120" />

# gamlss2: Infrastructure for Flexible Distributional Regression

## Overview

The primary purpose of this package is to facilitate the creation of
advanced infrastructures designed to enhance the GAMLSS modeling
framework. Notably, the `gamlss2` package represents a significant
overhaul of its predecessor,
[`gamlss`](https://cran.r-project.org/package=gamlss), with a key
emphasis on improving estimation speed and incorporating more flexible
infrastructures. These enhancements enable the seamless integration of
various algorithms into GAMLSS, including gradient boosting, Bayesian
estimation, regression trees, and forests, fostering a more versatile
and powerful modeling environment.

Moreover, the package expands its compatibility by supporting all model
terms from the base R [`mgcv`](https://cran.r-project.org/package=mgcv)
package. Additionally, the `gamlss2` package introduces the capability
to accommodate more than four parameter families. Essentially, this
means that users can now specify any type of model using these new
infrastructures, making the package highly flexible and accommodating to
a wide range of modeling requirements.

- The main model function is
  [`gamlss2()`](https://gamlss-dev.github.io/gamlss2/man/gamlss2.html).
- The default optimizer functions is
  [`RS()`](https://gamlss-dev.github.io/gamlss2/man/RS_CG.html).
  Optimizer functions can be exchanged.
- Most important methods: `summary()`,
  [`plot()`](https://gamlss-dev.github.io/gamlss2/man/plots.html),
  [`predict()`](https://gamlss-dev.github.io/gamlss2/man/predict.gamlss2.html).
- Easy development of new family objects, see
  [`?gamlss2,family`](https://gamlss-dev.github.io/gamlss2/man/gamlss2.family.html).
- User-specific “special” terms are possible, see
  [`?special_terms`](https://gamlss-dev.github.io/gamlss2/man/special_terms.html).

For examples, please visit the manual pages.

``` r
help(package = "gamlss2")
```

## Installation

The development version of `gamlss2` can be installed via

``` r
install.packages("gamlss2",
  repos = c("https://gamlss-dev.R-universe.dev",
            "https://cloud.R-project.org"))
```

## Licence

The package is available under the [General Public License version
3](https://www.gnu.org/licenses/gpl-3.0.html) or [version
2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)

## Vignettes

- The basic workflow is illustrated in the vignette [First
  Steps](https://gamlss-dev.github.io/gamlss2/vignettes/firststeps.html).
- To learn about distribution families in `gamlss2`, see the vignette
  [Families](https://gamlss-dev.github.io/gamlss2/vignettes/families.html),
  which provides an overview and examples on how to implement custom
  families.
- For information on using and setting up new special model terms (e.g.,
  neural networks), the vignette
  [Specials](https://gamlss-dev.github.io/gamlss2/vignettes/specials.html)
  provides a concise introduction.
- Model assessment and calibration using the `topmodels` package are
  discussed in the
  [Evaluation](https://gamlss-dev.github.io/gamlss2/vignettes/evaluation.html)
  vignette.

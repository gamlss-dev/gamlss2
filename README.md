# gamlss2: GAMLSS Modeling with Advanced Flexible Infrastructures

## Installation

The development version of _gamlss2_ can be installed via

```{r installation-github, eval=FALSE}
devtools::install_github("gamlss-dev/gamlss2")
```

## Overview 

The primary purpose of this package is to facilitate the creation of advanced infrastructures
designed to enhance the GAMLSS modeling framework. Notably, the _gamlss2_ package represents a
significant overhaul of its predecessor, _gamlss_, with a key emphasis on improving estimation
speed and incorporating more flexible infrastructures. These enhancements enable the seamless
integration of various algorithms into GAMLSS, including gradient boosting, Bayesian estimation,
regression trees, and forests, fostering a more versatile and powerful modeling environment.

Moreover, the package expands its compatibility by supporting all model terms from the base
R _mgcv_ package. Additionally, the _gamlss2_ package introduces the capability to
accommodate more than four parameter families. Essentially, this means that users can now
specify any type of model using these new infrastructures, making the package highly
flexible and accommodating to a wide range of modeling requirements.

* The main model function is `gamlss2()`.
* The default optimizer functions is `RS()`. Optimizer functions can be exchanged.
* Most important methods: `summary()`, `plot()`, `predict()`.
* Easy development of new family objects, see `?family.gamlss2`.
* User-specific "special" terms are possible, see `?special_terms`.

For examples, please visit the manual pages.

```{r installation-github, eval=FALSE}
help(package = "gamlss2")
```

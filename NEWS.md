# gamlss2 0.1-0

* Local standard errors for linear coefficients in `summary()` and `confint()`
  now account for smooth and random-effect coefficients within each
  distributional predictor. Previously, these coefficients were treated as
  known offsets, which could severely understate standard errors. The local
  method stores corrected linear covariance once after fitting and extracts
  it directly for summaries, confidence intervals, and linear-effect plots.
  RS fits use working curvature and omit cross-parameter covariance.

* Local linear covariance now uses selected solves for small systems and
  matrix-free, preconditioned solves for larger systems. Term designs stay
  separate, sparse and binned designs retain their compact representation,
  and observation chunks bound temporary memory. Failed convergence is
  reported explicitly without falling back to a large dense inverse.

* Covariance calculations for `fk()` terms now preserve their identifying
  zero coefficient. Inference remains conditional on the selected knots.

* First version of 'gamlss2' providing a fresh reimplementaton of the classic
  'gamlss' package while being more modular and facilitating the creation of
  advanced terms and models.



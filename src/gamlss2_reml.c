#define USE_FC_LEN_T

#include <math.h>
#include <float.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#ifndef FCONE
# define FCONE
#endif

/* Coupled local REML updates for multiple penalty matrices. */
SEXP calc_smooth_reml(
  SEXP X, SEXP z, SEXP w, SEXP XWX, SEXP XWz, SEXP penalties,
  SEXP lambda, SEXP index, SEXP penalty_basis
)
{
  if(!isReal(X) || !isMatrix(X) || !isReal(XWX) || !isMatrix(XWX) ||
      !isNewList(penalties) || !isNewList(penalty_basis))
    error("invalid local REML inputs");
  int nr = nrows(X), p = ncols(X), n = length(z), m = length(penalties);
  if(nr < 1 || p < 1 || m < 2 || nrows(XWX) != p || ncols(XWX) != p ||
      !isReal(z) || !isReal(w) || XLENGTH(w) != n ||
      !isReal(XWz) || XLENGTH(XWz) != p || !isReal(lambda) ||
      XLENGTH(lambda) != m || length(penalty_basis) != m)
    error("incompatible local REML inputs");

  int binned = !isNull(index);
  if(binned) {
    if(!isInteger(index) || XLENGTH(index) != n)
      error("invalid local REML bin index");
    for(int i = 0; i < n; i++)
      if(INTEGER(index)[i] < 1 || INTEGER(index)[i] > nr)
        error("invalid local REML bin index");
  } else if(nr != n) {
    error("incompatible local REML response");
  }

  const double **Sptr = (const double **) R_alloc((size_t) m, sizeof(double *));
  const double **Tptr = (const double **) R_alloc((size_t) m, sizeof(double *));
  int r = -1;
  for(int k = 0; k < m; k++) {
    SEXP S = VECTOR_ELT(penalties, k), T = VECTOR_ELT(penalty_basis, k);
    if(!isReal(S) || !isMatrix(S) || nrows(S) != p || ncols(S) != p ||
        !isReal(T) || !isMatrix(T) || nrows(T) < 1 || nrows(T) != ncols(T) ||
        (r >= 0 && nrows(T) != r))
      error("incompatible local REML penalties");
    if(r < 0) r = nrows(T);
    Sptr[k] = REAL(S);
    Tptr[k] = REAL(T);
  }

  SEXP coefficients, fitted, P, lambdas;
  PROTECT(coefficients = allocVector(REALSXP, p));
  PROTECT(fitted = allocVector(REALSXP, n));
  PROTECT(P = allocMatrix(REALSXP, p, p));
  PROTECT(lambdas = duplicate(lambda));
  double *b = REAL(coefficients), *fit = REAL(fitted), *Pptr = REAL(P);
  double *lptr = REAL(lambdas);
  double *work = (double *) R_alloc(
    (size_t) nr + (size_t) r * r + 3 * (size_t) m, sizeof(double)
  );
  double *xb = work, *G = xb + nr, *quadratic = G + (size_t) r * r;
  double *trace_penalty = quadratic + m, *next = trace_penalty + m;
  const double *B = REAL(XWX), *c = REAL(XWz);
  const double *weights = REAL(w), *response = REAL(z);
  int N = 0;
  long double zWz = 0.0L;
  for(int i = 0; i < n; i++) {
    if(weights[i] != 0.0) N++;
    zWz += (long double) weights[i] * response[i] * response[i];
  }
  for(int k = 0; k < m; k++)
    if(!R_FINITE(lptr[k]) || lptr[k] <= 0.0)
      error("invalid local REML smoothing parameter");

  const char upper = 'U', no_transpose = 'N';
  const int increment = 1, nrhs = 1;
  const double one = 1.0, zero = 0.0;
  int info = 0;
  double edf = 0.0;
  int final = 0;

  for(int it = 0; it <= 50; it++) {
    if(it == 50) final = 1;

    for(int j = 0; j < p; j++) {
      for(int i = 0; i <= j; i++) {
        R_xlen_t ij = i + (R_xlen_t) j * p;
        long double value = B[ij];
        for(int k = 0; k < m; k++)
          value += (long double) lptr[k] * Sptr[k][ij];
        Pptr[ij] = (double) value;
      }
      b[j] = c[j];
    }
    F77_CALL(dpotrf)(&upper, &p, Pptr, &p, &info FCONE);
    if(info != 0)
      error("local REML Cholesky factorization failed (info = %d)", info);
    F77_CALL(dpotrs)(&upper, &p, &nrhs, Pptr, &p, b, &p, &info FCONE);
    if(info != 0)
      error("local REML triangular solve failed (info = %d)", info);
    F77_CALL(dpotri)(&upper, &p, Pptr, &p, &info FCONE);
    if(info != 0)
      error("local REML inverse failed (info = %d)", info);

    long double edf_value = 0.0L;
    for(int j = 0; j < p; j++) {
      for(int i = 0; i < p; i++) {
        R_xlen_t ij = i + (R_xlen_t) j * p;
        R_xlen_t ji = j + (R_xlen_t) i * p;
        edf_value += (long double) B[ij] * Pptr[i <= j ? ij : ji];
      }
    }
    edf = (double) edf_value;

    if(final) {
      F77_CALL(dgemv)(
        &no_transpose, &nr, &p, &one, REAL(X), &nr, b, &increment,
        &zero, xb, &increment FCONE
      );
      for(int i = 0; i < n; i++)
        fit[i] = xb[binned ? INTEGER(index)[i] - 1 : i];
      break;
    }

    long double rss = zWz;
    for(int j = 0; j < p; j++) {
      rss -= 2.0L * (long double) b[j] * c[j];
      for(int i = 0; i < p; i++)
        rss += (long double) b[i] * B[i + (R_xlen_t) j * p] * b[j];
    }
    if(rss <= sqrt(DBL_EPSILON) * zWz) {
      F77_CALL(dgemv)(
        &no_transpose, &nr, &p, &one, REAL(X), &nr, b, &increment,
        &zero, xb, &increment FCONE
      );
      rss = 0.0L;
      for(int i = 0; i < n; i++) {
        double residual = response[i] - xb[binned ? INTEGER(index)[i] - 1 : i];
        rss += (long double) weights[i] * residual * residual;
      }
    }
    double df = N - edf;
    if(!R_FINITE(df) || df <= sqrt(DBL_EPSILON))
      error("invalid local REML residual degrees of freedom");
    double sig2 = (double) rss / df;
    if(!R_FINITE(sig2) || sig2 < 0.0)
      error("invalid local REML working variance");

    for(int j = 0; j < r; j++) {
      for(int i = 0; i <= j; i++) {
        R_xlen_t ij = i + (R_xlen_t) j * r;
        long double value = 0.0L;
        for(int k = 0; k < m; k++)
          value += (long double) lptr[k] * Tptr[k][ij];
        G[ij] = (double) value;
      }
    }
    F77_CALL(dpotrf)(&upper, &r, G, &r, &info FCONE);
    if(info != 0)
      error("local REML penalty factorization failed (info = %d)", info);
    F77_CALL(dpotri)(&upper, &r, G, &r, &info FCONE);
    if(info != 0)
      error("local REML penalty inverse failed (info = %d)", info);

    for(int k = 0; k < m; k++) {
      long double q = 0.0L, trP = 0.0L, trG = 0.0L;
      for(int j = 0; j < p; j++) {
        long double Sbj = 0.0L;
        for(int i = 0; i < p; i++) {
          R_xlen_t ij = i + (R_xlen_t) j * p;
          R_xlen_t ji = j + (R_xlen_t) i * p;
          Sbj += (long double) Sptr[k][ij] * b[i];
          trP += (long double) Sptr[k][ij] * Pptr[i <= j ? ij : ji];
        }
        q += (long double) b[j] * Sbj;
      }
      for(int j = 0; j < r; j++) {
        for(int i = 0; i < r; i++) {
          R_xlen_t ij = i + (R_xlen_t) j * r;
          R_xlen_t ji = j + (R_xlen_t) i * r;
          trG += (long double) Tptr[k][ij] * G[i <= j ? ij : ji];
        }
      }
      quadratic[k] = (double) q;
      trace_penalty[k] = (double) (trG - trP);
    }

    double change = 0.0;
    for(int k = 0; k < m; k++) {
      double q = quadratic[k], a = trace_penalty[k];
      if(q < 0.0 && q > -sqrt(DBL_EPSILON)) q = 0.0;
      if(a < 0.0 && a > -sqrt(DBL_EPSILON)) a = 0.0;
      if(q <= 0.0) {
        next[k] = a <= 0.0 ? lptr[k] : 1e7;
      } else if(a <= 0.0 || sig2 == 0.0) {
        next[k] = 1e-7;
      } else {
        next[k] = lptr[k] * sig2 * a / q;
      }
      if(!R_FINITE(next[k])) next[k] = 1e7;
      if(next[k] < 1e-7) next[k] = 1e-7;
      if(next[k] > 1e7) next[k] = 1e7;
      double step = fabs(log(next[k]) - log(lptr[k]));
      if(step > change) change = step;
    }
    for(int k = 0; k < m; k++) lptr[k] = next[k];
    if(change < 1e-7) final = 1;
  }

  for(int j = 0; j < p; j++)
    for(int i = j + 1; i < p; i++)
      Pptr[i + (R_xlen_t) j * p] = Pptr[j + (R_xlen_t) i * p];

  SEXP rval, names;
  PROTECT(rval = allocVector(VECSXP, 6));
  PROTECT(names = allocVector(STRSXP, 6));
  SET_VECTOR_ELT(rval, 0, coefficients);
  SET_VECTOR_ELT(rval, 1, fitted);
  SET_VECTOR_ELT(rval, 2, ScalarReal(edf));
  SET_VECTOR_ELT(rval, 3, lambdas);
  SET_VECTOR_ELT(rval, 4, P);
  SET_VECTOR_ELT(rval, 5, ScalarReal(n - edf));
  const char *labels[] = {"coefficients", "fitted.values", "edf", "lambdas", "vcov", "df"};
  for(int i = 0; i < 6; i++) SET_STRING_ELT(names, i, mkChar(labels[i]));
  setAttrib(rval, R_NamesSymbol, names);
  UNPROTECT(6);
  return rval;
}

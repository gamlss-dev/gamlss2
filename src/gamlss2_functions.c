#define USE_FC_LEN_T

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
// #include <omp.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rconfig.h>

#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Complex.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>

#ifndef FCONE
# define FCONE
#endif

/* Compute reduced weights and residuals. */
void calc_Xe(SEXP ind, SEXP weights, SEXP e, SEXP xweights, SEXP xrres, SEXP order)
{
  int i;
  int j = 0;
  int n = length(ind);
  int k = 0;

  PROTECT(xrres);
  PROTECT(xrres);

  double *weightsptr = REAL(weights);
  double *eptr = REAL(e);
  double *xweightsptr = REAL(xweights);
  double *xrresptr = REAL(xrres);
  int *indptr = INTEGER(ind);
  int *orderptr = INTEGER(order);

  xweightsptr[0] = 0.0;
  xrresptr[0] = 0.0;

  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {
      ++j;
      xweightsptr[j] = 0.0;
      xrresptr[j] = 0.0;
    }

    k = orderptr[i] - 1;

    xweightsptr[j] += weightsptr[k];
    xrresptr[j] += weightsptr[k] * eptr[k];
  }

  UNPROTECT(2);
}

/* Fast block diagonal crossproduct with weights. */
SEXP calc_XWX(SEXP x, SEXP w, SEXP index)
{
  int nr = nrows(x);
  int nc = ncols(x);
  int nc_index = ncols(index);
  int i, j, k;

  double *xptr = REAL(x);
  double *wptr = REAL(w);
  int *iptr = INTEGER(index);

  SEXP rval;
  PROTECT(rval = allocMatrix(REALSXP, nc, nc));
  double *rvalptr = REAL(rval);

  for(j = 0; j < nc; j++) {
    for(k = 0; k <= j; k++) {
      rvalptr[j + k * nc] = 0.0;
      rvalptr[k + j * nc] = 0.0;
    }
  }

  for(j = 0; j < nc_index; j++) {
    for(k = 0; k < nc_index; k++) {
      for(i = 0; i < nr; i++) {
        if((iptr[i + j * nr] < 0) || (iptr[i + k * nr] < 0))
          continue;
        rvalptr[iptr[i + j * nr] - 1 + (iptr[i + k * nr] - 1) * nc] += xptr[i + (iptr[i + j * nr] - 1) * nr] * (1.0 / wptr[i]) * xptr[i + (iptr[i + k * nr] - 1) * nr];
      }
    }
  }

  UNPROTECT(1);
  return rval;
}

/* Fused dense weighted crossproducts using symmetric BLAS updates. */
SEXP calc_XWXz(SEXP x, SEXP w, SEXP z)
{
  if(!isReal(x) || !isMatrix(x))
    error("'x' must be a numeric matrix");
  if(!isReal(w) || !isReal(z))
    error("'w' and 'z' must be numeric");

  int nr = nrows(x);
  int nc = ncols(x);
  if(nr < 1 || nc < 1)
    error("'x' must have positive dimensions");
  if(XLENGTH(w) != nr || XLENGTH(z) != nr)
    error("incompatible dimensions in weighted crossproducts");

  const double *xptr = REAL(x);
  const double *wptr = REAL(w);
  const double *zptr = REAL(z);

  SEXP XWX;
  PROTECT(XWX = allocMatrix(REALSXP, nc, nc));
  double *XWXptr = REAL(XWX);

  SEXP XWz;
  PROTECT(XWz = allocVector(REALSXP, nc));
  double *XWzptr = REAL(XWz);

  SEXP zWz;
  PROTECT(zWz = allocVector(REALSXP, 1));
  long double zWzvalue = 0.0;

  /* Keep the scaled work matrix near 8 MiB, while retaining sufficiently
     large BLAS calls for tall matrices. */
  const size_t max_work = 1048576;
  int block_rows = (int) (max_work / (size_t) nc);
  if(block_rows < 1)
    block_rows = 1;
  if(block_rows > nr)
    block_rows = nr;

  size_t matrix_work = (size_t) block_rows * (size_t) nc;
  double *work = (double *) R_alloc(
    matrix_work + 2 * (size_t) block_rows, sizeof(double)
  );
  double *sqrtw = work + matrix_work;
  double *zw = sqrtw + block_rows;

  const char upper = 'U';
  const char transpose = 'T';
  const double one = 1.0;
  const int increment = 1;
  int first = 1;

  for(int start = 0; start < nr; start += block_rows) {
    int nb = nr - start;
    if(nb > block_rows)
      nb = block_rows;

    for(int i = 0; i < nb; i++) {
      double wi = wptr[start + i];
      if(!R_FINITE(wi) || wi < 0.0)
        error("'w' must contain finite nonnegative values");
      sqrtw[i] = sqrt(wi);
      double zi = zptr[start + i];
      zw[i] = zi * sqrtw[i];
      double zi2 = zi * zi;
      zWzvalue += (long double) (wi * zi2);
    }

    for(int j = 0; j < nc; j++) {
      const double *xcol = xptr + (size_t) j * nr + start;
      double *workcol = work + (size_t) j * nb;
      for(int i = 0; i < nb; i++)
        workcol[i] = xcol[i] * sqrtw[i];
    }

    double beta = first ? 0.0 : 1.0;
    F77_CALL(dsyrk)(
      &upper, &transpose, &nc, &nb, &one, work, &nb,
      &beta, XWXptr, &nc FCONE FCONE
    );
    F77_CALL(dgemv)(
      &transpose, &nb, &nc, &one, work, &nb, zw,
      &increment, &beta, XWzptr, &increment FCONE
    );
    first = 0;
  }

  /* DSYRK writes one triangle only. */
  for(int j = 0; j < nc; j++)
    for(int i = j + 1; i < nc; i++)
      XWXptr[i + (size_t) j * nc] = XWXptr[j + (size_t) i * nc];

  REAL(zWz)[0] = (double) zWzvalue;

  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(rval, 0, XWX);
  SET_VECTOR_ELT(rval, 1, XWz);
  SET_VECTOR_ELT(rval, 2, zWz);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 3));
  SET_STRING_ELT(nrval, 0, mkChar("XWX"));
  SET_STRING_ELT(nrval, 1, mkChar("XWz"));
  SET_STRING_ELT(nrval, 2, mkChar("zWz"));
  setAttrib(rval, R_NamesSymbol, nrval);

  UNPROTECT(5);
  return rval;
}

/* Fused direct smooth-fit criterion and final-state kernel. */
SEXP calc_smooth_wfit(
  SEXP XWX, SEXP XWz, SEXP penalties, SEXP lambda,
  SEXP ridge, SEXP zWz, SEXP nobs, SEXP K,
  SEXP criterion, SEXP final
)
{
  if(!isReal(XWX) || !isMatrix(XWX))
    error("'XWX' must be a numeric matrix");
  if(!isReal(XWz))
    error("'XWz' must be numeric");
  if(!isNewList(penalties))
    error("'penalties' must be a list");
  if(!isReal(lambda))
    error("'lambda' must be numeric");

  int p = nrows(XWX);
  if(p < 1 || ncols(XWX) != p)
    error("'XWX' must be a non-empty square matrix");
  if(XLENGTH(XWz) != p)
    error("incompatible dimensions for 'XWz'");

  int m = length(penalties);
  if(XLENGTH(lambda) < m)
    error("not enough smoothing parameters");

  const double *XWXptr = REAL(XWX);
  const double *XWzptr = REAL(XWz);
  const double *lptr = REAL(lambda);
  const double **Sptr = (const double **) R_alloc(
    (size_t) m, sizeof(double *)
  );

  for(int k = 0; k < m; k++) {
    SEXP Sk = VECTOR_ELT(penalties, k);
    if(!isReal(Sk) || !isMatrix(Sk) ||
        nrows(Sk) != p || ncols(Sk) != p)
      error("invalid penalty matrix");
    Sptr[k] = REAL(Sk);
  }

  double ridge_value = asReal(ridge);
  double zWz_value = asReal(zWz);
  double n_value = asReal(nobs);
  double K_value = asReal(K);
  int criterion_value = asInteger(criterion);
  int return_final = asLogical(final) == TRUE;

  if(!R_FINITE(ridge_value) || ridge_value < 0.0)
    error("'ridge' must be finite and nonnegative");
  if(!R_FINITE(zWz_value) || !R_FINITE(n_value) || n_value <= 0.0)
    error("invalid criterion constants");
  if(criterion_value < 1 || criterion_value > 5)
    error("invalid smoothness criterion");

  SEXP P = R_NilValue;
  double *Pptr;
  if(return_final) {
    PROTECT(P = allocMatrix(REALSXP, p, p));
    Pptr = REAL(P);
  } else {
    Pptr = (double *) R_alloc(
      (size_t) p * (size_t) p, sizeof(double)
    );
  }

  /* DPOTRF reads one triangle only. Assemble the upper triangle directly,
     avoiding the intermediate Sl and A matrices created by the R path. */
  for(int j = 0; j < p; j++) {
    for(int i = 0; i <= j; i++) {
      R_xlen_t ij = i + (R_xlen_t) j * p;
      double value = i == j ? ridge_value : 0.0;
      for(int k = 0; k < m; k++)
        value += lptr[k] * Sptr[k][ij];
      Pptr[ij] = XWXptr[ij] + value;
    }
  }

  const char upper = 'U';
  int info = 0;
  F77_CALL(dpotrf)(&upper, &p, Pptr, &p, &info FCONE);
  if(info != 0)
    error("native smooth-fit Cholesky factorization failed (info = %d)", info);

  SEXP coefficients = R_NilValue;
  double *bptr;
  if(return_final) {
    PROTECT(coefficients = allocVector(REALSXP, p));
    bptr = REAL(coefficients);
  } else {
    bptr = (double *) R_alloc((size_t) p, sizeof(double));
  }
  for(int i = 0; i < p; i++)
    bptr[i] = XWzptr[i];

  const int nrhs = 1;
  F77_CALL(dpotrs)(
    &upper, &p, &nrhs, Pptr, &p, bptr, &p, &info FCONE
  );
  if(info != 0)
    error("native smooth-fit triangular solve failed (info = %d)", info);

  F77_CALL(dpotri)(&upper, &p, Pptr, &p, &info FCONE);
  if(info != 0)
    error("native smooth-fit inverse failed (info = %d)", info);

  long double edf_accumulator = 0.0L;
  for(int j = 0; j < p; j++) {
    for(int i = 0; i < p; i++) {
      R_xlen_t ij = i + (R_xlen_t) j * p;
      R_xlen_t ji = j + (R_xlen_t) i * p;
      edf_accumulator +=
        (long double) XWXptr[ij] *
        (long double) Pptr[i <= j ? ij : ji];
    }
  }

  long double cross_accumulator = 0.0L;
  for(int i = 0; i < p; i++)
    cross_accumulator +=
      (long double) bptr[i] * (long double) XWzptr[i];

  double *XWXb = (double *) R_alloc((size_t) p, sizeof(double));
  const char no_transpose = 'N';
  const int increment = 1;
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemv)(
    &no_transpose, &p, &p, &one, XWXptr, &p, bptr,
    &increment, &zero, XWXb, &increment FCONE
  );

  long double quadratic_accumulator = 0.0L;
  for(int i = 0; i < p; i++)
    quadratic_accumulator +=
      (long double) bptr[i] * (long double) XWXb[i];

  double edf = (double) edf_accumulator;
  double rss = zWz_value -
    2.0 * (double) cross_accumulator +
    (double) quadratic_accumulator;

  double value;
  switch(criterion_value) {
    case 1: {
      double denominator = n_value - edf;
      value = rss * n_value / (denominator * denominator);
      break;
    }
    case 2:
      value = rss + 2.0 * edf;
      break;
    case 3:
      value = rss + K_value * edf;
      break;
    case 4:
      value = rss + 2.0 * edf +
        (2.0 * edf * (edf + 1.0)) / (n_value - edf - 1.0);
      break;
    default:
      value = rss + log(n_value) * edf;
  }

  if(!return_final)
    return ScalarReal(value);

  /* DPOTRI writes the requested triangle only. Mirror it for the vcov
     matrix returned by the final fit. */
  for(int j = 0; j < p; j++)
    for(int i = j + 1; i < p; i++)
      Pptr[i + (R_xlen_t) j * p] = Pptr[j + (R_xlen_t) i * p];

  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(rval, 0, coefficients);
  SET_VECTOR_ELT(rval, 1, ScalarReal(edf));
  SET_VECTOR_ELT(rval, 2, P);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 3));
  SET_STRING_ELT(nrval, 0, mkChar("coefficients"));
  SET_STRING_ELT(nrval, 1, mkChar("edf"));
  SET_STRING_ELT(nrval, 2, mkChar("vcov"));
  setAttrib(rval, R_NamesSymbol, nrval);

  UNPROTECT(4);
  return rval;
}

/* Compute working response and weights for the Gaussian family. */
SEXP update_Gaussian(SEXP peta, SEXP y, SEXP eta, SEXP j)
{
  int n = length(eta);
  int i;
  const char *k = CHAR(STRING_ELT(j, 0));

  if(!isReal(y)) {
    if(isInteger(y)) {
      y = coerceVector(y, REALSXP);
    } else {
      error("Argument 'y' must be numeric or integer in update().");
    }
  }

  SEXP hess;
  PROTECT(hess = allocVector(REALSXP, n));

  SEXP z;
  PROTECT(z = allocVector(REALSXP, n));

  double *zptr = REAL(z);
  double *hessptr = REAL(hess);

  double *yptr = REAL(y);
  double *etaptr = REAL(eta);
  double *muptr = REAL(VECTOR_ELT(peta, 0));
  double *sigmaptr = REAL(VECTOR_ELT(peta, 1));

  int is_mu = strcmp(k, "mu") == 0;

  if(is_mu) {
    for(i = 0; i < n; i++) {
      double s2 = sigmaptr[i] * sigmaptr[i];
      double s2_inv = 1.0 / s2;
      double score = (yptr[i] - muptr[i]) * s2_inv;
      hessptr[i] = s2_inv;
      zptr[i] = etaptr[i] + score / hessptr[i];
    }
  } else {
    for(i = 0; i < n; i++) {
      double s2 = sigmaptr[i] * sigmaptr[i];
      double ymu = yptr[i] - muptr[i];
      double ymu2 = ymu * ymu;
      double score = (ymu2 / s2) - 1.0;
      hessptr[i] = 2.0;
      zptr[i] = etaptr[i] + score / hessptr[i];
    }
  }

  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 2));

  SET_VECTOR_ELT(rval, 0, z);
  SET_VECTOR_ELT(rval, 1, hess);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 2));

  SET_STRING_ELT(nrval, 0, mkChar("eta"));
  SET_STRING_ELT(nrval, 1, mkChar("weights"));

  setAttrib(rval, R_NamesSymbol, nrval);

  UNPROTECT(4);

  return rval;
}


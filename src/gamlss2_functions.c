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


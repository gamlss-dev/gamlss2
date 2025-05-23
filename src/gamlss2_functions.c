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

/* Compute working response and weights for the Gaussian family. */
SEXP z_weights_Gaussian(SEXP y, SEXP eta, SEXP peta, SEXP j)
{
  int n = length(eta);
  int i;
  const char *k = CHAR(STRING_ELT(j, 0));

  if(!isReal(y)) {
    if(isInteger(y)) {
      y = coerceVector(y, REALSXP);
    } else {
      error("Argument 'y' must be numeric or integer in z_weights().");
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

  SET_STRING_ELT(nrval, 0, mkChar("z"));
  SET_STRING_ELT(nrval, 1, mkChar("weights"));

  setAttrib(rval, R_NamesSymbol, nrval);

  UNPROTECT(4);

  return rval;
}


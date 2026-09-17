#define USE_FC_LEN_T

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
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

/* Fast NCV score for ordered, symmetric lag neighborhoods. */
SEXP calc_ncv_lag(SEXP Q, SEXP e, SEXP lag_, SEXP target_)
{
  if(!isReal(Q) || !isMatrix(Q) || !isReal(e) || !isInteger(lag_) ||
     length(lag_) != 1 || !isInteger(target_) || length(target_) != 1)
    error("invalid arguments to calc_ncv_lag");
  const int n = nrows(Q), p = ncols(Q), lag = INTEGER(lag_)[0];
  if(length(e) != n || lag < 0 || INTEGER(target_)[0] < 0)
    error("incompatible arguments to calc_ncv_lag");
  const double *q = REAL(Q), *ep = REAL(e);
  if(lag == 0) {
    double value = 0.0;
    for(int i = 0; i < n; i++) {
      double h = 0.0;
      for(int j = 0; j < p; j++) {
        const double v = q[i + (R_xlen_t)j * n];
        h += v * v;
      }
      const double den = 1.0 - h;
      if(!R_FINITE(ep[i]) || !R_FINITE(den) || den <= sqrt(DBL_EPSILON))
        return ScalarReal(R_PosInf);
      value += (ep[i] / den) * (ep[i] / den);
    }
    return ScalarReal(R_FINITE(value) ? value : R_PosInf);
  }
  double value = 0.0;
  const char upper = 'U';
  const int one = 1;
  int info;
  for(int i = 0; i < n; i++) {
    const int first = i > lag ? i - lag : 0;
    const int last = i + lag < n - 1 ? i + lag : n - 1;
    const int m = last - first + 1;
    const int target = i - first;
    double *a = (double *) R_alloc((size_t)m * (size_t)m, sizeof(double));
    double *rhs = (double *) R_alloc((size_t)m, sizeof(double));
    for(int col = 0; col < m; col++) {
      const int row_obs = first + col;
      for(int row = 0; row <= col; row++) {
        const int col_obs = first + row;
        double h = 0.0;
        for(int j = 0; j < p; j++)
          h += q[row_obs + (R_xlen_t)j * n] * q[col_obs + (R_xlen_t)j * n];
        a[row + (R_xlen_t)col * m] = (row == col ? 1.0 : 0.0) - h;
      }
    }
    for(int row = 0; row < m; row++) {
      rhs[row] = ep[first + row];
      if(!R_FINITE(rhs[row])) return ScalarReal(R_PosInf);
    }
    F77_CALL(dpotrf)(&upper, &m, a, &m, &info FCONE);
    if(info != 0 || !R_FINITE(a[0])) return ScalarReal(R_PosInf);
    F77_CALL(dpotrs)(&upper, &m, &one, a, &m, rhs, &m, &info FCONE);
    if(info != 0 || !R_FINITE(rhs[target])) return ScalarReal(R_PosInf);
    value += rhs[target] * rhs[target];
    if(!R_FINITE(value)) return ScalarReal(R_PosInf);
  }
  return ScalarReal(value);
}

/* NCV score and exact derivatives with respect to log(lambda). */
SEXP calc_ncv_lag_gradient(SEXP Q, SEXP e, SEXP Xw, SEXP b, SEXP P,
                           SEXP penalties, SEXP lag_)
{
  if(!isReal(Q) || !isMatrix(Q) || !isReal(e) || !isReal(Xw) ||
     !isMatrix(Xw) || !isReal(b) || !isReal(P) || !isMatrix(P) ||
     !isNewList(penalties) || !isInteger(lag_) || length(lag_) != 1)
    error("invalid arguments to calc_ncv_lag_gradient");
  const int n = nrows(Q), p = ncols(Q), lag = INTEGER(lag_)[0];
  if(n <= 0 || nrows(Xw) != n || ncols(Xw) != p || length(b) != p ||
     nrows(P) != p || ncols(P) != p || lag < 1)
    error("incompatible arguments to calc_ncv_lag_gradient");
  const int K = length(penalties);
  const double *q = REAL(Q), *ep = REAL(e), *xp = REAL(Xw);
  const double *bp = REAL(b), *pptr = REAL(P);
  for(int k = 0; k < K; k++) {
    SEXP sk = VECTOR_ELT(penalties, k);
    if(!isReal(sk) || !isMatrix(sk) || nrows(sk) != p || ncols(sk) != p)
      error("invalid penalty in calc_ncv_lag_gradient");
  }
  SEXP ans, grad;
  PROTECT(ans = allocVector(VECSXP, 2));
  PROTECT(grad = allocVector(REALSXP, K));
  double *gp = REAL(grad), value = 0.0;
  for(int k = 0; k < K; k++) gp[k] = 0.0;
  double **wks = (double **) R_alloc((size_t)K, sizeof(double *));
  double **vks = (double **) R_alloc((size_t)K, sizeof(double *));
  for(int k = 0; k < K; k++) {
    const double *sk = REAL(VECTOR_ELT(penalties, k));
    double *w = (double *) R_alloc((size_t)p * (size_t)p, sizeof(double));
    double *v = (double *) R_alloc((size_t)p, sizeof(double));
    wks[k] = w; vks[k] = v;
    for(int col = 0; col < p; col++) for(int row = 0; row < p; row++) {
      double z = 0.0;
      for(int a1 = 0; a1 < p; a1++) {
        double t = 0.0;
        for(int a2 = 0; a2 < p; a2++)
          t += sk[a1 + (R_xlen_t)a2 * p] * pptr[a2 + (R_xlen_t)col * p];
        z += pptr[row + (R_xlen_t)a1 * p] * t;
      }
      w[row + (R_xlen_t)col * p] = z;
    }
    for(int row = 0; row < p; row++) {
      double z = 0.0;
      for(int col = 0; col < p; col++) z += w[row + (R_xlen_t)col * p] * bp[col];
      v[row] = z;
    }
  }
  const char upper = 'U';
  const int one = 1;
  int info;
  for(int i = 0; i < n; i++) {
    const int first = i > lag ? i - lag : 0;
    const int last = i + lag < n - 1 ? i + lag : n - 1;
    const int m = last - first + 1, target = i - first;
    double *a = (double *) R_alloc((size_t)m * (size_t)m, sizeof(double));
    double *u = (double *) R_alloc((size_t)m, sizeof(double));
    for(int col = 0; col < m; col++) {
      const int row_obs = first + col;
      for(int row = 0; row <= col; row++) {
        const int col_obs = first + row;
        double h = 0.0;
        for(int j = 0; j < p; j++)
          h += q[row_obs + (R_xlen_t)j * n] * q[col_obs + (R_xlen_t)j * n];
        a[row + (R_xlen_t)col * m] = (row == col ? 1.0 : 0.0) - h;
      }
    }
    for(int row = 0; row < m; row++) {
      u[row] = ep[first + row];
      if(!R_FINITE(u[row])) { UNPROTECT(2); return ScalarReal(R_PosInf); }
    }
    F77_CALL(dpotrf)(&upper, &m, a, &m, &info FCONE);
    if(info != 0) { UNPROTECT(2); return ScalarReal(R_PosInf); }
    F77_CALL(dpotrs)(&upper, &m, &one, a, &m, u, &m, &info FCONE);
    if(info != 0 || !R_FINITE(u[target])) { UNPROTECT(2); return ScalarReal(R_PosInf); }
    value += u[target] * u[target];
    for(int k = 0; k < K; k++) {
      const double *w = wks[k], *v = vks[k];
      double *rhs = (double *) R_alloc((size_t)m, sizeof(double));
      for(int row = 0; row < m; row++) {
        const int obs = first + row;
        double de = 0.0, dh = 0.0;
        for(int col = 0; col < p; col++) {
          de += xp[obs + (R_xlen_t)col * n] * v[col];
          double qv = 0.0;
          for(int jj = 0; jj < m; jj++) {
            const int obs2 = first + jj;
            double wx = 0.0;
            for(int cc = 0; cc < p; cc++)
              wx += w[col + (R_xlen_t)cc * p] * xp[obs2 + (R_xlen_t)cc * n];
            qv += wx * u[jj];
          }
          dh -= xp[obs + (R_xlen_t)col * n] * qv;
        }
        rhs[row] = de + dh;
      }
      F77_CALL(dpotrs)(&upper, &m, &one, a, &m, rhs, &m, &info FCONE);
      if(info != 0 || !R_FINITE(rhs[target])) { UNPROTECT(2); return ScalarReal(R_PosInf); }
      gp[k] += 2.0 * u[target] * rhs[target];
    }
  }
  SET_VECTOR_ELT(ans, 0, ScalarReal(value));
  SET_VECTOR_ELT(ans, 1, grad);
  UNPROTECT(2);
  return ans;
}

/* Compute reduced weights and residuals. */
SEXP calc_Xe(SEXP ind, SEXP weights, SEXP e, SEXP xweights, SEXP xrres, SEXP order)
{
  int i;
  int j = 0;
  int n = length(ind);
  int k = 0;

  PROTECT(xweights);
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
  return R_NilValue;
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
SEXP calc_XWXz_cached(SEXP x, SEXP w, SEXP z, SEXP cached)
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

  int reuse = !isNull(cached);
  if(reuse && (!isReal(cached) || !isMatrix(cached) ||
      nrows(cached) != nc || ncols(cached) != nc))
    error("invalid cached weighted crossproduct");

  SEXP XWX;
  PROTECT(XWX = reuse ? cached : allocMatrix(REALSXP, nc, nc));
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

  size_t matrix_work = reuse ? 0 : (size_t) block_rows * (size_t) nc;
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
      if(!reuse) sqrtw[i] = sqrt(wi);
      double zi = zptr[start + i];
      zw[i] = zi * wi;
      double zi2 = zi * zi;
      zWzvalue += (long double) (wi * zi2);
    }

    if(!reuse) {
      for(int j = 0; j < nc; j++) {
        const double *xcol = xptr + (size_t) j * nr + start;
        double *workcol = work + (size_t) j * nb;
        for(int i = 0; i < nb; i++)
          workcol[i] = xcol[i] * sqrtw[i];
      }
    }

    double beta = first ? 0.0 : 1.0;
    if(!reuse) {
      F77_CALL(dsyrk)(
        &upper, &transpose, &nc, &nb, &one, work, &nb,
        &beta, XWXptr, &nc FCONE FCONE
      );
    }
    F77_CALL(dgemv)(
      &transpose, &nb, &nc, &one, xptr + start, &nr, zw,
      &increment, &beta, XWzptr, &increment FCONE
    );
    first = 0;
  }

  /* DSYRK writes one triangle only. */
  if(!reuse) {
    for(int j = 0; j < nc; j++)
      for(int i = j + 1; i < nc; i++)
        XWXptr[i + (size_t) j * nc] = XWXptr[j + (size_t) i * nc];
  }

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

/* Preserve the original three-argument entry point. */
SEXP calc_XWXz(SEXP x, SEXP w, SEXP z)
{
  return calc_XWXz_cached(x, w, z, R_NilValue);
}

/* Stable RSS and score residuals from the observation-space fit. */
SEXP calc_smooth_residual(SEXP X, SEXP coefficients, SEXP z, SEXP w, SEXP index, SEXP derivative)
{
  if(!isReal(X) || !isMatrix(X) || !isReal(coefficients) || !isReal(z) || !isReal(w))
    error("invalid smooth residual inputs");
  int nr = nrows(X), p = ncols(X), n = length(z);
  if(nr < 1 || p < 1 || XLENGTH(coefficients) != p || XLENGTH(w) != n)
    error("incompatible smooth residual inputs");
  int binned = !isNull(index);
  if(binned) {
    if(!isInteger(index) || XLENGTH(index) != n)
      error("invalid smooth residual bin index");
    for(int i = 0; i < n; i++)
      if(INTEGER(index)[i] < 1 || INTEGER(index)[i] > nr)
        error("invalid smooth residual bin index");
  } else if(nr != n) {
    error("incompatible smooth residual response");
  }
  if(!isNull(derivative) && (!isReal(derivative) || XLENGTH(derivative) != p))
    error("invalid smooth residual derivative");
  double *xb = (double *) R_alloc((size_t) 2 * nr, sizeof(double));
  double *wr = xb + nr;
  for(int i = 0; i < nr; i++) wr[i] = 0.0;
  const char no_transpose = 'N', transpose = 'T';
  const int increment = 1;
  const double one = 1.0, zero = 0.0;
  F77_CALL(dgemv)(
    &no_transpose, &nr, &p, &one, REAL(X), &nr, REAL(coefficients),
    &increment, &zero, xb, &increment FCONE
  );
  long double rss = 0.0L;
  for(int i = 0; i < n; i++) {
    int j = binned ? INTEGER(index)[i] - 1 : i;
    double residual = xb[j] - REAL(z)[i];
    wr[j] += REAL(w)[i] * residual;
    rss += (long double) (REAL(w)[i] * (residual * residual));
  }
  SEXP residual;
  PROTECT(residual = allocVector(REALSXP, p));
  F77_CALL(dgemv)(
    &transpose, &nr, &p, &one, REAL(X), &nr, wr, &increment,
    &zero, REAL(residual), &increment FCONE
  );
  long double rss1 = 0.0L;
  if(!isNull(derivative))
    for(int i = 0; i < p; i++) rss1 += 2.0L * REAL(residual)[i] * REAL(derivative)[i];
  SEXP rval, names;
  PROTECT(rval = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(rval, 0, ScalarReal((double) rss));
  SET_VECTOR_ELT(rval, 1, residual);
  SET_VECTOR_ELT(rval, 2, isNull(derivative) ? R_NilValue : ScalarReal((double) rss1));
  PROTECT(names = allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, mkChar("rss"));
  SET_STRING_ELT(names, 1, mkChar("residual"));
  SET_STRING_ELT(names, 2, mkChar("rss1"));
  setAttrib(rval, R_NamesSymbol, names);
  UNPROTECT(3);
  return rval;
}

/* Weighted Demmler-Reinsch basis from the generalized eigenproblem. */
SEXP calc_smooth_dr(SEXP XWX, SEXP XWz, SEXP penalty, SEXP ridge, SEXP rank)
{
  if(!isReal(XWX) || !isMatrix(XWX) || nrows(XWX) < 1 ||
      ncols(XWX) != nrows(XWX))
    error("'XWX' must be a non-empty numeric square matrix");
  int p = nrows(XWX);
  if(!isReal(XWz) || XLENGTH(XWz) != p || !isReal(penalty) ||
      !isMatrix(penalty) || nrows(penalty) != p || ncols(penalty) != p)
    error("incompatible Demmler-Reinsch inputs");
  double ridge_value = asReal(ridge);
  int rank_value = asInteger(rank);
  if(!R_FINITE(ridge_value) || ridge_value < 0.0 ||
      rank_value < -1 || rank_value > p)
    error("invalid Demmler-Reinsch constants");

  SEXP T, d, h, M, g;
  PROTECT(T = allocMatrix(REALSXP, p, p));
  PROTECT(d = allocVector(REALSXP, p));
  PROTECT(h = allocVector(REALSXP, p));
  PROTECT(M = allocMatrix(REALSXP, p, p));
  PROTECT(g = allocVector(REALSXP, p));
  double *Tptr = REAL(T);
  double *G = (double *) R_alloc((size_t) p * p, sizeof(double));
  const double *B = REAL(XWX);
  const double *S = REAL(penalty);
  for(int j = 0; j < p; j++) {
    for(int i = 0; i <= j; i++) {
      R_xlen_t ij = i + (R_xlen_t) j * p;
      R_xlen_t ji = j + (R_xlen_t) i * p;
      Tptr[ij] = 0.5 * (S[ij] + S[ji]);
      G[ij] = B[ij] + (i == j ? ridge_value : 0.0);
    }
  }

  /* DSYGVD transforms S with the Cholesky factor of B + ridge I and
     returns eigenvectors normalized in that metric. No explicit inverse
     or separate products with an inverse Cholesky factor are needed. */
  int itype = 1;
  const char vectors = 'V';
  const char upper = 'U';
  int lwork = -1, liwork = -1, info = 0;
  double work_query;
  int iwork_query;
  F77_CALL(dsygvd)(
    &itype, &vectors, &upper, &p, Tptr, &p, G, &p, REAL(d),
    &work_query, &lwork, &iwork_query, &liwork, &info FCONE FCONE
  );
  if(info != 0)
    error("Demmler-Reinsch workspace query failed (info = %d)", info);
  lwork = (int) work_query;
  liwork = iwork_query;
  double *work = (double *) R_alloc((size_t) lwork, sizeof(double));
  int *iwork = (int *) R_alloc((size_t) liwork, sizeof(int));
  F77_CALL(dsygvd)(
    &itype, &vectors, &upper, &p, Tptr, &p, G, &p, REAL(d),
    work, &lwork, iwork, &liwork, &info FCONE FCONE
  );
  if(info != 0)
    error("Demmler-Reinsch eigendecomposition failed (info = %d)", info);

  double scale = 1.0;
  for(int i = 0; i < p; i++)
    if(fabs(REAL(d)[i]) > scale) scale = fabs(REAL(d)[i]);
  double tolerance = sqrt(DBL_EPSILON) * scale;
  for(int i = 0; i < p; i++) {
    if(REAL(d)[i] < -tolerance)
      error("penalty matrix is not positive semi-definite");
    /* Keep the declared null space exactly unpenalized, including at
       lambda = 1e10 where tiny positive eigenvalues would matter. */
    if(REAL(d)[i] < 0.0 || (rank_value >= 0 && i < p - rank_value)) {
      if(fabs(REAL(d)[i]) > tolerance)
        error("penalty rank does not match its null space");
      REAL(d)[i] = 0.0;
    }
  }

  const char transpose = 'T';
  const int increment = 1;
  const double one = 1.0, zero = 0.0;
  F77_CALL(dsyrk)(
    &upper, &transpose, &p, &p, &one, Tptr, &p,
    &zero, REAL(M), &p FCONE FCONE
  );
  for(int j = 0; j < p; j++) {
    REAL(h)[j] = 1.0 - ridge_value * REAL(M)[j + (R_xlen_t) j * p];
    /* For a nearly zero B metric, subtraction from one loses precision. */
    if(REAL(h)[j] < sqrt(DBL_EPSILON)) {
      double *BT = work;
      const double *column = Tptr + (R_xlen_t) j * p;
      F77_CALL(dsymv)(
        &upper, &p, &one, B, &p, column, &increment,
        &zero, BT, &increment FCONE
      );
      long double value = 0.0L;
      for(int i = 0; i < p; i++) value += (long double) column[i] * BT[i];
      REAL(h)[j] = (double) value;
    }
    for(int i = j + 1; i < p; i++)
      REAL(M)[i + (R_xlen_t) j * p] = REAL(M)[j + (R_xlen_t) i * p];
  }
  F77_CALL(dgemv)(
    &transpose, &p, &p, &one, Tptr, &p, REAL(XWz), &increment,
    &zero, REAL(g), &increment FCONE
  );

  SEXP rval, names;
  PROTECT(rval = allocVector(VECSXP, 6));
  PROTECT(names = allocVector(STRSXP, 6));
  const char *labels[] = {"T", "d", "h", "M", "ridge", "g"};
  SEXP fields[] = {T, d, h, M, ridge, g};
  for(int i = 0; i < 6; i++) {
    SET_VECTOR_ELT(rval, i, fields[i]);
    SET_STRING_ELT(names, i, mkChar(labels[i]));
  }
  setAttrib(rval, R_NamesSymbol, names);
  UNPROTECT(7);
  return rval;
}

/* Fused DR criterion and log-lambda derivative. */
SEXP calc_smooth_dr_eval(
  SEXP decomposition, SEXP lambda, SEXP zWz, SEXP nobs, SEXP K, SEXP criterion
)
{
  if(!isNewList(decomposition) || length(decomposition) < 6)
    error("invalid Demmler-Reinsch decomposition");
  SEXP d = VECTOR_ELT(decomposition, 1);
  SEXP h = VECTOR_ELT(decomposition, 2);
  SEXP M = VECTOR_ELT(decomposition, 3);
  SEXP g = VECTOR_ELT(decomposition, 5);
  int p = length(d);
  if(!isReal(d) || !isReal(h) || !isReal(g) || !isReal(M) ||
      XLENGTH(h) != p || XLENGTH(g) != p || !isMatrix(M) ||
      nrows(M) != p || ncols(M) != p)
    error("incompatible Demmler-Reinsch decomposition");
  double l = asReal(lambda), n = asReal(nobs), k = asReal(K);
  double ridge = asReal(VECTOR_ELT(decomposition, 4));
  int code = asInteger(criterion);
  if(!R_FINITE(l) || l < 0.0 || !R_FINITE(n) || n <= 0.0 ||
      code < 1 || code > 5)
    error("invalid Demmler-Reinsch criterion constants");
  double *alpha = (double *) R_alloc((size_t) 3 * p, sizeof(double));
  double *alpha1 = alpha + p;
  double *Ma = alpha1 + p;
  long double rss = asReal(zWz), edf = 0.0L, edf1 = 0.0L;
  for(int i = 0; i < p; i++) {
    double t = l * REAL(d)[i], q = 1.0 / (1.0 + t);
    if(!R_FINITE(q) || q <= 0.0)
      error("invalid smoothing parameter in Demmler-Reinsch fit");
    alpha[i] = REAL(g)[i] * q;
    alpha1[i] = -alpha[i] * t * q;
    long double residual = (long double) REAL(g)[i] - alpha[i];
    rss += residual * residual - (long double) REAL(g)[i] * REAL(g)[i];
    edf += (long double) REAL(h)[i] * q;
    edf1 -= (long double) REAL(h)[i] * q * q * t;
  }
  const char upper = 'U';
  const int increment = 1;
  const double one = 1.0, zero = 0.0;
  F77_CALL(dsymv)(
    &upper, &p, &one, REAL(M), &p, alpha, &increment,
    &zero, Ma, &increment FCONE
  );
  long double rss1 = 0.0L;
  for(int i = 0; i < p; i++) {
    rss -= (long double) ridge * alpha[i] * Ma[i];
    rss1 += 2.0L * ((long double) alpha[i] - REAL(g)[i] -
      (long double) ridge * Ma[i]) * alpha1[i];
  }
  double r = (double) rss, e = (double) edf, value;
  double rss_multiplier = 1.0, edf_multiplier;
  switch(code) {
    case 1: {
      double den = n - e;
      value = r * n / (den * den);
      rss_multiplier = n / (den * den);
      edf_multiplier = 2.0 * r * n / (den * den * den);
      break;
    }
    case 2:
      value = r + 2.0 * e;
      edf_multiplier = 2.0;
      break;
    case 3:
      value = r + k * e;
      edf_multiplier = k;
      break;
    case 4: {
      double den = n - e - 1.0;
      value = r + 2.0 * e + 2.0 * e * (e + 1.0) / den;
      edf_multiplier = 2.0 +
        2.0 * ((2.0 * e + 1.0) * den + e * (e + 1.0)) / (den * den);
      break;
    }
    default:
      value = r + log(n) * e;
      edf_multiplier = log(n);
  }
  SEXP rval, names;
  PROTECT(rval = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(rval, 0, ScalarReal(value));
  SET_VECTOR_ELT(rval, 1, ScalarReal(
    rss_multiplier * (double) rss1 + edf_multiplier * (double) edf1
  ));
  PROTECT(names = allocVector(STRSXP, 4));
  SET_STRING_ELT(names, 0, mkChar("value"));
  SET_STRING_ELT(names, 1, mkChar("gradient"));
  SET_VECTOR_ELT(rval, 2, ScalarReal(r));
  SET_VECTOR_ELT(rval, 3, ScalarReal(e));
  SET_STRING_ELT(names, 2, mkChar("rss"));
  SET_STRING_ELT(names, 3, mkChar("edf"));
  setAttrib(rval, R_NamesSymbol, names);
  UNPROTECT(2);
  return rval;
}

/* Coefficients and covariance at lambda in the DR basis. */
static double smooth_dr_fit(
  int p, const double *T, const double *d, const double *h, const double *g,
  double lambda, double *b, double *P
)
{
  double *work = (double *) R_alloc((size_t) p + (size_t) p * p, sizeof(double));
  double *alpha = work, *scaled = alpha + p;
  long double edf = 0.0L;
  for(int j = 0; j < p; j++) {
    double q = 1.0 / (1.0 + lambda * d[j]);
    if(!R_FINITE(q) || q <= 0.0)
      error("invalid smoothing parameter in Demmler-Reinsch fit");
    alpha[j] = g[j] * q;
    edf += (long double) h[j] * q;
    double scale = sqrt(q);
    for(int i = 0; i < p; i++)
      scaled[i + (R_xlen_t) j * p] = T[i + (R_xlen_t) j * p] * scale;
  }
  const char no_transpose = 'N', upper = 'U';
  const int increment = 1;
  const double one = 1.0, zero = 0.0;
  F77_CALL(dgemv)(
    &no_transpose, &p, &p, &one, T, &p, alpha, &increment,
    &zero, b, &increment FCONE
  );
  F77_CALL(dsyrk)(
    &upper, &no_transpose, &p, &p, &one, scaled, &p,
    &zero, P, &p FCONE FCONE
  );
  for(int j = 0; j < p; j++)
    for(int i = j + 1; i < p; i++)
      P[i + (R_xlen_t) j * p] = P[j + (R_xlen_t) i * p];
  return (double) edf;
}

SEXP calc_smooth_dr_fit(SEXP decomposition, SEXP lambda)
{
  if(!isNewList(decomposition) || length(decomposition) < 6)
    error("invalid Demmler-Reinsch decomposition");
  SEXP T = VECTOR_ELT(decomposition, 0), d = VECTOR_ELT(decomposition, 1);
  SEXP h = VECTOR_ELT(decomposition, 2), g = VECTOR_ELT(decomposition, 5);
  int p = length(d);
  if(p < 1 || !isReal(T) || !isMatrix(T) || nrows(T) != p || ncols(T) != p ||
      !isReal(d) || !isReal(h) || XLENGTH(h) != p || !isReal(g) || XLENGTH(g) != p)
    error("incompatible Demmler-Reinsch fit inputs");
  double l = asReal(lambda);
  if(!R_FINITE(l) || l < 0.0) error("invalid smoothing parameter");
  SEXP b, P, rval, names;
  PROTECT(b = allocVector(REALSXP, p));
  PROTECT(P = allocMatrix(REALSXP, p, p));
  double edf = smooth_dr_fit(p, REAL(T), REAL(d), REAL(h), REAL(g), l, REAL(b), REAL(P));
  PROTECT(rval = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(rval, 0, b);
  SET_VECTOR_ELT(rval, 1, ScalarReal(edf));
  SET_VECTOR_ELT(rval, 2, P);
  PROTECT(names = allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, mkChar("coefficients"));
  SET_STRING_ELT(names, 1, mkChar("edf"));
  SET_STRING_ELT(names, 2, mkChar("vcov"));
  setAttrib(rval, R_NamesSymbol, names);
  UNPROTECT(4);
  return rval;
}

/* Local ML updates with shared matrix and residual workspaces. */
SEXP calc_smooth_ml(
  SEXP X, SEXP z, SEXP w, SEXP XWX, SEXP XWz, SEXP penalty,
  SEXP lambda, SEXP null_dim, SEXP index, SEXP decomposition, SEXP penalty_root
)
{
  if(!isReal(X) || !isMatrix(X) || !isReal(XWX) || !isMatrix(XWX))
    error("invalid local ML matrices");
  int nr = nrows(X), p = ncols(X), n = length(z);
  if(nr < 1 || p < 1 || nrows(XWX) != p || ncols(XWX) != p ||
      !isReal(z) || !isReal(w) || XLENGTH(w) != n ||
      !isReal(XWz) || XLENGTH(XWz) != p || !isReal(penalty) ||
      !isMatrix(penalty) || nrows(penalty) != p || ncols(penalty) != p)
    error("incompatible local ML inputs");
  int binned = !isNull(index);
  if(binned) {
    if(!isInteger(index) || XLENGTH(index) != n)
      error("invalid local ML bin index");
    for(int i = 0; i < n; i++)
      if(INTEGER(index)[i] < 1 || INTEGER(index)[i] > nr)
        error("invalid local ML bin index");
  } else if(nr != n) {
    error("incompatible local ML response");
  }
  double l = asReal(lambda), null = asReal(null_dim);
  if(!R_FINITE(l) || l < 0.0 || !R_FINITE(null) || null < 0.0 || null > p)
    error("invalid local ML constants");
  int use_dr = !isNull(decomposition);
  const double *Tptr = NULL, *dptr = NULL, *hptr = NULL, *gptr = NULL;
  if(use_dr) {
    if(!isNewList(decomposition) || length(decomposition) < 6)
      error("invalid local ML decomposition");
    SEXP T = VECTOR_ELT(decomposition, 0);
    SEXP d = VECTOR_ELT(decomposition, 1);
    SEXP h = VECTOR_ELT(decomposition, 2);
    SEXP g = VECTOR_ELT(decomposition, 5);
    if(!isReal(T) || !isMatrix(T) || nrows(T) != p || ncols(T) != p ||
        !isReal(d) || XLENGTH(d) != p || !isReal(h) || XLENGTH(h) != p ||
        !isReal(g) || XLENGTH(g) != p || asReal(VECTOR_ELT(decomposition, 4)) != 0.0)
      error("incompatible local ML decomposition");
    Tptr = REAL(T); dptr = REAL(d); hptr = REAL(h); gptr = REAL(g);
  }

  const double *U = NULL, *values = NULL;
  if(!use_dr) {
    if(!isNewList(penalty_root) || length(penalty_root) != 2)
      error("invalid local ML penalty basis");
    SEXP d = VECTOR_ELT(penalty_root, 0), vectors = VECTOR_ELT(penalty_root, 1);
    if(!isReal(d) || XLENGTH(d) != p || !isReal(vectors) ||
        !isMatrix(vectors) || nrows(vectors) != p || ncols(vectors) != p)
      error("incompatible local ML penalty basis");
    U = REAL(vectors); values = REAL(d);
  }

  SEXP coefficients, fitted, P;
  PROTECT(coefficients = allocVector(REALSXP, p));
  PROTECT(fitted = allocVector(REALSXP, n));
  PROTECT(P = allocMatrix(REALSXP, p, p));
  double *b = REAL(coefficients), *fit = REAL(fitted), *Pptr = REAL(P);
  double *work = (double *) R_alloc((size_t) nr + 2 * (size_t) p, sizeof(double));
  double *xb = work, *Sb = xb + nr, *alpha = Sb + p;
  const double *B = REAL(XWX), *c = REAL(XWz), *S = REAL(penalty);
  const double *weights = REAL(w), *response = REAL(z);
  int N = 0;
  long double zWz = 0.0L, residual_constant = 0.0L;
  for(int i = 0; i < n; i++) {
    if(weights[i] != 0.0) N++;
    zWz += (long double) (weights[i] * (response[i] * response[i]));
  }
  if(use_dr) {
    residual_constant = zWz;
    for(int i = 0; i < p; i++)
      residual_constant -= (long double) gptr[i] * gptr[i];
  }
  const char upper = 'U', no_transpose = 'N', transpose = 'T';
  const int increment = 1, nrhs = 1;
  const double one = 1.0, zero = 0.0;
  int info = 0;
  double edf = 0.0;

  /* The last pass computes the returned state at the updated lambda.
     The existing 50 Schall updates and stopping rule are unchanged. */
  int final = 0;
  for(int it = 0; it <= 50; it++) {
    if(it == 50) final = 1;
    if(use_dr && final) {
      edf = smooth_dr_fit(p, Tptr, dptr, hptr, gptr, l, b, Pptr);
    } else if(use_dr) {
      long double edf_value = 0.0L;
      for(int i = 0; i < p; i++) {
        double q = 1.0 / (1.0 + l * dptr[i]);
        alpha[i] = gptr[i] * q;
        edf_value += (long double) hptr[i] * q;
      }
      edf = (double) edf_value;
    } else {
      for(int j = 0; j < p; j++) {
        for(int i = 0; i <= j; i++) {
          R_xlen_t ij = i + (R_xlen_t) j * p;
          Pptr[ij] = B[ij] + l * S[ij];
        }
        b[j] = c[j];
      }
      F77_CALL(dpotrf)(&upper, &p, Pptr, &p, &info FCONE);
      if(info != 0) error("local ML Cholesky factorization failed (info = %d)", info);
      F77_CALL(dpotrs)(&upper, &p, &nrhs, Pptr, &p, b, &p, &info FCONE);
      if(info != 0) error("local ML triangular solve failed (info = %d)", info);
      F77_CALL(dpotri)(&upper, &p, Pptr, &p, &info FCONE);
      if(info != 0) error("local ML inverse failed (info = %d)", info);
      long double edf_value = 0.0L;
      for(int j = 0; j < p; j++) {
        for(int i = 0; i < p; i++) {
          R_xlen_t ij = i + (R_xlen_t) j * p;
          R_xlen_t ji = j + (R_xlen_t) i * p;
          edf_value += (long double) B[ij] * Pptr[i <= j ? ij : ji];
        }
      }
      edf = (double) edf_value;
    }
    /* With a zero-ridge DR metric, RSS and b'Sb need only the spectral
       coefficients. Form observation residuals if subtraction loses precision,
       and compute the fitted values once at the final smoothing parameter. */
    long double rss = residual_constant;
    if(use_dr && !final) {
      for(int i = 0; i < p; i++) {
        long double difference = (long double) gptr[i] - alpha[i];
        rss += difference * difference;
      }
    }
    if(!use_dr || final || rss <= sqrt(DBL_EPSILON) * zWz) {
      if(use_dr && !final) {
        F77_CALL(dgemv)(
          &no_transpose, &p, &p, &one, Tptr, &p, alpha, &increment,
          &zero, b, &increment FCONE
        );
      }
      F77_CALL(dgemv)(
        &no_transpose, &nr, &p, &one, REAL(X), &nr, b, &increment,
        &zero, xb, &increment FCONE
      );
      if(final) {
        for(int i = 0; i < n; i++) fit[i] = xb[binned ? INTEGER(index)[i] - 1 : i];
        break;
      }
      rss = 0.0L;
      for(int i = 0; i < n; i++) {
        double residual = response[i] - xb[binned ? INTEGER(index)[i] - 1 : i];
        rss += (long double) (weights[i] * (residual * residual));
      }
    }
    long double quadratic = 0.0L;
    if(use_dr) {
      for(int i = 0; i < p; i++)
        quadratic += (long double) dptr[i] * alpha[i] * alpha[i];
    } else {
      F77_CALL(dgemv)(
        &transpose, &p, &p, &one, U, &p, b, &increment,
        &zero, Sb, &increment FCONE
      );
      for(int i = 0; i < p; i++)
        quadratic += (long double) values[i] * Sb[i] * Sb[i];
    }
    double sig2 = (double) rss / (N - edf);
    double tau2 = (double) quadratic / (edf - null);
    if(tau2 < 1e-7) tau2 = 1e-7;
    double old = l;
    l = sig2 / tau2;
    if(l < 1e-7) l = 1e-7;
    if(l > 1e7) l = 1e7;
    if(fabs(l - old) < 1e-7 || l > 1e10) final = 1;
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
  SET_VECTOR_ELT(rval, 3, ScalarReal(l));
  SET_VECTOR_ELT(rval, 4, P);
  SET_VECTOR_ELT(rval, 5, ScalarReal(n - edf));
  const char *labels[] = {"coefficients", "fitted.values", "edf", "lambdas", "vcov", "df"};
  for(int i = 0; i < 6; i++) SET_STRING_ELT(names, i, mkChar(labels[i]));
  setAttrib(rval, R_NamesSymbol, names);
  UNPROTECT(5);
  return rval;
}

/* Derivatives of a quadratic smoothing criterion from an existing fit. */
static SEXP smooth_wfit_gradient(
  int p, int m, const double *XWXptr, const double *XWzptr,
  const double **Sptr, const double *lptr, const double *Pptr,
  const double *bptr, const double *XWXb, double ridge_value,
  double rss, double edf, long double edf_accumulator,
  double n_value, double K_value, int criterion_value, const double *root, const double *residual
)
{
  const char upper = 'U';
  const int increment = 1;
  const double one = 1.0;
  const double zero = 0.0;
  /* P is stored in its upper triangle. DSYMM avoids mirroring it and
     computes P B P once for all penalties, sharing the fitted state. */
  double *PB = (double *) R_alloc((size_t) p * p, sizeof(double));
  double *PBP = NULL;
  long double single_edf_derivative = 0.0L;
  const char left = 'L';
  const char right = 'R';
  if(m > 1 && root != NULL) {
    /* B = R'R, so P B P = (R P)'(R P). The cached factor replaces
       two dense symmetric products with a triangular product and SYRK. */
    for(int j = 0; j < p; j++) {
      for(int i = 0; i < p; i++) {
        R_xlen_t ij = i + (R_xlen_t) j * p;
        R_xlen_t ji = j + (R_xlen_t) i * p;
        PB[ij] = Pptr[i <= j ? ij : ji];
      }
    }
    const char no_transpose = 'N';
    const char nonunit = 'N';
    const char transpose = 'T';
    F77_CALL(dtrmm)(
      &left, &upper, &no_transpose, &nonunit, &p, &p, &one,
      root, &p, PB, &p FCONE FCONE FCONE FCONE
    );
    PBP = (double *) R_alloc((size_t) p * p, sizeof(double));
    F77_CALL(dsyrk)(
      &upper, &transpose, &p, &p, &one, PB, &p,
      &zero, PBP, &p FCONE FCONE
    );
  } else {
    F77_CALL(dsymm)(
      &left, &upper, &p, &p, &one, Pptr, &p, XWXptr, &p,
      &zero, PB, &p FCONE FCONE
    );
  }
  if(m == 1) {
    /* lambda S = P^-1 - B - ridge I. This identity needs only PB,
       rather than a second matrix product for the single-penalty trace.
       Accumulate the cancellation-prone difference in long double. */
    single_edf_derivative = -edf_accumulator;
    for(int j = 0; j < p; j++) {
      for(int i = 0; i < p; i++) {
        R_xlen_t ij = i + (R_xlen_t) j * p;
        R_xlen_t ji = j + (R_xlen_t) i * p;
        single_edf_derivative += (long double) PB[ij] * (
          (long double) PB[ji] + (long double) ridge_value *
          (long double) Pptr[i <= j ? ij : ji]
        );
      }
    }
  } else if(root == NULL) {
    PBP = (double *) R_alloc((size_t) p * p, sizeof(double));
    F77_CALL(dsymm)(
      &right, &upper, &p, &p, &one, Pptr, &p, PB, &p,
      &zero, PBP, &p FCONE FCONE
    );
  }

  double rss_multiplier = 1.0;
  double edf_multiplier;
  switch(criterion_value) {
    case 1: {
      double den = n_value - edf;
      rss_multiplier = n_value / (den * den);
      edf_multiplier = 2.0 * rss * n_value / (den * den * den);
      break;
    }
    case 2: edf_multiplier = 2.0; break;
    case 3: edf_multiplier = K_value; break;
    case 4: {
      double den = n_value - edf - 1.0;
      edf_multiplier = 2.0 +
        2.0 * ((2.0 * edf + 1.0) * den + edf * (edf + 1.0)) /
        (den * den);
      break;
    }
    default: edf_multiplier = log(n_value);
  }

  SEXP gradient;
  PROTECT(gradient = allocVector(REALSXP, m));
  double *Sb = (double *) R_alloc((size_t) p, sizeof(double));
  double *PSb = (double *) R_alloc((size_t) p, sizeof(double));
  for(int k = 0; k < m; k++) {
    F77_CALL(dsymv)(
      &upper, &p, &one, Sptr[k], &p, bptr, &increment,
      &zero, Sb, &increment FCONE
    );
    F77_CALL(dsymv)(
      &upper, &p, &one, Pptr, &p, Sb, &increment,
      &zero, PSb, &increment FCONE
    );
    long double rss_derivative = 0.0L;
    long double edf_derivative = 0.0L;
    for(int i = 0; i < p; i++)
      rss_derivative += (long double) PSb[i] *
        (residual == NULL ? (long double) XWXb[i] - (long double) XWzptr[i] :
          (long double) residual[i]);
    if(m == 1) {
      edf_derivative = single_edf_derivative;
    } else {
      for(int j = 0; j < p; j++) {
        for(int i = 0; i < p; i++) {
          R_xlen_t ij = i + (R_xlen_t) j * p;
          R_xlen_t ji = j + (R_xlen_t) i * p;
          edf_derivative -= (long double) lptr[k] *
            (long double) PBP[root != NULL && i > j ? ji : ij] *
            (long double) Sptr[k][i <= j ? ij : ji];
        }
      }
    }
    REAL(gradient)[k] = -2.0 * lptr[k] * rss_multiplier *
      (double) rss_derivative + edf_multiplier * (double) edf_derivative;
  }
  UNPROTECT(1);
  return gradient;
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
  int result_mode = asInteger(final);
  int return_final = result_mode == 1 || result_mode == 3 || result_mode == 4;
  int return_gradient = result_mode == 2;

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

  SEXP gradient = R_NilValue;
  if(return_gradient || result_mode == 4) {
    PROTECT(gradient = smooth_wfit_gradient(
      p, m, XWXptr, XWzptr, Sptr, lptr, Pptr, bptr, XWXb,
      ridge_value, rss, edf, edf_accumulator, n_value, K_value, criterion_value, NULL, NULL
    ));
  }
  if(return_gradient) {
    SEXP rval, names;
    PROTECT(rval = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(rval, 0, ScalarReal(value));
    SET_VECTOR_ELT(rval, 1, gradient);
    PROTECT(names = allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("value"));
    SET_STRING_ELT(names, 1, mkChar("gradient"));
    setAttrib(rval, R_NamesSymbol, names);
    UNPROTECT(3);
    return rval;
  }

  if(!return_final)
    return ScalarReal(value);

  /* DPOTRI writes the requested triangle only. Mirror it for the vcov
     matrix returned by the final fit. */
  for(int j = 0; j < p; j++)
    for(int i = j + 1; i < p; i++)
      Pptr[i + (R_xlen_t) j * p] = Pptr[j + (R_xlen_t) i * p];

  SEXP rval;
  int fields = result_mode == 4 ? 6 : result_mode == 3 ? 5 : 3;
  PROTECT(rval = allocVector(VECSXP, fields));
  SET_VECTOR_ELT(rval, 0, coefficients);
  SET_VECTOR_ELT(rval, 1, ScalarReal(edf));
  SET_VECTOR_ELT(rval, 2, P);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, fields));
  SET_STRING_ELT(nrval, 0, mkChar("coefficients"));
  SET_STRING_ELT(nrval, 1, mkChar("edf"));
  SET_STRING_ELT(nrval, 2, mkChar("vcov"));
  if(result_mode == 3 || result_mode == 4) {
    SET_VECTOR_ELT(rval, 3, ScalarReal(value));
    SET_VECTOR_ELT(rval, 4, ScalarReal(rss));
    SET_STRING_ELT(nrval, 3, mkChar("value"));
    SET_STRING_ELT(nrval, 4, mkChar("rss"));
  }
  if(result_mode == 4) {
    SET_VECTOR_ELT(rval, 5, gradient);
    SET_STRING_ELT(nrval, 5, mkChar("gradient"));
  }
  setAttrib(rval, R_NamesSymbol, nrval);

  UNPROTECT(result_mode == 4 ? 5 : 4);
  return rval;
}

/* A lazy gradient reuses the inverse and coefficients returned in mode 3. */
SEXP calc_smooth_wfit_gradient_root(
  SEXP XWX, SEXP XWz, SEXP penalties, SEXP lambda, SEXP state,
  SEXP ridge, SEXP nobs, SEXP K, SEXP criterion, SEXP root
)
{
  if(!isReal(XWX) || !isMatrix(XWX) ||
      nrows(XWX) < 1 || ncols(XWX) != nrows(XWX))
    error("'XWX' must be a non-empty numeric square matrix");
  int p = nrows(XWX);
  if(!isReal(XWz) || XLENGTH(XWz) != p ||
      !isNewList(penalties) || !isReal(lambda))
    error("invalid smoothing gradient inputs");
  int m = length(penalties);
  if(XLENGTH(lambda) < m || !isNewList(state) || length(state) < 5)
    error("invalid smoothing gradient state");
  SEXP coefficients = VECTOR_ELT(state, 0);
  SEXP P = VECTOR_ELT(state, 2);
  if(!isReal(coefficients) || XLENGTH(coefficients) != p ||
      !isReal(P) || !isMatrix(P) || nrows(P) != p || ncols(P) != p)
    error("incompatible smoothing gradient state");
  const double **Sptr = (const double **) R_alloc((size_t) m, sizeof(double *));
  for(int k = 0; k < m; k++) {
    SEXP Sk = VECTOR_ELT(penalties, k);
    if(!isReal(Sk) || !isMatrix(Sk) || nrows(Sk) != p || ncols(Sk) != p)
      error("invalid penalty matrix");
    Sptr[k] = REAL(Sk);
  }
  double ridge_value = asReal(ridge);
  double n_value = asReal(nobs);
  int criterion_value = asInteger(criterion);
  if(!R_FINITE(ridge_value) || ridge_value < 0.0 ||
      !R_FINITE(n_value) || n_value <= 0.0 ||
      criterion_value < 1 || criterion_value > 5)
    error("invalid smoothing gradient constants");

  if(!isNull(root) && (!isReal(root) || !isMatrix(root) ||
      nrows(root) != p || ncols(root) != p))
    error("invalid smoothing gradient factor");
  const double *residual = NULL;
  SEXP names = getAttrib(state, R_NamesSymbol);
  if(TYPEOF(names) == STRSXP && XLENGTH(names) == XLENGTH(state)) {
    for(int i = 0; i < length(state); i++) {
      if(strcmp(CHAR(STRING_ELT(names, i)), "residual") == 0) {
        SEXP value = VECTOR_ELT(state, i);
        if(!isReal(value) || XLENGTH(value) != p)
          error("invalid smoothing gradient residual");
        residual = REAL(value);
        break;
      }
    }
  }
  const double *XWXptr = REAL(XWX);
  const double *Pptr = REAL(P);
  long double edf_accumulator = 0.0L;
  for(int j = 0; j < p; j++) {
    for(int i = 0; i < p; i++) {
      R_xlen_t ij = i + (R_xlen_t) j * p;
      edf_accumulator += (long double) XWXptr[ij] * (long double) Pptr[ij];
    }
  }
  double *XWXb = (double *) R_alloc((size_t) p, sizeof(double));
  const char no_transpose = 'N';
  const int increment = 1;
  const double one = 1.0;
  const double zero = 0.0;
  F77_CALL(dgemv)(
    &no_transpose, &p, &p, &one, XWXptr, &p, REAL(coefficients),
    &increment, &zero, XWXb, &increment FCONE
  );
  SEXP gradient;
  PROTECT(gradient = smooth_wfit_gradient(
    p, m, XWXptr, REAL(XWz), Sptr, REAL(lambda), Pptr,
    REAL(coefficients), XWXb, ridge_value, asReal(VECTOR_ELT(state, 4)),
    asReal(VECTOR_ELT(state, 1)), edf_accumulator,
    n_value, asReal(K), criterion_value, isNull(root) ? NULL : REAL(root), residual
  ));
  UNPROTECT(1);
  return gradient;
}

/* Preserve the gradient entry point without a cached factor. */
SEXP calc_smooth_wfit_gradient(
  SEXP XWX, SEXP XWz, SEXP penalties, SEXP lambda, SEXP state,
  SEXP ridge, SEXP nobs, SEXP K, SEXP criterion
)
{
  return calc_smooth_wfit_gradient_root(
    XWX, XWz, penalties, lambda, state, ridge, nobs, K, criterion, R_NilValue
  );
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

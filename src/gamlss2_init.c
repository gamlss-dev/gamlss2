#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define USE_FC_LEN_T

SEXP calc_Xe(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP calc_XWX(SEXP, SEXP, SEXP);
SEXP calc_XWXz(SEXP, SEXP, SEXP);
SEXP calc_XWXz_cached(SEXP, SEXP, SEXP, SEXP);
SEXP calc_smooth_wfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP calc_smooth_wfit_gradient(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP calc_smooth_dr(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP calc_smooth_dr_eval(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP calc_smooth_ml(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP calc_smooth_wfit_gradient_root(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP calc_smooth_residual(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP calc_smooth_dr_fit(SEXP, SEXP);
SEXP update_Gaussian(SEXP, SEXP, SEXP, SEXP);
SEXP calc_ncv_lag(SEXP, SEXP, SEXP, SEXP);
SEXP calc_ncv_lag_gradient(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static R_CallMethodDef callMethods[] = {
  {"calc_Xe", (DL_FUNC) &calc_Xe, 6},
  {"calc_XWX", (DL_FUNC) &calc_XWX, 3},
  {"calc_XWXz", (DL_FUNC) &calc_XWXz, 3},
  {"calc_XWXz_cached", (DL_FUNC) &calc_XWXz_cached, 4},
  {"calc_smooth_wfit", (DL_FUNC) &calc_smooth_wfit, 10},
  {"calc_smooth_wfit_gradient", (DL_FUNC) &calc_smooth_wfit_gradient, 9},
  {"calc_smooth_dr", (DL_FUNC) &calc_smooth_dr, 5},
  {"calc_smooth_dr_eval", (DL_FUNC) &calc_smooth_dr_eval, 6},
  {"calc_smooth_ml", (DL_FUNC) &calc_smooth_ml, 11},
  {"calc_smooth_wfit_gradient_root", (DL_FUNC) &calc_smooth_wfit_gradient_root, 10},
  {"calc_smooth_residual", (DL_FUNC) &calc_smooth_residual, 6},
  {"calc_smooth_dr_fit", (DL_FUNC) &calc_smooth_dr_fit, 2},
  {"update_Gaussian", (DL_FUNC) &update_Gaussian, 4},
  {"calc_ncv_lag", (DL_FUNC) &calc_ncv_lag, 4},
  {"calc_ncv_lag_gradient", (DL_FUNC) &calc_ncv_lag_gradient, 7},
  {NULL, NULL, 0}
};

void R_init_gamlss2(DllInfo* info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

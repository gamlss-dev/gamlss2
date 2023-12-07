#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define USE_FC_LEN_T

SEXP calc_Xe(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP calc_XWX(SEXP, SEXP, SEXP);

static R_CallMethodDef callMethods[] = {
  {"calc_Xe", (DL_FUNC) &calc_Xe, 6},
  {"calc_XWX", (DL_FUNC) &calc_XWX, 3},
  {NULL, NULL, 0}
};

void R_init_sourcetools(DllInfo* info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}


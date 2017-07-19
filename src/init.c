#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "IDEAL_C.h"
#include "pi.h"

static R_NativePrimitiveArgType IDEAL_t[] = {
  INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, 
  REALSXP, REALSXP, REALSXP, 
  REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, LGLSXP, LGLSXP, STRSXP, LGLSXP,
  LGLSXP, INTSXP
};

static R_NativePrimitiveArgType simpi_t[] = {
  INTSXP, INTSXP
};

static const R_CMethodDef cMethods[] = {
  {"IDEAL", (DL_FUNC) &IDEAL, 23, IDEAL_t},
  {"simpi", (DL_FUNC) &simpi, 2, simpi_t},
  {NULL, NULL, 0, NULL}
};

void R_init_pscl(DllInfo *dll)
{
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

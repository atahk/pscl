#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "IDEAL_C.h"
#include "pi.h"

/* same signature as IDEAL */
void IDEAL_v2(int *n1, int *m1, int *d1, double *y1, int *maxiter1, int *thin1,
	      int *impute1, int *mda, double *xpriormeans1,
	      double *xpriorprec1, double *bpriormeans1, double *bpriorprec1,
	      double *xstart1, double *bstart1, double *xoutput, double *boutput,
	      int *burnin1, int *usefile, int *bsave, char **filename1, int *verbose1,
	      int *limitvoters1, int *usevoter);

void IDEAL_v3(int *n1, int *m1, int *d1, double *y1, int *maxiter1, int *thin1,
	      int *impute1, int *mda, double *xpriormeans1,
	      double *xpriorprec1, double *bpriormeans1, double *bpriorprec1,
	      double *xstart1, double *bstart1, double *xoutput, double *boutput,
	      int *burnin1, int *usefile, int *bsave, char **filename1, int *verbose1,
	      int *limitvoters1, int *usevoter);

void IDEAL_v4(int *n1, int *m1, int *d1, double *y1, int *maxiter1, int *thin1,
	      int *impute1, int *mda, double *xpriormeans1,
	      double *xpriorprec1, double *bpriormeans1, double *bpriorprec1,
	      double *xstart1, double *bstart1, double *xoutput, double *boutput,
	      int *burnin1, int *usefile, int *bsave, char **filename1, int *verbose1,
	      int *limitvoters1, int *usevoter);

void IDEAL_v5(int *n1, int *m1, int *d1, double *y1, int *maxiter1, int *thin1,
	      int *impute1, int *mda, double *xpriormeans1,
	      double *xpriorprec1, double *bpriormeans1, double *bpriorprec1,
	      double *xstart1, double *bstart1, double *xoutput, double *boutput,
	      int *burnin1, int *usefile, int *bsave, char **filename1, int *verbose1,
	      int *limitvoters1, int *usevoter, int *nthreads1);

static R_NativePrimitiveArgType IDEAL_t[] = {
  INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, LGLSXP, LGLSXP, STRSXP, LGLSXP,
  LGLSXP, INTSXP
};

static R_NativePrimitiveArgType IDEAL_v5_t[] = {
  INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, LGLSXP, LGLSXP, STRSXP, LGLSXP,
  LGLSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType simpi_t[] = {
  INTSXP, INTSXP
};

static const R_CMethodDef cMethods[] = {
  {"IDEAL",    (DL_FUNC) &IDEAL,    23, IDEAL_t},
  {"IDEAL_v2", (DL_FUNC) &IDEAL_v2, 23, IDEAL_t},
  {"IDEAL_v3", (DL_FUNC) &IDEAL_v3, 23, IDEAL_t},
  {"IDEAL_v4", (DL_FUNC) &IDEAL_v4, 23, IDEAL_t},
  {"IDEAL_v5", (DL_FUNC) &IDEAL_v5, 24, IDEAL_v5_t},
  {"simpi",    (DL_FUNC) &simpi,     2, simpi_t},
  {NULL, NULL, 0, NULL}
};

void R_init_pscl(DllInfo *dll)
{
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

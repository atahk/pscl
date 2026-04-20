/****************************************************************
 ** update ideal points, conditional on ystar and beta (v2 backend)
 **
 ** Same logic as updatex(), but the posterior-mean solve and the
 ** MVN sample are done in a single Cholesky of the precision matrix
 ** via bayesreg_chol() (replacing bayesreg + rmvnorm).
 ****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>
#include "util.h"
#include "ideal.h"

void updatex_v2(double **ystar, int **ok, double **beta,
		double **x, double **xp, double **xpv,
		int n, int m, int d,
		int impute)
{
  int i, j, k, l;
  extern double **bpb, *xprior, **xpriormat, *xbar, *bpw, **w;
  extern double **v2_Lx, *v2_dx, *v2_zx, *v2_rx;

  /* LHS for x-regression: ystar plus the negative intercept */
  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      w[i][j] = ystar[i][j] + beta[j][d];
    }
  }

  if (impute == 0) {                      /* stingy filtering of missing data */
    for (i = 0; i < n; i++) {
      for (k = 0; k < d; k++) {
	bpw[k] = 0.0;
	xbar[k] = 0.0;
	xprior[k] = 0.0;
	for (l = 0; l < d; l++) {
	  bpb[k][l] = 0.0;
	  xpriormat[k][l] = 0.0;
	}
      }
      for (k = 0; k < d; k++) {
	xprior[k] = xp[i][k];
	xpriormat[k][k] = xpv[i][k];
      }
      crosscheckx(beta, w, ok, m, d, i, bpb, bpw);
      bayesreg_chol(bpb, bpw, xprior, xpriormat, xbar, x[i], d,
		    v2_Lx, v2_dx, v2_zx, v2_rx);
    }
  }

  if (impute == 1) {                      /* allow propagation of missing data */
    crossprod(beta, m, d, bpb);           /* get B'B once and reuse */
    for (i = 0; i < n; i++) {
      for (k = 0; k < d; k++) {
	bpw[k] = 0.0;
	xbar[k] = 0.0;
	xprior[k] = 0.0;
	for (l = 0; l < d; l++) {
	  xpriormat[k][l] = 0.0;
	}
      }
      for (k = 0; k < d; k++) {
	xprior[k] = xp[i][k];
	xpriormat[k][k] = xpv[i][k];
      }
      crossxyi(beta, w, m, d, i, bpw);
      bayesreg_chol(bpb, bpw, xprior, xpriormat, xbar, x[i], d,
		    v2_Lx, v2_dx, v2_zx, v2_rx);
    }
  }

  return;
}

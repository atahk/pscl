/****************************************************************
 ** update item parameters, conditional on x and ystar (v2 backend)
 **
 ** Same logic as updateb()/updatebusevoter(), but the posterior-mean
 ** solve and the MVN sample are done in a single Cholesky of the
 ** precision matrix via bayesreg_chol().
 ****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include "util.h"
#include "ideal.h"

void updateb_v2(double **ystar, int **ok, double **beta, double **xreg,
		double **bp, double **bpv,
		int n, int m, int d,
		int impute)
{
  int j, k, q;
  extern double **xpx, **bpriormat, *bprior, *bbar, *xpy;
  extern double **v2_Lb, *v2_db, *v2_zb, *v2_rb;

  q = d + 1;

  for (j = 0; j < q; j++) {
    xpy[j] = 0.0;
    for (k = 0; k < q; k++) {
      xpx[j][k] = 0.0;
      bpriormat[j][k] = 0.0;
    }
  }

  if (impute == 0) {
    for (j = 0; j < m; j++) {
      for (k = 0; k < q; k++) {
	bpriormat[k][k] = bpv[j][k];
	bprior[k] = bp[j][k];
      }
      crosscheck(xreg, ystar, ok, n, q, j, xpx, xpy);
      bayesreg_chol(xpx, xpy, bprior, bpriormat, bbar, beta[j], q,
		    v2_Lb, v2_db, v2_zb, v2_rb);
    }
  }

  if (impute == 1) {
    crossprod(xreg, n, q, xpx);
    for (j = 0; j < m; j++) {
      for (k = 0; k < q; k++) {
	bpriormat[k][k] = bpv[j][k];
	bprior[k] = bp[j][k];
      }
      crossxyj(xreg, ystar, n, q, j, xpy);
      bayesreg_chol(xpx, xpy, bprior, bpriormat, bbar, beta[j], q,
		    v2_Lb, v2_db, v2_zb, v2_rb);
    }
  }

  return;
}

void updateb_v2_usevoter(double **ystar, int **ok, double **beta, double **xreg,
			 double **bp, double **bpv,
			 int n, int m, int d,
			 int impute, int *usevoter)
{
  int j, k, q;
  extern double **xpx, **bpriormat, *bprior, *bbar, *xpy;
  extern double **v2_Lb, *v2_db, *v2_zb, *v2_rb;

  q = d + 1;

  for (j = 0; j < q; j++) {
    xpy[j] = 0.0;
    for (k = 0; k < q; k++) {
      xpx[j][k] = 0.0;
      bpriormat[j][k] = 0.0;
    }
  }

  if (impute == 0) {
    for (j = 0; j < m; j++) {
      for (k = 0; k < q; k++) {
	bpriormat[k][k] = bpv[j][k];
	bprior[k] = bp[j][k];
      }
      crosscheckusevoter(xreg, ystar, ok, n, q, j, xpx, xpy, usevoter);
      bayesreg_chol(xpx, xpy, bprior, bpriormat, bbar, beta[j], q,
		    v2_Lb, v2_db, v2_zb, v2_rb);
    }
  }

  if (impute == 1) {
    crossprodusevoter(xreg, n, q, xpx, usevoter);
    for (j = 0; j < m; j++) {
      for (k = 0; k < q; k++) {
	bpriormat[k][k] = bpv[j][k];
	bprior[k] = bp[j][k];
      }
      crossxyjusevoter(xreg, ystar, n, q, j, xpy, usevoter);
      bayesreg_chol(xpx, xpy, bprior, bpriormat, bbar, beta[j], q,
		    v2_Lb, v2_db, v2_zb, v2_rb);
    }
  }

  return;
}

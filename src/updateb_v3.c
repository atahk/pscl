/****************************************************************
 ** update item parameters, conditional on x and ystar (v3 backend)
 **
 ** v3 = v2 (fused Cholesky) + missing-data subtract trick.
 **
 ** Stingy (impute=0) rebuilds X'X for each item in v1/v2. Here we
 ** precompute the full X'X once per call and subtract contributions
 ** from missing voters for each item j. Work on the q x q part drops
 ** from O(n * m * q^2) to O(n * q^2 + N_miss * q^2).
 ****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include "util.h"
#include "ideal.h"

void updateb_v3(double **ystar, int **ok, double **beta, double **xreg,
		double **bp, double **bpv,
		int n, int m, int d,
		int impute,
		int *n_miss_j, int **miss_i_for_j,
		double **xpx_total)
{
  int i, j, k, l, q, idx, nmiss;
  double xk, *xrow, ystarij;
  extern double **xpx, **bpriormat, *bprior, *bbar, *xpy;
  extern double **v2_Lb, *v2_db, *v2_zb, *v2_rb;

  q = d + 1;

  /* Zero the persistent prior-precision matrix once; only its diagonal
     is written inside the loop (see bpriormat[k][k] = bpv[j][k]). */
  for (k = 0; k < q; k++) {
    for (l = 0; l < q; l++) bpriormat[k][l] = 0.0;
  }

  if (impute == 0) {
    crossprod(xreg, n, q, xpx_total);

    for (j = 0; j < m; j++) {
      /* xpx <- full; then subtract missing-voter contributions. */
      for (k = 0; k < q; k++) {
	for (l = 0; l < q; l++) xpx[k][l] = xpx_total[k][l];
      }
      nmiss = n_miss_j[j];
      for (idx = 0; idx < nmiss; idx++) {
	i = miss_i_for_j[j][idx];
	xrow = xreg[i];
	for (k = 0; k < q; k++) {
	  xk = xrow[k];
	  for (l = 0; l < q; l++) xpx[k][l] -= xk * xrow[l];
	}
      }

      /* xpy: observed-only sum over i. */
      for (k = 0; k < q; k++) xpy[k] = 0.0;
      for (i = 0; i < n; i++) {
	if (ok[i][j]) {
	  xrow = xreg[i];
	  ystarij = ystar[i][j];
	  for (k = 0; k < q; k++) xpy[k] += xrow[k] * ystarij;
	}
      }

      for (k = 0; k < q; k++) {
	bpriormat[k][k] = bpv[j][k];
	bprior[k] = bp[j][k];
      }

      bayesreg_chol(xpx, xpy, bprior, bpriormat, bbar, beta[j], q,
		    v2_Lb, v2_db, v2_zb, v2_rb);
    }
  }

  if (impute == 1) {
    crossprod(xreg, n, q, xpx);
    for (j = 0; j < m; j++) {
      for (k = 0; k < q; k++) xpy[k] = 0.0;
      crossxyj(xreg, ystar, n, q, j, xpy);
      for (k = 0; k < q; k++) {
	bpriormat[k][k] = bpv[j][k];
	bprior[k] = bp[j][k];
      }
      bayesreg_chol(xpx, xpy, bprior, bpriormat, bbar, beta[j], q,
		    v2_Lb, v2_db, v2_zb, v2_rb);
    }
  }

  return;
}

void updateb_v3_usevoter(double **ystar, int **ok, double **beta, double **xreg,
			 double **bp, double **bpv,
			 int n, int m, int d,
			 int impute, int *usevoter,
			 int *n_miss_j, int **miss_i_for_j,
			 double **xpx_total)
{
  int i, j, k, l, q, idx, nmiss;
  double xk, *xrow, ystarij;
  extern double **xpx, **bpriormat, *bprior, *bbar, *xpy;
  extern double **v2_Lb, *v2_db, *v2_zb, *v2_rb;

  q = d + 1;

  for (k = 0; k < q; k++) {
    for (l = 0; l < q; l++) bpriormat[k][l] = 0.0;
  }

  if (impute == 0) {
    crossprodusevoter(xreg, n, q, xpx_total, usevoter);

    for (j = 0; j < m; j++) {
      for (k = 0; k < q; k++) {
	for (l = 0; l < q; l++) xpx[k][l] = xpx_total[k][l];
      }
      nmiss = n_miss_j[j];
      for (idx = 0; idx < nmiss; idx++) {
	i = miss_i_for_j[j][idx];
	if (usevoter[i] > 0) {         /* only subtract what was added */
	  xrow = xreg[i];
	  for (k = 0; k < q; k++) {
	    xk = xrow[k];
	    for (l = 0; l < q; l++) xpx[k][l] -= xk * xrow[l];
	  }
	}
      }

      for (k = 0; k < q; k++) xpy[k] = 0.0;
      for (i = 0; i < n; i++) {
	if (ok[i][j] && (usevoter[i] > 0)) {
	  xrow = xreg[i];
	  ystarij = ystar[i][j];
	  for (k = 0; k < q; k++) xpy[k] += xrow[k] * ystarij;
	}
      }

      for (k = 0; k < q; k++) {
	bpriormat[k][k] = bpv[j][k];
	bprior[k] = bp[j][k];
      }

      bayesreg_chol(xpx, xpy, bprior, bpriormat, bbar, beta[j], q,
		    v2_Lb, v2_db, v2_zb, v2_rb);
    }
  }

  if (impute == 1) {
    crossprodusevoter(xreg, n, q, xpx, usevoter);
    for (j = 0; j < m; j++) {
      for (k = 0; k < q; k++) xpy[k] = 0.0;
      crossxyjusevoter(xreg, ystar, n, q, j, xpy, usevoter);
      for (k = 0; k < q; k++) {
	bpriormat[k][k] = bpv[j][k];
	bprior[k] = bp[j][k];
      }
      bayesreg_chol(xpx, xpy, bprior, bpriormat, bbar, beta[j], q,
		    v2_Lb, v2_db, v2_zb, v2_rb);
    }
  }

  return;
}

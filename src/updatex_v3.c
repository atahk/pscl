/****************************************************************
 ** update ideal points, conditional on ystar and beta (v3 backend)
 **
 ** v3 = v2 (fused Cholesky) + missing-data subtract trick.
 **
 ** The stingy branch (impute=0) rebuilds B'B for each legislator in
 ** the v1/v2 backends. Here we precompute the full sum sum_j beta_j
 ** beta_j' once per call (shared across all n legislators) and, for
 ** each legislator i, subtract the contributions from the j's where
 ** y[i,j] is missing. Complexity on the d x d part drops from
 ** O(n * m * d^2) to O(m * d^2 + N_miss * d^2), where N_miss is the
 ** total count of missing (i,j) cells.
 **
 ** The B'w side (per legislator) still iterates over observed j's
 ** because w[i,j] varies with i; no free lunch there.
 ****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>
#include "util.h"
#include "ideal.h"

void updatex_v3(double **ystar, int **ok, double **beta,
		double **x, double **xp, double **xpv,
		int n, int m, int d,
		int impute,
		int *n_miss_i, int **miss_j_for_i,
		double **bpb_total)
{
  int i, j, k, l, idx, nmiss;
  double bk, *betarow, wij;
  extern double **bpb, *xprior, **xpriormat, *xbar, *bpw, **w;
  extern double **v2_Lx, *v2_dx, *v2_zx, *v2_rx;

  /* LHS for x-regression: ystar plus the negative intercept */
  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      w[i][j] = ystar[i][j] + beta[j][d];
    }
  }

  if (impute == 0) {
    /* Full B'B once: shared by all legislators. */
    crossprod(beta, m, d, bpb_total);

    for (i = 0; i < n; i++) {
      /* bpb <- full; then subtract missing contributions for this i. */
      for (k = 0; k < d; k++) {
	for (l = 0; l < d; l++) {
	  bpb[k][l] = bpb_total[k][l];
	}
      }
      nmiss = n_miss_i[i];
      for (idx = 0; idx < nmiss; idx++) {
	j = miss_j_for_i[i][idx];
	betarow = beta[j];
	for (k = 0; k < d; k++) {
	  bk = betarow[k];
	  for (l = 0; l < d; l++) {
	    bpb[k][l] -= bk * betarow[l];
	  }
	}
      }

      /* bpw: observed-only sum over j. */
      for (k = 0; k < d; k++) bpw[k] = 0.0;
      for (j = 0; j < m; j++) {
	if (ok[i][j]) {
	  betarow = beta[j];
	  wij = w[i][j];
	  for (k = 0; k < d; k++) bpw[k] += betarow[k] * wij;
	}
      }

      /* Prior setup. */
      for (k = 0; k < d; k++) {
	xbar[k] = 0.0;
	xprior[k] = xp[i][k];
	for (l = 0; l < d; l++) xpriormat[k][l] = 0.0;
	xpriormat[k][k] = xpv[i][k];
      }

      bayesreg_chol(bpb, bpw, xprior, xpriormat, xbar, x[i], d,
		    v2_Lx, v2_dx, v2_zx, v2_rx);
    }
  }

  if (impute == 1) {
    /* No missing to subtract in the impute=1 regime. Same as v2. */
    crossprod(beta, m, d, bpb);
    for (i = 0; i < n; i++) {
      for (k = 0; k < d; k++) {
	bpw[k] = 0.0;
	xbar[k] = 0.0;
	xprior[k] = 0.0;
	for (l = 0; l < d; l++) xpriormat[k][l] = 0.0;
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

/***********************************************************************
** Fused Bayesian-regression + MVN sample via a single Cholesky of the
** posterior precision.  Replaces bayesreg() + rmvnorm() in the v2 backend.
**
** Given data X'X, X'y, prior mean b0 and prior precision P0:
**   P        = X'X + P0                           (posterior precision)
**   P = L L'                                      (Cholesky)
**   L u = (X'y + P0 b0),  L' m = u                (posterior mean m)
**   L' w = z,  z ~ N(0,I)                         (Cov(w) = P^{-1})
**   sample   = m + w
**
** choldc() stores L's strict-lower triangle in L[i][j] (i>j) and its
** diagonal in d[i] (see chol.c), so triangular solves read L[i][j]/d[i].
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <Rmath.h>
#include "util.h"
#include "ideal.h"
#include "chol.h"

void bayesreg_chol(double **xpx, double *xpy,
		   double *bp, double **priormat,
		   double *bpost, double *sample,
		   int p,
		   double **L, double *d, double *z, double *r)
{
  int i, j;
  double s;

  /* L <- P = X'X + P0 (full p x p; choldc only needs upper triangle) */
  for (i = 0; i < p; i++) {
    for (j = 0; j < p; j++) {
      L[i][j] = xpx[i][j] + priormat[i][j];
    }
  }

  /* r <- X'y + P0 b0 */
  for (i = 0; i < p; i++) {
    r[i] = xpy[i];
    for (j = 0; j < p; j++) {
      r[i] += priormat[i][j] * bp[j];
    }
  }

  /* Cholesky: P = L L' (diag in d, strict-lower in L[i][j], i>j) */
  choldc(L, p, d);

  /* Forward solve L u = r; store u in bpost */
  for (i = 0; i < p; i++) {
    s = r[i];
    for (j = 0; j < i; j++) s -= L[i][j] * bpost[j];
    bpost[i] = s / d[i];
  }

  /* Back solve L' mu = u; bpost now holds the posterior mean */
  for (i = p - 1; i >= 0; i--) {
    s = bpost[i];
    for (j = i + 1; j < p; j++) s -= L[j][i] * bpost[j];
    bpost[i] = s / d[i];
  }

  /* z ~ N(0,I) */
  for (i = 0; i < p; i++) z[i] = norm_rand();

  /* Back solve L' w = z; sample <- w has Cov = P^{-1} */
  for (i = p - 1; i >= 0; i--) {
    s = z[i];
    for (j = i + 1; j < p; j++) s -= L[j][i] * sample[j];
    sample[i] = s / d[i];
  }

  /* sample <- m + w */
  for (i = 0; i < p; i++) sample[i] += bpost[i];
}

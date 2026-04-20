/***********************************************************************
 ** bayesreg_chol_omp.c -- thread-safe variant of bayesreg_chol.
 **
 ** Identical math to bayesreg_chol() in bayesreg_chol.c; the only
 ** change is that the N(0,I) draw is taken from a per-thread PCG
 ** stream (pcg_norm) instead of R's global norm_rand(), so this can
 ** be called concurrently from an OpenMP parallel region.
 ************************************************************************/

#include "pcg.h"
#include "chol.h"
#include "bayesreg_chol_omp.h"

void bayesreg_chol_omp(double **xpx, double *xpy,
                       double *bp, double **priormat,
                       double *bpost, double *sample,
                       int p,
                       double **L, double *d, double *z, double *r,
                       pcg_t *rng)
{
  int i, j;
  double s;

  for (i = 0; i < p; i++) {
    for (j = 0; j < p; j++) {
      L[i][j] = xpx[i][j] + priormat[i][j];
    }
  }

  for (i = 0; i < p; i++) {
    r[i] = xpy[i];
    for (j = 0; j < p; j++) {
      r[i] += priormat[i][j] * bp[j];
    }
  }

  choldc(L, p, d);

  /* Forward solve L u = r; store u in bpost */
  for (i = 0; i < p; i++) {
    s = r[i];
    for (j = 0; j < i; j++) s -= L[i][j] * bpost[j];
    bpost[i] = s / d[i];
  }

  /* Back solve L' mu = u */
  for (i = p - 1; i >= 0; i--) {
    s = bpost[i];
    for (j = i + 1; j < p; j++) s -= L[j][i] * bpost[j];
    bpost[i] = s / d[i];
  }

  /* z ~ N(0,I) from per-thread PCG */
  for (i = 0; i < p; i++) z[i] = pcg_norm(rng);

  /* Back solve L' w = z */
  for (i = p - 1; i >= 0; i--) {
    s = z[i];
    for (j = i + 1; j < p; j++) s -= L[j][i] * sample[j];
    sample[i] = s / d[i];
  }

  for (i = 0; i < p; i++) sample[i] += bpost[i];
}

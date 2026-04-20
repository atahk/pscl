/****************************************************************
 ** updatex_v3_omp.c -- x-update for IDEAL_v5, called inside a
 ** caller-owned OpenMP parallel region.
 **
 ** Same math as updatex_v3(): full B'B once, per-legislator
 ** subtract of missing rows, observed-only B'w, then fused
 ** Cholesky/MVN sample via bayesreg_chol_omp(). The crossprod step
 ** is wrapped in `#pragma omp single`; the legislator loop is an
 ** orphaned `#pragma omp for`.
 **
 ** w = ystar + beta[:,d] is NOT materialised to memory. It is
 ** recomputed on the fly here inside the B'w reduction, where
 ** beta[j] is already in register (we read betarow[0..d-1] for the
 ** reduction, so betarow[d] is one contiguous element further). The
 ** add is free, and we save an entire n*m store in updatey plus an
 ** n*m load here.
 ****************************************************************/

#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ideal.h"
#include "pcg.h"
#include "scratch_omp.h"
#include "bayesreg_chol_omp.h"
#include "updatex_v3_omp.h"

void updatex_v3_omp(double **ystar, int **ok, double **beta,
                    double **x, double **xp, double **xpv,
                    int n, int m, int d,
                    int impute,
                    int *n_miss_i, int **miss_j_for_i,
                    double **bpb_total,
                    scratch_t *scr, pcg_t *rng, int n_threads)
{
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  scratch_t *s = &scr[tid];
  pcg_t     *r = &rng[tid];
  int i, j, k, l, idx, nmiss;
  double bk, *betarow, wij;
  (void)n_threads;

  /* Full B'B: O(m*d^2), shared read-only across the legislator loop.
     Small enough to keep serial. */
  #pragma omp single
  {
    crossprod(beta, m, d, bpb_total);
  }

  if (impute == 0) {
    #pragma omp for schedule(static)
    for (i = 0; i < n; i++) {
      for (k = 0; k < d; k++) {
        for (l = 0; l < d; l++) s->bpb[k][l] = bpb_total[k][l];
      }
      nmiss = n_miss_i[i];
      for (idx = 0; idx < nmiss; idx++) {
        j = miss_j_for_i[i][idx];
        betarow = beta[j];
        for (k = 0; k < d; k++) {
          bk = betarow[k];
          for (l = 0; l < d; l++) s->bpb[k][l] -= bk * betarow[l];
        }
      }

      /* Branchless B'w. updatey_omp writes ystar for every cell
         (observed and missing), so multiplying by ok[i][j] (0 or 1)
         is mathematically equivalent to the conditional sum, and
         lets the compiler vectorise the FMA loop -- the
         if (ok[i][j]) form had a ~30% mispredict rate on H117. */
      for (k = 0; k < d; k++) s->bpw[k] = 0.0;
      for (j = 0; j < m; j++) {
        betarow = beta[j];
        double mask = (double)ok[i][j];
        wij = mask * (ystar[i][j] + betarow[d]);
        for (k = 0; k < d; k++) s->bpw[k] += betarow[k] * wij;
      }

      for (k = 0; k < d; k++) {
        s->xbar[k]   = 0.0;
        s->xprior[k] = xp[i][k];
        for (l = 0; l < d; l++) s->xpriormat[k][l] = 0.0;
        s->xpriormat[k][k] = xpv[i][k];
      }

      bayesreg_chol_omp(s->bpb, s->bpw, s->xprior, s->xpriormat,
                        s->xbar, x[i], d,
                        s->Lx, s->dx, s->zx, s->rx, r);
    }
  }

  if (impute == 1) {
    #pragma omp for schedule(static)
    for (i = 0; i < n; i++) {
      for (k = 0; k < d; k++) {
        for (l = 0; l < d; l++) s->bpb[k][l] = bpb_total[k][l];
      }
      for (k = 0; k < d; k++) s->bpw[k] = 0.0;
      for (j = 0; j < m; j++) {
        betarow = beta[j];
        wij = ystar[i][j] + betarow[d];
        for (k = 0; k < d; k++) s->bpw[k] += betarow[k] * wij;
      }
      for (k = 0; k < d; k++) {
        s->xbar[k]   = 0.0;
        s->xprior[k] = xp[i][k];
        for (l = 0; l < d; l++) s->xpriormat[k][l] = 0.0;
        s->xpriormat[k][k] = xpv[i][k];
      }
      bayesreg_chol_omp(s->bpb, s->bpw, s->xprior, s->xpriormat,
                        s->xbar, x[i], d,
                        s->Lx, s->dx, s->zx, s->rx, r);
    }
  }
}

/****************************************************************
 ** updateb_v3_omp.c -- b-update for IDEAL_v5, called inside a
 ** caller-owned OpenMP parallel region.
 **
 ** Mirrors updateb_v3() / updateb_v3_usevoter(): full X'X once,
 ** per-item subtract of missing voters, observed-only X'ystar, then
 ** fused Cholesky/MVN sample. The `crossprod` step is wrapped in
 ** `#pragma omp single`; the item loop is an orphaned `#pragma omp
 ** for`.
 ****************************************************************/

#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ideal.h"
#include "pcg.h"
#include "scratch_omp.h"
#include "bayesreg_chol_omp.h"
#include "updateb_v3_omp.h"

void updateb_v3_omp(double **ystar, int **ok, double **beta, double **xreg,
                    double **bp, double **bpv,
                    int n, int m, int d,
                    int impute,
                    int *n_miss_j, int **miss_i_for_j,
                    double **xpx_total,
                    scratch_t *scr, pcg_t *rng, int n_threads)
{
  int q = d + 1;
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  scratch_t *s = &scr[tid];
  pcg_t     *r = &rng[tid];
  int i, j, k, l, idx, nmiss;
  double xk, *xrow, ystarij;
  (void)n_threads;

  #pragma omp single
  {
    crossprod(xreg, n, q, xpx_total);
  }

  if (impute == 0) {
    #pragma omp for schedule(static)
    for (j = 0; j < m; j++) {
      for (k = 0; k < q; k++) {
        for (l = 0; l < q; l++) s->xpx[k][l] = xpx_total[k][l];
      }
      nmiss = n_miss_j[j];
      for (idx = 0; idx < nmiss; idx++) {
        i = miss_i_for_j[j][idx];
        xrow = xreg[i];
        for (k = 0; k < q; k++) {
          xk = xrow[k];
          for (l = 0; l < q; l++) s->xpx[k][l] -= xk * xrow[l];
        }
      }

      /* Branchless X'y. ystar[i][j] is written for every cell by
         updatey_omp, so multiplying by ok[i][j] (0/1) drops the
         missing-cell contribution and lets clang vectorise the
         inner FMA loop. The if-form had a ~30% mispredict rate. */
      for (k = 0; k < q; k++) s->xpy[k] = 0.0;
      for (i = 0; i < n; i++) {
        xrow = xreg[i];
        ystarij = ystar[i][j] * (double)ok[i][j];
        for (k = 0; k < q; k++) s->xpy[k] += xrow[k] * ystarij;
      }

      for (k = 0; k < q; k++) {
        for (l = 0; l < q; l++) s->bpriormat[k][l] = 0.0;
        s->bpriormat[k][k] = bpv[j][k];
        s->bprior[k] = bp[j][k];
      }

      bayesreg_chol_omp(s->xpx, s->xpy, s->bprior, s->bpriormat,
                        s->bbar, beta[j], q,
                        s->Lb, s->db, s->zb, s->rb, r);
    }
  }

  if (impute == 1) {
    #pragma omp for schedule(static)
    for (j = 0; j < m; j++) {
      for (k = 0; k < q; k++) {
        for (l = 0; l < q; l++) s->xpx[k][l] = xpx_total[k][l];
      }
      for (k = 0; k < q; k++) s->xpy[k] = 0.0;
      for (i = 0; i < n; i++) {
        xrow = xreg[i];
        ystarij = ystar[i][j];
        for (k = 0; k < q; k++) s->xpy[k] += xrow[k] * ystarij;
      }
      for (k = 0; k < q; k++) {
        for (l = 0; l < q; l++) s->bpriormat[k][l] = 0.0;
        s->bpriormat[k][k] = bpv[j][k];
        s->bprior[k] = bp[j][k];
      }
      bayesreg_chol_omp(s->xpx, s->xpy, s->bprior, s->bpriormat,
                        s->bbar, beta[j], q,
                        s->Lb, s->db, s->zb, s->rb, r);
    }
  }
}

void updateb_v3_omp_usevoter(double **ystar, int **ok, double **beta, double **xreg,
                             double **bp, double **bpv,
                             int n, int m, int d,
                             int impute, int *usevoter,
                             int *n_miss_j, int **miss_i_for_j,
                             double **xpx_total,
                             scratch_t *scr, pcg_t *rng, int n_threads)
{
  int q = d + 1;
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  scratch_t *s = &scr[tid];
  pcg_t     *r = &rng[tid];
  int i, j, k, l, idx, nmiss;
  double xk, *xrow, ystarij;
  (void)n_threads;

  #pragma omp single
  {
    crossprodusevoter(xreg, n, q, xpx_total, usevoter);
  }

  if (impute == 0) {
    #pragma omp for schedule(static)
    for (j = 0; j < m; j++) {
      for (k = 0; k < q; k++) {
        for (l = 0; l < q; l++) s->xpx[k][l] = xpx_total[k][l];
      }
      nmiss = n_miss_j[j];
      for (idx = 0; idx < nmiss; idx++) {
        i = miss_i_for_j[j][idx];
        if (usevoter[i] > 0) {
          xrow = xreg[i];
          for (k = 0; k < q; k++) {
            xk = xrow[k];
            for (l = 0; l < q; l++) s->xpx[k][l] -= xk * xrow[l];
          }
        }
      }

      /* Branchless variant: mask = ok * (usevoter > 0). */
      for (k = 0; k < q; k++) s->xpy[k] = 0.0;
      for (i = 0; i < n; i++) {
        xrow = xreg[i];
        double mask = (double)(ok[i][j] & (usevoter[i] > 0));
        ystarij = ystar[i][j] * mask;
        for (k = 0; k < q; k++) s->xpy[k] += xrow[k] * ystarij;
      }

      for (k = 0; k < q; k++) {
        for (l = 0; l < q; l++) s->bpriormat[k][l] = 0.0;
        s->bpriormat[k][k] = bpv[j][k];
        s->bprior[k] = bp[j][k];
      }

      bayesreg_chol_omp(s->xpx, s->xpy, s->bprior, s->bpriormat,
                        s->bbar, beta[j], q,
                        s->Lb, s->db, s->zb, s->rb, r);
    }
  }

  if (impute == 1) {
    #pragma omp for schedule(static)
    for (j = 0; j < m; j++) {
      for (k = 0; k < q; k++) {
        for (l = 0; l < q; l++) s->xpx[k][l] = xpx_total[k][l];
      }
      for (k = 0; k < q; k++) s->xpy[k] = 0.0;
      for (i = 0; i < n; i++) {
        xrow = xreg[i];
        ystarij = ystar[i][j] * (double)(usevoter[i] > 0);
        for (k = 0; k < q; k++) s->xpy[k] += xrow[k] * ystarij;
      }
      for (k = 0; k < q; k++) {
        for (l = 0; l < q; l++) s->bpriormat[k][l] = 0.0;
        s->bpriormat[k][k] = bpv[j][k];
        s->bprior[k] = bp[j][k];
      }
      bayesreg_chol_omp(s->xpx, s->xpy, s->bprior, s->bpriormat,
                        s->bbar, beta[j], q,
                        s->Lb, s->db, s->zb, s->rb, r);
    }
  }
}

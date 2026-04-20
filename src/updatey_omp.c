/*****************************************************************
 ** updatey_omp.c -- ystar update for IDEAL_v5, called inside a
 ** caller-owned OpenMP parallel region.
 **
 ** Row-parallel over the n legislators via an orphaned
 ** `#pragma omp for`. Each thread uses its own pcg_t (rng[tid]) so
 ** there is no shared RNG on the hot path. When called outside a
 ** parallel region (e.g. with OpenMP disabled), the for-loop
 ** executes in a single-thread team -- i.e. serially on rng[0].
 **
 ** y is an int8 encoding {0, 1, 9=missing} rather than double: the
 ** hot inner loop is memory-bandwidth bound, so cutting y read
 ** traffic 8x matters. The LHS for the x-regression,
 **   w[i][j] = ystar[i][j] + beta[j][d],
 ** is NOT materialised to memory -- updatex_v3_omp recomputes it on
 ** the fly inside its B'w reduction, where beta[j] is already in
 ** register, so the extra add is effectively free and we save a
 ** full n*m store + another n*m load later.
 **
 ** IDEAL_v5 opens one parallel region per MCMC iteration and calls
 ** all three _omp functions inside it, so the fork/join cost is
 ** amortised across updatey + updatex + updateb instead of being
 ** paid three times per iteration.
 *****************************************************************/

#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "pcg.h"
#include "dtnorm_omp.h"
#include "updatey_omp.h"

static inline double pcg_rnorm(pcg_t *rng, double mu) {
  return mu + pcg_norm(rng);
}

void updatey_omp(double **ystar, signed char **y,
                 double **x, double **beta,
                 int n, int m, int d, pcg_t *rng, int n_threads)
{
#ifdef _OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  pcg_t *r = &rng[tid];
  int i, j, k;
  double *xrow, *brow, mu;
  signed char *yrow;
  signed char yij;
  (void)n_threads;

  #pragma omp for schedule(static)
  for (i = 0; i < n; i++) {
    xrow = x[i];
    yrow = y[i];
    for (j = 0; j < m; j++) {
      brow = beta[j];
      mu = -brow[d];
      for (k = 0; k < d; k++) mu += brow[k] * xrow[k];
      yij = yrow[j];
      if (yij == 9)
        ystar[i][j] = pcg_rnorm(r, mu);
      else
        ystar[i][j] = dtnorm_omp(r, mu, (int)yij);
    }
  }
}

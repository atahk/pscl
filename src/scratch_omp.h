/***************************************************************
 ** scratch_omp.h -- per-thread scratch buffers for IDEAL_v5's
 ** parallel updatex / updateb steps.
 **
 ** Each thread gets its own copy of the small workspace matrices
 ** and vectors previously held as file-static globals in the v1-v4
 ** backends. Buffers are tiny (sized O(d^2) and O(q^2) where q=d+1,
 ** and d is typically 1 or 2), so per-thread allocation is cheap.
 ***************************************************************/

#ifndef SCRATCH_OMP_H
#define SCRATCH_OMP_H

typedef struct {
  /* x-regression scratch (size d) */
  double **bpb, *bpw, *xbar, *xprior, **xpriormat;
  double **Lx, *dx, *zx, *rx;

  /* b-regression scratch (size q = d+1) */
  double **xpx, *xpy, *bprior, *bbar, **bpriormat;
  double **Lb, *db, *zb, *rb;
} scratch_t;

scratch_t *scratch_alloc(int n_threads, int d);
void       scratch_free (scratch_t *scr, int n_threads);

#endif

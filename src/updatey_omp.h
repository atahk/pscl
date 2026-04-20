#ifndef UPDATEY_OMP_H
#define UPDATEY_OMP_H

#include "pcg.h"

/* v5 ystar update. `y` is an int8 encoding of the original vote matrix
 * (0 = negative, 1 = positive, 9 = missing) so the inner loop reads 1 byte
 * per cell instead of 8. w = ystar + beta[:,d] is NOT materialised here;
 * updatex_v3_omp recomputes it on-the-fly inside its B'w reduction, where
 * beta[j] is already in register. */
void updatey_omp(double **ystar, signed char **y,
                 double **x, double **beta,
                 int n, int m, int d, pcg_t *rng, int n_threads);

#endif

#ifndef UPDATEX_V3_OMP_H
#define UPDATEX_V3_OMP_H

#include "pcg.h"
#include "scratch_omp.h"

void updatex_v3_omp(double **ystar, int **ok, double **beta,
                    double **x, double **xp, double **xpv,
                    int n, int m, int d,
                    int impute,
                    int *n_miss_i, int **miss_j_for_i,
                    double **bpb_total,
                    scratch_t *scr, pcg_t *rng, int n_threads);

#endif

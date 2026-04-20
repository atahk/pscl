#ifndef UPDATEB_V3_OMP_H
#define UPDATEB_V3_OMP_H

#include "pcg.h"
#include "scratch_omp.h"

void updateb_v3_omp(double **ystar, int **ok, double **beta, double **xreg,
                    double **bp, double **bpv,
                    int n, int m, int d,
                    int impute,
                    int *n_miss_j, int **miss_i_for_j,
                    double **xpx_total,
                    scratch_t *scr, pcg_t *rng, int n_threads);

void updateb_v3_omp_usevoter(double **ystar, int **ok, double **beta, double **xreg,
                             double **bp, double **bpv,
                             int n, int m, int d,
                             int impute, int *usevoter,
                             int *n_miss_j, int **miss_i_for_j,
                             double **xpx_total,
                             scratch_t *scr, pcg_t *rng, int n_threads);

#endif

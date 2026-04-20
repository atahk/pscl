#ifndef BAYESREG_CHOL_OMP_H
#define BAYESREG_CHOL_OMP_H

#include "pcg.h"

void bayesreg_chol_omp(double **xpx, double *xpy,
                       double *bp, double **priormat,
                       double *bpost, double *sample,
                       int p,
                       double **L, double *d, double *z, double *r,
                       pcg_t *rng);

#endif

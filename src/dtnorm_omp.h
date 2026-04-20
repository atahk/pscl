#ifndef DTNORM_OMP_H
#define DTNORM_OMP_H

#include "pcg.h"

/* Standard truncated-normal sample: draws Z ~ N(0,1) conditional on
 * Z > lower_bound. Three-regime Devroye/Robert (1995) rejection. */
double dtnorm_std_omp(pcg_t *rng, const double lower_bound);

/* Draws X ~ N(mu, 1) truncated to (0, inf) if `positive`, else to
 * (-inf, 0). sd is hardwired to 1: in IDEAL_v5 ystar always has unit
 * variance, so we drop the sd parameter rather than pay its multiply
 * hundreds of millions of times.
 *
 * Inlined so the call overhead disappears at the ~300M call/run rate.
 * The rejection loop itself stays out-of-line inside dtnorm_std_omp. */
static inline double dtnorm_omp(pcg_t *rng, double mu, int positive) {
  if (!positive)
    return mu - dtnorm_std_omp(rng,  mu);
  else
    return mu + dtnorm_std_omp(rng, -mu);
}

#endif

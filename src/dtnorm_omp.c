/*************************************************************************
 ** dtnorm_omp.c -- thread-safe truncated-normal sampler for IDEAL_v5.
 **
 ** Mirrors dtnorm.c's three-regime Devroye/Robert (1995) rejection
 ** scheme, but draws from a per-thread pcg_t state rather than R's
 ** global RNG, so many threads can run it concurrently.
 **
 ** dtnorm_omp() itself is `static inline` in the header; only the
 ** standard-form rejection loop lives here (too bulky to profitably
 ** inline, and the cost is dominated by the loop body's work anyway).
 *************************************************************************/

#include <math.h>
#include "pcg.h"
#include "dtnorm_omp.h"

double dtnorm_std_omp(pcg_t *rng, const double lower_bound) {
  double y, e;
  if (lower_bound < 0.0) {
    do {
      y = pcg_norm(rng);
    } while (y <= lower_bound);
    return y;
  } else if (lower_bound < 0.75) {
    do {
      y = fabs(pcg_norm(rng));
    } while (y <= lower_bound);
    return y;
  } else {
    do {
      y = pcg_exp(rng);
      e = pcg_exp(rng);
    } while (e * lower_bound * lower_bound <= 0.5 * y * y);
    return y / lower_bound + lower_bound;
  }
}

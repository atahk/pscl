/***************************************************************
 ** pcg.c -- PCG32 seeding + Ziggurat table init.
 **
 ** The hot-path primitives (pcg32_next / pcg_unif / pcg_norm /
 ** pcg_exp) live as `static inline` in pcg.h so callers across
 ** translation units can inline them without needing LTO. Only
 ** the cold seeding routine and the ziggurat tables live here.
 **
 ** PCG32 reference: Melissa O'Neill, "PCG: A Family of Simple Fast
 ** Space-Efficient Statistically Good Algorithms for Random Number
 ** Generation" (2014), XSH-RR variant with 64-bit state / 32-bit
 ** output.
 **
 ** Ziggurat reference: Marsaglia & Tsang, "The Ziggurat Method for
 ** Generating Random Variables", J. Stat. Software 5(8), 2000.
 ** Table-construction follows their original recurrence:
 **   x[N-1] = R
 **   x[i]   = sqrt(-2 ln(V/x[i+1] + f(x[i+1])))     for i = N-2..1
 **   x[0]   = V / f(R)        (covers the rectangle + tail)
 ** with R = 3.6541528853610088, V = 9.91256303526217e-3 (volume of
 ** each layer). The Ziggurat tables are filled lazily on the first
 ** pcg_seed() call; in IDEAL_v5 that loop is serial, so no race.
 ***************************************************************/

#include <math.h>
#include "pcg.h"

uint32_t pcg_zig_ki[128];
double   pcg_zig_wi[128];
double   pcg_zig_fi[128];

static int pcg_zig_inited = 0;

static void pcg_zig_init(void) {
  /* Marsaglia constants for N=128 layers. */
  const double m1 = 2147483648.0;          /* 2^31 */
  const double r  = PCG_ZIG_R;             /* base-strip right edge */
  const double v  = 9.91256303526217e-3;   /* layer volume */
  double dn = r;
  double tn = r;
  double q  = v / exp(-0.5 * dn * dn);
  int i;

  pcg_zig_ki[0]   = (uint32_t)((dn / q) * m1);
  pcg_zig_ki[1]   = 0;                      /* always rejects fast-path */
  pcg_zig_wi[0]   = q / m1;
  pcg_zig_wi[127] = dn / m1;
  pcg_zig_fi[0]   = 1.0;
  pcg_zig_fi[127] = exp(-0.5 * dn * dn);

  for (i = 126; i >= 1; i--) {
    dn = sqrt(-2.0 * log(v / dn + exp(-0.5 * dn * dn)));
    pcg_zig_ki[i + 1] = (uint32_t)((dn / tn) * m1);
    tn = dn;
    pcg_zig_fi[i] = exp(-0.5 * dn * dn);
    pcg_zig_wi[i] = dn / m1;
  }
}

void pcg_seed(pcg_t *rng, uint64_t seed, uint64_t stream) {
  if (!pcg_zig_inited) {
    pcg_zig_init();
    pcg_zig_inited = 1;
  }
  rng->state = 0u;
  rng->inc   = (stream << 1u) | 1u;   /* must be odd */
  (void)pcg32_next(rng);
  rng->state += seed;
  (void)pcg32_next(rng);
  rng->cached_norm = 0.0;
  rng->has_cached  = 0;
}

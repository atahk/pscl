/***************************************************************
 ** pcg.h -- PCG32 RNG for thread-local streams in IDEAL_v5.
 **
 ** Inline, self-contained. Seeded deterministically from R's
 ** RNG at IDEAL_v5 entry so set.seed() on the R side still
 ** drives reproducibility (given a fixed thread count).
 **
 ** pcg_norm() uses the Marsaglia & Tsang (2000) Ziggurat method
 ** with 128 layers. The fast path -- ~99% of draws -- is a layer
 ** lookup, scaled multiply, and integer compare; no log/sin/cos/
 ** sqrt are touched. The slow paths (tail and wedge, together <1%)
 ** still touch libm but are rare enough not to matter. Tables are
 ** built lazily on the first pcg_seed() call.
 **
 ** The hot-path primitives (pcg32_next, pcg_unif, pcg_norm,
 ** pcg_exp) are defined `static inline` here so callers in other
 ** translation units can inline them without needing LTO. Seeding
 ** (pcg_seed) and table init stay out-of-line in pcg.c.
 ***************************************************************/

#ifndef PCG_H
#define PCG_H

#include <stdint.h>
#include <math.h>

typedef struct {
  uint64_t state;
  uint64_t inc;         /* stream selector, must be odd */
  double   cached_norm; /* unused since switch to Ziggurat; kept for ABI */
  int      has_cached;
} pcg_t;

/* Cold: seed the state from a 64-bit seed and a 64-bit stream id.
 * Also initialises the ziggurat tables on first call. */
void pcg_seed(pcg_t *rng, uint64_t seed, uint64_t stream);

/* Ziggurat tables, populated by pcg_seed(). 128 layers + r constant. */
extern uint32_t pcg_zig_ki[128];
extern double   pcg_zig_wi[128];
extern double   pcg_zig_fi[128];
#define PCG_ZIG_R 3.6541528853610088   /* right edge of base strip */

/* Core 32-bit PCG advance: update state, return permuted output. */
static inline uint32_t pcg32_next(pcg_t *rng) {
  uint64_t oldstate = rng->state;
  rng->state = oldstate * 6364136223846793005ULL + rng->inc;
  uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
  uint32_t rot = (uint32_t)(oldstate >> 59u);
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

/* Uniform(0,1) with 32-bit precision, strictly in the open interval.
 * (u + 0.5) / 2^32 maps {0..2^32-1} -> (0, 1) exclusive on both ends, so
 * log(pcg_unif(...)) is always safe. */
static inline double pcg_unif(pcg_t *rng) {
  uint32_t u = pcg32_next(rng);
  return ((double)u + 0.5) * (1.0 / 4294967296.0);
}

/* Standard normal via Ziggurat (Marsaglia & Tsang 2000, 128 layers).
 * Fast path: one pcg32_next, one table lookup, one multiply, one
 * unsigned compare. ~99% of draws return here. The remaining ~1%
 * fall into either the tail (iz==0; exp-rejection) or a wedge
 * (compare against the true density). */
static inline double pcg_norm(pcg_t *rng) {
  for (;;) {
    uint32_t u = pcg32_next(rng);
    int32_t  hz = (int32_t)u;             /* signed; sign bit drives result sign */
    uint32_t iz = u & 127u;                /* layer index (0..127) */
    /* Branchless |hz| as uint32_t. Avoids UB at hz == INT32_MIN
     * (negation of unsigned is well-defined modulo 2^32). */
    uint32_t az = (uint32_t)hz;
    if (hz < 0) az = -az;
    if (az < pcg_zig_ki[iz]) {
      return (double)hz * pcg_zig_wi[iz];
    }
    if (iz == 0) {
      /* Tail: sample N(0,1) | |x| > R via Marsaglia exp-rejection. */
      double x, y;
      do {
        x = -log(pcg_unif(rng)) * (1.0 / PCG_ZIG_R);
        y = -log(pcg_unif(rng));
      } while (y + y < x * x);
      return (hz > 0) ? (PCG_ZIG_R + x) : -(PCG_ZIG_R + x);
    }
    /* Wedge between layer iz-1 and the density curve. */
    double x = (double)hz * pcg_zig_wi[iz];
    if (pcg_zig_fi[iz] + pcg_unif(rng) * (pcg_zig_fi[iz - 1] - pcg_zig_fi[iz])
        < exp(-0.5 * x * x))
      return x;
    /* else: rejected, loop. */
  }
}

/* Standard exponential, rate 1: -log(U). pcg_unif is in (0,1) so log is
 * safe (argument never hits 0). */
static inline double pcg_exp(pcg_t *rng) {
  return -log(pcg_unif(rng));
}

#endif

/****************************************************************
 ** IDEAL_v5: experimental backend for ideal() with OpenMP-level
 ** parallelism inside a single chain.
 **
 ** v5 = v4 (fused Cholesky + missing-data subtract + contiguous
 ** storage) with all three Gibbs updates parallelised across
 ** legislators/items via OpenMP. Each thread carries its own PCG32
 ** RNG state, seeded deterministically from R's RNG at entry, so
 ** set.seed() still drives reproducibility *for a fixed thread
 ** count*. Different thread counts will produce different streams
 ** (but equally valid posterior draws).
 **
 ** One `#pragma omp parallel` region is opened per MCMC iteration
 ** and the three _omp updates share it, so fork/join is paid once
 ** per iter rather than per update.
 **
 ** The y vote matrix is converted to an int8 encoding (0/1/9) at
 ** entry to cut read bandwidth 8x in updatey's hot loop, which is
 ** otherwise memory-bandwidth bound. w = ystar + beta[:,d] is NOT
 ** materialised -- updatex_v3_omp folds the `+ beta[j][d]` into its
 ** B'w reduction where beta[j] is already in register.
 ****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/* omp.h must precede Rinternals.h: Rinternals.h #defines `match` as a macro,
 * which clashes with the `match` clause keyword in OpenMP 5.1's omp.h. This
 * file does not use SEXP/Rinternals facilities, so we simply do not include
 * Rinternals.h at all. */
#ifdef _OPENMP
#include <omp.h>
#endif
#include <R.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "ideal.h"
#include "pcg.h"
#include "scratch_omp.h"
#include "updatey_omp.h"
#include "updatex_v3_omp.h"
#include "updateb_v3_omp.h"

static double **v5_bpb_total;
static double **v5_xpx_total;
static int  *v5_n_miss_i;
static int **v5_miss_j_for_i;
static int  *v5_n_miss_j;
static int **v5_miss_i_for_j;

static void v5_build_missing_index(int **ok, int n, int m)
{
  int i, j, cnt, idx;

  v5_n_miss_i = ivector(n);
  v5_n_miss_j = ivector(m);
  for (i = 0; i < n; i++) v5_n_miss_i[i] = 0;
  for (j = 0; j < m; j++) v5_n_miss_j[j] = 0;

  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      if (!ok[i][j]) {
        v5_n_miss_i[i]++;
        v5_n_miss_j[j]++;
      }
    }
  }

  v5_miss_j_for_i = (int **) calloc(n, sizeof(int *));
  for (i = 0; i < n; i++) {
    cnt = v5_n_miss_i[i];
    if (cnt > 0) {
      v5_miss_j_for_i[i] = ivector(cnt);
      idx = 0;
      for (j = 0; j < m; j++) {
        if (!ok[i][j]) v5_miss_j_for_i[i][idx++] = j;
      }
    } else {
      v5_miss_j_for_i[i] = NULL;
    }
  }

  v5_miss_i_for_j = (int **) calloc(m, sizeof(int *));
  for (j = 0; j < m; j++) {
    cnt = v5_n_miss_j[j];
    if (cnt > 0) {
      v5_miss_i_for_j[j] = ivector(cnt);
      idx = 0;
      for (i = 0; i < n; i++) {
        if (!ok[i][j]) v5_miss_i_for_j[j][idx++] = i;
      }
    } else {
      v5_miss_i_for_j[j] = NULL;
    }
  }
}

static void v5_free_missing_index(int n, int m)
{
  int i, j;
  for (i = 0; i < n; i++) {
    if (v5_miss_j_for_i[i] != NULL) free(v5_miss_j_for_i[i]);
  }
  free(v5_miss_j_for_i);
  for (j = 0; j < m; j++) {
    if (v5_miss_i_for_j[j] != NULL) free(v5_miss_i_for_j[j]);
  }
  free(v5_miss_i_for_j);
  free(v5_n_miss_i);
  free(v5_n_miss_j);
}

/* Draw a uint64_t seed for each thread from R's RNG. Using two 32-bit
 * draws keeps the seed well-mixed across set.seed() values. */
static pcg_t *v5_seed_rng_array(int n_threads)
{
  pcg_t *rng = (pcg_t *) calloc((size_t)n_threads, sizeof(pcg_t));
  int t;
  uint64_t seed, stream;
  for (t = 0; t < n_threads; t++) {
    uint64_t hi = (uint64_t)(unif_rand() * 4294967296.0);
    uint64_t lo = (uint64_t)(unif_rand() * 4294967296.0);
    seed = (hi << 32) ^ lo;
    stream = (uint64_t)t + 1ULL;
    pcg_seed(&rng[t], seed, stream);
  }
  return rng;
}

void IDEAL_v5(int *n1, int *m1, int *d1, double *y1, int *maxiter1, int *thin1,
              int *impute1, int *mda, double *xpriormeans1,
              double *xpriorprec1, double *bpriormeans1, double *bpriorprec1,
              double *xstart1, double *bstart1, double *xoutput, double *boutput,
              int *burnin1, int *usefile, int *bsave, char **filename1, int *verbose1,
              int *limitvoters1, int *usevoter, int *nthreads1)
{
  int e, xocursor, bocursor, xlength, blength, q, iter;
  int inloop, **ok, burnin, n, m, d, nm, maxiter, thin, impute, verbose, n_threads;
  double **ystar, **x, **xreg, **y, **beta, **bp, **bpv, nm_doub, iterPerCent;
  double **xp, **xpv, *xtemp, *btemp;
  signed char **y8;
  signed char *y8_data;
  pcg_t *rng;
  scratch_t *scr;
  FILE *ofp;
  int i, j;

  n = *n1;
  m = *m1;
  d = *d1;
  maxiter = *maxiter1;
  thin = *thin1;
  impute = *impute1;
  verbose = *verbose1;
  burnin = *burnin1;
  n_threads = *nthreads1;
  if (n_threads < 1) n_threads = 1;

  iter = 0;
  nm = n * m;
  q = d + 1;

  y     = dmatrix_contig(n, m);
  ystar = dmatrix_contig(n, m);

  beta = dmatrix_contig(m, q);
  bp   = dmatrix_contig(m, q);
  bpv  = dmatrix_contig(m, q);

  x    = dmatrix_contig(n, d);
  xreg = dmatrix_contig(n, q);
  xp   = dmatrix_contig(n, d);
  xpv  = dmatrix_contig(n, d);

  ok = imatrix_contig(n, m);

  if (*usefile == 1) {
    ofp = fopen(R_ExpandFileName(*filename1), "a");
    if (ofp == NULL) calcerror("Can't open outfile file!\n");
  }

  GetRNGstate();

  dvecTOdmat(y1, y, n, m);
  dvecTOdmat(bpriormeans1, bp, m, q);
  dvecTOdmat(bpriorprec1, bpv, m, q);
  dvecTOdmat(xpriormeans1, xp, n, d);
  dvecTOdmat(xpriorprec1, xpv, n, d);
  dvecTOdmat(xstart1, x, n, d);
  dvecTOdmat(bstart1, beta, m, q);

  /* Seed per-thread RNG from R's RNG before we stop using it in the
   * parallel region. Serial sections below (updatex_v3, updateb_v3)
   * still use R's RNG, which is fine. */
  rng = v5_seed_rng_array(n_threads);

  xtemp = dvector(n * d);
  xlength = n * d;
  if (burnin == 0) {
    if (*usefile == 1) {
    } else {
      xocursor = n * d - 1;
      dmatTOdvec(xoutput, x, n, d);
    }
  } else {
    xocursor = -1;
  }

  btemp = dvector(m * q);
  blength = m * q;
  if (burnin == 0) {
    if (*bsave == 1) {
      if (*usefile == 1) {
      } else {
        bocursor = m * q - 1;
        dmatTOdvec(boutput, beta, m, q);
      }
    }
  } else {
    bocursor = -1;
  }

  nm_doub = check(y, ok, n, m);
  (void) nm_doub;
  (void) nm;

  v5_build_missing_index(ok, n, m);

  /* Build int8 encoding of y for updatey's hot loop: 0 = negative,
   * 1 = positive, 9 = missing. Cuts the y read traffic 8x vs the
   * double matrix. After this, the double y is no longer needed. */
  y8_data = (signed char *) calloc((size_t)n * (size_t)m, sizeof(signed char));
  y8      = (signed char **) calloc((size_t)n, sizeof(signed char *));
  for (i = 0; i < n; i++) y8[i] = y8_data + (size_t)i * (size_t)m;
  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      double v = y[i][j];
      if (v == 9.0)     y8[i][j] = 9;
      else if (v > 0.0) y8[i][j] = 1;
      else              y8[i][j] = 0;
    }
  }
  free_dmatrix_contig(y);
  y = NULL;

  v5_bpb_total = dmatrix_contig(d, d);
  v5_xpx_total = dmatrix_contig(q, q);

  scr = scratch_alloc(n_threads, d);

  /* Optional per-step profiling timers. Enabled via env var
     PSCL_PROFILE=1; otherwise a no-op. Report at end of run. */
  double t_y = 0.0, t_x = 0.0, t_xreg = 0.0, t_b = 0.0;
  int profile = 0;
  {
    const char *pf = getenv("PSCL_PROFILE");
    if (pf != NULL && pf[0] == '1') profile = 1;
  }

  while (iter < maxiter) {

    for (inloop = 0; inloop < thin; inloop++) {
      iter++;

      if (verbose) {
        iterPerCent = iter / (maxiter * 1.0) * 20.0;
        if (iterPerCent == round(iterPerCent)) {
          Rprintf("\nCurrent Iteration: %d (%.0lf%% of %d iterations requested)",
                  iter,
                  round(iterPerCent * 5.0),
                  maxiter);
        }
      }

      if (iter > maxiter)
        break;

      /* One parallel region per iteration: fork/join paid once,
         amortised across updatey + updatex + makexreg + updateb.
         Internal `#pragma omp for` and `#pragma omp single` steps
         are separated by implicit barriers so data dependencies
         are preserved without further regions. */
#ifdef _OPENMP
      double ts0 = profile ? omp_get_wtime() : 0.0;
      double ts1 = 0.0, ts2 = 0.0, ts3 = 0.0, ts4 = 0.0;
#endif
      #pragma omp parallel num_threads(n_threads)
      {
        updatey_omp(ystar, y8, x, beta, n, m, d, rng, n_threads);
#ifdef _OPENMP
        #pragma omp master
        { if (profile) ts1 = omp_get_wtime(); }
#endif

        updatex_v3_omp(ystar, ok, beta, x, xp, xpv, n, m, d, impute,
                       v5_n_miss_i, v5_miss_j_for_i, v5_bpb_total,
                       scr, rng, n_threads);
#ifdef _OPENMP
        #pragma omp master
        { if (profile) ts2 = omp_get_wtime(); }
#endif

        #pragma omp single
        {
          makexreg(xreg, x, n, d, q);
        }
#ifdef _OPENMP
        #pragma omp master
        { if (profile) ts3 = omp_get_wtime(); }
#endif

        if (*limitvoters1 > 0)
          updateb_v3_omp_usevoter(ystar, ok, beta, xreg, bp, bpv, n, m, d,
                                  impute, usevoter,
                                  v5_n_miss_j, v5_miss_i_for_j, v5_xpx_total,
                                  scr, rng, n_threads);
        else
          updateb_v3_omp(ystar, ok, beta, xreg, bp, bpv, n, m, d, impute,
                         v5_n_miss_j, v5_miss_i_for_j, v5_xpx_total,
                         scr, rng, n_threads);
#ifdef _OPENMP
        #pragma omp master
        { if (profile) ts4 = omp_get_wtime(); }
#endif
      }

#ifdef _OPENMP
      if (profile) {
        t_y    += ts1 - ts0;
        t_x    += ts2 - ts1;
        t_xreg += ts3 - ts2;
        t_b    += ts4 - ts3;
      }
#endif

      R_CheckUserInterrupt();
    }

    if (iter >= burnin) {
      if (*usefile == 1) {
        dmatTOdvec(xtemp, x, n, d);
        fprintf(ofp, "%d", iter);
        for (e = 0; e < xlength; e++) fprintf(ofp, ",%f", xtemp[e]);
        if (*bsave != 1) fprintf(ofp, "\n");
      } else {
        dmatTOdvec(xtemp, x, n, d);
        for (e = 0; e < xlength; e++) {
          xocursor++;
          xoutput[xocursor] = xtemp[e];
        }
      }

      if (*bsave == 1) {
        if (*usefile == 1) {
          dmatTOdvec(btemp, beta, m, q);
          for (e = 0; e < blength; e++) fprintf(ofp, ",%f", btemp[e]);
          fprintf(ofp, "\n");
        } else {
          dmatTOdvec(btemp, beta, m, q);
          for (e = 0; e < blength; e++) {
            bocursor++;
            boutput[bocursor] = btemp[e];
          }
        }
      }
    }
  }

  PutRNGstate();

  if (*usefile == 1) fclose(ofp);

  /* double y was freed right after the int8 conversion; free y8 here. */
  free(y8);
  free(y8_data);
  free_dmatrix_contig(ystar);
  free_dmatrix_contig(beta);
  free_dmatrix_contig(bp);
  free_dmatrix_contig(bpv);
  free_dmatrix_contig(x);
  free_dmatrix_contig(xreg);
  free_dmatrix_contig(xp);
  free_dmatrix_contig(xpv);
  free_imatrix_contig(ok);

  free(xtemp);
  free(btemp);

  scratch_free(scr, n_threads);

  free_dmatrix_contig(v5_bpb_total);
  free_dmatrix_contig(v5_xpx_total);
  v5_free_missing_index(n, m);

  free(rng);

  if (profile) {
    double tot = t_y + t_x + t_xreg + t_b;
    Rprintf("\n[PSCL_PROFILE] n_threads=%d   totals (s):\n", n_threads);
    Rprintf("  updatey  %.3f  (%.1f%%)\n", t_y,    100.0 * t_y    / tot);
    Rprintf("  updatex  %.3f  (%.1f%%)\n", t_x,    100.0 * t_x    / tot);
    Rprintf("  makexreg %.3f  (%.1f%%)\n", t_xreg, 100.0 * t_xreg / tot);
    Rprintf("  updateb  %.3f  (%.1f%%)\n", t_b,    100.0 * t_b    / tot);
    Rprintf("  total    %.3f\n", tot);
  }

  return;
}

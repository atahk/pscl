/****************************************************************
 ** IDEAL_v4: experimental backend for ideal().
 **
 ** v4 = v3 (fused Cholesky + missing-data subtract) with every
 ** matrix allocated in contiguous storage (dmatrix_contig /
 ** imatrix_contig). The public (double **, int **) API is unchanged,
 ** so every inner helper (crossprod, crosscheck*, updatey,
 ** updatex_v3, updateb_v3, bayesreg_chol, etc.) works unmodified.
 **
 ** Rationale: the v1/v2/v3 dmatrix allocates each row as its own
 ** heap block, which defeats the hardware prefetcher and prevents
 ** the compiler from autovectorising the n x m inner loops in
 ** updatey and the cross-product accumulators. Contiguous storage
 ** restores sequential row-major access, which is the hot pattern
 ** throughout the sampler.
 ****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include <R_ext/Utils.h>
#include "util.h"
#include "ideal.h"

extern double **bpb, *xprior, **xpriormat, *xbar, *bpw, **w;
extern double **xpx, **bpriormat, *bprior, *bbar, *xpy;
extern double **v2_Lx, *v2_dx, *v2_zx, *v2_rx;
extern double **v2_Lb, *v2_db, *v2_zb, *v2_rb;

static double **v4_bpb_total;
static double **v4_xpx_total;
static int  *v4_n_miss_i;
static int **v4_miss_j_for_i;
static int  *v4_n_miss_j;
static int **v4_miss_i_for_j;

static void v4_build_missing_index(int **ok, int n, int m)
{
  int i, j, cnt, idx;

  v4_n_miss_i = ivector(n);
  v4_n_miss_j = ivector(m);
  for (i = 0; i < n; i++) v4_n_miss_i[i] = 0;
  for (j = 0; j < m; j++) v4_n_miss_j[j] = 0;

  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      if (!ok[i][j]) {
	v4_n_miss_i[i]++;
	v4_n_miss_j[j]++;
      }
    }
  }

  v4_miss_j_for_i = (int **) calloc(n, sizeof(int *));
  for (i = 0; i < n; i++) {
    cnt = v4_n_miss_i[i];
    if (cnt > 0) {
      v4_miss_j_for_i[i] = ivector(cnt);
      idx = 0;
      for (j = 0; j < m; j++) {
	if (!ok[i][j]) v4_miss_j_for_i[i][idx++] = j;
      }
    } else {
      v4_miss_j_for_i[i] = NULL;
    }
  }

  v4_miss_i_for_j = (int **) calloc(m, sizeof(int *));
  for (j = 0; j < m; j++) {
    cnt = v4_n_miss_j[j];
    if (cnt > 0) {
      v4_miss_i_for_j[j] = ivector(cnt);
      idx = 0;
      for (i = 0; i < n; i++) {
	if (!ok[i][j]) v4_miss_i_for_j[j][idx++] = i;
      }
    } else {
      v4_miss_i_for_j[j] = NULL;
    }
  }
}

static void v4_free_missing_index(int n, int m)
{
  int i, j;
  for (i = 0; i < n; i++) {
    if (v4_miss_j_for_i[i] != NULL) free(v4_miss_j_for_i[i]);
  }
  free(v4_miss_j_for_i);
  for (j = 0; j < m; j++) {
    if (v4_miss_i_for_j[j] != NULL) free(v4_miss_i_for_j[j]);
  }
  free(v4_miss_i_for_j);
  free(v4_n_miss_i);
  free(v4_n_miss_j);
}

void IDEAL_v4(int *n1, int *m1, int *d1, double *y1, int *maxiter1, int *thin1,
	      int *impute1, int *mda, double *xpriormeans1,
	      double *xpriorprec1, double *bpriormeans1, double *bpriorprec1,
	      double *xstart1, double *bstart1, double *xoutput, double *boutput,
	      int *burnin1, int *usefile, int *bsave, char **filename1, int *verbose1,
	      int *limitvoters1, int *usevoter)
{
  int e, xocursor, bocursor, xlength, blength, q, iter;
  int inloop, **ok, burnin, n, m, d, nm, maxiter, thin, impute, verbose;
  double **ystar, **x, **xreg, **y, **beta, **bp, **bpv, nm_doub, iterPerCent;
  double **xp, **xpv, *xtemp, *btemp;
  FILE *ofp;

  n = *n1;
  m = *m1;
  d = *d1;
  maxiter = *maxiter1;
  thin = *thin1;
  impute = *impute1;
  verbose = *verbose1;
  burnin = *burnin1;

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

  v4_build_missing_index(ok, n, m);

  /* scratch shared with v1 path (same extern names, contig storage) */
  bpb       = dmatrix_contig(d, d);
  bpw       = dvector(d);
  xbar      = dvector(d);
  xprior    = dvector(d);
  xpriormat = dmatrix_contig(d, d);
  w         = dmatrix_contig(n, m);

  xpy       = dvector(q);
  xpx       = dmatrix_contig(q, q);
  bbar      = dvector(q);
  bprior    = dvector(q);
  bpriormat = dmatrix_contig(q, q);

  v2_Lx = dmatrix_contig(d, d);
  v2_dx = dvector(d);
  v2_zx = dvector(d);
  v2_rx = dvector(d);

  v2_Lb = dmatrix_contig(q, q);
  v2_db = dvector(q);
  v2_zb = dvector(q);
  v2_rb = dvector(q);

  v4_bpb_total = dmatrix_contig(d, d);
  v4_xpx_total = dmatrix_contig(q, q);

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

      updatey(ystar, y, x, beta, n, m, d, iter);
      updatex_v3(ystar, ok, beta, x, xp, xpv, n, m, d, impute,
		 v4_n_miss_i, v4_miss_j_for_i, v4_bpb_total);
      makexreg(xreg, x, n, d, q);

      if (*limitvoters1 > 0)
	updateb_v3_usevoter(ystar, ok, beta, xreg, bp, bpv, n, m, d, impute, usevoter,
			    v4_n_miss_j, v4_miss_i_for_j, v4_xpx_total);
      else
	updateb_v3(ystar, ok, beta, xreg, bp, bpv, n, m, d, impute,
		   v4_n_miss_j, v4_miss_i_for_j, v4_xpx_total);

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

  free_dmatrix_contig(y);
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

  free_dmatrix_contig(bpb);
  free(bpw);
  free(xbar);
  free(xprior);
  free_dmatrix_contig(xpriormat);
  free_dmatrix_contig(w);
  free(xpy);
  free_dmatrix_contig(xpx);
  free(bbar);
  free(bprior);
  free_dmatrix_contig(bpriormat);

  free_dmatrix_contig(v2_Lx);
  free(v2_dx);
  free(v2_zx);
  free(v2_rx);
  free_dmatrix_contig(v2_Lb);
  free(v2_db);
  free(v2_zb);
  free(v2_rb);

  free_dmatrix_contig(v4_bpb_total);
  free_dmatrix_contig(v4_xpx_total);
  v4_free_missing_index(n, m);

  return;
}

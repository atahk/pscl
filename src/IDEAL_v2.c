/****************************************************************
 ** IDEAL_v2: experimental backend for ideal().
 **
 ** Identical control flow to IDEAL(), but the Gibbs updates for x
 ** and beta use a single Cholesky of the posterior precision (see
 ** bayesreg_chol.c) instead of Gauss-Jordan inversion of the
 ** covariance followed by a separate Cholesky for sampling.
 **
 ** The classic IDEAL path in IDEAL.c is left byte-identical so the
 ** two backends can be benchmarked side-by-side from R via the
 ** .backend argument to ideal().
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

/* shared globals defined in IDEAL.c (re-used scratch) */
extern double **bpb, *xprior, **xpriormat, *xbar, *bpw, **w;
extern double **xpx, **bpriormat, *bprior, *bbar, *xpy;

/* v2-only Cholesky workspaces */
double **v2_Lx, *v2_dx, *v2_zx, *v2_rx;
double **v2_Lb, *v2_db, *v2_zb, *v2_rb;

void IDEAL_v2(int *n1, int *m1, int *d1, double *y1, int *maxiter1, int *thin1,
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

  y = dmatrix(n, m);
  ystar = dmatrix(n, m);

  beta = dmatrix(m, q);
  bp  = dmatrix(m, q);
  bpv = dmatrix(m, q);

  x = dmatrix(n, d);
  xreg = dmatrix(n, q);
  xp  = dmatrix(n, d);
  xpv = dmatrix(n, d);

  ok = imatrix(n, m);

  if (*usefile == 1) {
    ofp = fopen(R_ExpandFileName(*filename1), "a");
    if (ofp == NULL) {
      calcerror("Can't open outfile file!\n");
    }
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

  /* scratch shared with v1 path */
  bpb = dmatrix(d, d);
  bpw = dvector(d);
  xbar = dvector(d);
  xprior = dvector(d);
  xpriormat = dmatrix(d, d);
  w = dmatrix(n, m);

  xpy = dvector(q);
  xpx = dmatrix(q, q);
  bbar = dvector(q);
  bprior = dvector(q);
  bpriormat = dmatrix(q, q);

  /* v2 Cholesky workspaces */
  v2_Lx = dmatrix(d, d);
  v2_dx = dvector(d);
  v2_zx = dvector(d);
  v2_rx = dvector(d);

  v2_Lb = dmatrix(q, q);
  v2_db = dvector(q);
  v2_zb = dvector(q);
  v2_rb = dvector(q);

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
      updatex_v2(ystar, ok, beta, x, xp, xpv, n, m, d, impute);
      makexreg(xreg, x, n, d, q);

      if (*limitvoters1 > 0)
	updateb_v2_usevoter(ystar, ok, beta, xreg, bp, bpv, n, m, d, impute, usevoter);
      else
	updateb_v2(ystar, ok, beta, xreg, bp, bpv, n, m, d, impute);

      R_CheckUserInterrupt();
    }

    if (iter >= burnin) {
      if (*usefile == 1) {
	dmatTOdvec(xtemp, x, n, d);
	fprintf(ofp, "%d", iter);
	for (e = 0; e < xlength; e++) {
	  fprintf(ofp, ",%f", xtemp[e]);
	}
	if (*bsave != 1) {
	  fprintf(ofp, "\n");
	}
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
	  for (e = 0; e < blength; e++) {
	    fprintf(ofp, ",%f", btemp[e]);
	  }
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

  if (*usefile == 1) {
    fclose(ofp);
  }

  free_dmatrix(y, n);
  free_dmatrix(ystar, n);
  free_dmatrix(beta, m);
  free_dmatrix(bp, m);
  free_dmatrix(bpv, m);
  free_dmatrix(x, n);
  free_dmatrix(xreg, n);
  free_dmatrix(xp, n);
  free_dmatrix(xpv, n);
  free_imatrix(ok, n);

  free(xtemp);
  free(btemp);

  free_dmatrix(bpb, d);
  free(bpw);
  free(xbar);
  free(xprior);
  free_dmatrix(xpriormat, d);
  free_dmatrix(w, n);
  free(xpy);
  free_dmatrix(xpx, q);
  free(bbar);
  free(bprior);
  free_dmatrix(bpriormat, q);

  free_dmatrix(v2_Lx, d);
  free(v2_dx);
  free(v2_zx);
  free(v2_rx);

  free_dmatrix(v2_Lb, q);
  free(v2_db);
  free(v2_zb);
  free(v2_rb);

  return;
}

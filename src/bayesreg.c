/***********************************************************************
** given data xpx (X'X, a p by p matrix), xpy (X'y, a p by 1 vector),
** priors bp (a p by 1 vector), and bpv (a p by p precision matrix)
** return the posterior mean and variance for the coefficients in the
** regression of y on X
**
** n.b., this is specific to probit, where var(e)=sigma^2=1
** by assumption (for identification)
**
** (C) Simon Jackman, Stanford University, 2001
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "util.h"
#include "ideal.h"

void bayesreg(double **xpx, double *xpy, 
	      double *bp, double **priormat,
	      double *bpost, double **vpost,
	      int p)
{
  int j,k;
  double *bpb;

  bpb = dvector(p);

  for(j=0;j<p;j++){
    bpost[j]=0.0;
    for(k=0;k<p;k++){
      vpost[j][k] = xpx[j][k] + priormat[j][k];  /* sum precisions */
    }
  }

  for(j=0;j<p;j++){
    bpb[j]=0.0;
    for(k=0;k<p;k++){
      bpb[j] += priormat[j][k]*bp[k]; /* prior, weighted by precision */
    }
    bpost[j] = xpy[j] + bpb[j];       /* add precision-weighted prior */
  }

  gaussj(vpost,p,bpost,1);       /* vpost inverted, bpost is posterior mean */


  free(bpb);


  return;
}

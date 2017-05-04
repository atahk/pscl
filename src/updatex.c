/****************************************************************
 ** update ideal points, conditional on ystar and beta
 **
 ** simon jackman, dept of political science, stanford university
 ** sep 2001
 ****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>
#include "util.h"
#include "ideal.h"

void updatex(double **ystar, int **ok, double **beta, 
	     double **x, double **xp, double **xpv,
	     int n, int m, int d, 
	     int impute)
{
  int i, j, k, l;
  extern double *xxprod, **xxchol, *xz, *xxp, **xxa;
  extern double **bpb, *xprior, **xpriormat, *xbar, **xvpost, *bpw, **w;

  /*
  Rprintf("xp 1: %d\n",xp[1][1]);
  Rprintf("xp 4: %d\n",xp[4][1]);
  Rprintf("xp 5: %d\n",xp[5][1]);
  Rprintf("xp 6: %d\n",xp[6][1]);
  Rprintf("xp 9: %d\n",xp[9][1]);
  printmat(xp,n,d);
  */
  /* form dependent variable */
  //w = dmatrix(n,m);
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      w[i][j] = ystar[i][j] + beta[j][d]; /*d+1*, add negative intercept to LHS */
    }
  }
  
  /* INITIALIZE */
  /* bpb = dmatrix(d,d);
     bpw  = dvector(d);
     xbar = dvector(d);
     xvpost = dmatrix(d,d);
     xprior = dvector(d);
     xpriormat = dmatrix(d,d); */

  //Rprintf("Point 2\n");

  if(impute==0){                     /* stingy filtering of missing data */
    for(i=0;i<n;i++){               /* loop over legislators */
      for(k=0;k<d;k++){             /* zero out */
	bpw[k]=0.0; 
	xbar[k]=0.0; 
	xprior[k]=0.0;
	for(l=0;l<d;l++){
	  bpb[k][l]=0.0;
	  xpriormat[k][l] = 0.0;
	  bpb[k][l]=0.0;
	  xvpost[k][l]=0.0;
	}
      }
      for (k=0;k<d;k++){            /* re-initializations */
	xprior[k] = xp[i][k];        /* copy prior mean */
	xpriormat[k][k] = xpv[i][k];  /* copy (diagonal) prior precision */
      }
      //Rprintf("\nupdatex: calling crosscheckx\n");
      crosscheckx(beta,w,ok,m,d,i,bpb,bpw); /* xprods, filter missings */
      //Rprintf("\nupdatex: calling bayesreg\n");
      bayesreg(bpb,bpw,xprior,xpriormat,xbar,xvpost,d);
      // Rprintf("\nupdatex: calling rmvnorm\n");
      // rmvnorm(x[i],xbar,xvpost,d);    /* do the sampling */
      rmvnorm(x[i],xbar,xvpost,d, xxprod, xxchol, xz, xxp, xxa); /* do the sampling */  
    }
  }
  
  if(impute==1){                     /* allow propagation of missing data */
    crossprod(beta,m,d,bpb);         /* gete bpb once and store */  
    for(i=0;i<n;i++){               /* loop over legislators */
      for(k=0;k<d;k++){             /* zero out */
	bpw[k]=0.0; 
	xbar[k]=0.0; 
	xprior[k]=0.0;
	for(l=0;l<d;l++){
	  bpb[k][l]=0.0;
	  xpriormat[k][l] = 0.0;
	  bpb[k][l]=0.0;
	  xvpost[k][l]=0.0;
	}
      }
      for (k=0;k<d;k++){            /* re-initializations */
	xprior[k] = xp[i][k];        /* copy prior mean */
	xpriormat[k][k] = xpv[i][k];  /* copy (diagonal) prior precision */
      }
      //Rprintf("\nupdatex: calling crosscheckx\n");
      crossxyi(beta,w,m,d,i,bpw);    /* bpw */
      //Rprintf("\nupdatex: calling bayesreg\n");
      bayesreg(bpb,bpw,xprior,xpriormat,xbar,xvpost,d);
      // Rprintf("\nupdatex: calling rmvnorm\n");
      rmvnorm(x[i],xbar,xvpost,d, xxprod, xxchol, xz, xxp, xxa);    /* do the sampling */
    }
  }


  /* free_dmatrix(bpb,d);
     free_dmatrix(w,n);
     free(bpw);
     free(xbar);
     free(xprior);
     free_dmatrix(xvpost,d);
     free_dmatrix(xpriormat,d); */

  return;
}

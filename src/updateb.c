/****************************************************************
 ** update proposal parameters, conditional on x and ystar
 ** 
 ** simon jackman, dept of political science, stanford university
 ** sep 2001
 ****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include "util.h"
#include "ideal.h"

void updateb(double **ystar, int **ok, double **beta, double **xreg,
	     double **bp, double **bpv,
	     int n, int m, int d, 
	     int impute)
{
  int j,k,q;
  extern double *bxprod, **bchol, *bz, *bbp, **bba,
    **xpx, **bvpost, **bpriormat, *bprior, *bbar, *xpy;

  q = d + 1;

  /*   xpy = dvector(q); */
  /*   xpx = dmatrix(q,q); */
  /*   bbar = dvector(q); */
  /*   bprior = dvector(q); */
  /*   bvpost = dmatrix(q,q); */
  /*   bpriormat = dmatrix(q,q); */
  
  for (j=0;j<q;j++){               /* initialize */
    xpy[j] = 0.0;
    for (k=0;k<q;k++){
      xpx[j][k] = 0.0;
      bvpost[j][k] = 0.0;
      bpriormat[j][k] = 0.0;
    }
  }

  if(impute==0){                      /* stingy filtering of missing data */
    for (j=0;j<m;j++){               /* loop over roll calls */
      for (k=0;k<q;k++){             /* initializations */
	bpriormat[k][k] = bpv[j][k];   /* diagonal prior precision */
	bprior[k]=bp[j][k];           /* copy prior mean */  
      }                               /* end initializations */
      //Rprintf("\nupdateb: calling crosscheck\n");
      crosscheck(xreg,ystar,ok,n,q,j,xpx,xpy); /* screen out missing data */
      //Rprintf("\nupdateb: calling bayesreg\n");
      bayesreg(xpx,xpy,bprior,bpriormat,bbar,bvpost,q);
      // Rprintf("\nupdateb: calling rmvnorm...");
      rmvnorm(beta[j],bbar,bvpost,q, bxprod, bchol, bz, bbp, bba);    
      // Rprintf("done\n");
    }
  }

  if(impute==1){                      /* allow propagation of missing data */
    crossprod(xreg,n,q,xpx);          /* get xpx once and store */
    for(j=0;j<m;j++){                /* loop over roll calls */
      for (k=0;k<q;k++){             /* initializations */
	bpriormat[k][k] = bpv[j][k];   /* diagonal prior precision */
	bprior[k]=bp[j][k];           /* copy prior mean */
      }                               /* end initializations */
      //Rprintf("\nupdateb: calling crossxyj\n");
      crossxyj(xreg,ystar,n,q,j,xpy); /* xpy */
      //Rprintf("\nupdateb: calling bayesreg\n");
      bayesreg(xpx,xpy,bprior,bpriormat,bbar,bvpost,q);

      // Rprintf("\nupdateb: calling rmvnorm...");
      rmvnorm(beta[j],bbar,bvpost,q, bxprod, bchol, bz, bbp, bba);   
      // Rprintf("done\n");
    }
  }


/*   free_dmatrix(xpx,q); */
/*   free_dmatrix(bvpost,q); */
/*   free_dmatrix(bpriormat,q); */
/*   free(xpy); */
/*   free(bbar); */
/*   free(bprior); */

  return;
}


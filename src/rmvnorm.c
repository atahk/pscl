/**************************************************
 * multivariate normal sampling
 *
 *
 * simon jackman, dept of political science, stanford university
 ***************************************************************/

#include <stdio.h>
#include <math.h>
#include <Rmath.h>
#include <R_ext/Print.h>
#include "util.h"
#include "ideal.h"

void rmvnorm(double *theta, double *mu, double **sigma, int k,
	     double *xprod, double **chol, double *z,
	     double *p, double **a)
{
  int i,j;
/*   double *xprod, **chol; */
/*   double *z; */
/*   long idum=(-13); */

/*   z = dvector(k); */
/*   xprod = dvector(k); */
/*   chol = dmatrix(k,k); */
  //Rprintf("Ready for decomposition.\n");
  xchol(sigma,chol,k,p,a);
  
  for(i=0;i<k;i++){
    xprod[i] = 0.0;                  /* initialize xprod */
    z[i] = norm_rand();              /* sample from univariate normals */
  }
  
  for(i=0;i<k;i++){
    for(j=0;j<k;j++){
      xprod[i] += chol[i][j]*z[j];   /* multiply by Choleski of sigma */
    }
  }

  for(i=0;i<k;i++){
    theta[i] = mu[i] + xprod[i];     /* add in means */
  }
  
/*   free(z); */
/*   free(xprod); */
/*   free_dmatrix(chol,k); */

  return;
}

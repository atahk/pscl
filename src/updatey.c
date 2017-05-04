/*****************************************************************
 **
 ** update ystar
 **
 ** simon jackman, dept of political science, stanford university
 ** mar 2001
 *****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <Rmath.h>
#include "util.h"
#include "ideal.h"

/* putting stuff in y star */

void updatey(double **ystar, double **y, double **x, double **beta,
	     int n, int m, int d, int iter)
{
  int i,j,k;
  double *xrow, *brow, mu, sd;
  sd = 1.0;
  //  float z;

  for(i=0;i<n;i++){               /* loop over legislators */
    xrow = x[i];

    for(j=0;j<m;j++){             /* loop over proposals */
      brow = beta[j];
      
      mu = -1.0*brow[d];          /* minus intercept */ /*d+1*/
     
      for(k=0;k<d;k++){
	mu += brow[k]*xrow[k];
      }
      if (y[i][j]==9.0){         /* sample untruncated, missing responses */
	ystar[i][j] = rnorm(mu,1.0);
      }
      else{                      /* sample from truncated normals */
	ystar[i][j] = dtnorm(&mu,&sd,&y[i][j]);  /* try two */
      }
    }
  }
  return;
}

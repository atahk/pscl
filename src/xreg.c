/* form regressors for update of beta 
 *
 * simon jackman, dept of political science,
 * stanford university */

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "ideal.h"

void makexreg(double **xreg, double **x, int n, int d, int q)
{

  int i,j;
  
  for (i=0;i<n;i++){               /* initialize */
    xreg[i][d] = -1.0;             /* add negative intercept on end */
    for (j=0;j<d;j++)
      xreg[i][j] = x[i][j];        /* copy over xreg */
  }             

}

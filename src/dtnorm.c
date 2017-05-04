/*************************************************************************
 ** truncated Normal sampling
 **
 ** simon jackman, dept of political science, stanford university
 ** feb 2000
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include "util.h"

static double zero = 0.0;
static double arg = 0.5;
static double iarg = 0.5;
static double pupper = 0.5;
static double one = 1.0;


double dtnorm(double *mu, double *sd, double *y)
{
  double f, z;
  double norm=0;
  
  if (*y==0.0){
    z = *mu/(*sd);
    if(z<1.6){
      /* try rejection sampling */
      do {
	norm = rnorm(*mu, *sd);
      } while (norm >= 0.0);
    }
    else{
      /* otherwise use inverse-uniform method, z is always positive */
      /* work with natural logarithms to avoid underflows */
      f = -exp_rand();
      
      pupper = pnorm(z,zero,one,0,1);
      arg = f + pupper;
      iarg = qnorm(arg,zero,one,1,1);

      norm = *mu + (*sd)*iarg;
    }
  }
  else{  /* Y=1 */
    /* try rejection sampling */
    z = *mu/(*sd);
    if(z>-1.6){
	    do {
		    norm = rnorm(*mu, *sd);
	    } while (norm <= 0.0);
    }
    else{
      /* otherwise use inverse-uniform method, n.b., z is always neg */
      /* work with natural logarithms to avoid underflows */
      f = -exp_rand();

      pupper = pnorm(z,zero,one,1,1);
      arg = f + pupper;
      iarg = qnorm(arg,zero,one,0,1);

      norm = *mu + (*sd)*iarg;
    }
  }
  return(norm);
}

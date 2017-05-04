#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

void simpi(int *n, int *z)
{
  int i;
  double d;

  GetRNGstate();

  for(i=0;i<*n;i++){
    d = hypot(unif_rand(),unif_rand());
    if(d<1.0) 
      (*z)++;
  }

  PutRNGstate();
  return;
}

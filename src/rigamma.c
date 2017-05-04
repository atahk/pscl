#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double r_sd(double s, double df)
{
  double root, r, g;
  
  g = rchisq(df);
  r = s/g;
  root = sqrt(r);

  // Rprintf("r_sd: s,g, s/g = %14.4lf %14.4lf %14.4lf\n",s,g,r);

  return(root);
}

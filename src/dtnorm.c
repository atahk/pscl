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

double dtnorm_std(const double lower_bound)
{
  double y;
  if (lower_bound < 0.0)
    {
      do {
	y = norm_rand();
      } while (y <= lower_bound);
      return y;
    }
  else if (lower_bound < 0.75)
    {
      do {
	y = fabs(norm_rand());
      } while (y <= lower_bound);
      return y;
    }
  else
    {
      do {
	y = exp_rand();
      } while (exp_rand() * lower_bound * lower_bound <= 0.5 * y * y);
      return y / lower_bound + lower_bound;
    }
}

double dtnorm(const double mu, const double sd, const double y)
{
  if (y <= 0.0)
    return mu - sd * dtnorm_std(mu / sd);
  else
    return mu + sd * dtnorm_std(- mu / sd);
}

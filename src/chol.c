#include <stdio.h>
#include <math.h>
#include "util.h"

/* from NR, p97:
 * Given a positive-definite symmetric matrix a[1..n][1..n], this
 * routine constructs its Cholesky decomposition, A = LL'.  On input,
 * only the upper triangle of A need be given; it is not modified.
 * The Cholesky factor L is returned in the lower triangle of A,
 * except for its diagonal elements which are returned in p[1..n]
 */
void choldc(double **a, int n, double p[])
{
  void calcerror(char error_text[]);
  int i,j,k;
  double sum;

  for (i=0;i<n;i++){
    for (j=0;j<n;j++){
      for (sum=a[i][j],k=i-1;k>=0;k--)
	sum -= a[i][k]*a[j][k];
      if (i == j) {
	if (sum <= 0.0)
	  calcerror("Cholesky decomposition error: Matrix is not positive definite\n");
	p[i]=sqrt(sum);
      }
      else a[j][i]=sum/p[i];
    }
  }
}



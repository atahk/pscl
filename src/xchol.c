#include <stdio.h>
#include "util.h"
#include "chol.h"

void xchol(double **aorig, double **chol, int n, double *p, double **a)
{
  int i,j;
  //double **a, *p;

  // p = dvector(n);
  // a = dmatrix(n,n);

  /* printf("xchol: n = %d\n",n); */
  /* printf("xchol: starting reassignments\n"); */
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      /* printf("%12.7lf",aorig[i][j]); */
      a[i][j] = aorig[i][j];
      chol[i][j] = 0.0;
    }
    /* fprintf(stdout,"\n"); */
  }
  /* printf("xchol: calling cholesky routine\n"); */
  choldc(a,n,p);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      chol[i][j]=((i > j) ? a[i][j] : ( i == j ? p[i] : 0.0));
      if (i > j) 
	chol[i][j]=a[i][j];
      else 
	chol[i][j]=(i ==j ? p[i] : 0.0);
    }
  }
  // free(p);
  // free_dmatrix(a,n);
}

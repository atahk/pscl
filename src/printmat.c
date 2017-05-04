#include <stdio.h>
#include "util.h"
#include <R_ext/Print.h>

void printmat(double **mat, int nr, int nc)
{
  int i,j;

  for(i=0;i<nr;i++){
    for(j=0;j<nc;j++){
      Rprintf("mat[%d][%d]=%2.3lf ",
	      i,j,mat[i][j]);
    }
    Rprintf("\n");
  }
  return;
}

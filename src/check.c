/****************************************************************
 ** report some simple summary statistics for the roll call data
 **
 ** simon jackman, dept of political science, stanford university
 ** feb 2000
 ****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include "util.h"

double check(double **data, int **ok, int n, int m)
{
  int i, j;
  double *yeas, *nummiss, *inummiss, *x, nok;
  
  yeas = dvector(m);
  x = dvector(n);
  nummiss = dvector(m);
  inummiss = dvector(n);

  for(i=0;i<n;i++){
    x[i]=0.0;
    inummiss[i]=0.0;
  }
  for(j=0;j<m;j++){
    yeas[j]=0.0;
    nummiss[j]=0.0;
  }
  
  nok = 0.0;

  /* end preliminaries */

  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if (data[i][j]==1.0){
	x[i]++;                     /* yeas for legislator i */
	yeas[j]++;                  /* yeas for bill j */
      }
      if (data[i][j]==9.0){
	inummiss[i]++;              /* missing for legislator i */
	nummiss[j]++;               /* missing for bill j */
      }
      else{                         /* if not missing */
	nok++;                      /* total number ok */
	ok[i][j]=1;                 /* indicator, i,j-th decision ok? */
      }
    }
  }

  /*
  for(j=0;j<m;j++){
    Rprintf("Vote %4d: Yea: %5.2lf, Nay: %5.2lf NA: %5.2lf\n",
	    j,yeas[j],n-yeas[j]-nummiss[j],nummiss[j]);
  }

  for(i=0;i<n;i++){
    Rprintf("Legislator %4d: Percent Yea: %5.2lf%, Missing Data %5.2lf%\n",
	    i,
	    x[i]/m*100.0,
	    inummiss[i]/m*100.0);
  }

  */

  free(yeas);
  free(x);
  free(nummiss);
  free(inummiss);

  return(nok);
}

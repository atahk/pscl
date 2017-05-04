#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "util.h"
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

/*****************************************************************************
 ** Linear equation solution by Gauss-Jordan elimination.
 ** Given the system Ax=b, this routine returns x in b and
 ** A^{-1} in A.  Source: _NR_ 2.1, p39. 
 ****************************************************************************/

void gaussj(double **a, int n, double *b, int m)
{
  int *indxc, *indxr, *ipiv;
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv,temp;
  
  indxc=ivector(n);
  indxr=ivector(n);
  ipiv=ivector(n);

  for (j=0;j<n;j++) ipiv[j]=0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++) 
      if (ipiv[j] != 1)
        for (k=0;k<n;k++) {
          if (ipiv[k] == 0) {
            if (fabs(a[j][k]) >= big) {
              big=fabs(a[j][k]);
              irow=j;
              icol=k;
            }
          } else if (ipiv[k] > 1) calcerror("Error in Gauss-Jordan elimination: Singular Matrix\n");
        }
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
      SWAP(b[irow],b[icol])
    }
    indxr[i]=irow;
    indxc[i]=icol;

    if (a[icol][icol] == 0.0) calcerror("Error in Gauss-Jordan elimination: Singular Matrix\n");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<n;l++) a[icol][l] *= pivinv;
    b[icol] *= pivinv;

    for (ll=0;ll<n;ll++)
      if (ll != icol) {
        dum=a[ll][icol];
        a[ll][icol]=0.0;
        for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
        b[ll] -= b[icol]*dum;
      }
  }
  for (l=n-1; l>=0;l--) {
    if (indxr[l] != indxc[l])
      for (k=0;k<n;k++)
        SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }

  free(ipiv);
  free(indxr);
  free(indxc);

  return;
}      

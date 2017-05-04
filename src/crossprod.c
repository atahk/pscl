/*****************************************************************
 **
 ** cross-products routine
 **
 ** given a n-by-k matrix x, write x'x (k-by-k matrix)
 **
 ** simon jackman, dept of political science, stanford university
 ** feb 2000
 *****************************************************************/

void crossprod(double **x, int n, int p, double **xpx)
{
  int i, j, k;
  double *xrow;
  
  for (j=0;j<p;j++) {           /* initialization */
    for (k=0;k<p;k++) {
      xpx[j][k]=0.0;
    }
  }

  for (i=0;i<n;i++){
    xrow = x[i];
    for(j=0;j<p;j++){
      for(k=0;k<p;k++){
        xpx[j][k] += xrow[j]*xrow[k];
      }
    }
  }

  return;
}


void crossprodslow(double **x, int n, int p, double **xpx)
{
  int i, j, k;

  for (j=0;j<p;j++) {           /* initialization */
    for (k=0;k<p;k++) {
      xpx[j][k]=0.0;
    }
  }

  for (i=0;i<n;i++){
    for(j=0;j<p;j++){
      for(k=0;k<p;k++){
        xpx[j][k] += x[i][j]*x[i][k];
      }
    }
  }

  return;
}




void crossxy(double **x, double *y, int n, int k, double *xpy)
{
  int i,j;
  double *xrow, yi;

  /* xrow = dvector(k); */

  for (j=0;j<k;j++){
    xpy[j]=0.0;                /* initialization */
  }

  for (i=0;i<n;i++){          /* loop over observations */
    xrow = x[i];
    yi = y[i];
    for (j=0;j<k;j++){
      xpy[j] += xrow[j]*y[i];
    }
  }

  return;
}

void crossxyj(double **x, double **y, int n, int k, int p, double *xpy)
{
  int i,j;
  double *xrow, yip;

  for (j=0;j<k;j++){    /* initialize */
    xpy[j]=0.0;
  }

  for (i=0;i<n;i++){    /* loop over observations */
    xrow = x[i];
    yip = y[i][p];
    for(j=0;j<k;j++){
      xpy[j] += xrow[j]*yip;
    }
  }
  return;
}

/* use this for forming beta'w in updatex */
void crossxyi(double **beta, double **w, int m, int d, int p, double *bpw)
{
  int j,k;
  double *betarow, wpj;

  for (k=0;k<d;k++){   /* initialize */
    bpw[k]=0.0;
  }

  for (j=0;j<m;j++){   /* loop over roll calls */
    betarow = beta[j];
    wpj = w[p][j];
    for (k=0;k<d;k++){ 
      bpw[k] += betarow[k]*wpj;
    }
  }
  return;
}

void crossxyd(double **x, double *y, int n, int k, double *xpy)
{
  int i,j;
  double *xrow, yi;

  for (j=0;j<k;j++){   /* initialize */
    xpy[j]=0.0;
  }

  for (i=0;i<n;i++){   /* loop over obs */
    xrow = x[i];
    yi = y[i];
    for (j=0;j<k;j++){ /* loop over cols */
      xpy[j] += xrow[j]*yi;
    }
  }
  return;
}


/* used in updating beta_j, without checks for missing data */
void crossall(double **x, double **ystar, int n, int q, int j, 
	      double **xpx, double *xpy)
{
  int i,k,l;
  double *xrow, ystarij;

  for(i=0;i<n;i++){                   /* loop over obs */
    xrow = x[i]; 
    ystarij=ystar[i][j];
    for(k=0;k<q;k++){                 /* loop over cols */
      xpy[k] += xrow[k]*ystarij;     
      for(l=0;l<q;l++){
	xpx[k][l] += xrow[k]*xrow[l]; 
      }
    }
  }
  return;
}

/* used in updating beta_j, with checks for missing y_{ij} */
void crosscheck(double **x, double **ystar, int **ok,
		int n, int q, int j, 
		double **xpx, double *xpy)
{
  int i,k,l, okij;
  double *xrow, ystarij, xk, xl;
  for(k=0;k<q;k++){        /* initializations */
    xpy[k]=0.0;
    for(l=0;l<q;l++){
      xpx[k][l]=0.0;
    }
  }

  for(i=0;i<n;i++){        /* loop over legislators */
    okij = ok[i][j];
    if(okij){
      xrow = x[i];
      ystarij = ystar[i][j];
      for(k=0;k<q;k++){                    /* loop over cols */
	xk = xrow[k];
	xpy[k] += xk*ystarij;            /* X'y contribution */
	for(l=0;l<q;l++){
	  xl = xrow[l];
	  xpx[k][l] += xk*xl;           /* X'X contribution */
	}
      }
    }
  }
  return;
}

/* used in updating x_i, with checks for missing y_{ij} */
void crosscheckx(double **beta, double **w, int **ok,
		 int m, int d, int i, 
		 double **bpb, double *bpw)
{
  int j,k,l, okij=1;
  double *betarow, wij, bk;

  for(k=0;k<d;k++){
    bpw[k]=0.0;
    for(l=0;l<d;l++){
      bpb[k][l]=0.0;
    }
  }

  for(j=0;j<m;j++){        /* loop over roll calls/items */
    wij = w[i][j];
    okij = ok[i][j];
    if(okij){
      betarow = beta[j];
      for(k=0;k<d;k++){
	bk = betarow[k];
	bpw[k] += bk*wij;  /* bpw contribution */
	for(l=0;l<d;l++){
	  bpb[k][l] += bk*betarow[l]; /* bpb contribution */
	}
      }
    }
  }

  return;
}

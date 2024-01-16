#include <stdio.h>
#include <stdlib.h>
#include <R_ext/Print.h>
#include <R_ext/Error.h>
#include "util.h"

void memallocerror(void)
{
  error("Memory allocation error.\n");
}

void calcerror(char error_text[])
{
  error("%s",error_text);
}

int *ivector(long n)
/* allocate an int vector with subscript range v[nl..nh] */
{

	int *v = malloc(n*sizeof(int));
	if (!v) memallocerror();
	return v;
}
double *dvector(long n)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;
  v = (double *) calloc(n, sizeof(double));
  if (!v) memallocerror();
  return v;
}

double **dmatrix(long nr, long nc)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i;
  double **m;
  m = (double **) calloc(nr,sizeof(double*));
  if (!m) memallocerror();
  for(i=0; i<nr; i++) {
    m[i]= (double *) calloc(nc,sizeof(double));
  }
  return m;
}

int **imatrix(long nr, long nc)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i;
  int **m;

  m = (int **) calloc(nr,sizeof(int*));
  if (!m) memallocerror();
  for(i=0; i<nr; i++) {
    m[i] = (int *) calloc(nc,sizeof(int));
  }
  return m;
}

void free_dmatrix(double **m, long nr)
/* free a double matrix allocated by dmatrix() */
{
  long i;

  for(i=0; i<nr; i++) {
    free(m[i]);
  }
  free(m);
}

void free_imatrix(int **m, long nr)
/* free an int matrix allocated by imatrix() */
{
  long i;

  for(i=0; i<nr; i++) {
    free(m[i]);
  }
  free(m);
}

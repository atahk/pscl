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

/* Contiguous-storage variants: one big data block plus a row-pointer
   array pointing into it. Element access m[i][j] is unchanged, so any
   existing helper taking (double **, int **) works without modification. */

double **dmatrix_contig(long nr, long nc)
{
  double **m;
  double *data;
  long i;

  if (nr <= 0 || nc <= 0) return NULL;

  m = (double **) calloc(nr, sizeof(double *));
  if (!m) memallocerror();
  data = (double *) calloc((size_t) nr * (size_t) nc, sizeof(double));
  if (!data) memallocerror();

  for (i = 0; i < nr; i++) m[i] = data + i * nc;
  return m;
}

int **imatrix_contig(long nr, long nc)
{
  int **m;
  int *data;
  long i;

  if (nr <= 0 || nc <= 0) return NULL;

  m = (int **) calloc(nr, sizeof(int *));
  if (!m) memallocerror();
  data = (int *) calloc((size_t) nr * (size_t) nc, sizeof(int));
  if (!data) memallocerror();

  for (i = 0; i < nr; i++) m[i] = data + i * nc;
  return m;
}

void free_dmatrix_contig(double **m)
{
  if (m == NULL) return;
  free(m[0]);    /* the contiguous data block */
  free(m);       /* the row-pointer array */
}

void free_imatrix_contig(int **m)
{
  if (m == NULL) return;
  free(m[0]);
  free(m);
}

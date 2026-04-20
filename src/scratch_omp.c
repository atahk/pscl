/***************************************************************
 ** scratch_omp.c -- allocator/deallocator for per-thread scratch.
 ***************************************************************/

#include <stdlib.h>
#include "util.h"
#include "scratch_omp.h"

scratch_t *scratch_alloc(int n_threads, int d)
{
  int q = d + 1;
  int t;
  scratch_t *scr = (scratch_t *) calloc((size_t)n_threads, sizeof(scratch_t));
  for (t = 0; t < n_threads; t++) {
    scratch_t *s = &scr[t];
    s->bpb       = dmatrix_contig(d, d);
    s->bpw       = dvector(d);
    s->xbar      = dvector(d);
    s->xprior    = dvector(d);
    s->xpriormat = dmatrix_contig(d, d);
    s->Lx        = dmatrix_contig(d, d);
    s->dx        = dvector(d);
    s->zx        = dvector(d);
    s->rx        = dvector(d);

    s->xpx       = dmatrix_contig(q, q);
    s->xpy       = dvector(q);
    s->bprior    = dvector(q);
    s->bbar      = dvector(q);
    s->bpriormat = dmatrix_contig(q, q);
    s->Lb        = dmatrix_contig(q, q);
    s->db        = dvector(q);
    s->zb        = dvector(q);
    s->rb        = dvector(q);
  }
  return scr;
}

void scratch_free(scratch_t *scr, int n_threads)
{
  int t;
  for (t = 0; t < n_threads; t++) {
    scratch_t *s = &scr[t];
    free_dmatrix_contig(s->bpb);
    free(s->bpw);
    free(s->xbar);
    free(s->xprior);
    free_dmatrix_contig(s->xpriormat);
    free_dmatrix_contig(s->Lx);
    free(s->dx);
    free(s->zx);
    free(s->rx);

    free_dmatrix_contig(s->xpx);
    free(s->xpy);
    free(s->bprior);
    free(s->bbar);
    free_dmatrix_contig(s->bpriormat);
    free_dmatrix_contig(s->Lb);
    free(s->db);
    free(s->zb);
    free(s->rb);
  }
  free(scr);
}

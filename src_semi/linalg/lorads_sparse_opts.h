#ifndef LORADS_SPARSE_OPT_H
#define LORADS_SPARSE_OPT_H

#include <stdio.h>
#include "lorads.h"

extern lorads_int csp_nnz_cols ( lorads_int n, lorads_int *Ap );
extern void spMul(lorads_int n, lorads_int nnz, double a,
                  lorads_int *Ai, lorads_int *Aj, double *Ax, double *x,
                  lorads_int k, double *y);
extern void tsp_decompress( lorads_int n, lorads_int nnz,
                            lorads_int *Ci, double *Cx,
                            lorads_int *Ai, lorads_int *Aj, double *Ax );
#endif

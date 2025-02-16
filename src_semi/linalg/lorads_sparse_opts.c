#include <math.h>
#include <assert.h>
#include "lorads_sparse_opts.h"
#include "lorads_vec_opts.h"

#ifdef MEMDEBUG
#include "memwatch.h"
#endif


extern lorads_int csp_nnz_cols ( lorads_int n, lorads_int *Ap ) {
    
    lorads_int nzcols = 0;
    
    for ( lorads_int i = 0; i < n; ++i ) {
        nzcols += ( Ap[i + 1] - Ap[i] > 0 );
    }
    return nzcols;
}

extern void spMul(lorads_int n, lorads_int nnz, double a,
                  lorads_int *Ai, lorads_int *Aj, double *Ax, double *x,
                  lorads_int k, double *y){
    if (a == 0.0){
        return;
    }
    lorads_int row = 0;
    lorads_int col = 0;
    for (int i = 0; i < nnz; ++i){
        row = Ai[i];
        col = Aj[i];
        axpy(&k, &Ax[i], &x[col], &n, &y[row], &n);
    }
}

/* Decompress a column to triplet matrix */
extern void tsp_decompress( lorads_int n, lorads_int nnz, lorads_int *Ci, double *Cx, lorads_int *Ai, lorads_int *Aj, double *Ax ) {

    lorads_int j = 0, idthresh = n;

    for ( lorads_int k = 0; k < nnz; ++k ) {
        while ( Ci[k] >= idthresh ) {
            j += 1;
            idthresh += n - j;
        }
        Ai[k] = Ci[k] - idthresh + n;
        Aj[k] = j;
        Ax[k] = Cx[k];
    }

    return;
}

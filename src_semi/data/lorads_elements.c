

#include "def_lorads_elements.h"
#include "lorads_utils.h"
#include "lorads_vec_opts.h"

extern void addDense(double *alpha, void *constrVal, double *vec)
{
    dense_vec *dense = (dense_vec *)constrVal;
    lorads_int incx = 1;
    axpy(&dense->nnz, alpha, dense->val, &incx, vec, &incx);
}

extern void addSparse(double *alpha, void *constrVal, double *vec)
{
    sparse_vec *sparse = (sparse_vec *)constrVal;
    int idx = 0;
    for (lorads_int i = 0; i < sparse->nnz; ++i)
    {
        idx = sparse->nnzIdx[i];
        vec[idx] += alpha[0] * sparse->val[i];
    }
}

extern void zeroDense(void *constrVal)
{
    dense_vec *dense = (dense_vec *)constrVal;
    LORADS_ZERO(dense->val, double, dense->nnz);
}

extern void zeroSparse(void *constrVal)
{
    sparse_vec *sparse = (sparse_vec *)constrVal;
    for (lorads_int i = 0; i < sparse->nnz; ++i){
        sparse->val[i] = 0;
    }
}

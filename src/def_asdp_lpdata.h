#ifndef DEF_ASDP_LPDATA_H
#define DEF_ASDP_LPDATA_H

#ifdef HEADERPATH
#include "interface/asdp.h"
#include "src/asdp_linsolver.h"
#include "src/asdp_rk_mat.h"
#else
#include "asdp.h"
#include "def_asdp_rk_mat.h"
#endif
/* Implementations of the LP coefficient matrix
In ASDP, we implement three data structure for LP coefficient matrix:

1. Zero matrix: all entries are zero
2. Sparse matrix: only non-zero entries are stored
3. Dense matrix: all entries are stored

*/

typedef enum{
    LP_COEFF_ZERO,
    LP_COEFF_SPARSE,
    LP_COEFF_DENSE
}lp_coeff_type;

typedef struct {
    int nRows;
    lp_coeff_type dataType;
    void *dataMat;
    
    asdp_retcode (*create)              (void **, int, int, int *, double *);
    void         (*destroy)             (void **);
    void         (*mul_inner_rk_double) (void *, double *, double *);
    void         (*weight_sum)          (void *, double *, double *);
    void         (*scaleData)           (void *, double );
}lp_coeff;

typedef struct {
    int nRows;
}lp_coeff_zero;

typedef struct {
    int nRows;
    int nnz;
    int *rowPtr;
    double *val;
}lp_coeff_sparse;

typedef struct {
    int nRows;
    double *val;
}lp_coeff_dense;


#endif // DEF_ASDP_LPDATA_H

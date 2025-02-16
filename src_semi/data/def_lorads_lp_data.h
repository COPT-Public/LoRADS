#ifndef DEF_LORADS_LP_DATA
#define DEF_LORADS_LP_DATA

#include "lorads.h"
typedef enum{
    LP_COEFF_ZERO,
    LP_COEFF_SPARSE,
    LP_COEFF_DENSE
}lp_coeff_type;

typedef struct {
    lorads_int nRows;
    lp_coeff_type dataType;
    void *dataMat;

    void  (*create)              (void **, lorads_int, lorads_int, lorads_int *, double *);
    void  (*destroy)             (void **);
    void  (*mul_inner_rk_double) (void *, double *, double *);
    void  (*weight_sum)          (void *, double *, double *);
//    void         (*scaleData)           (void *, double );
}lp_coeff;

typedef struct {
    lorads_int nRows;
}lp_coeff_zero;

typedef struct {
    lorads_int nRows;
    lorads_int nnz;
    lorads_int *rowPtr;
    double *val;
}lp_coeff_sparse;

typedef struct {
    lorads_int nRows;
    double *val;
}lp_coeff_dense;



#endif

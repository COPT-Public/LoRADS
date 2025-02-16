#ifndef DEF_LORADS_ELEMENTS
#define DEF_LORADS_ELEMENTS

#include "lorads.h"

typedef enum{
    LORADS_DENSE_VEC,
    LORADS_SPARSE_VEC
}vecType;

typedef struct{
    vecType type;
    void *data;
    void (*add) (double *, void*, double *);
    void (*zero) (void *);
}lorads_vec;

typedef struct{
    lorads_int nnz;
    lorads_int *nnzIdx;
    double *val;
}sparse_vec;

typedef struct{
    lorads_int nnz;
    double *val;
}dense_vec;

typedef struct{
    lorads_int rank;
    lorads_int nRows;
    double *matElem;
}lorads_sdp_dense;


typedef struct{
    lorads_int nCols;
    double *matElem;
}lorads_lp_dense;

#endif
#ifndef DEF_ASDP_RK_MAT_H
#define DEF_ASDP_RK_MAT_H

typedef struct {
    // structure for optimizer
    int rank; // rank of the matrix, i.e., column number of the matrix
    int nRows; // row number of the matrix
    double *matElem; // element of the matrix
} asdp_rk_mat_dense;

typedef struct {
    int nLPCols;
    double *matElem;
}asdp_rk_mat_lp;


typedef enum {
    ASDP_DENSE,
    ASDP_SPARSE
} constrType;

typedef struct {
    void *constrVal;
    constrType type;
    void (*add) (double *, void*, double *);
    void (*zero) (void *);
}constrValStruct;

typedef struct {
    int nnz;
    int *nnzIdx;
    double *val;
}constrValSparse;

typedef struct {
    int nnz;
    double *val;
}constrValDense;


#endif // DEF_ASDP_RK_MAT_H

#ifndef sparse_opts_h
#define sparse_opts_h

#include "asdp_sdpdata.h"
#ifdef __cplusplus
extern "C" {
#endif

extern void csp_Axpby( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y );
extern double csp_sum_abs( int n, int *Ap, int *Ai, double *Ax );
extern double csp_fro_norm( int n, int *Ap, int *Ai, double *Ax );
extern void csp_aApB( int n, int nnz, double a, int *Al, double *Ax, double *Bx );
extern int csp_nnz_cols ( int n, int *Ap );
extern void csp_dump( int n, int *Ap, int *Ai, double *Ax, double *v );

extern void tsp_decompress( int n, int nnz, int *Ci, double *Cx, int *Ai, int *Aj, double *Ax );
extern int tsp_r1_extract( int n, int nnz, int *Ai, int *Aj, double *Ax, double *sgn, double *a );
extern void tsp_scal( double a, int nnz, double *Ax );
extern double tsp_sum_abs( int nnz, int *Ai, int *Aj, double *Ax );
extern double tsp_fro_norm( int nnz, int *Ai, int *Aj, double *Ax );
extern void tsp_dump( int n, int nnz, int *Ai, int *Aj, double *Ax, double *v );
extern double tsp_quadform( int n, int nnz, int *Ai, int *Aj, double *Ax, double *v );
extern void tsp_Axpby( int n, int nnz, double a, int *Ai, int *Aj, double *Ax, double *x, double *y );
extern void tsp_ASemiLowby( int n, int nnz, double a, int *Ai, int *Aj, double *Ax, double *x, double *y );
extern void tsp_ASemiLowbyMul(int n, int nnz, double a, int *Ai, int *Aj, double *Ax, double *x, int k, double *y);
extern void ATspx(int n, int m, double *A, int nnz, int *idx, double *x, double *res);


struct  sparse_matrix;
typedef struct sparse_matrix *sparse_matrix_t;
/* status of the routines */
typedef enum
{
    SPARSE_STATUS_SUCCESS           = 0,    /* the operation was successful */
    SPARSE_STATUS_NOT_INITIALIZED   = 1,    /* empty handle or matrix arrays */
    SPARSE_STATUS_ALLOC_FAILED      = 2,    /* internal error: memory allocation failed */
    SPARSE_STATUS_INVALID_VALUE     = 3,    /* invalid input value */
    SPARSE_STATUS_EXECUTION_FAILED  = 4,    /* e.g. 0-diagonal element for triangular solver, etc. */
    SPARSE_STATUS_INTERNAL_ERROR    = 5,    /* internal error */
    SPARSE_STATUS_NOT_SUPPORTED     = 6     /* e.g. operation for double precision doesn't support other types */
} sparse_status_t;

typedef enum
{
    SPARSE_INDEX_BASE_ZERO  = 0,           /* C-style */
    SPARSE_INDEX_BASE_ONE   = 1            /* Fortran-style */
} sparse_index_base_t;
typedef enum
{
    SPARSE_DIAG_NON_UNIT    = 50,           /* triangular matrix with non-unit diagonal */
    SPARSE_DIAG_UNIT        = 51            /* triangular matrix with unit diagonal */
} sparse_diag_type_t;
typedef enum
{
    SPARSE_FILL_MODE_LOWER  = 40,           /* lower triangular part of the matrix is stored */
    SPARSE_FILL_MODE_UPPER  = 41,            /* upper triangular part of the matrix is stored */
    SPARSE_FILL_MODE_FULL   = 42            /* upper triangular part of the matrix is stored */
} sparse_fill_mode_t;
typedef enum
{
    SPARSE_MATRIX_TYPE_GENERAL            = 20,   /*    General case                    */
    SPARSE_MATRIX_TYPE_SYMMETRIC          = 21,   /*    Triangular part of              */
    SPARSE_MATRIX_TYPE_HERMITIAN          = 22,   /*    the matrix is to be processed   */
    SPARSE_MATRIX_TYPE_TRIANGULAR         = 23,
    SPARSE_MATRIX_TYPE_DIAGONAL           = 24,   /* diagonal matrix; only diagonal elements will be processed */
    SPARSE_MATRIX_TYPE_BLOCK_TRIANGULAR   = 25,
    SPARSE_MATRIX_TYPE_BLOCK_DIAGONAL     = 26    /* block-diagonal matrix; only diagonal blocks will be processed */
} sparse_matrix_type_t;

struct matrix_descr {
    sparse_matrix_type_t  type;       /* matrix type: general, diagonal or triangular / symmetric / hermitian */
    sparse_fill_mode_t    mode;       /* upper or lower triangular part of the matrix ( for triangular / symmetric / hermitian case) */
    sparse_diag_type_t    diag;       /* unit or non-unit diagonal ( for triangular / symmetric / hermitian case) */
};

/* sparse_status_t mkl_sparse_d_create_csr(       sparse_matrix_t     *A,
                                         const sparse_index_base_t indexing,
                                         const int             rows,
                                         const int             cols,
                                               int             *rows_start,
                                               int             *rows_end,
                                               int             *col_indx,
                                               double              *values );

sparse_status_t mkl_sparse_ee_init (int* pm);
sparse_status_t mkl_sparse_d_ev (char *which, int *pm, sparse_matrix_t A, struct matrix_descr descrA, int k0, int *k, double *E, double *X, double *res);
sparse_status_t mkl_sparse_destroy( sparse_matrix_t  A );

extern void sparse_Min_Eig(sdp_coeff *slackVar, double *minEig);*/
#ifdef __cplusplus
}
#endif

#endif /* sparse_opts_h */

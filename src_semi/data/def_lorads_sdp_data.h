#ifndef DEF_LORADS_SDP_DATA
#define DEF_LORADS_SDP_DATA

#include "lorads.h"
#include "def_lorads_elements.h"


typedef struct {
    lorads_int row;
    lorads_int col;
    lorads_int index;
} SparseElement;

typedef struct DictNode {
    SparseElement element;
    struct DictNode *next;
} DictNode;

typedef struct {
    DictNode **table;
    lorads_int size;
} Dict;

typedef enum {
    SDP_COEFF_ZERO,
    SDP_COEFF_SPARSE,
    SDP_COEFF_DENSE,
} sdp_coeff_type;

typedef struct{
    lorads_int        nSDPCol;

    sdp_coeff_type dataType;
    void      *dataMat;

    void (*create)                      ( void **, lorads_int, lorads_int, lorads_int *, double * );
    void (*scal)                        ( void *, double );
    lorads_int (*getnnz)                ( void * );
    void (*getmatnz)                    ( void *, lorads_int * );
    void (*add2buffer)                  ( void *, double, lorads_int *, double *);
    void (*destroy)                     ( void ** );
    void (*view)                        ( void * );

    void (*mul_rk)                      (void *, lorads_sdp_dense *, double *);
    void (*mv)                          (void *, double *, double *, lorads_int);
    void (*mul_inner_rk_double)         (void *, lorads_sdp_dense *, lorads_sdp_dense *, double *, void *, sdp_coeff_type);
    void (*zeros)                       (void *);
    void (*add_sdp_coeff)               (void *, void *, double, sdp_coeff_type);
    void (*nrm1)                        (void *, double *);
    void (*nrm2Square)                  (void *, double *);
    void (*nrmInf)                      (void *, double *);
    void (*statNnz)                     (lorads_int *);
    void (*scaleData)                   (void *, double);
    void (*collectNnzPos)               (void *, lorads_int *, lorads_int *, lorads_int *);
    void (*reConstructIndex)            (void *, Dict *);
}sdp_coeff;


typedef struct {

    /* In a zero matrix. There is nothing but an integer recording matrix dimension */

    lorads_int nSDPCol;

} sdp_coeff_zero;


typedef struct {

    /* In a sparse matrix, we adopt the triplet format i, j, x */
    lorads_int      nSDPCol;
    lorads_int      nTriMatElem;
    lorads_int      *triMatCol;
    lorads_int      *triMatRow;
    double          *triMatElem;
//    lorads_int      **rowCol2NnzIdx;
    lorads_int      *nnzIdx2ResIdx;
} sdp_coeff_sparse;


typedef struct {

    /* In a dense matrix, we store an n * (n + 1) / 2 array in packed format */
    lorads_int     nSDPCol;
    double        *dsMatElem;
//    lorads_int    **rowCol2NnzIdx;
    double         *fullMat;   // UVt full matrix, not dataMat full version
} sdp_coeff_dense;




#endif
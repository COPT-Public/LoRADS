#ifndef DEF_ASDP_SDPDATA_H
#define DEF_ASDP_SDPDATA_H

#ifdef HEADERPATH
#include "interface/asdp.h"
#include "src/asdp_linsolver.h"
#include "src/asdp_rk_mat.h"
#else
#include "asdp.h"
#include "def_asdp_rk_mat.h"
#endif

/* Implementations of the SDP coefficient matrix
 In ASDP, we implement five data structures for SDP coefficient matrix
 
 1. Zero matrix
 2. Sparse matrix
 3. Dense matrix
 4. Sparse rank one matrix
 5. Denser rank one matrix
 
 Each type of coefficient will be asociated with a data type, its eigen-decomposition
 and several methods that make it work in the ASDP solver.
 
 */
typedef enum {
    SDP_COEFF_ZERO,
    SDP_COEFF_SPARSE,
    SDP_COEFF_DENSE,
    SDP_COEFF_SPR1,
    SDP_COEFF_DSR1,
} sdp_coeff_type;

struct sdp_coeff_s {
    
    int        nSDPCol;
    
    sdp_coeff_type dataType;
    void      *dataMat;
    
    int eigRank;
    double *eigVals;
    double *eigVecs;
    
    asdp_retcode (*create) ( void **, int, int, int *, double * );
    void (*scal)          ( void *, double );
    double (*norm)        ( void *, int );
    asdp_retcode (*eig)  ( void *, int *, double *, double **, double ** );
    int (*getnnz)         ( void * );
    void (*dump)          ( void *, double * );
    void (*getmatnz)      ( void *, int * );
    void (*add2buffer)    ( void *, double, int *, double *);
    void (*destroy)       ( void ** );
    void (*view)          ( void * );
    int  (*iseye)         ( void *, double * );
    int  (*isunitcol)     ( void *, int * );
    // mul_rk(void *A, asdp_rk_mat *X, asdp_rk_mat *AX)
    asdp_retcode (*mul_rk) (void *, asdp_rk_mat_dense *, double *);
    asdp_retcode (*mul_inner_rk_double) (void *, asdp_rk_mat_dense *, asdp_rk_mat_dense *, double *, void *, sdp_coeff_type, int);
    void         (*zeros)        (void *);
    void         (*add_sdp_coeff) (void *, void *, double, sdp_coeff_type);
    void         (*nrm1) (void *, double *);
    void         (*nrm2Square) (void *, double *);
    void         (*nrmInf) (void *, double *);
    void         (*statNnz)(int *);
    void         (*addPreprocessRankOneConeDetect) (void *, int *);
    void         (*scaleData) (void *, double);
//    void         (*chol) (void *, double, int *);// just for slack variable 
};

typedef struct sdp_coeff_s sdp_coeff;

typedef struct {
    
    /* In a zero matrix. There is nothing but an integer recording matrix dimension */
    
    int nSDPCol;
    
} sdp_coeff_zero;

typedef struct {
    
    /* In a sparse matrix, we adopt the triplet format i, j, x */
    int     nSDPCol;
    int     nTriMatElem;
    int    *triMatCol;
    int    *triMatRow;
    double *triMatElem;
    int    **rowCol2NnzIdx;
} sdp_coeff_sparse;


typedef struct {
    
    /* In a dense matrix, we store the matrix in full col-major order, but only the lower triangle is maintained */
    int     nSDPCol;
    double *dsMatElem;
    int    **rowCol2NnzIdx;
} sdp_coeff_dense;


typedef struct {
    
    /* In the sparse rank 1 structure, we store the number of
       SDP columns as well as the coefficient sign, so that
     
       A = s * a * a',
        s: spR1FactorSign
     
       where a is represented using a triplet format.
       Also a full vector is used to save the effort to expand the sparse vector
     
     */
    
    int     nSDPCol; // dimension
    double  spR1FactorSign; ///< Include scale, may not equal to +1.0 or -1.0
    int     nSpR1FactorElem; // nonzero element number
    int    *spR1MatIdx; // indices
    double *spR1MatElem; // element
    double *spR1MatFactor; // record all factor
    double *nnzVal;
} sdp_coeff_spr1;


typedef struct {
    
    /* A dense rank 1 matrix is stored in a similar way to the sparse version
       There is no sparse representation of the factor
     */

    int     nSDPCol;
    double  r1FactorSign; ///< Include scale, may not equal to +1.0 or -1.0
    double *r1MatFactor;
    double *fullLowTriElem;
} sdp_coeff_dsr1;

extern asdp_retcode dataMatDenseMultiRkMat(void *A, asdp_rk_mat_dense *X, double *AX);
extern void dataMatDenseZeros(void *A);
extern asdp_retcode dataMatSparseMultiRkMat(void *A, asdp_rk_mat_dense *X, double *AX);
extern void dataMatSparseZeros(void *A);
extern void rkMatZeroDestroy(void **rk_mat);
extern void rkMatRowSparseDestroy(void **rk_mat);
extern void rkMatRowDenseDestroy(void **rk_mat);
#endif /* def_asdp_sdpdata_h */

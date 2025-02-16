#ifndef DEF_LORADS_CGS
#define DEF_LORADS_CGS

#include "def_lorads_elements.h"

typedef enum {

    CG_ITER_STATUS_OK,
    CG_ITER_STATUS_NUMERICAL,
    CG_ITER_STATUS_MAXITER,
    CG_ITER_STATUS_FAILED

} iter_status_lorads_cg;



typedef struct {

    lorads_int nr; // dimension * r
    lorads_int nnzRowOneCol; // non zero row number
    lorads_int nnzRow; // non zero row number * r
    lorads_int m; // non zero constraint number
    lorads_int *rowIdxInv; // when row sparse, it used
    lorads_int *rowIdx; // when row sparse, it used

    lorads_sdp_dense **a; // data, not allocated when created
    double *A; // use when row sparse

    double *rIter; // r
    double *rIterNew; // r_new
    double *pIter; // p
    double *qIter; // q
    double *qIterNew; // q_new
    double *QIter; // Q
    double *iterVec; // solution x

    /* Pre-conditioner */
    lorads_int useJacobi;
//    double *JacobiPrecond;


    // Statistics
    double solveTime;
    double cgDuration;
    double resiNorm;
    lorads_int    iter;

    iter_status_lorads_cg solStatus;

    // parameters
    lorads_int nRestartFreq;

    void  *MMat;
    void (*Mvec)  (void *, double *, double *);

} lorads_cg_linsys;




#endif
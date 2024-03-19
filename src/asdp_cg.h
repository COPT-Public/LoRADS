#ifndef ASDP_CG_H
#define ASDP_CG_H


#ifdef HEADERPATH
#include "linalg/def_asdp_linsolver.h"
#include "interface/asdp.h"
#else
#include "def_asdp_linsolver.h"
#include "asdp.h"
#include "def_asdp_sdpdata.h"
#include "def_asdp_conic.h"
#endif

typedef enum {
    
    CG_ITER_STATUS_OK,
    CG_ITER_STATUS_NUMERICAL,
    CG_ITER_STATUS_MAXITER,
    CG_ITER_STATUS_FAILED
    
} iter_status_asdp_cg;



typedef struct {

    int nr; // dimension * r
    int nnzRowOneCol; // non zero row number
    int nnzRow; // non zero row number * r
    int m; // non zero constraint number
    int *rowIdxInv; // when row sparse, it used
    int *rowIdx; // when row sparse, it used

    asdp_rk_mat_dense **a; // data, not allocated when created
    double *A; // use when row sparse

    double *rIter; // r
    double *rIterNew; // r_new
    double *pIter; // p
    double *qIter; // q
    double *qIterNew; // q_new
    double *QIter; // Q
    double *iterVec; // solution x

    /* Pre-conditioner */
    int useJacobi;
    double *JacobiPrecond;
    

    // Statistics
    double solveTime;
    double cgDuration;
    double resiNorm;
    int    iter;

    iter_status_asdp_cg solStatus;

    // parameters
    int maxIter;
    int nRestartFreq;
    double absTol;
    double relTol;
    
    void  *MMat;
    void (*Mvec)  (void *, double *, double *);

} asdp_cg_linsys;



extern asdp_retcode ASDPCGSolverCreate(asdp_cg_linsys **pCGSolver, int blkDim, int rank, int nConstr);
extern void ASDPCGSolverReCreate(asdp_cg_linsys **pCGSolver, int blkDim, int rank, int nConstr);
extern asdp_retcode CGSolve(void *linSys, double *x, double *b, double constr_vio);
extern void CGSetData (asdp_cg_linsys *cg, void *MMat, void (*Mvec) (void *, double *, double *));
extern void CGSolverClear(void *pCGSolver);

#endif //ASDP_CG_H

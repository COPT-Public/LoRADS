/** @file def\_asdp\_linsolver.h
 *  @brief ASDP linear system solver
 */
#ifndef DEF_ASDP_LINSOLVER_H
#define DEF_ASDP_LINSOLVER_H

#ifdef HEADERPATH
#include "interface/asdp.h"
#else
#include "asdp.h"
#endif

/* In ASDP, there are two cases where positve definite matrices need to
   be factorized: when inverting the dual matrix and when solving the Schur complement system.
   According to the sparsity of the matrix, we choose differnt data structures to do the job
 
   Dual matrix: sparse/dense Cholesky
   Schur matrix: sparse/dense Cholesky + Pre-conditioned conjugate gradient/residual
 
 */

typedef enum {
    /* Direct solver supports both Schur complement and matrix variable */
    ASDP_LINSYS_DENSE_DIRECT,
    ASDP_LINSYS_SMALL_DIRECT, /* Tailored for extremely small cones of size <= 3 */
    ASDP_LINSYS_SPARSE_DIRECT,
    
    /* Iterative solver is only used for Schur complement */
    ASDP_LINSYS_DENSE_ITERATIVE,
    ASDP_LINSYS_DENSE_INDEFINITE
    
} linsys_type;

typedef struct {
    
    int nCol;
    void *chol;
    linsys_type LinType;
    
    asdp_retcode (*cholCreate) ( void **, int );
    void (*cholSetParam) ( void *, void * );
    asdp_retcode (*cholSymbolic) ( void *, int *, int * );
    asdp_retcode (*cholNumeric) ( void *, int *, int *, double * );
    asdp_retcode (*cholPsdCheck) ( void *, int *, int *, double *, int * );
    
    void (*cholFSolve) ( void *, int, double *, double * );
    void (*cholBSolve) ( void *, int, double *, double * );
    asdp_retcode (*cholSolve) ( void *, int, double *, double * );
    asdp_retcode (*cholGetDiag) ( void *, double * );
    void (*cholInvert) ( void *, double *, double * );
    
    void (*cholDestroy) ( void ** );
    
    int nSolves;
    int nFactorizes;
    
} asdp_linsys_fp;

/* Sparse direct */
typedef struct {
    
    int nCol;
    int *colMatBeg;
    int *colMatIdx;
    double *colMatElem;
    
    double *dWork;
    
    void *pt[64];
    int iparm[64];
    
} pardiso_linsys;

/* Dense direct */
typedef struct {
    
    int nCol;
    double *dFullMatElem;
    
} lapack_flinsys;

typedef struct {
    
    int nCol;
    double *dFullMatElem;
    int *iAuxiIPIV;
    int iLwork;
    double *dWork;
    
} lapack_indef_flinsys;

typedef struct {
    
    int nCol;
    double dSmallMatElem[9];
    
} small_linsys;

typedef struct {
    
    int useJacobi;
    int maxIter;
    int nRestartFreq;
    double absTol;
    double relTol;
    
} iterative_params;

typedef enum {
    
    ITERATIVE_STATUS_OK,
    ITERATIVE_STATUS_NUMERICAL,
    ITERATIVE_STATUS_MAXITER,
    ITERATIVE_STATUS_FAILED
    
} iter_status;

/* Dense iterative */
typedef struct {
    
    int nCol;
    
    double *fullMatElem;
    double *iterResi;
    double *iterResiNew;
    double *iterDirection;
    double *preInvResi;
    double *MTimesDirection;
    double *iterVec;
    double *rhsBuffer;
    
    /* Pre-conditioner */
    int useJacobi;
    double *JacobiPrecond;
    lapack_flinsys *lap;
    
    /* Statistics */
    double iterResiNorm;
    double solveTime;
    
    int nIters;
    iter_status solStatus;
    int nSolves;
    
    iterative_params params;
    
} iterative_linsys;

/* Sparse direct */
#define PARDISO_RET_OK          ( 0)
#define PARDISO_RET_INDEFINITE  (-4)
#define PARDISO_SYM_POSDEFINITE ( 2)
#define PARDISO_SYM_INDEFINITE  (-2)
#define PARDISO_PHASE_SYM       (11)      // Pardiso symbolic analysis
#define PARDISO_PHASE_SYM_FAC   (12)      // Pardiso symbolic analysis
#define PARDISO_PHASE_FAC       (22)      // Pardiso numerical factorization
#define PARDISO_PHASE_FORWARD  (331)
#define PARDISO_PHASE_BACKWARD (333)
#define PARDISO_PHASE_SOLVE     (33)      // Solve linear system
#define PARDISO_PHASE_FREE      (-1)      // Free internal data structure

#define PARDISO_PARAM_NONDEFAULT    (0)
#define PARDISO_PARAM_SYMBOLIC      (1)
#define PARDISO_PARAM_SYMBOLIC_MMD  (0)
#define PARDISO_PARAM_SYMBOLIC_ND   (2)
#define PARDISO_PARAM_REFINEMENT    (7)
#define PARDISO_PARAM_INPLACE       (5)
#define PARDISO_PARAM_PERTURBATION  (9)
#define PARDISO_PARAM_SCALING      (10)
#define PARDISO_PARAM_MATCHING     (12)
#define PARDISO_PARAM_FACNNZ       (17)
#define PARDISO_PARAM_FACFLOP      (18)
#define PARDISO_PARAM_THREADS      (33)
#define PARDISO_PARAM_INDEX        (34)
#define PARDISO_PARAM_INDEX_C       (1)
#define PARDISO_PARAM_DIAGONAL     (55)
#define PARDISO_PARAM_DIAGONAL_ON   (1)

#define set_pardiso_param(iparm, param, val) iparm[param] = val
#define get_pardiso_output(iparm, param) iparm[param]

extern void pardisoinit ( void *, int *, int * );
extern void pardiso     ( void     *, int    *, int *, int *, int *, int *,
                          double   *, int    *, int *, int *, int *, int *,
                          int *, double      *, double   *, int * );
extern void pardiso_getdiag ( const void *, void *, void *, const int *, int * );




#endif /* DEF_ASDP_LINSOLVER_H */

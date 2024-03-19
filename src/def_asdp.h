#ifndef DEF_ASDP_H
#define DEF_ASDP_H

#ifdef HEADERPATH
#include "interface/asdp.h"
#include "interface/asdp_conic.h"
#include "interface/asdp_utils.h"
#include "interface/asdp_user_data.h"
#else
#include "asdp.h"
#include "asdp_conic.h"
#include "asdp_utils.h"
#include "asdp_user_data.h"
#include "asdp_cg.h"
#include "def_asdp_lanczos.h"
#include "asdp_lanczos.h"
#endif

#define NUM_INT_PARAM 20
#define NUM_DBL_PARAM 20

#define INT_FEATURE_I_NULLOBJ 0
#define INT_FEATURE_I_MANYCONES 1
#define INT_FEATURE_I_NOPINTERIOR 2
#define INT_FEATURE_I_NODINTERIOR 3
#define INT_FEATURE_I_VERYDENSE 4
#define INT_FEATURE_I_IMPTRACE 5
#define INT_FEATURE_I_IMPYBOUND 6
#define INT_FEATURE_N_SUMCONEDIMS 7
#define INT_FEATURE_N_MAXCONEDIM 8
#define INT_FEATURE_N_CONES 9
#define INT_FEATURE_N_ROWS 10
#define INT_FEATURE_N_SPSDPCONES 11
#define INT_FEATURE_N_DSSDPCONES 12
#define INT_FEATURE_N_LPCONES 13
#define INT_FEATURE_N_BNDCONES 14
#define INT_FEATURE_N_ZEORMATS 15
#define INT_FEATURE_N_SPMATS 16
#define INT_FEATURE_N_DSMATS 17
#define INT_FEATURE_N_SPR1MATS 18
#define INT_FEATURE_N_DSR1MATS 19

#define DBL_FEATURE_OBJFRONORM 0
#define DBL_FEATURE_OBJONENORM 1
#define DBL_FEATURE_RHSFRONORM 2
#define DBL_FEATURE_RHSONENORM 3
#define DBL_FEATURE_RHSINFNORM 4
#define DBL_FEATURE_OBJSCALING 5
#define DBL_FEATURE_RHSSCALING 6
#define DBL_FEATURE_DATAFRONORM 7
#define DBL_FEATURE_DATAONENORM 8
#define DBL_FEATURE_IMPYBOUNDUP 9
#define DBL_FEATURE_IMPYBOUNDLOW 10
#define DBL_FEATURE_IMPTRACEX 11

#define ASDP_EPSILON 1e-8

#define BMMethod 1
#define ADMMMethod 0

struct lbfgs_node_internal
{
    int allElem;

    double *s;    // difference of R, Rk - Rk-1
    double *y;    // difference of grad
    double beta;  // 1 / <y, s>
    double alpha; // beta <s, q>
    struct lbfgs_node_internal *next;
    struct lbfgs_node_internal *prev;
};

typedef struct lbfgs_node_internal lbfgs_node;

struct asdp_solver_internal
{

    /* User data */
    int nRows;      // constraint number
    double *rowRHS; // b of Ax = b

    /* Cones */
    int nCones; // sdp cones block number
    asdp_cone **ACones;
    asdp_cg_linsys **CGLinsys;
    asdp_lanczos **LanczoSys;
    double **LanczosStart;

    /* Auxiliary variable SDPCone */
    constrValStruct **constrVal; // constraint violation [iCone][iConstr]
    double *constrValSum;        // constraint violation [iConstr]
    double *ARDSum;              // q1 in line search of BM
    double *ADDSum;              // q2 in line search of BM
    double **bLinSys;            // for solving linear system, bInLinSys[iCone]
    int *rankElem;               // all rank
    double *M1temp;              // M1 for solving linear system
    double *bestDualVar;
    asdp_rk_mat_dense **M2temp; // M2 for solving linear system

    /* Variables SDPCone */
    asdp_rk_mat_dense **U;    // admm variable, and lbfgs descent direction D
    asdp_rk_mat_dense **V;    // admm variable only
    asdp_rk_mat_dense **R;    // average variable for storage and BM variable
    asdp_rk_mat_dense **Grad; // grad of R

    /* Auxiliary variable LPCone (SDPCone but rank is 1, dim is 1) */
    asdp_rk_mat_lp *rLp;
    asdp_rk_mat_lp *uLp;
    asdp_rk_mat_lp *vLp;
    asdp_rk_mat_lp *vlagLp;
    asdp_rk_mat_lp *lagLp;
    asdp_rk_mat_lp *gradLp;
    asdpLPCone *lpCone;
    int nLpCols;
    constrValStruct **constrValLP;

    /* BM lbfgs and ADMM Variables */
    double *dualVar;
    asdp_rk_mat_dense **Vlag; // ADMM dual variable and lbfgs gradient
    asdp_rk_mat_dense **lag;

    /* BM lbfgs Variables */
    int hisRecT;          // all node number is hisRecT + 1
    lbfgs_node *lbfgsHis; // record difference of primal variable history and gradient history, lbfgsHis[iCone]
    int maxBMInIter;      // inner iteration
    int maxBMOutIter;     // outer interation, increase penalty

    /* ADMM Parameters */
    int maxInter;

    /* Monitor */
    int whichMethod;
    int nIterCount;
    double cgTime;
    int cgIter;
    int checkSolTimes;
    double traceSum;

    /* Convergence criterion */
    double pObjVal;
    double dObjVal;
    double pInfeas;
    double dInfeas;
    double cObjNrm1;
    double cObjNrm2;
    double cObjNrmInf;
    double bRHSNrm1;
    double bRHSNrmInf;
    double bRHSNrm2;
    double *dimacError;
    double *constrVio;
    double **UsubV;
    double *uSubvLp;
    double *negLambd;

    asdp_status AStatus;

    /* Starting time */
    double dTimeBegin;

    /* Parameters */
    double rho; // BM and ADMM
    double rhoMax;
    int strategy;

    // scale
    double cScaleFactor;
    double bScaleFactor;

    // check exit bm
    int *rank_max;
    double *sparsitySDPCoeff;
    double overallSparse;
    int    nnzSDPCoeffSum;
    int    SDPCoeffSum;
};

#endif /* def_asdp_h */

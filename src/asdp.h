#ifndef ASDP_H
#define ASDP_H

/* Macro for MATLAB */
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define calloc mxCalloc
#define printf mexPrintf
#define free mxFree
#else
#endif

// #define ASDP_LANCZOS_DEBUG
// #define ASDP_CONIC_DEBUG
// #define ASDP_KKT_DEBUG
// #define ASDP_LINSYS_DEBUG
// #define ASDP_CONJGRAD_DEBUG
// #define ASDP_ALGO_DEBUG
#define ASDP_SPARSE_CONE_THRESHOLD (0.3)
#define USE_CG
// #define USE_DI

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

typedef enum
{

    ASDP_RETCODE_OK,
    ASDP_RETCODE_FAILED,
    ASDP_RETCODE_MEMORY,
    ASDP_RETCODE_EXIT,
    ASDP_RETCODE_RANK,
    ASDP_RECODE_SPLIT

} asdp_retcode;

typedef enum
{

    ASDP_UNKNOWN,
    ASDP_DUAL_FEASIBLE,
    ASDP_DUAL_OPTIMAL,
    ASDP_PRIMAL_DUAL_OPTIMAL,
    ASDP_PRIMAL_OPTIMAL,
    ASDP_MAXITER,
    ASDP_SUSPECT_INFEAS_OR_UNBOUNDED,
    ASDP_INFEAS_OR_UNBOUNDED,
    ASDP_TIMELIMIT,
    ASDP_USER_INTERRUPT,
    ASDP_INTERNAL_ERROR,
    ASDP_NUMERICAL

} asdp_status;

typedef struct asdp_solver_internal asdp;
#include "def_asdp_rk_mat.h"

#define ASDP_INFINITY 1e+30

// Integer Parameters
#define INT_PARAM_MAXITER 0
#define INT_PARAM_CORRECTORA 1
#define INT_PARAM_CORRECTORB 2
#define INT_PARAM_THREADS 3

#define FLAG_U 1
#define FLAG_V -1
#define FLAG_S 0
#define FLAG_R 2
#define FLAG_INI 3
#define FLAG_OBJ 4
#define FLAG_WSUM 5
#define FLAG_UVt 6

// BM status
#define EASY ('e')
#define MEDIUM ('m')
#define HARD ('h')
#define SUPER ('s')

// Strategy option
#define STRATEGY_MIN_BISECTION 1
#define STRATEGY_MAX_CUT 2
#define STRATEGY_SNL 3
#define STRATEGY_OPF 4
#define STRATEGY_GENERAL 5
#define STRATEGY_GAMMA_DYNAMIC 6

// Version information
#define VERSION_MAJOR 1
#define VERSION_MINOR 0
#define VERSION_TECHNICAL 0

// Build date
#define BUILD_DATE_YEAR 2023
#define BUILD_DATE_MONTH 09
#define BUILD_DATE_DAY 27

// Termination Criteria
#define ASDP_TERMINATION_TOL (1e-5)
#define TERMINATION_SCHEME_BM
// #define TERMINATION_SCHEME_OFFICIAL

#ifdef TERMINATION_SCHEME_OFFICIAL
#undef TERMINATION_SCHEME_BM
#else
#ifndef TERMINATION_SCHEME_BM
#define TERMINATION_SCHEME_BM
#endif
#endif

// #define NO_PENALTY_METHOD
#ifdef __cplusplus
extern "C"
{
#endif

    extern asdp_retcode ASDPOptimize(asdp *ASolver, int rhoFreq, double rhoFactor, int rhoStrategy, double tau, double gamma, double rhoMin, double ori_start, double timeLimit);
    extern asdp_retcode ASDPDualOptimize(asdp *ASolver);
    extern asdp_retcode ASDP_BMOptimize(asdp *ASolver, double endBMTol, double endBMTol_pd, double endTauTol, double endBMALSub, double ori_start, int is_rank_max, int *pre_mainiter, int *pre_miniter, double timeLimit);
    extern int BMWarmStart(char *filename, double **RPointer, double **dualVarRPointer, double *rhoRPointer);
    extern void ASDP_BMtoADMM(asdp *ASolver, double heuristic);
    extern void BMWarmStart2ADMM(asdp *ASolver, double *R, double *dualVar, double rho);
    //    extern asdp_retcode endBMwarmStart(asdp *ASolver);
    extern void ASDPCheckDimacErr(asdp *ASolver);
    extern void ASDPCheckDimacErrBMCriteria(asdp *ASolver);
    //    extern asdp_retcode ASDPCheckSolverStatus(asdp *ASolver);
    extern asdp_retcode ASDPCheckSolverStatusBM(asdp *ASolver);
    extern void ASDP_SCALE(asdp *ASolver);
    extern void ASDPEndProgram(asdp *ASolver);
    extern void detectMaxCutProb(asdp *ASolver, int *blkDim, int *maxCut);
    extern void detectSparsitySDPCoeff(asdp *ASolver);
    extern void freeDetectSparsitySDPCoeff(asdp *ASolver);
    extern asdp_retcode ASDPPreprocess(asdp *ASolver, int *BlkDims);
    extern void destroyPreprocess(asdp *ASolver);

#ifdef __cplusplus
}
#endif

#endif /* asdp_h */

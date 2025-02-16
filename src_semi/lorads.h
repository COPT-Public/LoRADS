#ifndef LORADS_H
#define LORADS_H

/* Macro for MATLAB */
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define calloc mxCalloc
#define printf mexPrintf
#define free mxFree
#else
#endif

#define LORADS_SPARSE_CONE_THRESHOLD  (0.3)
#define USE_CG

#ifdef MEMDEBUG
#include "memwatch.h"
#endif


// Type define
#ifdef MAC_INT64
#define lorads_int     int64_t
#endif

#ifdef UNIX_INT64
#define lorads_int     int64_t
#endif

#ifdef INT32
#define lorads_int     int
#endif




typedef enum {
    LORADS_RETCODE_OK,
    LORADS_RETCODE_FAILED,
    LORADS_RETCODE_MEMORY,
    LORADS_RETCODE_EXIT,
    LORADS_RETCODE_RANK,
} lorads_retcode;

typedef enum{
    LORADS_UNKNOWN,
    LORADS_PRIMAL_DUAL_OPTIMAL,
    LORADS_PRIMAL_OPTIMAL,
    LORADS_MAXITER,
    LORADS_TIME_LIMIT,
}lorads_status;


#define LORADS_INFINITY          1e+30

// ALM status
#define EASY                   ('e')
#define MEDIUM                 ('m')
#define HARD                   ('h')
#define SUPER                  ('s')

#define RET_CODE_OK            (0)
#define RET_CODE_TIME_OUT      (1)
#define RET_CODE_NUM_ERR       (4)
#define RET_CODE_BAD_ITER      (8)
extern int MAX_ALM_SUB_ITER;


// Build date
#define BUILD_DATE_YEAR         (2024)
#define BUILD_DATE_MONTH        (08)
#define BUILD_DATE_DAY          (27)


// some constants
#define LORADS_TERMINATION_TOL (1e-5)


#include <getopt.h>
#include <stdbool.h>

typedef struct{
    char *fname;
    double initRho;
    double rhoMax;
    double rhoCellingALM;
    double rhoCellingADMM;
    lorads_int maxALMIter;
    lorads_int maxADMMIter;
    double timesLogRank;
    lorads_int rhoFreq;
    double rhoFactor;
    double ALMRhoFactor;
    double phase1Tol;
    double phase2Tol;
    double timeSecLimit;
    double heuristicFactor;
    lorads_int lbfgsListLength;
    double endTauTol;
    double endALMSubTol;
    bool l2Rescaling;
    lorads_int reoptLevel;
    lorads_int dyrankLevel;
    bool highAccMode;
} lorads_params;




#define NO_PENALTY_METHOD
#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif /* lorads_h */

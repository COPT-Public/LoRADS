#ifndef LORADS_UTILS_H
#define LORADS_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "lorads.h"

/* Define macros */
#define LORADS_ERROR_TRACE printf("File [%30s] Line [%d]\n", __FILE__, __LINE__)
#define LORADS_FREE(vars)    \
    do                    \
    {                     \
        if (vars)          \
        {                 \
            free((vars));  \
            (vars) = NULL; \
        }                 \
    } while (0)
#define LORADS_INIT(var, type, size) (var) = (type *)calloc(size, sizeof(type))
#define LORADS_REALLOC(var, type, size) (var) = (type *)realloc(var, sizeof(type) * (size))
#define LORADS_MEMCPY(dst, src, type, size) memcpy(dst, src, sizeof(type) * (size))
#define LORADS_ZERO(var, type, size) memset(var, 0, sizeof(type) * (size))
#define LORADS_NULLCHECK(var)            \
    if (!(var))                        \
    {                                  \
        LORADS_ERROR_TRACE;             \
    }
#define LORADS_MEMCHECK(var)             \
    if (!(var))                        \
    {                                    \
        LORADS_ERROR_TRACE;                                 \
    }
#define LORADS_MAX(x, y) ((x) > (y) ? (x) : (y))
#define LORADS_MIN(x, y) ((x) < (y) ? (x) : (y))
#define LORADS_ABS(x) fabs(x)
#define LORADS_DIMAC_ERROR_CONSTRVIO_L1 (0)
#define LORADS_DIMAC_ERROR_PDGAP (1)
#define LORADS_DIMAC_ERROR_DUALFEASIBLE_L1 (2)

#define PI (3.1415926)

#define PACK_NNZ(n) ((n) * ((n) + 1) / 2)
#define PACK_IDX(n, i, j) (lorads_int)((2 * (n) - (j)-1) * (j) / 2) + (i)
#define FULL_IDX(n, i, j) ((j) * (n) + (i))
#define PACK_ENTRY(A, n, i, j) (A[(lorads_int)((2 * (n) - (j)-1) * (j) / 2) + (i)])
#define FULL_ENTRY(A, n, i, j) (A[(j) * (n) + (i)])

#ifdef __cplusplus
extern "C"
{
#endif

extern double LUtilGetTimeStamp(void);
extern void LUtilMatSymmetrize(lorads_int n, double *v);
extern lorads_int LUtilCheckIfAscending(lorads_int n, lorads_int *idx);
extern void LUtilDescendSortIntByInt(lorads_int *data, lorads_int *ref, lorads_int low, lorads_int up);
extern void LUtilSortIntbyDbl(lorads_int *data, double *ref, lorads_int low, lorads_int up);
extern void LUtilAscendSortDblByInt(double *data, lorads_int *ref, lorads_int low, lorads_int up);
extern void LUtilAscendSortDblByInt(double *data, lorads_int *ref, lorads_int low, lorads_int up);

extern void LUtilStartCtrlCCheck(void);
extern lorads_int LUtilCheckCtrlC(void);
extern void LUtilResetCtrl(void);

extern lorads_int LUtilGetGlobalMKLThreads(void);
extern void LUtilSetGlobalMKLThreads(lorads_int nTargetThreads);

extern lorads_int LUtilUpdateCheckEma(double *current_ema, double *old_ema, double new_value, double alpha, double threshold, lorads_int update_interval, lorads_int *counter);
extern void REALLOC(double **data, lorads_int nOld, lorads_int nNew);

#ifdef __cplusplus
}
#endif

#endif /* LORADS_UTILS_H */

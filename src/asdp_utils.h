#ifndef ASDP_UTILS_H
#define ASDP_UTILS_H

#ifdef HEADERPATH
#include "interface/asdp.h"
#else
#include "asdp.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* Define macros */
#ifdef SILENT_SOLVER
#define ASDP_printf
#else
#define asdp_printf printf
#endif

#define ASDP_FREE(var)    \
    do                    \
    {                     \
        if (var)          \
        {                 \
            free((var));  \
            (var) = NULL; \
        }                 \
    } while (0)
#define ASDP_INIT(var, type, size) (var) = (type *)calloc(size, sizeof(type))
#define ASDP_REALLOC(var, type, size) (var) = (type *)realloc(var, sizeof(type) * (size))
#define ASDP_MEMCPY(dst, src, type, size) memcpy(dst, src, sizeof(type) * (size))
#define ASDP_ZERO(var, type, size) memset(var, 0, sizeof(type) * (size))
#define ASDP_NULLCHECK(var)            \
    if (!(var))                        \
    {                                  \
        retcode = ASDP_RETCODE_FAILED; \
        goto exit_cleanup;             \
    }
#define ASDP_MEMCHECK(var)             \
    if (!(var))                        \
    {                                  \
        retcode = ASDP_RETCODE_MEMORY; \
        goto exit_cleanup;             \
    }
#define ASDP_ERROR_TRACE printf("File [%30s] Line [%d]\n", __FILE__, __LINE__)
#define ASDP_CALL(func)             \
    retcode = (func);               \
    if (retcode != ASDP_RETCODE_OK) \
    {                               \
        goto exit_cleanup;          \
    }
#define ASDP_STATUS_CHECK(func)       \
    retcode = (func);                 \
    if (retcode == ASDP_RETCODE_EXIT) \
    {                                 \
        goto exit_cleanup;            \
    }
#define ASDP_MAX(x, y) ((x) > (y) ? (x) : (y))
#define ASDP_MIN(x, y) ((x) < (y) ? (x) : (y))
#define ASDP_ABS(x) fabs(x)

#define PACK_NNZ(n) ((n) * ((n) + 1) / 2)
#define PACK_IDX(n, i, j) (int)((2 * (n) - (j)-1) * (j) / 2) + (i)
#define FULL_IDX(n, i, j) ((j) * (n) + (i))
#define PACK_ENTRY(A, n, i, j) (A[(int)((2 * (n) - (j)-1) * (j) / 2) + (i)])
#define FULL_ENTRY(A, n, i, j) (A[(j) * (n) + (i)])

#define ASDP_PROFILER(func, ntest)                                \
    double HLocalTProfiler = AUtilGetTimeStamp();                 \
    for (int iTest = 0; iTest < (ntest); ++iTest)                 \
    {                                                             \
        (func);                                                   \
        printf("Run %d. Elapsed time: %f\n",                      \
               iTest + 1, AUtilGetTimeStamp() - HLocalTProfiler); \
    }                                                             \
    printf("Function Profiler: Line %d of %s by %d runs. "        \
           "Average running time: %fs\n",                         \
           __LINE__, __FILE__,                                    \
           (ntest), (AUtilGetTimeStamp() - HLocalTProfiler) / ntest);

#define ASDP_CODE_PROFILER_START double tASDPStart = AUtilGetTimeStamp()
#define ASDP_CODE_PROFILER_END             \
    printf("Code Profiler Line %d of %s. " \
           "Running time: %fs\n",          \
           __LINE__, __FILE__,             \
           (AUtilGetTimeStamp() - tASDPStart))

#define set_func_pointer(A, B) (A = (typeof(A))B)

#ifdef __cplusplus
extern "C"
{
#endif

    extern double AUtilGetTimeStamp(void);
    extern void AUtilMatSymmetrize(int n, double *v);
    extern int AUtilCheckIfAscending(int n, int *idx);
    extern void AUtilDescendSortIntByInt(int *data, int *ref, int low, int up);
    extern void AUtilSortIntbyDbl(int *data, double *ref, int low, int up);
    extern void AUtilAscendSortDblByInt(double *data, int *ref, int low, int up);
    extern void AUtilAscendSortDblByInt(double *data, int *ref, int low, int up);
    extern void AUtilPrintDblContent(int n, double *d);
    extern void AUtilPrintIntContent(int n, int *d);
    extern double AUtilPrintDblSum(int n, double *d);
    extern double AUtilPrintDblAbsSum(int n, double *d);
    extern int AUtilCheckUserInterrupt(void);

    extern void AUtilStartCtrlCCheck(void);
    extern int AUtilCheckCtrlC(void);
    extern void AUtilResetCtrl(void);

    extern int AUtilGetGlobalMKLThreads(void);
    extern void AUtilSetGlobalMKLThreads(int nTargetThreads);

    extern void ASDP_ONE(double *var, int size);
    extern int AUtilUpdateCheckEma(double *current_ema, double *old_ema, double new_value, double alpha, double threshold, int update_interval, int *counter);
extern void REALLOC(double **data, int nOld, int nNew);
#ifdef __cplusplus
}
#endif

#endif /* ASDP_UTILS_H */

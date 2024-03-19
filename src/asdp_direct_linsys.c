#ifdef HEADERPATH
#include "interface/asdp.h"
#include "interface/asdp_utils.h"
#include "src/def_asdp_linsolver.h"
#include "src/asdp_linsolver.h"
#include "src/dense_opts.h"
#include "src/vec_opts.h"
#else
#include "asdp.h"
#include "asdp_utils.h"
#include "dense_opts.h"
#include "vec_opts.h"
#include "asdp_cg.h"
#include "asdp_debug.h"
#include "asdp_direct_linsys.h"
#endif


#ifdef MEMDEBUG
#include "memwatch.h"
#endif

#include <math.h>
extern asdp_retcode ASDPDirectLinSysCreate(asdp_direct_linsys **linsys, int blkDim, int rank, int nConstr){
    asdp_retcode retcode = ASDP_RETCODE_OK;
    asdp_direct_linsys *linsysPtr;    
    ASDP_INIT(linsysPtr, asdp_direct_linsys, 1);
    ASDP_MEMCHECK(linsys);
    linsysPtr->uplo = 'L';
    linsysPtr->n = blkDim * rank;
    linsysPtr->nConstr = nConstr;
    linsysPtr->blkDim = blkDim;
    linsysPtr->rank = rank;

    linsysPtr->lda = linsysPtr->n;
    linsysPtr->ldb = linsysPtr->n;
    linsysPtr->info = 0;

    *linsys = linsysPtr;
exit_cleanup:
    return retcode;
}
extern void FullMatrixDebug(double *A, int n){
    for (int row = 0; row < n; ++row){
        for (int col = 0; col < n; ++col){
            printf("row:%d, col: %d, val:%f\n", row, col, A[col * n + row]);
            if (col > 5){
                break;
            }
        }
        if (row > 5){
            break;
        }
    }
}

extern void reconstructFullMatrix(double *a, int nr, int m, double *A){
    for (int r = 0; r < m; ++r){
        for (int i = 0; i < nr; ++i){
            for (int j = 0; j < nr; ++j){
                A[i * nr + j] += (a[r * nr + i] * a[r * nr + j]);
            }
        }
    }
    for (int i = 0; i < nr; ++i){
        A[i * nr + i] += 1.0;
    }
}



extern asdp_retcode DirectLinSysSolve(void *linSysPtr, double *a, double *b, double *res){
    asdp_direct_linsys *linSys = (asdp_direct_linsys *) linSysPtr;
    asdp_retcode retcode = ASDP_RETCODE_OK;
    int nr = linSys->n;
    int m = linSys->nConstr;
    double *A;
    ASDP_INIT(A, double, nr * nr);
    ASDP_MEMCHECK(A);
    ASDP_ZERO(A, double, nr * nr);
    reconstructFullMatrix(a, nr, m, A);
    ASDP_MEMCPY(res, b, double, nr);
#ifdef UNDER_BLAS
    dpotrf_(&(linSys->uplo), &(nr), A, &(linSys->lda), &(linSys->info));
    dpotrs_(&(linSys->uplo), &(nr), &(AIntConstantOne), A, &(linSys->lda), res, &(linSys->ldb), &(linSys->info));
#else
    dpotrf(&(linSys->uplo), &(nr), A, &(linSys->lda), &(linSys->info));
    dpotrs(&(linSys->uplo), &(nr), &(AIntConstantOne), A, &(linSys->lda), res, &(linSys->ldb), &(linSys->info));
#endif
exit_cleanup:
    ASDP_FREE(A);
    return retcode;
}

extern void directLinSysClear(void *linSysPtr){
    if (!linSysPtr){
        return;
    }
    asdp_direct_linsys *linSys = (asdp_direct_linsys *)linSysPtr;
    ASDP_FREE(linSys);
}

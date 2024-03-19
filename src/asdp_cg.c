#ifdef HEADERPATH
#include "interface/asdp.h"
#include "interface/asdp_utils.h"
#include "src/def_asdp_linsolver.h"
#include "src/dense_opts.h"
#include "src/vec_opts.h"
#else
#include "asdp.h"
#include "asdp_utils.h"
#include "def_asdp_linsolver.h"
#include "dense_opts.h"
#include "vec_opts.h"
#include "asdp_cg.h"
#include "asdp_debug.h"
#include "def_asdp_rk_mat.h"
#include "def_asdp_conic.h"
#include "asdp_algo.h"
#endif

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

#include <math.h>

extern asdp_retcode ASDPCGSolverCreate(asdp_cg_linsys **pCGSolver, int blkDim, int rank, int nConstr){
    /*  create a cg solver for one cone
        blkDim is the dimension of cone: n
        rank is the rank of cone: r
        nConstr is the number of constraints: m
     */
    asdp_retcode retcode = ASDP_RETCODE_OK;
    ASDP_NULLCHECK(pCGSolver);

    asdp_cg_linsys *cg = NULL;
    ASDP_INIT(cg, asdp_cg_linsys, 1);
    ASDP_MEMCHECK(cg);

    int nr = blkDim * rank; // dim = n * r

    cg->nr = nr; // set dimension of liear system
    cg->m = nConstr;
    cg->iter = 0;

    ASDP_INIT(cg->rIter, double, nr);
    ASDP_INIT(cg->rIterNew, double, nr);
    ASDP_INIT(cg->pIter, double, nr);
    ASDP_INIT(cg->qIter, double, nr);
    ASDP_INIT(cg->qIterNew, double, nr);
    ASDP_INIT(cg->QIter, double, nr);

    cg->useJacobi = 0;

    ASDP_INIT(cg->JacobiPrecond, double, nr);

    /* Set parameters */
    cg->absTol = 1e-8;
    cg->relTol = 1e-10;
    cg->maxIter = ASDP_MAX(100, cg->nr / 10);
    cg->nRestartFreq = 20;
    cg->solveTime = 0.0;

    cg->a = NULL;
    *pCGSolver = cg;
exit_cleanup:
    return retcode;
}


extern void ASDPCGSolverReCreate(asdp_cg_linsys **pCGSolver, int blkDim, int rank, int nConstr){
    asdp_cg_linsys *cg = *pCGSolver;
    int nr = blkDim * rank;
    cg->nr = nr; // set dimension of liear system
    cg->m = nConstr;
    cg->iter = 0;
    ASDP_REALLOC(cg->rIter, double, nr);
    ASDP_REALLOC(cg->rIterNew, double, nr);
    ASDP_REALLOC(cg->pIter, double, nr);
    ASDP_REALLOC(cg->qIter, double, nr);
    ASDP_REALLOC(cg->qIterNew, double, nr);
    ASDP_REALLOC(cg->QIter, double, nr);
    ASDP_REALLOC(cg->JacobiPrecond, double, nr);
}


extern void CGSetData (asdp_cg_linsys *cg, void *MMat, void (*Mvec) (void *, double *, double *)){
//    if (cg->MMat || cg->Mvec){
//        return;
//    }
    cg->MMat = MMat;
    cg->Mvec = Mvec;
}


extern asdp_retcode CGSolve(void *linSys, double *x, double *b, double constr_vio){

    asdp_cg_linsys *cg = (asdp_cg_linsys *)linSys;
    asdp_retcode retcode = ASDP_RETCODE_OK;
    cg->solStatus = CG_ITER_STATUS_MAXITER; // Initialize to max iter
    int cgMaxIter = cg->maxIter;

    /* CG working scalers */
    double qTr = 0.0; // for alpha, beta
    double qTrNew = 0.0; // for beta 
    double pTQ = 0.0; // for alpha
    double alpha = 0.0;
    double beta = 0.0;
    double resiNorm = 0.0; // for stopping 
    double bNorm = 0.0; // for rescale

    double *r = cg->rIter;
    double *p = cg->pIter;
    double *q = cg->qIter;
    double *qNew = cg->qIterNew;
    double *Q= cg->QIter;
    int nr = cg->nr;
    int m = cg->m;
    cg->absTol = ASDP_MAX((constr_vio * 1e-3), 1e-8);
    cg->relTol = ASDP_MAX((constr_vio * 1e-5), 1e-10);

    /* Get restart frequency */
    int nRestartFreq = ASDP_MAX(cg->nRestartFreq, 20);
    // set constant
    int incx = 1;
    double one = 1.0;
    double minusOne = -1.0;

    // start solving
    double cgStartTime = AUtilGetTimeStamp();
    double cgDuration = 0.0;

    /* Compute initial b norm */
    bNorm = nrm1(&nr, b, &incx);
    //    if ( bNorm > 1e-8 ) {
    //        scaleFactor = one/bNorm;
    //        // b = b / bNorm;
    //        scal(&nr, &scaleFactor, b, &incx);
    //        // a = a / bNorm;
    //        double sqrtScale = sqrt(scaleFactor);
    //        // int dim = nr * m;
    //        // scal(&dim, &sqrtScale, a, &incx);
    //        for (int i = 0; i < cg->m; ++i){
    //            a[i]->scal(a[i]->rk_mat, sqrtScale);
    //        }
    //    }
    
#ifdef ASDP_CG_DEBUG
    int dim = nr * m;
    double aNorm;
    aNorm = nrm2(&dim, a, &incx);
    asdp_printf("a norm after rescale %f\n", aNorm);
#endif
    
    // change a and initial solution
    // double *x= updateVar->matElem;
    
#ifdef ASDP_CG_DEBUG
    int Adim = nr * nr;
    double ANorm;
    ANorm = nrm2(&Adim, a, &incx);
    asdp_printf("a norm after solver Init %f\n", ANorm);
#endif
//    asdp_rk_mat_dense *shell;
//    ASDP_INIT(shell, asdp_rk_mat_dense, 1);
    /* Compute initial residual */
//    linSysProduct(ACone, weight, noUpdateVar, shell, updateVar->matElem, r);
    cg->Mvec(cg->MMat, x, r);
    // r = b - r
    axpy(&nr, &minusOne, b, &incx, r, &incx);
    scal(&nr, &minusOne, r, &incx);
    
    
    /* Compute initial residual norm */
    resiNorm = nrm2(&nr, r, &incx);
    if ( resiNorm < cg->absTol || resiNorm / bNorm < cg->relTol) {
        cg->solStatus = CG_ITER_STATUS_OK;
        goto exit_cleanup;
    }
    
    double zero = 0.0;
    // p = r;
    // ASDP_MEMCPY(p, r, double, nr);
    scal(&nr, &zero, p, &incx);
    axpy(&nr, &one, r, &incx, p, &incx);
    // q = r;
    // ASDP_MEMCPY(q, r, double, nr);
    scal(&nr, &zero, q, &incx);
    axpy(&nr, &one, r, &incx, q, &incx);
    // Initialize qTr
    qTr = dot(&nr, q, &incx, r, &incx);
    cg->iter = 0;
    for (int k = 0; k < cgMaxIter; ++k){
        cg->iter += 1;
#ifdef ASDP_CG_DEBUG
        asdp_printf("CGIter %d, normResi: %f\n", k, resiNorm);
#endif
//        linSysProduct(ACone, weight, noUpdateVar, shell, p, Q);
        cg->Mvec(cg->MMat, p, Q);
        qTr = dot(&nr, q, &incx, r, &incx);
        pTQ = dot(&nr, p, &incx, Q, &incx);
        alpha = qTr / pTQ;
        // x = x + alpha * p;
        axpy(&nr, &alpha, p, &incx, x, &incx);
        // r = r - alpha * Q;
        double negAlpha = -alpha;
        axpy(&nr, &negAlpha, Q, &incx, r, &incx);
        resiNorm = nrm2(&nr, r, &incx);
        cg->resiNorm = resiNorm;
        if ( resiNorm / bNorm < cg->relTol ) {
            cg->solStatus = CG_ITER_STATUS_OK;
            goto exit_cleanup;
        }
        if (k%nRestartFreq == 0){
            // r = b - a * x;
            // linSysProduct(ACone, weight, noUpdateVar, shell, updateVar->matElem, r);
            cg->Mvec(cg->MMat, x, r);
            axpy(&nr, &minusOne, b, &incx, r, &incx);
            scal(&nr, &minusOne, r, &incx);
            // p = r;
            // ASDP_MEMCPY(p, r, double, nr);
            scal(&nr, &zero, p, &incx);
            axpy(&nr, &one, r, &incx, p, &incx);
            // q = r;
            // ASDP_MEMCPY(q, r, double, nr);
            scal(&nr, &zero, q, &incx);
            axpy(&nr, &one, r, &incx, q, &incx);
            // Initialize qTr
            qTr = dot(&nr, q, &incx, r, &incx);
        }
        // qNew = r
        // ASDP_MEMCPY(qNew, r, double, nr);
        scal(&nr, &zero, qNew, &incx);
        axpy(&nr, &one, r, &incx, qNew, &incx);
        
        qTrNew = dot(&nr, qNew, &incx, r, &incx);
        beta = qTrNew / qTr;
        // p = rNew + beta * p;
        scal(&nr, &beta, p, &incx);
        axpy(&nr, &one, r, &incx, p, &incx);

        // update qTr
        qTr = qTrNew;
        // update q
        // ASDP_MEMCPY(q, qNew, double, nr);
        scal(&nr, &zero, q, &incx);
        axpy(&nr, &one, qNew, &incx, q, &incx);

        if ( resiNorm != resiNorm ) {
            cg->solStatus = CG_ITER_STATUS_NUMERICAL;
            retcode = ASDP_RETCODE_FAILED;
            goto exit_cleanup;
        }
    }

exit_cleanup:
    /* Collect solution time statistics */
    cg->cgDuration = AUtilGetTimeStamp() - cgStartTime;
    cg->solveTime += cg->cgDuration;
    return retcode;
}

extern void CGSolverClear(void *pCGSolver){
    if ( !pCGSolver ) {
        return;
    }
    asdp_cg_linsys *cg = (asdp_cg_linsys *) pCGSolver;
    ASDP_FREE(cg->rIter);
    ASDP_FREE(cg->rIterNew);
    ASDP_FREE(cg->pIter);
    ASDP_FREE(cg->qIter);
    ASDP_FREE(cg->qIterNew);
    ASDP_FREE(cg->QIter);
    ASDP_FREE(cg->JacobiPrecond);
    ASDP_FREE(cg);
}

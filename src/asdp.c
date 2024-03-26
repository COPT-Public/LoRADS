#ifdef HEADERPATH
#include "interface/def_asdp.h"
#include "interface/asdp.h"
#include "interface/asdp_conic.h"
#include "interface/asdp_utils.h"
#include "interface/asdp_user_data.h"
#include "interface/asdp_algo.h"
#include "src/asdp_sdpdata.h"
#include "src/asdp_debug.h"
#include "src/vec_opts.h"
#else
#include "def_asdp.h"
#include <stdio.h>
#include "asdp.h"
#include "asdp_conic.h"
#include "asdp_utils.h"
#include "def_asdp_user_data.h"
#include "asdp_user_data.h"
#include "asdp_algo.h"
#include "asdp_sdpdata.h"
#include "asdp_debug.h"
#include "vec_opts.h"
#include "dense_opts.h"
#include "def_asdp_func_set.h"
#endif

#include <math.h>

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

extern void ASDPPrintLog(asdp *ASolver, int Iter, double *dErr, double time)
{
    asdp_printf("Iter:%d objVal:%5.5e dualObj:%5.5e ConstrVio(1):%5.5e ConstrVio(Inf):%5.5e PDGap:%5.5e rho:%3.2f cgIter:%d trace:%3.2f Time:%3.2f\n", Iter, ASolver->pObjVal, ASolver->dObjVal, dErr[ASDP_DIMAC_ERROR_CONSTRVIO], dErr[ASDP_DIMAC_ERROR_CONSTRVIO] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf), dErr[ASDP_DIMAC_ERROR_PDGAP], ASolver->rho, ASolver->cgIter, ASolver->traceSum, time);
}

extern void BMPrintLog(asdp *ASolver, int Iter, int minIter, double *dErr, double time)
{
    asdp_printf("Iter:%d objVal:%5.5e dualObj:%5.5e ConstrVio(1):%5.5e ConstrVio(Inf):%5.5e PDGap:%5.5e rho:%3.2f minIter:%d trace:%3.2f Time:%3.2f\n", Iter, ASolver->pObjVal, ASolver->dObjVal, dErr[ASDP_DIMAC_ERROR_CONSTRVIO], dErr[ASDP_DIMAC_ERROR_CONSTRVIO] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf), dErr[ASDP_DIMAC_ERROR_PDGAP], ASolver->rho, minIter, ASolver->traceSum, time);
}

static void ASDPIGetStatistics(asdp *ASolver)
{
    /* Collect statistics */
    asdp_printf("  Collecting statistcs \n");

    /* Get sum of conic dimensions. Used to convert the identity matrix to Frobenius norm */
    int sumConeDims = 0;
    int maxConeDim = 0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        maxConeDim = ASDP_MAX(maxConeDim, AConeGetDim(ASolver->ACones[iCone]));
        sumConeDims += AConeGetDim(ASolver->ACones[iCone]);
    }

    /* Get norms */
    double objOneNorm = 0.0;
    double objFroNorm = 0.0;
    double dataOneNorm = 0.0;
    double dataFroNorm = 0.0;
    double dFroNormTmp = 0.0;
    double rhsOneNorm = 0.0;
    double rhsFroNorm = 0.0;
    double rhsInfNorm = 0.0;
    return;
}

extern void detectMaxCutProb(asdp *ASolver, int *BlkDims, int *maxCut)
{
    maxCut[0] = -1;
    int stat = 0;
    asdp_cone *cone = ASolver->ACones[0];

    if (cone->cone == ASDP_CONETYPE_DENSE_SDP)
    {
        asdp_cone_sdp_dense *denseCone = (asdp_cone_sdp_dense *)cone->coneData;
        //        if (ASolver->nRows != denseCone->nRow){
        //            return;
        //        }
        for (int iRow = 0; iRow < denseCone->nRow; ++iRow)
        {
            sdp_coeff *sdpCoeff = denseCone->sdpRow[iRow];
            if (sdpCoeff->dataType == SDP_COEFF_SPARSE)
            {
                sdp_coeff_sparse *sdpSparse = sdpCoeff->dataMat;
                if (sdpSparse->nTriMatElem == 1)
                {
                    if (sdpSparse->triMatRow[0] == sdpSparse->triMatCol[0])
                    {
                        stat += 1;
                    }
                }
            }
        }
    }
    if (stat >= ASolver->nRows - 1)
    {
        maxCut[0] = 1;
    }
}

extern void detectSparsitySDPCoeff(asdp *ASolver)
{
    ASDP_INIT(ASolver->sparsitySDPCoeff, double, ASolver->nCones * ASolver->nRows);
    ASolver->nnzSDPCoeffSum = 0;
    ASolver->SDPCoeffSum = 0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_cone *ACone = ASolver->ACones[iCone];
        ACone->nnzStatCoeff(ACone->coneData, &ASolver->sparsitySDPCoeff[iCone * ASolver->nRows], &ASolver->nnzSDPCoeffSum, &ASolver->SDPCoeffSum);
    }
    ASolver->overallSparse = (double)ASolver->nnzSDPCoeffSum / (double)ASolver->SDPCoeffSum / (double)ASolver->nRows;
}

extern void freeDetectSparsitySDPCoeff(asdp *ASolver)
{
    ASDP_FREE(ASolver->sparsitySDPCoeff);
}

extern asdp_retcode ASDPPreprocess(asdp *ASolver, int *BlkDims)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;

    /* Start optimization */
    asdp_printf("\nASDP: software for semi-definite programming \n\n");
    asdp_printf("---------------------------------------------\n");
    ASolver->dTimeBegin = AUtilGetTimeStamp();

    /* Process conic data */
    asdp_printf("Pre-solver starts \n");
    asdp_printf("  Processing the cones \n");
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASDP_CALL(AConeProcData(ASolver->ACones[iCone]));                     // set data
        ASDP_CALL(AConePresolveData(ASolver->ACones[iCone], BlkDims[iCone])); // detect rank structure
    }

#ifdef ASDP_CONIC_DEBUG
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_printf("view:%d", iCone);
        AConeView(ASolver->ACones[iCone]);
    }
#endif

    asdp_printf("  End preprocess \n");

    // Detect `sum of sdp coeff` and `sum of sdp coeff and obj` strcuture
    ASDPSumSDPData(ASolver);
    ASDPDetectSparsityOfSumSDP(ASolver);
    // Calculate obj nrm1 of obj
    ASDPNrm1Obj(ASolver);
    ASDPNrm2Obj(ASolver);

    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        if (ASolver->ACones[iCone]->cone == ASDP_CONETYPE_DENSE_SDP)
        {
            constrValDense *dense;
            ASDP_INIT(dense, constrValDense, 1);
            dense->nnz = ASolver->nRows;
            ASDP_INIT(dense->val, double, ASolver->nRows);
            ASolver->constrVal[iCone]->constrVal = (void *)dense;
            ASolver->constrVal[iCone]->add = addDense;
            ASolver->constrVal[iCone]->zero = zeroDense;
            ASolver->constrVal[iCone]->type = ASDP_DENSE;
        }
        else if (ASolver->ACones[iCone]->cone == ASDP_CONETYPE_SPARSE_SDP)
        {
            constrValSparse *sparse;
            ASDP_INIT(sparse, constrValSparse, 1);
            asdp_cone_sdp_sparse *cone = (asdp_cone_sdp_sparse *)ASolver->ACones[iCone]->coneData;
            ASDP_INIT(sparse->val, double, cone->nRowElem);
            sparse->nnz = cone->nRowElem;
            ASDP_INIT(sparse->nnzIdx, int, cone->nRowElem);
            ASDP_MEMCPY(sparse->nnzIdx, cone->rowIdx, int, cone->nRowElem);
            ASolver->constrVal[iCone]->constrVal = (void *)sparse;
            ASolver->constrVal[iCone]->add = addSparse;
            ASolver->constrVal[iCone]->zero = zeroSparse;
            ASolver->constrVal[iCone]->type = ASDP_SPARSE;
        }
    }

    if (ASolver->nLpCols > 0)
    {
        asdpLPCone *cone = (asdpLPCone *)ASolver->lpCone;
        asdp_cone_lp *lpCone = (asdp_cone_lp *)cone->coneData;
        for (int iCol = 0; iCol < ASolver->nLpCols; ++iCol)
        {
            if (lpCone->lpCol[iCol]->dataType == LP_COEFF_DENSE)
            {
                constrValDense *dense;
                ASDP_INIT(dense, constrValDense, 1);
                dense->nnz = ASolver->nRows;
                ASDP_INIT(dense->val, double, ASolver->nRows);
                ASolver->constrValLP[iCol]->constrVal = (void *)dense;
                ASolver->constrValLP[iCol]->add = addDense;
                ASolver->constrValLP[iCol]->zero = zeroDense;
                ASolver->constrValLP[iCol]->type = ASDP_DENSE;
            }
            else if (lpCone->lpCol[iCol]->dataType == LP_COEFF_SPARSE)
            {
                constrValSparse *sparse;
                ASDP_INIT(sparse, constrValSparse, 1);
                lp_coeff_sparse *lpCol = (lp_coeff_sparse *)lpCone->lpCol[iCol]->dataMat;
                ASDP_INIT(sparse->val, double, lpCol->nnz);
                sparse->nnz = lpCol->nnz;
                ASDP_INIT(sparse->nnzIdx, int, lpCol->nnz);
                ASDP_MEMCPY(sparse->nnzIdx, lpCol->rowPtr, int, lpCol->nnz);
                ASolver->constrValLP[iCol]->constrVal = (void *)sparse;
                ASolver->constrValLP[iCol]->add = addSparse;
                ASolver->constrValLP[iCol]->zero = zeroSparse;
                ASolver->constrValLP[iCol]->type = ASDP_SPARSE;
            }
            else if (lpCone->lpCol[iCol]->dataType == LP_COEFF_ZERO)
            {
                ASDP_ERROR_TRACE;
            }
        }
    }

exit_cleanup:
    return retcode;
}

extern void destroyPreprocess(asdp *ASolver)
{
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        AConeDestroyPresolveData(ASolver->ACones[iCone]);
        AConeDestroyProcData(ASolver->ACones[iCone]);
    }
    ASDP_FREE(ASolver->constrVio);
}

extern void ASDP_BMtoADMM(asdp *ASolver, double heuristic)
{
    ASolver->rho *= heuristic;
    // copy V -> U
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_rk_mat_dense *U = ASolver->U[iCone];
        asdp_rk_mat_dense *V = ASolver->V[iCone];
        asdp_rk_mat_dense *Vlag = ASolver->Vlag[iCone];
        ASDP_MEMCPY(U->matElem, V->matElem, double, V->nRows * V->rank);
        int n = U->nRows * U->rank;
#ifndef NO_PENALTY_METHOD
        ASDP_ZERO(Vlag->matElem, double, n);
#endif
    }
    // double sa = 0.5;
    // int incx = 1;
    // scal(&ASolver->nRows, &sa, ASolver->dualVar, &incx);
    if (ASolver->nLpCols > 0)
    {
        asdp_rk_mat_lp *u = ASolver->uLp;
        asdp_rk_mat_lp *v = ASolver->vLp;
        asdp_rk_mat_lp *vlag = ASolver->vlagLp;
        ASDP_MEMCPY(u->matElem, v->matElem, double, v->nLPCols);
#ifndef NO_PENALTY_METHOD
        ASDP_ZERO(vlag->matElem, double, v->nLPCols);
#endif
    }
}

extern void BMWarmStart2ADMM(asdp *ASolver, double *R, double *dualVar, double rho)
{
    int nElem = 0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_rk_mat_dense *U = ASolver->U[iCone];
        asdp_rk_mat_dense *V = ASolver->V[iCone];
        asdp_rk_mat_dense *Vlag = ASolver->Vlag[iCone];
        nElem = V->nRows * V->rank;
        ASDP_MEMCPY(U->matElem, &R[iCone * nElem], double, nElem);
        ASDP_MEMCPY(V->matElem, &R[iCone * nElem], double, nElem);

        for (int i = 0; i < nElem; ++i)
        {
            Vlag->matElem[i] = normalRandom();
        }
        int incx = 1;
        double froNorm = nrm2(&nElem, Vlag->matElem, &incx);
        rscl(&nElem, &froNorm, Vlag->matElem, &incx);
    }
    ASDP_MEMCPY(ASolver->dualVar, dualVar, double, ASolver->nRows);
    ASolver->rho = rho;
}

// extern asdp_retcode endBMwarmStart(asdp *ASolver)
//{
//     asdp_retcode retcode = ASDP_RETCODE_OK;
//     ASDPCalObj(ASolver, FLAG_V);
//     ASDPCalDualObj(ASolver);
//     ASDPInitConstrValAll(ASolver, ASolver->uLp, ASolver->vLp, ASolver->U, ASolver->V);
//     ASDPInitConstrValSum(ASolver);
//     ASDPUpdateDimacError(ASolver);
//     asdp_printf("**Complete ALM+BM warm start\n");
//     asdp_printf("objVal:%5.5e dualObj:%5.5e ConstrVio:%5.5e Assym:%5.5e DualInfe:%5.5e PDGap:%5.5e rho:%3.2f\n", ASolver->pObjVal, ASolver->dObjVal, ASolver->dimacError[ASDP_DIMAC_ERROR_CONSTRVIO], 0.0, ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE], ASolver->dimacError[ASDP_DIMAC_ERROR_PDGAP], ASolver->rho);
//     ASDPCheckDimacErr(ASolver);
//     ASDP_STATUS_CHECK(ASDPCheckSolverStatus(ASolver));
//
// exit_cleanup:
//     return retcode;
// }

extern void ASDP_SCALE(asdp *ASolver)
{
    // Calculate obj nrmInf of obj
    ASDPNrmInfObj(ASolver);
    int incx = 1;
    ASolver->bRHSNrm1 = nrm1(&ASolver->nRows, ASolver->rowRHS, &incx);
    ASolver->bRHSNrmInf = fabs(ASolver->rowRHS[idamax(&ASolver->nRows, ASolver->rowRHS, &incx)]);
    ASolver->bRHSNrm2 = nrm2(&ASolver->nRows, ASolver->rowRHS, &incx);
    // scale b
    if (ASolver->bRHSNrmInf > 1e+5 || ASolver->cObjNrmInf > 1e+5)
    {
        double bScaleFactor = 1 / ASDP_MAX(ASolver->bRHSNrmInf, 1);
        double cScaleFactor = 1 / ASDP_MAX(ASolver->cObjNrmInf, 1);
        if (ASolver->nLpCols > 0)
        {
            ASolver->lpCone->scaleData(ASolver->lpCone->coneData, cScaleFactor, bScaleFactor);
        }

        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            ASolver->ACones[iCone]->dataScale(ASolver->ACones[iCone]->coneData, cScaleFactor, bScaleFactor);
        }

        vvscl(&ASolver->nRows, &bScaleFactor, ASolver->rowRHS);
        ASolver->cScaleFactor = cScaleFactor;
        ASolver->bScaleFactor = bScaleFactor;
    }
    else
    {
        ASolver->cScaleFactor = 1.0;
        ASolver->bScaleFactor = 1.0;
    }
}

extern asdp_retcode ASDP_BMOptimize(asdp *ASolver, double endBMTol, double endBMTol_pd, double endTauTol, double endBMALSub, double ori_start, int is_rank_max, int *pre_mainiter, int *pre_miniter, double timeLimit)
{
    asdp_func *aFunc;
    ASDPInitFuncSet(&aFunc, ASolver->nLpCols);
    asdp_retcode retcode = ASDP_RETCODE_OK;
    ASolver->whichMethod = BMMethod;
    double start = AUtilGetTimeStamp();
    // Calculate obj nrmInf of obj
    ASDPNrmInfObj(ASolver);

    int nConstr = ASolver->nRows;
    int incx = 1;
    ASolver->bRHSNrm1 = nrm1(&nConstr, ASolver->rowRHS, &incx);
    ASolver->bRHSNrmInf = fabs(ASolver->rowRHS[idamax(&nConstr, ASolver->rowRHS, &incx)]);

    double rho_certicate = 1.0e-1;
    double rho_certicate_tol = rho_certicate / ASolver->rho;
    double rho_creticate_val;
    double constrVioTol = 1e-5;
    double minusOne = -1.0;
    double tau = 0.0;
    double lagNormSquare = 0.0;
    double lagNorm2 = 0.0;
    int minIter = 0 + *pre_miniter;
    double bestInfe = 1e+30;
    aFunc->InitConstrValAll(ASolver, ASolver->rLp, ASolver->rLp, ASolver->R, ASolver->R);
    aFunc->InitConstrValSum(ASolver);

    aFunc->BMCalGrad(ASolver, ASolver->rLp, ASolver->gradLp, ASolver->R, ASolver->Grad, &lagNormSquare);
    rho_creticate_val = sqrt(lagNormSquare) / (1 + ASolver->cObjNrmInf);
    char difficulty = HARD;
    int localIter = 0;
    double constrVioCrite = 0.0;
    int clearLBFGS = 0;
    double innerProduct = 0.0;
    int rank_flag = 0;
    int mainiter_start = *pre_mainiter;

    for (int k = 0 + *pre_mainiter; k < ASolver->maxBMOutIter; ++k)
    {
        int rhoInnerIter = 0;
        double alpha = 0.1;      // The EMA smoothing factor
        double threshold = 0.05; // The threshold for significant proportion change
        double current_ema = 0;  // Initialize the EMA value
        double old_ema = 0;      // Initialize the old EMA value
        int update_interval = 5; // Specify the update interval
        int counter = 1;         // Initialize the counter
        *pre_mainiter = k;
        while (difficulty != EASY)
        {
            localIter = 0;
            if ((!AUtilUpdateCheckEma(&current_ema, &old_ema, rho_creticate_val, alpha, threshold, update_interval, &counter) || rank_flag >= 15) && !is_rank_max)
            {
                break;
                // goto UpdateRho;
            }
            if (rho_creticate_val <= rho_certicate_tol)
            {
                break;
            }

            // Solving AL subproblem
            while (rho_creticate_val - rho_certicate_tol > endBMALSub)
            {
                if (AUtilGetTimeStamp() - ori_start >= timeLimit)
                {
                    goto END_BM_WarmStart;
                }
                if (localIter % 300 == 0)
                {
                    clearLBFGS = 0;
                }
                aFunc->LBFGSDirection(ASolver, ASolver->lbfgsHis, ASolver->gradLp, ASolver->uLp, ASolver->Grad, ASolver->U, clearLBFGS);
                aFunc->LBFGSDirUseGrad(ASolver, ASolver->uLp, ASolver->gradLp, ASolver->U, ASolver->Grad);
                double *q0 = ASolver->M1temp;
                // q0 = b - A(RRt)
                ASDP_MEMCPY(q0, ASolver->rowRHS, double, ASolver->nRows);
                axpy(&(ASolver->nRows), &minusOne, ASolver->constrValSum, &incx, q0, &incx);
                constrVioCrite = nrm2(&(ASolver->nRows), q0, &incx) / (1 + ASolver->bRHSNrmInf);
                if (bestInfe > constrVioCrite)
                {
                    bestInfe = constrVioCrite;
                    aFunc->copyRtoV(ASolver->rLp, ASolver->vLp, ASolver->R, ASolver->V, ASolver->nCones);
                    ASDP_MEMCPY(ASolver->bestDualVar, ASolver->dualVar, double, ASolver->nRows);
                }
                if ((bestInfe <= endBMTol || ASolver->dimacError[ASDP_DIMAC_ERROR_PDGAP] <= endBMTol_pd) && k >= 1)
                {
                    goto END_BM_WarmStart;
                }

                double p12[2];
                aFunc->BMCalq12p12(ASolver, ASolver->rLp, ASolver->uLp, ASolver->R, ASolver->U, ASolver->ARDSum, ASolver->ADDSum, p12);
                int rootNum = BMLineSearch(ASolver->rho, ASolver->nRows, ASolver->dualVar, p12[0], p12[1], q0, ASolver->ARDSum, ASolver->ADDSum, &tau);
                if (fabs(tau) < endTauTol)
                {
                    printf("update rho.\n");
                    goto UpdateRho;
                }
                if (rootNum == 0)
                {
                    asdp_printf("*Numerical Fail in BM minter %d, we use the best feasible solution as warm start\n", minIter);
                    goto END_BM_WarmStart;
                }
                // y = (-grad)
                // set lbfgs one
                aFunc->setAsNegGrad(ASolver, ASolver->gradLp, ASolver->Grad);
                // update R
                aFunc->BMupdateVar(ASolver, ASolver->rLp, ASolver->uLp, ASolver->R, ASolver->U, tau);
                // update gradient
                lagNormSquare = 0.0;
                // update constrValSum first
                double tauSquare = tau * tau;
                axpy(&(ASolver->nRows), &tau, ASolver->ARDSum, &incx, ASolver->constrValSum, &incx);
                axpy(&(ASolver->nRows), &(tauSquare), ASolver->ADDSum, &incx, ASolver->constrValSum, &incx);
                aFunc->BMCalGrad(ASolver, ASolver->rLp, ASolver->gradLp, ASolver->R, ASolver->Grad, &lagNormSquare);
                rho_creticate_val = sqrt(lagNormSquare) / (1 + ASolver->cObjNrmInf);
                aFunc->setlbfgsHisTwo(ASolver, ASolver->gradLp, ASolver->uLp, ASolver->Grad, ASolver->U, tau);
                minIter++;
                *pre_miniter = minIter;
                localIter++;
                clearLBFGS++;
                if (localIter > 800)
                {
                    break;
                }
            }
            ASDPUpdateDualVar(ASolver);
            aFunc->BMCalGrad(ASolver, ASolver->rLp, ASolver->gradLp, ASolver->R, ASolver->Grad, &lagNormSquare);
            rho_creticate_val = sqrt(lagNormSquare) / (1 + ASolver->cObjNrmInf);
            if (localIter <= 20)
            {
                difficulty = EASY;
            }
            else if (20 < localIter && localIter <= 100)
            {
                difficulty = MEDIUM;
                rank_flag += 2;
            }
            else if (100 < localIter)
            {
                difficulty = HARD;
                rank_flag += 3;
            }
            else if (400 < localIter)
            {
                difficulty = SUPER;
                rank_flag += 4;
            }

            rhoInnerIter++;
        }
        if (difficulty == EASY)
        {
            rank_flag = 0;
        }

    UpdateRho:
        do
        {
            ASolver->rho *= 2;
            aFunc->BMCalGrad(ASolver, ASolver->rLp, ASolver->gradLp, ASolver->R, ASolver->Grad, &lagNormSquare);
            rho_creticate_val = sqrt(lagNormSquare) / (1 + ASolver->cObjNrmInf);
            rho_certicate_tol = rho_certicate / ASolver->rho;
        } while (rho_certicate_tol >= rho_creticate_val);
    END_CURRENT_MAJOR_ITER:

        // Refresh some parameter
        difficulty = HARD;
        clearLBFGS = 0;
        if (k % 1 == 0)
        {
            aFunc->calObj(ASolver, FLAG_R);
            ASDPCalDualObj(ASolver);
            aFunc->updateDimac(ASolver);
            aFunc->copyRtoV(ASolver->rLp, ASolver->lagLp, ASolver->R, ASolver->lag, ASolver->nCones);
            ASDPNuclearNorm(ASolver);
            BMPrintLog(ASolver, k, minIter, ASolver->dimacError, AUtilGetTimeStamp() - ori_start);
#ifdef TERMINATION_SCHEME_OFFICIAL
            ASDPCheckDimacErr(ASolver);
#else
            ASDPCheckDimacErrBMCriteria(ASolver);
#endif
        }
#ifdef TERMINATION_SCHEME_OFFICIAL
        ASDP_STATUS_CHECK(ASDPCheckSolverStatus(ASolver));
#else
        ASDP_STATUS_CHECK(ASDPCheckSolverStatusBM(ASolver));
#endif

        if (rank_flag >= 15 && !is_rank_max)
        {
            rank_flag = 0;
            if (k - mainiter_start >= 2)
            {
                retcode = ASDP_RETCODE_RANK;
                *pre_mainiter += 1;
                goto exit_cleanup;
            }
        }
    }

END_BM_WarmStart:
    ASDP_MEMCPY(ASolver->dualVar, ASolver->bestDualVar, double, ASolver->nRows);
    aFunc->copyRtoV(ASolver->vLp, ASolver->rLp, ASolver->V, ASolver->R, ASolver->nCones);
    aFunc->copyRtoV(ASolver->vLp, ASolver->uLp, ASolver->V, ASolver->U, ASolver->nCones);
    aFunc->calObj(ASolver, FLAG_V);
    ASDPCalDualObj(ASolver);
    aFunc->updateDimac(ASolver);
    ASDPNuclearNorm(ASolver);
    asdp_printf("**Complete ALM+BM warm start\n");
    asdp_printf("objVal:%5.5e dualObj:%5.5e ConstrVio:%5.5e Assym:%5.5e DualInfe:%5.5e PDGap:%5.5e rho:%3.2f minIter:%d Time:%3.2f\n", ASolver->pObjVal, ASolver->dObjVal, ASolver->dimacError[ASDP_DIMAC_ERROR_CONSTRVIO], 0.0, ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE], ASolver->dimacError[ASDP_DIMAC_ERROR_PDGAP], ASolver->rho, minIter, AUtilGetTimeStamp() - ori_start);
#ifdef TERMINATION_SCHEME_OFFICIAL
    ASDPCheckDimacErr(ASolver);
    ASDP_STATUS_CHECK(ASDPCheckSolverStatus(ASolver));
#else
    ASDPCheckDimacErrBMCriteria(ASolver);
    ASDP_STATUS_CHECK(ASDPCheckSolverStatusBM(ASolver));
#endif
exit_cleanup:
    // printf("%d\n", retcode);
    return retcode;
}

extern asdp_retcode ASDPOptimize(asdp *ASolver, int rhoFreq, double rhoFactor, int rhoStrategy, double tau, double gamma, double rhoMin, double ori_start, double timeLimit)
{
    asdp_func *aFunc;
    ASDPInitFuncSet(&aFunc, ASolver->nLpCols);
    asdp_retcode retcode = ASDP_RETCODE_OK;
    ASolver->whichMethod = ADMMMethod;
    ASolver->cgIter = 0;
    asdp_printf("**Change method into ADMM Split method\n");
    // Initialize info
    // Initialize Constraint function value <A, UV^T>
    ASDPInitConstrValAll(ASolver, ASolver->uLp, ASolver->vLp, ASolver->U, ASolver->V);
    ASDPInitConstrValSum(ASolver);
    // Initialize Objective function value
    aFunc->calObj(ASolver, FLAG_U);
    ASolver->cgTime = 0.0;
    double start = AUtilGetTimeStamp();
    int k;
    double rho_0 = ASolver->rho;
    if (rho_0 >= ASolver->rhoMax)
    {
        rho_0 = ASolver->rhoMax * 0.75;
    }
    double iter_cross = rhoFreq * ASDP_MAX(ceil(log(ASolver->rhoMax / rho_0) / log(rhoFactor)), 1);
    double c = ASDP_MAX((ASolver->rhoMax - rho_0), 0) / log(iter_cross + 1);
    for (k = 0; k < ASolver->maxInter; ++k)
    {
        if (AUtilGetTimeStamp() - ori_start >= timeLimit)
        {
            goto exit_cleanup;
        }
        aFunc->admmUpdateVar(ASolver);
        // update dual var lambda
        ASDPUpdateDualVar(ASolver);
        if (k % rhoFreq == 0)
        {
            if (k >= iter_cross && rhoStrategy == 1)
            {
                if (c == 0)
                {
                    ASolver->rho = k == 0 ? ASolver->rhoMax : c * log(k + 1) + ASolver->rhoMax;
                }
                else
                {
                    ASolver->rho = k == 0 ? rho_0 : c * log(k + 1) + rho_0;
                }
            }
            else
            {
                ASolver->rho *= rhoFactor;
            }
        }
        if (k % 1 == 0)
        {
            double gradNorm = 0.0;
            aFunc->calObj(ASolver, FLAG_U);
            ASDPCalDualObj(ASolver);
            aFunc->updateDimac(ASolver);
            aFunc->copyRtoV(ASolver->rLp, ASolver->lagLp, ASolver->R, ASolver->lag, ASolver->nCones);
            ASDPNuclearNorm(ASolver);
            aFunc->BMCalGrad(ASolver, ASolver->rLp, ASolver->gradLp, ASolver->R, ASolver->Grad, &gradNorm);
            if (rhoStrategy == STRATEGY_GAMMA_DYNAMIC)
            {
                if (ASolver->dimacError[ASDP_DIMAC_ERROR_CONSTRVIO] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrm2) > tau * gradNorm)
                {
                    ASolver->rho *= gamma;
                }
                else
                {
                    ASolver->rho = ASDP_MAX(ASolver->rho / gamma, rhoMin);
                }
            }
            ASDPPrintLog(ASolver, k, ASolver->dimacError, AUtilGetTimeStamp() - start);
#ifdef TERMINATION_SCHEME_OFFICIAL
            ASDPCheckDimacErr(ASolver);
#endif
        }
#ifndef TERMINATION_SCHEME_OFFICIAL
        int nConstr = ASolver->nRows;
        aFunc->InitConstrValAll(ASolver, ASolver->rLp, ASolver->rLp, ASolver->R, ASolver->R);
        aFunc->InitConstrValSum(ASolver);
        double *constrVio = ASolver->constrVio;
        double one = 1.0;
        double minusOne = -1.0;
        axpbyAddition(&nConstr, &(one), ASolver->rowRHS, &(minusOne), ASolver->constrValSum, constrVio);
        int incx = 1;
        ASolver->dimacError[ASDP_DIMAC_ERROR_CONSTRVIO] = nrm2(&nConstr, constrVio, &incx) / (1 + ASolver->bRHSNrm1);
        ASDPCheckDimacErrBMCriteria(ASolver);
#endif
#ifdef TERMINATION_SCHEME_OFFICIAL
        ASDP_STATUS_CHECK(ASDPCheckSolverStatus(ASolver));
#else
        ASDP_STATUS_CHECK(ASDPCheckSolverStatusBM(ASolver));
#endif
    }
exit_cleanup:
    if (k == ASolver->maxInter)
    {
        ASolver->AStatus = ASDP_MAXITER;
        retcode = ASDP_RETCODE_EXIT;
    }
    else if (AUtilGetTimeStamp() - ori_start >= timeLimit)
    {
        ASolver->AStatus = ASDP_TIMELIMIT;
        retcode = ASDP_RETCODE_EXIT;
    }
    ASDP_FREE(aFunc);
    return retcode;
}


extern void ASDPCheckDimacErr(asdp *ASolver)
{

    double dMaxDimacsErr = 0.0;
    double *dErrs = ASolver->dimacError;
    for (int iElem = 0; iElem < 4; ++iElem)
    {
        if (iElem != ASDP_DIMAC_ERROR_ASSYMMETRY)
        {
            // Since the objective val is calculated by average of U and V
            dMaxDimacsErr = ASDP_MAX(dMaxDimacsErr, dErrs[iElem]);
        }
    }
    if (ASolver->strategy == STRATEGY_MIN_BISECTION)
    {
        if (dErrs[ASDP_DIMAC_ERROR_CONSTRVIO] < 1e-5 && dErrs[ASDP_DIMAC_ERROR_PDGAP] < 1e-3 && dErrs[ASDP_DIMAC_ERROR_DUALFEASIBLE] / dErrs[ASDP_DIMAC_ERROR_CONSTRVIO] > 10)
        {
            ASolver->rho /= 2;
        }
    }

    if (dMaxDimacsErr < ASDP_TERMINATION_TOL)
    {
        ASolver->AStatus = ASDP_PRIMAL_DUAL_OPTIMAL;
    }
    return;
}

extern void ASDPCheckDimacErrBMCriteria(asdp *ASolver)
{
    double *dErrs = ASolver->dimacError;
    double BMCriter = dErrs[ASDP_DIMAC_ERROR_CONSTRVIO] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf);
    //    asdp_printf("BMCriter: %5.2e \n", BMCriter);
    //    printf("diff gap:%f\n",  ASolver->dimacError[ASDP_ERROR_DIFF_GAP]);
    if (BMCriter < ASDP_TERMINATION_TOL)
    {
        ASolver->AStatus = ASDP_PRIMAL_OPTIMAL;
    }
    return;
}

static void printRes(double pObj, double dObj, double constrVio, double dualInfe, double pdgap, double constrVioInf, double dualInfeInf)
{
    asdp_printf("-----------------------------------------------------------------------\n");
    asdp_printf("Objective function Value are:\n");
    asdp_printf("\t 1.Primal Objective:            : %5.2e\n", pObj);
    asdp_printf("\t 2.Dual Objective:              : %5.2e\n", dObj);
    asdp_printf("Dimacs Error are:\n");
    asdp_printf("\t 1.Constraint Violation(1)      : %5.2e\n", constrVio);
    asdp_printf("\t 2.Dual Infeasibility(1)        : %5.2e\n", dualInfe);
    asdp_printf("\t 3.Primal Dual Gap              : %5.2e\n", pdgap);
    asdp_printf("\t 4.Primal Variable Semidefinite : %5.2e\n", 0.0);
    asdp_printf("\t 5.Constraint Violation(Inf)    : %5.2e\n", constrVioInf);
    asdp_printf("\t 6.Dual Infeasibility(Inf)      : %5.2e\n", dualInfeInf);
    asdp_printf("-----------------------------------------------------------------------\n");
}

extern asdp_retcode ASDPCheckSolverStatusBM(asdp *ASolver)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    if (ASolver->AStatus == ASDP_PRIMAL_OPTIMAL)
    {
        retcode = ASDP_RETCODE_EXIT;
    }
    return retcode;
}

extern void ASDPEndProgram(asdp *ASolver)
{
    asdp_printf("-----------------------------------------------------------------------\n");
    if (ASolver->AStatus == ASDP_MAXITER)
    {
        asdp_printf("End Program due to reaching `the maximum number of iterations`:\n");
    }
    else if (ASolver->AStatus == ASDP_PRIMAL_DUAL_OPTIMAL)
    {
        asdp_printf("End Program due to reaching `Official terminate criteria`:\n");
    }
    else if (ASolver->AStatus == ASDP_PRIMAL_OPTIMAL)
    {
        asdp_printf("End Program due to reaching `final terminate criteria`:\n");
    }
    else if (ASolver->AStatus == ASDP_TIMELIMIT)
    {
        asdp_printf("End Program due to reaching `Time limit`\n");
    }
    else if (ASolver->AStatus == ASDP_UNKNOWN)
    {
        asdp_printf("End Program but the status is unknown, please notify the authors\n");
    }

    double *negLambd = ASolver->negLambd;
    ASDP_ZERO(negLambd, double, ASolver->nRows);
    int incx = 1;
    double minusOne = -1.0;
    axpy(&(ASolver->nRows), &minusOne, ASolver->dualVar, &incx, negLambd, &incx);
    ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] = 0.0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_rk_mat_dense *U = ASolver->U[iCone];
        int nRows = U->nRows;

        asdp_cone *Cone = ASolver->ACones[iCone];
        Cone->sdp_slack_var->zeros(Cone->sdp_slack_var->dataMat);
        Cone->addObjCoeff(Cone->coneData, Cone->sdp_slack_var);
        Cone->sdpDataWSum(Cone->coneData, negLambd, Cone->sdp_slack_var);

        double minEig = 0.0;
        fds_syev_Min_eig(Cone->sdp_slack_var, &minEig);
        double dualVio = ASDP_ABS(ASDP_MIN(minEig, 0));
        ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] += dualVio;
    }
    ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] /= (ASolver->cObjNrm1 + 1);

    printRes(ASolver->pObjVal, ASolver->dObjVal, ASolver->dimacError[ASDP_DIMAC_ERROR_CONSTRVIO], ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE], ASolver->dimacError[ASDP_DIMAC_ERROR_PDGAP], ASolver->dimacError[ASDP_DIMAC_ERROR_CONSTRVIO] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf), ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] * (1 + ASolver->cObjNrm1) / (1 + ASolver->cObjNrmInf));
}

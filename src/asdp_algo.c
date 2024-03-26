#ifdef HEADERPATH
#include "interface/def_asdp.h"
#include "interface/asdp.h"
#include "interface/asdp_utils.h"
#include "interface/asdp_user_data.h"
#include "interface/asdp_file_io.h"
#include "interface/asdp_conic.h"
#include "interface/asdp_schur.h"
#include "interface/asdp_algo.h"
#include "interface/def_asdp_user_data.h"
#else
#include "def_asdp.h"
#include "asdp.h"
#include <stdio.h>
#include "asdp_utils.h"
#include "asdp_user_data.h"
#include "def_asdp_user_data.h"
#include "asdp_file_io.h"
#include "asdp_conic.h"
#include "asdp_algo.h"
#include "asdp_sdpdata.h"
#include "asdp_debug.h"
#include "vec_opts.h"
#include "dense_opts.h"
#include "sparse_opts.h"
#include "asdp_direct_linsys.h"
#include "asdp_lbfgs.h"
#include "asdp_chol.h"
#endif

#include <math.h>

#ifdef MEMDEBUG
#include "memwatch.h"
#endif
// #define ASDP_ONE_INIT
// #define ASDP_NORMAL_RANDOM_INIT

/* Define ASDP Solver interface */
extern asdp_retcode ASDPSolverCreate(asdp **pASolver)
{

    asdp_retcode retcode = ASDP_RETCODE_OK;

    if (!pASolver)
    {
        retcode = ASDP_RETCODE_FAILED;
        return retcode;
    }

    asdp *ASolver;
    ASDP_INIT(ASolver, asdp, 1);
    ASDP_MEMCHECK(ASolver);

    *pASolver = ASolver;

    AUtilStartCtrlCCheck();

exit_cleanup:
    return retcode;
}

extern void ASDPDestroy(asdp **pASolver)
{

    if (!pASolver)
    {
        return;
    }

    ASDP_FREE(*pASolver);

    asdp_printf("ASDP ends. Exiting \n");

    return;
}

extern double normalRandom()
{
    double u1 = (double)rand() / RAND_MAX; // uniform [0, 1]
    double u2 = (double)rand() / RAND_MAX;

    // Box-Muller transform
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);

    return z;
}

static asdp_retcode ASDPSetLpCone(asdpLPCone *lpCone, int nRows, int nLpCols, int *LpMatBeg, int *LpMatIdx, double *LpMatElem)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;

    lpCone->nCol = nLpCols;

    LPSetConeFuncPointer(lpCone);

    lpCone->coneCreate(&lpCone->coneData);

    lpCone->coneProcData(lpCone->coneData, nRows, nLpCols, LpMatBeg, LpMatIdx, LpMatElem);

    lpCone->conePresolveData(lpCone->coneData);
    // lpCone->coneView(lpCone->coneData);
exit_cleanup:
    return retcode;
}

extern asdp_retcode ASDPInitConeData(asdp *ASolver, user_data **SDPDatas, double **coneMatElem, int **coneMatBeg, int **coneMatIdx, int *BlkDims, int nConstrs, int nBlks, int nLpCols, int *LpMatBeg, int *LpMatIdx, double *LpMatElem)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;

    if (nLpCols > 0)
    {
        ASDPSetLpCone(ASolver->lpCone, ASolver->nRows, nLpCols, LpMatBeg, LpMatIdx, LpMatElem);
    }
    user_data *SDPData = NULL;
    for (int iCone = 0; iCone < nBlks; ++iCone)
    {
        ASDP_CALL(AUserDataCreate(&SDPData)); // initialize, set all equals to 0
        AUserDataSetConeData(SDPData, ASDP_CONETYPE_DENSE_SDP, nConstrs, BlkDims[iCone],
                             coneMatBeg[iCone], coneMatIdx[iCone], coneMatElem[iCone]);
        ASDP_CALL(ASDPSetCone(ASolver, iCone, SDPData));
        SDPDatas[iCone] = SDPData;
        SDPData = NULL;
    }
exit_cleanup:
    return retcode;
}
extern void ASDPDestroyConeData(asdp *ASolver)
{
    // destroy for function ASDPInitConeData_and_UVS
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        // SDPdata pointer is set as usrData here
        ASDP_FREE(ASolver->ACones[iCone]->usrData);
    }
    if (ASolver->nLpCols > 0)
    {
        ASolver->lpCone->coneDestroyData(&ASolver->lpCone->coneData);
    }
    ASDP_FREE(ASolver->lpCone);
}

extern void lpRandom(double *data, int n)
{
    for (int i = 0; i < n; ++i)
    {
        data[i] = (double)rand() / RAND_MAX;
        data[i] -= (double)rand() / RAND_MAX;
    }
}

extern void lpRandomDiag(double *data, int n, int nRows, int nCols)
{
    if (n == 0)
    {
        return;
    }
    int r = ASDP_MIN(nRows, nCols);
    ASDP_ZERO(data, double, n);
    for (int i = 0; i < r; ++i)
    {
        data[i * nRows + i] = 1 / sqrt(r);
    }
}

extern asdp_retcode ASDPInitBMVars(asdp *ASolver, int *BlkDims, int nBlks, int nLpCols)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    ASolver->hisRecT = 3;
    ASolver->maxBMInIter = 200;
    ASolver->maxBMOutIter = 200;

    ASDP_INIT(ASolver->rLp, asdp_rk_mat_lp, 1);
    ASDP_MEMCHECK(ASolver->rLp);
    ASolver->rLp->nLPCols = nLpCols;
    ASDP_INIT(ASolver->gradLp, asdp_rk_mat_lp, 1);
    ASDP_MEMCHECK(ASolver->gradLp);
    ASolver->gradLp->nLPCols = nLpCols;

    srand(925);
    ASDP_INIT(ASolver->ADDSum, double, ASolver->nRows);
    ASDP_MEMCHECK(ASolver->ADDSum);
    ASDP_INIT(ASolver->ARDSum, double, ASolver->nRows);
    ASDP_MEMCHECK(ASolver->ARDSum);
    ASDP_INIT(ASolver->M1temp, double, ASolver->nRows);
    ASDP_MEMCHECK(ASolver->M1temp);
    // R is variable of BMALM
    // R is also the average variable for ADMM split method
    ASDP_INIT(ASolver->R, asdp_rk_mat_dense *, nBlks);
    ASDP_MEMCHECK(ASolver->R);
    ASDP_INIT(ASolver->Grad, asdp_rk_mat_dense *, nBlks);
    ASDP_MEMCHECK(ASolver->Grad);
    int allElem = 0;
    if (nLpCols > 0)
    {
        allElem += nLpCols;
    }
    for (int iCone = 0; iCone < nBlks; ++iCone)
    {
        asdp_cone *ACone = ASolver->ACones[iCone];
        // set variable R
        asdp_rk_mat_dense *R;
        ASDP_INIT(R, asdp_rk_mat_dense, 1);
        ASDP_MEMCHECK(R);
        R->rank = ASolver->rankElem[iCone];
        R->nRows = BlkDims[iCone];
        ASDP_RANDOM_rk_MAT(R);
        //        ASDP_ONE_rk_MAT(R);
        ASDPSetVarUVS(ASolver, iCone, R, FLAG_R);
        allElem += R->nRows * R->rank;

        asdp_rk_mat_dense *Grad;
        ASDP_INIT(Grad, asdp_rk_mat_dense, 1);
        ASDP_MEMCHECK(Grad);
        Grad->rank = ASolver->rankElem[iCone];
        Grad->nRows = BlkDims[iCone];
        int n = Grad->nRows * Grad->rank;
        ASDP_INIT(Grad->matElem, double, n);
        ASolver->Grad[iCone] = Grad;
    }
    if (nLpCols > 0)
    {
        ASDP_INIT(ASolver->rLp->matElem, double, nLpCols);
        ASDP_MEMCHECK(ASolver->rLp->matElem);
        lpRandom(ASolver->rLp->matElem, nLpCols);
        ASDP_INIT(ASolver->gradLp->matElem, double, nLpCols);
        ASDP_MEMCHECK(ASolver->gradLp->matElem);
    }
    ASDP_INIT(ASolver->lbfgsHis, lbfgs_node, 1);
    ASDP_MEMCHECK(ASolver->lbfgsHis);
    // ALL  (hisRecT + 1) node
    lbfgs_node *head = ASolver->lbfgsHis;
    ASDP_INIT(head->s, double, allElem);
    ASDP_INIT(head->y, double, allElem);
    head->allElem = allElem;
    lbfgs_node *node = head;
    for (int nodeNum = 0; nodeNum < ASolver->hisRecT; ++nodeNum)
    {
        lbfgs_node *nextNode;
        ASDP_INIT(nextNode, lbfgs_node, 1);
        ASDP_INIT(nextNode->s, double, allElem);
        ASDP_MEMCHECK(nextNode->s);
        ASDP_INIT(nextNode->y, double, allElem);
        ASDP_MEMCHECK(nextNode->y);
        nextNode->allElem = allElem;
        node->next = nextNode;
        nextNode->prev = node;
        node = node->next;
        if (nodeNum == ASolver->hisRecT - 1)
        {
            node->next = head;
            head->prev = node;
        }
    }
    ASolver->lbfgsHis = head;

exit_cleanup:
    return retcode;
}

extern void ASDPDestroyBMVars(asdp *ASolver)
{
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASDP_FREE(ASolver->R[iCone]->matElem);
        ASDP_FREE(ASolver->R[iCone]);
        ASDP_FREE(ASolver->Grad[iCone]->matElem);
        ASDP_FREE(ASolver->Grad[iCone]);
    }
    lbfgs_node *node = ASolver->lbfgsHis;
    lbfgs_node *nextNode;
    for (int nodeNum = 0; nodeNum < ASolver->hisRecT; ++nodeNum)
    {
        ASDP_FREE(node->s);
        ASDP_FREE(node->y);
        nextNode = node->next;
        ASDP_FREE(node);
        node = nextNode;
    }
    ASDP_FREE(nextNode->s);
    ASDP_FREE(nextNode->y);
    ASDP_FREE(nextNode);
    ASDP_FREE(ASolver->R);
    ASDP_FREE(ASolver->Grad);
    if (ASolver->nLpCols > 0)
    {
        ASDP_FREE(ASolver->rLp->matElem);
    }
    ASDP_FREE(ASolver->rLp);
    ASDP_FREE(ASolver->gradLp);
    ASDP_FREE(ASolver->M1temp);
    ASDP_FREE(ASolver->ADDSum);
    ASDP_FREE(ASolver->ARDSum);
}

extern asdp_retcode ASDPInitADMMVars(asdp *ASolver, int *BlkDims, int nBlks, int nLpCols)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;

    ASDP_INIT(ASolver->uLp, asdp_rk_mat_lp, 1);
    ASDP_MEMCHECK(ASolver->uLp);
    ASolver->uLp->nLPCols = nLpCols;
    ASDP_INIT(ASolver->vLp, asdp_rk_mat_lp, 1);
    ASDP_MEMCHECK(ASolver->vLp);
    ASolver->vLp->nLPCols = nLpCols;
    ASDP_INIT(ASolver->vlagLp, asdp_rk_mat_lp, 1);
    ASDP_MEMCHECK(ASolver->vlagLp);
    ASDP_INIT(ASolver->lagLp, asdp_rk_mat_lp, 1);
    ASDP_MEMCHECK(ASolver->lagLp);
    ASolver->vlagLp->nLPCols = nLpCols;
    if (nLpCols > 0)
    {
        ASDP_INIT(ASolver->uLp->matElem, double, nLpCols);
        ASDP_MEMCHECK(ASolver->uLp->matElem);
        lpRandom(ASolver->uLp->matElem, nLpCols);

        ASDP_INIT(ASolver->vLp->matElem, double, nLpCols);
        ASDP_MEMCHECK(ASolver->vLp->matElem);
        lpRandom(ASolver->vLp->matElem, nLpCols);

        ASDP_INIT(ASolver->vlagLp->matElem, double, nLpCols);
        ASDP_MEMCHECK(ASolver->vlagLp->matElem);
        
        ASDP_INIT(ASolver->lagLp->matElem, double, nLpCols);
        ASDP_MEMCHECK(ASolver->lagLp->matElem);
        // lpRandom(ASolver->vlagLp->matElem, nLpCols);
        for (int i = 0; i < nLpCols; ++i)
        {
            ASolver->vlagLp->matElem[i] = 1.0;
            ASolver->lagLp->matElem[i] = 1.0;
        }
    }

    ASDP_INIT(ASolver->bLinSys, double *, nBlks);
    ASDP_MEMCHECK(ASolver->bLinSys);
    for (int iCone = 0; iCone < nBlks; ++iCone)
    {
        ASDP_INIT(ASolver->bLinSys[iCone], double, BlkDims[iCone] * ASolver->rankElem[iCone]);
        ASDP_MEMCHECK(ASolver->bLinSys[iCone]);
        ASDP_ZERO(ASolver->bLinSys[iCone], double, BlkDims[iCone] * ASolver->rankElem[iCone]);
    }
    ASDP_INIT(ASolver->V, asdp_rk_mat_dense *, nBlks);
    ASDP_MEMCHECK(ASolver->V);
    ASDP_INIT(ASolver->U, asdp_rk_mat_dense *, nBlks);
    ASDP_MEMCHECK(ASolver->U);
    ASDP_INIT(ASolver->Vlag, asdp_rk_mat_dense *, nBlks);
    ASDP_MEMCHECK(ASolver->Vlag);
    ASDP_INIT(ASolver->lag, asdp_rk_mat_dense *, nBlks);
    ASDP_MEMCHECK(ASolver->lag);

    user_data *SDPData = NULL;
    for (int iCone = 0; iCone < nBlks; ++iCone)
    {
        asdp_cone *ACone = ASolver->ACones[iCone];
        // set variable
        asdp_rk_mat_dense *U;
        ASDP_INIT(U, asdp_rk_mat_dense, 1);
        ASDP_MEMCHECK(U);
        U->rank = ASolver->rankElem[iCone];
        U->nRows = BlkDims[iCone];

        asdp_rk_mat_dense *V;
        ASDP_INIT(V, asdp_rk_mat_dense, 1);
        ASDP_MEMCHECK(V);
        V->rank = ASolver->rankElem[iCone];
        V->nRows = BlkDims[iCone];
        
        asdp_rk_mat_dense *lag;
        ASDP_INIT(lag, asdp_rk_mat_dense, 1);
        ASDP_MEMCHECK(lag);
        lag->rank = ASolver->rankElem[iCone];
        lag->nRows = BlkDims[iCone];
        ASDP_RANDOM_rk_MAT(lag);

        asdp_rk_mat_dense *S;
        ASDP_INIT(S, asdp_rk_mat_dense, 1);
        ASDP_MEMCHECK(S);
        S->rank = ASolver->rankElem[iCone];
        S->nRows = BlkDims[iCone];
        ASDP_RANDOM_rk_MAT(U);
        ASDP_RANDOM_rk_MAT(V);
        int n = S->nRows * S->rank;
        ASDP_INIT(S->matElem, double, n);
        ASDP_MEMCHECK(S->matElem);
        for (int i = 0; i < n; ++i)
        {
            S->matElem[i] = normalRandom();
        }
        int incx = 1;
        double froNorm = nrm2(&n, S->matElem, &incx);
        rscl(&n, &froNorm, S->matElem, &incx);
        //        ASDP_ONE_rk_MAT(U);
        //        ASDP_ONE_rk_MAT(V);
        //        ASDP_ONE_rk_MAT(S);
        ASDPSetVarUVS(ASolver, iCone, U, FLAG_U);
        ASDPSetVarUVS(ASolver, iCone, V, FLAG_V);
        ASDPSetVarUVS(ASolver, iCone, S, FLAG_S);
        ASDPSetVarUVS(ASolver, iCone, lag, -10);
    }
#ifdef USE_CG
    ASDPCGCreate(ASolver);
#endif
    ASDPLanczosCreate(ASolver);
    ASDPInitM1M2Temp(ASolver);
    ASDPInitAuxiCri(ASolver);
exit_cleanup:
    return retcode;
}

extern int AUG_RANK(asdp *ASolver, int *BlkDims, int nBlks, double aug_factor)
{
    // double rank
    asdp_rk_mat_dense **U = ASolver->U;
    asdp_rk_mat_dense **V = ASolver->V;
    asdp_rk_mat_dense **R = ASolver->R;
    asdp_rk_mat_dense **S = ASolver->Vlag;
    asdp_rk_mat_dense **Grad = ASolver->Grad;
    asdp_rk_mat_dense **M2temp = ASolver->M2temp;
    int allElem = 0;
    int ret = 1;
    if (ASolver->nLpCols > 0)
    {
        allElem += ASolver->nLpCols;
    }
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        int new_rank = ASDP_MIN(ceil(U[iCone]->rank * aug_factor), ASolver->rank_max[iCone]);
        int aug_rank = new_rank - U[iCone]->rank;
        int old_rank = U[iCone]->rank;
        // printf("old_rank: %d, new_rank: %d, max_rank: %d, aug_rank: %d\n", U[iCone]->rank, new_rank, ASolver->rank_max[iCone], aug_rank);
        if (new_rank >= ASolver->rank_max[iCone])
        {
            printf("**Rank truncated to sqrt(2m) on SDP Cone No.%d.\n", iCone);
        }

        REALLOC(&U[iCone]->matElem, old_rank * U[iCone]->nRows, new_rank * U[iCone]->nRows);
        // lpRandom(&U[iCone]->matElem[U[iCone]->nRows * old_rank], U[iCone]->nRows * aug_rank);

        lpRandomDiag(&U[iCone]->matElem[U[iCone]->nRows * old_rank], U[iCone]->nRows * aug_rank, U[iCone]->nRows, aug_rank);
        U[iCone]->rank = new_rank;

        REALLOC(&V[iCone]->matElem, old_rank * U[iCone]->nRows, new_rank * U[iCone]->nRows);
        // lpRandom(&V[iCone]->matElem[U[iCone]->nRows * old_rank], U[iCone]->nRows * aug_rank);

        lpRandomDiag(&V[iCone]->matElem[U[iCone]->nRows * old_rank], U[iCone]->nRows * aug_rank, U[iCone]->nRows, aug_rank);
        V[iCone]->rank = new_rank;

        REALLOC(&R[iCone]->matElem, old_rank * U[iCone]->nRows, new_rank * U[iCone]->nRows);
        // lpRandom(&R[iCone]->matElem[U[iCone]->nRows * old_rank], U[iCone]->nRows * aug_rank);
        lpRandomDiag(&R[iCone]->matElem[U[iCone]->nRows * old_rank], U[iCone]->nRows * aug_rank, U[iCone]->nRows, aug_rank);
        R[iCone]->rank = new_rank;

        REALLOC(&Grad[iCone]->matElem, old_rank * U[iCone]->nRows, new_rank * U[iCone]->nRows);
        // lpRandom(&Grad[iCone]->matElem[U[iCone]->nRows * old_rank], U[iCone]->nRows * aug_rank);
        lpRandomDiag(&Grad[iCone]->matElem[U[iCone]->nRows * old_rank], U[iCone]->nRows * aug_rank, U[iCone]->nRows, aug_rank);
        Grad[iCone]->rank = new_rank;

        REALLOC(&S[iCone]->matElem, old_rank * U[iCone]->nRows, new_rank * U[iCone]->nRows);
        // lpRandom(&S[iCone]->matElem[U[iCone]->nRows * old_rank], U[iCone]->nRows * aug_rank);
        lpRandomDiag(&S[iCone]->matElem[U[iCone]->nRows * old_rank], U[iCone]->nRows * aug_rank, U[iCone]->nRows, aug_rank);
        S[iCone]->rank = new_rank;

        REALLOC(&M2temp[iCone]->matElem, old_rank * U[iCone]->nRows, new_rank * U[iCone]->nRows);
        M2temp[iCone]->rank = new_rank;

        allElem += V[iCone]->rank * V[iCone]->nRows;

        REALLOC(&ASolver->UsubV[iCone], old_rank * U[iCone]->nRows, new_rank * U[iCone]->nRows);

        REALLOC(&ASolver->bLinSys[iCone], old_rank * ASolver->U[iCone]->nRows, new_rank * ASolver->U[iCone]->nRows);
    }
    ASDPCGReAllocate(ASolver);

    lbfgs_node *head = ASolver->lbfgsHis;
    head->allElem = allElem;
    ASDP_REALLOC(head->s, double, allElem);
    ASDP_REALLOC(head->y, double, allElem);
    lbfgs_node *node = head;
    for (int nodeNum = 0; nodeNum < ASolver->hisRecT; ++nodeNum)
    {
        lbfgs_node *nextNode = node->next;
        nextNode->allElem = allElem;
        ASDP_REALLOC(nextNode->s, double, allElem);
        ASDP_REALLOC(nextNode->y, double, allElem);
        node = node->next;
    }

    return CheckAllRankMax(ASolver, aug_factor);
}

extern void ASDPDestroyADMMVars(asdp *ASolver)
{
    // destroy for function ASDPInitConeData_and_UVS
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASDP_FREE(ASolver->UsubV[iCone]);
        ASDP_FREE(ASolver->M2temp[iCone]->matElem);
        ASDP_FREE(ASolver->M2temp[iCone]);
        ASDP_FREE(ASolver->U[iCone]->matElem);
        ASDP_FREE(ASolver->U[iCone]);
        ASDP_FREE(ASolver->V[iCone]->matElem);
        ASDP_FREE(ASolver->V[iCone]);
        ASDP_FREE(ASolver->Vlag[iCone]->matElem);
        ASDP_FREE(ASolver->Vlag[iCone]);
        ASDP_FREE(ASolver->lag[iCone]->matElem);
        ASDP_FREE(ASolver->lag[iCone]);
        ASDP_FREE(ASolver->bLinSys[iCone]);
    }
    if (ASolver->nLpCols > 0)
    {
        ASDP_FREE(ASolver->uLp->matElem);
        ASDP_FREE(ASolver->vLp->matElem);
        ASDP_FREE(ASolver->vlagLp->matElem);
        ASDP_FREE(ASolver->lagLp->matElem);
    }
    ASDP_FREE(ASolver->uLp);
    ASDP_FREE(ASolver->vlagLp);
    ASDP_FREE(ASolver->vLp);
    ASDP_FREE(ASolver->U);
    ASDP_FREE(ASolver->V);
    ASDP_FREE(ASolver->Vlag);
    ASDP_FREE(ASolver->lag);
    ASDP_FREE(ASolver->bLinSys);
    ASDP_FREE(ASolver->M2temp);
    ASDP_FREE(ASolver->UsubV);
    ASDP_FREE(ASolver->negLambd);
    if (ASolver->nLpCols > 0)
    {
        ASDP_FREE(ASolver->uSubvLp);
    }
    ASDPCGDetroy(ASolver);
    ASDPLanczosDetroy(ASolver);
}

extern asdp_retcode ASDPInitSolver(asdp *ASolver, int nRows, int nCones, int *blkDims, int nLpCols, double rho, double rhoMax, int maxIter, int strategy)
{
    ASolver->nLpCols = nLpCols;
    ASolver->strategy = strategy;
    // nRows: number of constraints
    // nCones: number of blks
    asdp_retcode retcode = ASDP_RETCODE_OK;

    ASolver->nRows = nRows;   // constraint number
    ASolver->nCones = nCones; // blk number

    ASDP_INIT(ASolver->rowRHS, double, nRows);
    ASDP_MEMCHECK(ASolver->rowRHS);

    // set dual variable
    ASDP_INIT(ASolver->dualVar, double, nRows);
    ASDP_ZERO(ASolver->dualVar, double, nRows);
    ASDP_MEMCHECK(ASolver->dualVar);

    ASDP_INIT(ASolver->bestDualVar, double, nRows);
    ASDP_ZERO(ASolver->bestDualVar, double, nRows);
    ASDP_MEMCHECK(ASolver->bestDualVar);

    // set Auxiliary variable
    ASDP_INIT(ASolver->constrVal, constrValStruct *, nCones);
    ASDP_MEMCHECK(ASolver->constrVal);
    for (int iCone = 0; iCone < nCones; ++iCone)
    {
        ASDP_INIT(ASolver->constrVal[iCone], constrValStruct, 1);
        ASDP_MEMCHECK(ASolver->constrVal[iCone]);
    }

    ASDP_INIT(ASolver->constrValSum, double, nRows);
    ASDP_MEMCHECK(ASolver->constrValSum);

    ASDP_INIT(ASolver->ACones, asdp_cone *, nCones); // allocate n cones
    ASDP_MEMCHECK(ASolver->ACones);

    for (int iCone = 0; iCone < nCones; ++iCone)
    {
        ASDP_CALL(AConeCreate(&ASolver->ACones[iCone])); // allocate n cone initialize zero
    }
    double dim = (double)nLpCols;
    for (int iCone = 0; iCone < nCones; ++iCone)
    {
        dim += (double)blkDims[iCone];
    }
    if (fabs(rho) < 1e-10)
    {
        ASolver->rho = 1.0 / dim;
    }
    else
    {
        ASolver->rho = rho;
    }

    ASolver->rhoMax = rhoMax;
    ASolver->maxInter = maxIter;

    ASDP_INIT(ASolver->dimacError, double, 5);
    ASDP_MEMCHECK(ASolver->dimacError);
    ASDP_ZERO(ASolver->dimacError, double, 4);

    if (ASolver->nLpCols > 0)
    {
        ASDP_INIT(ASolver->constrValLP, constrValStruct *, ASolver->nLpCols);
        ASDP_MEMCHECK(ASolver->constrValLP);
        for (int iLpCol = 0; iLpCol < ASolver->nLpCols; ++iLpCol)
        {
            ASDP_INIT(ASolver->constrValLP[iLpCol], constrValStruct, 1);
            ASDP_MEMCHECK(ASolver->constrValLP[iLpCol]);
        }
    }
    ASDP_INIT(ASolver->lpCone, asdpLPCone, 1); // allocate n cones
    ASDP_MEMCHECK(ASolver->lpCone);

exit_cleanup:
    return retcode;
}

extern asdp_retcode ASDPDetermineRank(asdp *ASolver, int *blkDims, double timesRank, int rankSpecify)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    int nCones = ASolver->nCones;
    ASDP_INIT(ASolver->rankElem, int, nCones);
    ASDP_MEMCHECK(ASolver->rankElem);
    ASDP_INIT(ASolver->rank_max, int, nCones);
    ASDP_MEMCHECK(ASolver->rank_max);
    for (int iCone = 0; iCone < nCones; ++iCone)
    {
        asdp_cone *ACone = ASolver->ACones[iCone];
        int nnzRows = 0;
        ACone->nnzStat(ACone->coneData, &nnzRows);
        if (rankSpecify == 0)
        {
            if (timesRank <= 1e-6)
            {
                ASolver->rankElem[iCone] = ASDP_MIN((int)sqrt(2 * nnzRows) + 1, blkDims[iCone]);
            }
            else
            {
                ASolver->rankElem[iCone] = ASDP_MIN(ceil(timesRank * log(blkDims[iCone])), ASDP_MIN((int)sqrt(2 * nnzRows) + 1, blkDims[iCone]));
            }
            ASolver->rankElem[iCone] = ASDP_MAX(1, ASolver->rankElem[iCone]);
        }
        else
        {
            ASolver->rankElem[iCone] = rankSpecify;
        }
        ASolver->rank_max[iCone] = ASDP_MIN((int)sqrt(2 * nnzRows) + 1, blkDims[iCone]);
    }
exit_cleanup:
    return retcode;
}

extern asdp_retcode ASDPCGCreate(asdp *ASolver)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    ASDP_INIT(ASolver->CGLinsys, asdp_cg_linsys *, ASolver->nCones);
    ASDP_MEMCHECK(ASolver->CGLinsys);
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASDPCGSolverCreate(&ASolver->CGLinsys[iCone], ASolver->U[iCone]->nRows, ASolver->U[iCone]->rank, ASolver->nRows);
    }
exit_cleanup:
    return retcode;
}

extern void ASDPCGReAllocate(asdp *ASolver)
{
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASDPCGSolverReCreate(&ASolver->CGLinsys[iCone], ASolver->U[iCone]->nRows, ASolver->U[iCone]->rank, ASolver->nRows);
    }
}

extern asdp_retcode ASDPLanczosCreate(asdp *ASolver)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    ASDP_INIT(ASolver->LanczoSys, asdp_lanczos *, ASolver->nCones);
    ASDP_MEMCHECK(ASolver->LanczoSys);
    ASDP_INIT(ASolver->LanczosStart, double *, ASolver->nCones);
    ASDP_MEMCHECK(ASolver->LanczosStart);
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ALanczosCreate(&ASolver->LanczoSys[iCone]);
        ALanczosInit(ASolver->LanczoSys[iCone], ASolver->U[iCone]->nRows, ASDP_MIN(ASolver->U[iCone]->nRows * 10, 100));
        ASDP_INIT(ASolver->LanczosStart[iCone], double, ASolver->U[iCone]->nRows);
        ASDP_MEMCHECK(ASolver->LanczosStart[iCone]);
        ASolver->LanczosStart[iCone][0] = 1.0;
    }
exit_cleanup:
    return retcode;
}

extern void ASDPCGDetroy(asdp *ASolver)
{
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        CGSolverClear(ASolver->CGLinsys[iCone]);
    }
    ASDP_FREE(ASolver->CGLinsys);
}

extern void ASDPLanczosDetroy(asdp *ASolver)
{
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ALanczosClear(ASolver->LanczoSys[iCone]);
        ALanczosDestroy(&ASolver->LanczoSys[iCone]);
        ASDP_FREE(ASolver->LanczosStart[iCone]);
    }
    ASDP_FREE(ASolver->LanczoSys);
    ASDP_FREE(ASolver->LanczosStart);
}

extern void ASDPDestroySolver(asdp *ASolver)
{
    // corresponding to `ASDPInitSolver`
    ASDP_FREE(ASolver->rowRHS);
    ASDP_FREE(ASolver->dualVar);
    ASDP_FREE(ASolver->bestDualVar);
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        AConeDestroy(&ASolver->ACones[iCone]);
        if (ASolver->constrVal[iCone]->type == ASDP_DENSE)
        {
            constrValDense *dense = (constrValDense *)ASolver->constrVal[iCone]->constrVal;
            ASDP_FREE(dense->val);
            ASDP_FREE(ASolver->constrVal[iCone]);
        }
        else if (ASolver->constrVal[iCone]->type == ASDP_SPARSE)
        {
            constrValSparse *sparse = (constrValSparse *)ASolver->constrVal[iCone]->constrVal;
            ASDP_FREE(sparse->nnzIdx);
            ASDP_FREE(sparse->val);
            ASDP_FREE(ASolver->constrVal[iCone]);
        }
    }
    ASDP_FREE(ASolver->constrVal);
    ASDP_FREE(ASolver->constrValSum);
    ASDP_FREE(ASolver->ACones);
    ASDP_FREE(ASolver->dimacError);
    if (ASolver->nLpCols > 0)
    {
        for (int iLPCol = 0; iLPCol < ASolver->nLpCols; ++iLPCol)
        {
            if (ASolver->constrValLP[iLPCol]->type == ASDP_DENSE)
            {
                constrValDense *dense = (constrValDense *)ASolver->constrValLP[iLPCol]->constrVal;
                ASDP_FREE(dense->val);
                ASDP_FREE(ASolver->constrValLP[iLPCol]);
            }
            else if (ASolver->constrValLP[iLPCol]->type == ASDP_SPARSE)
            {
                constrValSparse *sparse = (constrValSparse *)ASolver->constrValLP[iLPCol]->constrVal;
                ASDP_FREE(sparse->nnzIdx);
                ASDP_FREE(sparse->val);
                ASDP_FREE(ASolver->constrValLP[iLPCol]);
            }
            ASDP_FREE(ASolver->constrVal);
        }
        ASDP_FREE(ASolver->constrValLP);
        ALpConeDestroy(&ASolver->lpCone);
        ASDP_FREE(ASolver->lpCone);
    }
}

extern void ASDPDestroyRankElem(asdp *ASolver)
{
    ASDP_FREE(ASolver->rankElem);
}

extern asdp_retcode ASDPSetCone(asdp *ASolver, int iCone, void *userCone)
{

    asdp_retcode retcode = ASDP_RETCODE_OK;

    ASDP_CALL(AConeSetData(ASolver->ACones[iCone], userCone));

exit_cleanup:
    return retcode;
}

extern asdp_retcode ASDPSetVarUVS(asdp *ASolver, int iCone, void *dataUV, int FLAG_UV)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    asdp_rk_mat_dense *data = (asdp_rk_mat_dense *)dataUV;
    if (FLAG_UV == FLAG_U)
    {
        ASolver->U[iCone] = data;
    }
    else if (FLAG_UV == FLAG_V)
    {
        ASolver->V[iCone] = data;
    }
    else if (FLAG_UV == FLAG_S)
    {
        ASolver->Vlag[iCone] = data;
    }
    else if (FLAG_UV == FLAG_R)
    {
        ASolver->R[iCone] = data;
    }else{
        ASolver->lag[iCone] = data;
    }
exit_cleanup:
    return retcode;
}

extern asdp_retcode ASDP_ONE_rk_MAT(asdp_rk_mat_dense *U)
{
    // initialize rk_mat to all one
    asdp_retcode retcode = ASDP_RETCODE_OK;
    ASDP_INIT(U->matElem, double, U->nRows * U->rank);
    ASDP_MEMCHECK(U->matElem);
    ASDP_ONE(U->matElem, U->nRows * U->rank);
exit_cleanup:
    return retcode;
}

extern asdp_retcode ASDP_RANDOM_rk_MAT(asdp_rk_mat_dense *U)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    int n = U->nRows * U->rank;
    ASDP_INIT(U->matElem, double, n);
    ASDP_MEMCHECK(U->matElem);
    for (int i = 0; i < n; ++i)
    {
        U->matElem[i] = (double)rand() / RAND_MAX;
        U->matElem[i] -= (double)rand() / RAND_MAX;
    }
exit_cleanup:
    return retcode;
}

extern asdp_retcode ASDP_BMRANDOM_rk_MAT(asdp_rk_mat_dense *U)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    int n = U->nRows * U->rank;
    ASDP_INIT(U->matElem, double, n);
    ASDP_MEMCHECK(U->matElem);
    for (int i = 0; i < n; ++i)
    {
        U->matElem[i] = (double)rand() / RAND_MAX;
        U->matElem[i] -= (double)rand() / RAND_MAX;
    }
    int incx = 1;
    double froNorm = nrm2(&n, U->matElem, &incx);
    rscl(&n, &froNorm, U->matElem, &incx);
exit_cleanup:
    return retcode;
}

extern void ASDPSetDualObjective(asdp *ASolver, double *dObj)
{
    ASDP_MEMCPY(ASolver->rowRHS, dObj, double, ASolver->nRows);
    return;
}

extern void ASDPRkMatSub(asdp_rk_mat_dense *A, asdp_rk_mat_dense *B, double alpha, asdp_rk_mat_dense *S, int flag_UV)
{
    // A = A - alpha * B + S
    double negAlpha = -1 * alpha;
    int n = A->nRows * A->rank;
    int incx = 1;
    axpy(&n, &negAlpha, B->matElem, &incx, A->matElem, &incx);
    if (flag_UV == FLAG_U)
    {
        double one = 1.0;
        axpy(&n, &one, S->matElem, &incx, A->matElem, &incx);
    }
    else if (flag_UV == FLAG_V)
    {
        double minusOne = -1.0;
        axpy(&n, &minusOne, S->matElem, &incx, A->matElem, &incx);
    }
}

extern void ASDPUVt(sdp_coeff *UVt_w_sum, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V)
{
    /* sdp_coeff_w_sum is the sum of all sdp data in the cone,
     It has two data type
     - sparse: means most all sdp data coeff is sparse
     - dense:  there exist a sdp data coeff is dense

     calculate (UVt + VUt)
     only store lower semi triangular part
     */
    if (UVt_w_sum->dataType == SDP_COEFF_SPARSE)
    {
        sdp_coeff_sparse *sparse = (sdp_coeff_sparse *)UVt_w_sum->dataMat;
        if (sparse->nTriMatElem / pow(sparse->nSDPCol, 2.0) < 0.08)
        {
            // method1:
            int row = 0;
            int col = 0;
            int incx = 1;
            int incy = 1;
            for (int i = 0; i < sparse->nTriMatElem; ++i)
            {
                row = sparse->triMatRow[i];
                col = sparse->triMatCol[i];
                incx = U->nRows;
                incy = V->nRows;
                if (row != col)
                {
                    sparse->triMatElem[i] = dot(&(U->rank), &U->matElem[row], &incx, &V->matElem[col], &incy);
                    sparse->triMatElem[i] += dot(&(U->rank), &U->matElem[col], &incx, &V->matElem[row], &incy);
                }
                else
                {
                    sparse->triMatElem[i] = 2 * dot(&(U->rank), &U->matElem[row], &incx, &V->matElem[col], &incy);
                }
            }
        }
        else
        {
            // method2:
            double *fullDataMat;
            ASDP_INIT(fullDataMat, double, sparse->nSDPCol * sparse->nSDPCol);
            ASDP_ZERO(fullDataMat, double, sparse->nSDPCol * sparse->nSDPCol);
            fds_syr2k(ACharConstantUploLow, 'N', U->nRows, U->rank, 1.0, U->matElem, V->matElem, 0.0, fullDataMat);
            int idx = 0;
            int row = 0;
            int col = 0;
            for (int i = 0; i < sparse->nTriMatElem; ++i)
            {
                row = sparse->triMatRow[i];
                col = sparse->triMatCol[i];
                sparse->triMatElem[i] = fullDataMat[col * sparse->nSDPCol + row];
            }
            ASDP_FREE(fullDataMat);
        }
    }
    else if (UVt_w_sum->dataType == SDP_COEFF_DENSE)
    {
        sdp_coeff_dense *dense = (sdp_coeff_dense *)UVt_w_sum->dataMat;
        double *fullDataMat;
        ASDP_INIT(fullDataMat, double, dense->nSDPCol * dense->nSDPCol);
        ASDP_ZERO(fullDataMat, double, dense->nSDPCol * dense->nSDPCol);
        // alpha = 1.0, beta = 0.0;
        fds_syr2k(ACharConstantUploLow, 'N', U->nRows, U->rank, 1.0, U->matElem, V->matElem, 0.0, fullDataMat);
        int idx = 0;
        int row = 0;
        for (int col = 0; col < dense->nSDPCol; ++col)
        {
            ASDP_MEMCPY(&dense->dsMatElem[idx], &fullDataMat[dense->nSDPCol * col + row], double, dense->nSDPCol - col);
            row++;
            idx += (dense->nSDPCol - col);
        }
        ASDP_FREE(fullDataMat);
    }
}

extern void valRes(void *constrVal, constrType type, double **res)
{
    if (type == ASDP_DENSE)
    {
        constrValDense *dense = (constrValDense *)constrVal;
        *res = dense->val;
    }
    else if (type == ASDP_SPARSE)
    {
        constrValSparse *sparse = (constrValSparse *)constrVal;
        *res = sparse->val;
    }
}

extern void addDense(double *alpha, void *constrVal, double *vec)
{
    constrValDense *dense = (constrValDense *)constrVal;
    int incx = 1;
    axpy(&dense->nnz, alpha, dense->val, &incx, vec, &incx);
}

extern void addSparse(double *alpha, void *constrVal, double *vec)
{
    constrValSparse *sparse = (constrValSparse *)constrVal;
    int idx = 0;
    for (int i = 0; i < sparse->nnz; ++i)
    {
        idx = sparse->nnzIdx[i];
        vec[idx] += alpha[0] * sparse->val[i];
    }
}

extern void zeroDense(void *constrVal)
{
    constrValDense *dense = (constrValDense *)constrVal;
    int incx = 1;
    ASDP_ZERO(dense->val, double, dense->nnz);
}

extern void zeroSparse(void *constrVal)
{
    constrValSparse *sparse = (constrValSparse *)constrVal;
    for (int i = 0; i < sparse->nnz; ++i)
    {
        sparse->val[i] = 0;
    }
}

extern void ASDPInitConstrVal(asdp_cone *ACone, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *constrVal, int FLAG_UV)
{
    ASDPUVt(ACone->sdp_coeff_w_sum, U, V);
    // Calculate Constraint Value for one cone
    ACone->coneAUV(ACone->coneData, U, V, constrVal, ACone->sdp_coeff_w_sum, FLAG_UV);
}

extern void ASDPInitConstrValAll(asdp *ASolver, asdp_rk_mat_lp *uLpDummy, asdp_rk_mat_lp *vLpDummy, asdp_rk_mat_dense **U, asdp_rk_mat_dense **V)
{
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        double *constrVal;
        valRes(ASolver->constrVal[iCone]->constrVal, ASolver->constrVal[iCone]->type, &constrVal);
        ASDPInitConstrVal(ASolver->ACones[iCone], U[iCone], V[iCone], constrVal, FLAG_INI);
    }
}

extern void ASDPInitConstrValAllLP(asdp *ASolver, asdp_rk_mat_lp *uLp, asdp_rk_mat_lp *vLp, asdp_rk_mat_dense **U, asdp_rk_mat_dense **V)
{
    asdpLPCone *lpCone = ASolver->lpCone;
    ASDPInitConstrValAll(ASolver, uLp, vLp, U, V);
    for (int iCol = 0; iCol < lpCone->nCol; ++iCol)
    {
        double *constrVal;
        valRes(ASolver->constrValLP[iCol]->constrVal, ASolver->constrValLP[iCol]->type, &constrVal);
        lpCone->coneAUV(lpCone->coneData, uLp, vLp, constrVal, iCol);
    }
}

extern void ASDPInitConstrValObjVal(asdp_cone *ACone, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *constrVal, double *UVobjVal, int FLAG_UV)
{
    ASDPUVt(ACone->sdp_obj_sum, U, V);
    // Calculate Constraint Value for one cone
    ACone->objAUV(ACone->coneData, U, V, UVobjVal, ACone->sdp_obj_sum, FLAG_UV);
    ACone->coneAUV(ACone->coneData, U, V, constrVal, ACone->sdp_obj_sum, FLAG_UV);
}

extern void ASDPObjConstrValAll(asdp *ASolver, asdp_rk_mat_dense **U, asdp_rk_mat_dense **V, double *UVobjVal)
{
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASolver->constrVal[iCone]->zero(ASolver->constrVal[iCone]->constrVal);
        double *constrVal;
        valRes(ASolver->constrVal[iCone]->constrVal, ASolver->constrVal[iCone]->type, &constrVal);
        ASDPInitConstrValObjVal(ASolver->ACones[iCone], U[iCone], V[iCone], constrVal, UVobjVal, FLAG_INI);
    }
}

extern void ASDPObjConstrValAllLP(asdp *ASolver, asdp_rk_mat_lp *uLp, asdp_rk_mat_lp *vLp, asdp_rk_mat_dense **U, asdp_rk_mat_dense **V, double *UVobjVal)
{
    ASolver->lpCone->objAUV(ASolver->lpCone->coneData, uLp, vLp, UVobjVal);
    for (int iCol = 0; iCol < uLp->nLPCols; ++iCol)
    {
        double *constrVal;
        valRes(ASolver->constrValLP[iCol]->constrVal, ASolver->constrValLP[iCol]->type, &constrVal);
        ASolver->lpCone->coneAUV(ASolver->lpCone->coneData, uLp, vLp, constrVal, iCol);
    }
    ASDPObjConstrValAll(ASolver, U, V, UVobjVal);
}

extern void ASDPUpdateConstrVal(asdp_cone *ACone, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, constrValStruct *constrVal, int FLAG_UV)
{
    // constrVal = Acal(UVt)
    double *constrValRes;
    valRes(constrVal->constrVal, constrVal->type, &constrValRes);
    ASDPInitConstrVal(ACone, U, V, constrValRes, FLAG_UV);
}

extern void ASDPUpdateConstrValCG(asdp_cone *ACone, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *weight, int FLAG_UV)
{
    ASDPUVt(ACone->sdp_coeff_w_sum, U, V);
    ACone->coneAUV(ACone->coneData, U, V, weight, ACone->sdp_coeff_w_sum, FLAG_UV);
    if (ACone->cone == ASDP_CONETYPE_SPARSE_SDP)
    {
        double temp = 0.0;
        int idx = 0;
        asdp_cone_sdp_sparse *sparse = (asdp_cone_sdp_sparse *)ACone->coneData;
        for (int i = sparse->nRowElem - 1; i >= 0; i--)
        {
            temp = weight[i];
            weight[i] = 0.0;
            idx = sparse->rowIdx[i];
            weight[idx] = temp;
        }
    }
}

extern void ASDPUpdateConstrValLP(asdpLPCone *lp_cone, asdp_rk_mat_lp *uLp, asdp_rk_mat_lp *vLp, constrValStruct *constrVal, int iCol)
{
    double *constrValRes;
    valRes(constrVal->constrVal, constrVal->type, &constrValRes);
    lp_cone->coneAUV(lp_cone->coneData, uLp, vLp, constrValRes, iCol);
}

extern void ASDPInitConstrValSum(asdp *ASolver)
{
    ASDP_ZERO(ASolver->constrValSum, double, ASolver->nRows);
    double alpha = 1.0;
    int incx = 1;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASolver->constrVal[iCone]->add(&alpha, ASolver->constrVal[iCone]->constrVal, ASolver->constrValSum);
    }
}

extern void ASDPInitConstrValSumLP(asdp *ASolver)
{
    ASDP_ZERO(ASolver->constrValSum, double, ASolver->nRows);
    double alpha = 1.0;
    for (int iCol = 0; iCol < ASolver->nLpCols; ++iCol)
    {
        ASolver->constrValLP[iCol]->add(&alpha, ASolver->constrValLP[iCol]->constrVal, ASolver->constrValSum);
    }

    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASolver->constrVal[iCone]->add(&alpha, ASolver->constrVal[iCone]->constrVal, ASolver->constrValSum);
    }
}

extern void ASDPConstrValSumBMtemp(asdp *ASolver, double *q1)
{
    double alpha = 1.0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASolver->constrVal[iCone]->add(&alpha, ASolver->constrVal[iCone]->constrVal, q1);
    }
}

extern void ASDPConstrValSumBMtempLP(asdp *ASolver, double *q1)
{
    double alpha = 1.0;
    int incx = 1;
    for (int iCol = 0; iCol < ASolver->nLpCols; ++iCol)
    {
        ASolver->constrValLP[iCol]->add(&alpha, ASolver->constrValLP[iCol]->constrVal, q1);
    }
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASolver->constrVal[iCone]->add(&alpha, ASolver->constrVal[iCone]->constrVal, q1);
    }
}

extern asdp_retcode ASDPInitM1M2Temp(asdp *ASolver)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    ASDP_INIT(ASolver->M2temp, asdp_rk_mat_dense *, ASolver->nCones);
    ASDP_MEMCHECK(ASolver->M2temp);
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASDP_INIT(ASolver->M2temp[iCone], asdp_rk_mat_dense, 1);
        ASDP_MEMCHECK(ASolver->M2temp[iCone]);
        ASolver->M2temp[iCone]->rank = ASolver->U[iCone]->rank;
        ASolver->M2temp[iCone]->nRows = ASolver->U[iCone]->nRows;
        ASDP_INIT(ASolver->M2temp[iCone]->matElem, double, ASolver->M2temp[iCone]->rank * ASolver->M2temp[iCone]->nRows);
        ASDP_MEMCHECK(ASolver->M2temp[iCone]->matElem);
    }
exit_cleanup:
    return retcode;
}

extern asdp_retcode ASDPInitAuxiCri(asdp *ASolver)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    int nConstr = ASolver->nRows;
    ASDP_INIT(ASolver->constrVio, double, nConstr);
    ASDP_MEMCHECK(ASolver->constrVio);
    ASDP_INIT(ASolver->UsubV, double *, ASolver->nCones);
    ASDP_MEMCHECK(ASolver->UsubV);
    if (ASolver->nLpCols > 0)
    {
        ASDP_INIT(ASolver->uSubvLp, double, ASolver->nLpCols);
        ASDP_MEMCHECK(ASolver->uSubvLp);
    }
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASDP_INIT(ASolver->UsubV[iCone], double, ASolver->U[iCone]->rank * ASolver->U[iCone]->nRows);
        ASDP_MEMCHECK(ASolver->UsubV[iCone]);
    }
    ASDP_INIT(ASolver->negLambd, double, ASolver->nRows);
    ASDP_MEMCHECK(ASolver->negLambd);
exit_cleanup:
    return retcode;
}

extern void averageUV(asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, asdp_rk_mat_dense *UVavg)
{
    int n = U->rank * U->nRows;
    for (int i = 0; i < n; ++i)
    {
        UVavg->matElem[i] = (U->matElem[i] + V->matElem[i]) / 2;
    }
}

extern void ASDPCalObj(asdp *ASolver, int FLAG_BM_USEV)
{
    ASolver->pObjVal = 0.0;

    if (FLAG_BM_USEV == FLAG_U)
    {
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            asdp_rk_mat_dense *R = ASolver->R[iCone];
            asdp_cone *ACone = ASolver->ACones[iCone];
            sdp_coeff *dummy;
            asdp_rk_mat_dense *U = ASolver->U[iCone];
            asdp_rk_mat_dense *V = ASolver->V[iCone];
            averageUV(U, V, R);
            ASDPUVt(ACone->sdp_obj_sum, R, R);
            ACone->objAUV(ACone->coneData, R, R, &ASolver->pObjVal, ACone->sdp_obj_sum, FLAG_INI);
        }
    }
    else if (FLAG_BM_USEV == FLAG_R)
    {
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            asdp_rk_mat_dense *R = ASolver->R[iCone];
            asdp_cone *ACone = ASolver->ACones[iCone];
            sdp_coeff *dummy;
            ASDPUVt(ACone->sdp_obj_sum, R, R);
            ACone->objAUV(ACone->coneData, R, R, &ASolver->pObjVal, ACone->sdp_obj_sum, FLAG_INI);
        }
    }
    else if (FLAG_BM_USEV == FLAG_V)
    {
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            asdp_rk_mat_dense *R = ASolver->R[iCone];
            asdp_cone *ACone = ASolver->ACones[iCone];
            sdp_coeff *dummy;
            asdp_rk_mat_dense *V = ASolver->V[iCone];
            ASDPUVt(ACone->sdp_obj_sum, V, V);
            ACone->objAUV(ACone->coneData, V, V, &ASolver->pObjVal, ACone->sdp_obj_sum, FLAG_INI);
        }
    }
    ASolver->pObjVal /= ASolver->cScaleFactor;
}

extern void ASDPCalObjLP(asdp *ASolver, int FLAG_BM_USEV)
{
    ASolver->pObjVal = 0.0;
    if (FLAG_BM_USEV == FLAG_U)
    {
        double alpha = 0.5;
        int incx = 1;
        double zero = 0.0;
        scal(&ASolver->rLp->nLPCols, &zero, ASolver->rLp->matElem, &incx);
        axpy(&(ASolver->rLp->nLPCols), &alpha, ASolver->uLp->matElem, &incx, ASolver->rLp->matElem, &incx);
        axpy(&(ASolver->rLp->nLPCols), &alpha, ASolver->vLp->matElem, &incx, ASolver->rLp->matElem, &incx);
        ASolver->lpCone->objAUV(ASolver->lpCone->coneData, ASolver->rLp, ASolver->rLp, &ASolver->pObjVal);
    }
    else if (FLAG_BM_USEV == FLAG_R)
    {
        ASolver->lpCone->objAUV(ASolver->lpCone->coneData, ASolver->rLp, ASolver->rLp, &ASolver->pObjVal);
    }
    else if (FLAG_BM_USEV == FLAG_V)
    {
        ASolver->lpCone->objAUV(ASolver->lpCone->coneData, ASolver->vLp, ASolver->vLp, &ASolver->pObjVal);
    }
    if (FLAG_BM_USEV == FLAG_U)
    {
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            asdp_rk_mat_dense *R = ASolver->R[iCone];
            asdp_cone *ACone = ASolver->ACones[iCone];
            sdp_coeff *dummy;
            asdp_rk_mat_dense *U = ASolver->U[iCone];
            asdp_rk_mat_dense *V = ASolver->V[iCone];
            averageUV(U, V, R);
            ACone->objAUV(ACone->coneData, R, R, &ASolver->pObjVal, dummy, FLAG_OBJ);
        }
    }
    else if (FLAG_BM_USEV == FLAG_R)
    {
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            asdp_rk_mat_dense *R = ASolver->R[iCone];
            asdp_cone *ACone = ASolver->ACones[iCone];
            sdp_coeff *dummy;
            ACone->objAUV(ACone->coneData, R, R, &ASolver->pObjVal, dummy, FLAG_OBJ);
        }
    }
    else if (FLAG_BM_USEV == FLAG_V)
    {
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            asdp_rk_mat_dense *R = ASolver->R[iCone];
            asdp_cone *ACone = ASolver->ACones[iCone];
            sdp_coeff *dummy;
            asdp_rk_mat_dense *V = ASolver->V[iCone];
            ACone->objAUV(ACone->coneData, V, V, &ASolver->pObjVal, dummy, FLAG_OBJ);
        }
    }
}

extern void BMCalq12p12(asdp *ASolver, asdp_rk_mat_lp *rLpDummy, asdp_rk_mat_lp *dDummy, asdp_rk_mat_dense **R, asdp_rk_mat_dense **D, double *q1, double *q2, double *p12)
{
    int n = ASolver->nRows;
    double one = 1.0;
    double zero = 0.0;
    int incx = 1;
    // ASDP_ZERO(q1, double, n);
    // ASDP_ZERO(q2, double, n);
    scal(&n, &zero, q1, &incx);
    scal(&n, &zero, q2, &incx);
    p12[0] = 0.0;
    p12[1] = 0.0;
    ASDPObjConstrValAll(ASolver, R, D, &p12[0]);
    ASDPConstrValSumBMtemp(ASolver, q1);
    double Two = 2.0;
    scal(&n, &Two, q1, &incx);
    p12[0] *= 2;

    ASDPObjConstrValAll(ASolver, D, D, &p12[1]);
    ASDPConstrValSumBMtemp(ASolver, q2);
}

extern void BMCalq12p12LP(asdp *ASolver, asdp_rk_mat_lp *rLp, asdp_rk_mat_lp *d, asdp_rk_mat_dense **R, asdp_rk_mat_dense **D, double *q1, double *q2, double *p12)
{
    int n = ASolver->nRows;
    double one = 1.0;
    int incx = 1;
    ASDP_ZERO(q1, double, n);
    ASDP_ZERO(q2, double, n);
    p12[0] = 0.0;
    p12[1] = 0.0;
    ASDPObjConstrValAllLP(ASolver, rLp, d, R, D, &p12[0]);
    ASDPConstrValSumBMtempLP(ASolver, q1);
    double Two = 2.0;
    scal(&n, &Two, q1, &incx);
    p12[0] *= 2;

    ASDPObjConstrValAllLP(ASolver, d, d, D, D, &p12[1]);
    ASDPConstrValSumBMtempLP(ASolver, q2);
}
static double functionVal(double a, double b, double c, double d, double x)
{
    return a * pow(x, 4) + b * pow(x, 3) + c * pow(x, 2) + d * x;
}

extern int BMLineSearch(double rho, int n, double *lambd, double p1, double p2, double *q0, double *q1, double *q2, double *tau)
{
    int incx = 1;
    double q2nrm2 = nrm2(&n, q2, &incx);
    double a = rho * q2nrm2 * q2nrm2 / 2;
    double b = rho * dot(&n, q1, &incx, q2, &incx);
    // q0 = (1 / rho * lambd + q0)
    double rhoInv = 1 / rho;
    axpy(&n, &rhoInv, lambd, &incx, q0, &incx);
    double q1nrm = nrm2(&n, q1, &incx);
    double c = p2 - rho * dot(&n, q0, &incx, q2, &incx) + rho * q1nrm * q1nrm / 2;
    double d = p1 - rho * dot(&n, q0, &incx, q1, &incx);
    double roots[3] = {0.0, 0.0, 0.0};
    int rootNum = ASDPcubic_equation(4 * a, 3 * b, 2 * c, d, roots);
    double f0 = 0.0;
    double tauMax = 1.0;
    double f1 = functionVal(a, b, c, d, 1.0);
    double fRoot1 = 1e+30;
    double fRoot2 = 1e+30;
    double fRoot3 = 1e+30;
    if (rootNum >= 1)
    {
        if (roots[0] > 1e-20 && roots[0] <= tauMax)
        {
            fRoot1 = functionVal(a, b, c, d, roots[0]);
        }
    }
    if (rootNum >= 2)
    {
        if (roots[1] > 1e-20 && roots[1] <= tauMax)
        {
            fRoot2 = functionVal(a, b, c, d, roots[1]);
        }
    }
    if (rootNum == 3)
    {
        if (roots[2] > 1e-20 && roots[2] <= tauMax)
        {
            fRoot3 = functionVal(a, b, c, d, roots[2]);
        }
    }
    else if (rootNum == 0)
    {
        ;
    }
    double minFval = ASDP_MIN(ASDP_MIN(ASDP_MIN(ASDP_MIN(f0, f1), fRoot1), fRoot2), fRoot3);
    if (fabs(minFval - f0) < 1e-10)
    {
        tau[0] = 0.0;
    }
    if (fabs(minFval - f1) < 1e-10)
    {
        tau[0] = 1.0;
    }
    if (fabs(minFval - fRoot1) < 1e-10)
    {
        tau[0] = roots[0];
    }
    if (fabs(minFval - fRoot2) < 1e-10)
    {
        tau[0] = roots[1];
    }
    if (fabs(minFval - fRoot3) < 1e-10)
    {
        tau[0] = roots[2];
    }
    return rootNum;
}

extern void ASDPCalDualObj(asdp *ASolver)
{
    int n = ASolver->nRows;
    int one = 1;
    ASolver->dObjVal = dot(&n, ASolver->rowRHS, &one, ASolver->dualVar, &one);
    ASolver->dObjVal /= ASolver->bScaleFactor;
}

extern void AConeDetectDenseRowSparsity(sdp_coeff *sdp_coeff_w_sum, int *nNnzRowSum)
{
    sdp_coeff_dense *dense = (sdp_coeff_dense *)sdp_coeff_w_sum->dataMat;
    int n = dense->nSDPCol;
    int row = 0;
    int col = 0;
    int *nNnzRowSumTemp;
    ASDP_INIT(nNnzRowSumTemp, int, n);
    ASDP_ZERO(nNnzRowSumTemp, int, n);

    nNnzRowSum[0] = 0;
    for (int i = 0; i < (n + 1) * n / 2; ++i)
    {
        if (fabs(dense->dsMatElem[i]) > 0.0 && nNnzRowSumTemp[row] == 0)
        {
            nNnzRowSumTemp[row] = 1;
            nNnzRowSum[0] += nNnzRowSumTemp[row];
        }
        row++;
        if (row == n)
        {
            col++;
            row = col;
        }
    }
    ASDP_FREE(nNnzRowSumTemp);
}

extern void AConeDenseDetectSparsity(sdp_coeff **sdp_coeff_w_sum_pointer)
{
    // add modify sdp_coeff_w_sum if spRatio is smaller than 0.1
    sdp_coeff *sdp_coeff_w_sum = *sdp_coeff_w_sum_pointer;
    sdp_coeff_dense *dense = (sdp_coeff_dense *)sdp_coeff_w_sum->dataMat;
    int n = dense->nSDPCol;
    int row = 0;
    int col = 0;
    int nnz = 0;
    for (int i = 0; i < (n + 1) * n / 2; ++i)
    {
        if (fabs(dense->dsMatElem[i]) > 0.0)
        {
            nnz += 1;
        }
        row++;
        if (row == n)
        {
            col++;
            row = col;
        }
    }
    double spRatio = (double)nnz / (double)((n + 1) * n / 2);
    row = 0;
    col = 0;
    if (spRatio <= 0.1)
    {
        sdp_coeff_w_sum->dataType = SDP_COEFF_SPARSE;
        sdp_coeff_sparse *sparse;
        ASDP_INIT(sparse, sdp_coeff_sparse, 1);
        ASDP_INIT(sparse->triMatRow, int, nnz);
        ASDP_INIT(sparse->triMatCol, int, nnz);
        ASDP_INIT(sparse->triMatElem, double, nnz);
        sparse->nTriMatElem = nnz;
        sparse->nSDPCol = n;
        ASDP_INIT(sparse->rowCol2NnzIdx, int *, n);
        for (int row = 0; row < n; row++)
        {
            ASDP_INIT(sparse->rowCol2NnzIdx[row], int, row + 1);
            for (int i = 0; i < row + 1; ++i)
            {
                sparse->rowCol2NnzIdx[row][i] = -1;
            }
        }
        int count = 0;
        for (int i = 0; i < (n + 1) * n / 2; ++i)
        {
            if (fabs(dense->dsMatElem[i]) > 0.0)
            {
                sparse->triMatElem[count] = dense->dsMatElem[i];
                sparse->triMatRow[count] = row;
                sparse->triMatCol[count] = col;
                sparse->rowCol2NnzIdx[row][col] = count;
                count++;
            }
            row++;
            if (row == n)
            {
                col++;
                row = col;
            }
        }
        ASDP_FREE(dense->dsMatElem);
        ASDP_FREE(dense);
        sdp_coeff_w_sum->dataMat = (void *)sparse;
        sdp_coeff_w_sum->mul_rk = dataMatSparseMultiRkMat;
        sdp_coeff_w_sum->destroy = destroyForAuxiSparse;
        sdp_coeff_w_sum->zeros = dataMatSparseZeros;
    }
    else
    {
        sdp_coeff_w_sum->dataType = SDP_COEFF_DENSE;
        ASDP_INIT(dense->rowCol2NnzIdx, int *, dense->nSDPCol);
        for (int row = 0; row < dense->nSDPCol; ++row)
        {
            ASDP_INIT(dense->rowCol2NnzIdx[row], int, row + 1);
            for (int i = 0; i < row + 1; ++i)
            {
                dense->rowCol2NnzIdx[row][i] = -1;
            }
        }
        int count = 0;
        for (int col = 0; col < dense->nSDPCol; ++col)
        {
            for (int row = col; row < dense->nSDPCol; ++row)
            {
                dense->rowCol2NnzIdx[row][col] = count;
                count++;
            }
        }
        sdp_coeff_w_sum->destroy = destroyForAuxiDense;
    }
    *sdp_coeff_w_sum_pointer = sdp_coeff_w_sum;
}

static void modifySlackVar(sdp_coeff *w_sum, sdp_coeff **slackVarPointer)
{
    if (w_sum->dataType == SDP_COEFF_DENSE)
    {
        sdp_coeff *slackVar = *slackVarPointer;
        sdp_coeff_dense *denseSlack = (sdp_coeff_dense *)slackVar->dataMat;
        int n = w_sum->nSDPCol;
        ASDP_INIT(denseSlack->rowCol2NnzIdx, int *, n);
        sdp_coeff_dense *w_sumData = (sdp_coeff_dense *)w_sum->dataMat;
        for (int i = 0; i < n; ++i)
        {
            ASDP_INIT(denseSlack->rowCol2NnzIdx[i], int, i + 1);
            ASDP_MEMCPY(denseSlack->rowCol2NnzIdx[i], w_sumData->rowCol2NnzIdx[i], int, i + 1);
        }
        slackVar->destroy = destroyForAuxiDense;
//        slackVar->chol = ASDPPosDefDetermineDense;
        return;
    }
    else if (w_sum->dataType == SDP_COEFF_SPARSE)
    {
        sdp_coeff *slackVar = *slackVarPointer;
        slackVar->dataType = SDP_COEFF_SPARSE;
        sdp_coeff_dense *dense = (sdp_coeff_dense *)slackVar->dataMat;
        sdp_coeff_sparse *sparseSrc = (sdp_coeff_sparse *)w_sum->dataMat;
        sdp_coeff_sparse *sparse;
        int nnz = sparseSrc->nTriMatElem;
        ASDP_INIT(sparse, sdp_coeff_sparse, 1);
        sparse->nTriMatElem = nnz;
        ASDP_INIT(sparse->triMatRow, int, nnz);
        ASDP_INIT(sparse->triMatCol, int, nnz);
        ASDP_INIT(sparse->triMatElem, double, nnz);
        ASDP_MEMCPY(sparse->triMatRow, sparseSrc->triMatRow, int, nnz);
        ASDP_MEMCPY(sparse->triMatCol, sparseSrc->triMatCol, int, nnz);
        ASDP_ZERO(sparse->triMatElem, double, nnz);
        sparse->nSDPCol = sparseSrc->nSDPCol;
        ASDP_INIT(sparse->rowCol2NnzIdx, int *, sparse->nSDPCol);
        for (int row = 0; row < sparse->nSDPCol; ++row)
        {
            ASDP_INIT(sparse->rowCol2NnzIdx[row], int, row + 1);
            ASDP_MEMCPY(sparse->rowCol2NnzIdx[row], sparseSrc->rowCol2NnzIdx[row], int, row + 1);
        }
        slackVar->dataMat = (void *)sparse;
        slackVar->mul_rk = dataMatSparseMultiRkMat;
        slackVar->destroy = destroyForAuxiSparse;
        slackVar->zeros = dataMatSparseZeros;
//        slackVar->chol = ASDPPosDefDetermineSparse;
        *slackVarPointer = slackVar;
        ASDP_FREE(dense->dsMatElem);
        ASDP_FREE(dense);
    }
}

static void AConeFullSDPCoeff_W_sum(sdp_coeff **wSumPointer)
{
    // since w_sum will receive UVt, which is not symmetry
    // Hence, we need a full storage
    sdp_coeff *wSum = *wSumPointer;
    if (wSum->dataType == SDP_COEFF_DENSE)
    {
        sdp_coeff_dense *dense = (sdp_coeff_dense *)wSum->dataMat;
        int nSDPCol = dense->nSDPCol;
        ASDP_FREE(dense->dsMatElem);
        ASDP_INIT(dense->dsMatElem, double, nSDPCol *nSDPCol);
    }
    else if (wSum->dataType == SDP_COEFF_SPARSE)
    {
        sdp_coeff_sparse *sparse = (sdp_coeff_sparse *)wSum->dataMat;
        int nnz = 0;
        int row = 0;
        int col = 0;
        int *stat; // how many nnz in i-th column
        ASDP_INIT(stat, int, sparse->nSDPCol);
        ASDP_ZERO(stat, int, sparse->nSDPCol);
        for (int i = 0; i < sparse->nTriMatElem; ++i)
        {
            row = sparse->triMatRow[i];
            col = sparse->triMatCol[i];
            nnz += 2;
            stat[row] += 1;
            stat[col] += 1;
            if (row == col)
            {
                nnz -= 1;
                stat[row] -= 1;
            }
        }
        int *statCulSum;
        ASDP_INIT(statCulSum, int, sparse->nSDPCol);
        for (int i = 0; i < sparse->nSDPCol - 1; ++i)
        {
            statCulSum[i + 1] = stat[i] + statCulSum[i];
        }

        int *triMatColNew;
        ASDP_INIT(triMatColNew, int, nnz);
        int *triMatRowNew;
        ASDP_INIT(triMatRowNew, int, nnz);
        double *triMatElemNew;
        ASDP_INIT(triMatElemNew, double, nnz);
        double val;
        for (int i = 0; i < sparse->nTriMatElem; ++i)
        {
            row = sparse->triMatRow[i];
            col = sparse->triMatCol[i];
            val = sparse->triMatElem[i];
            triMatColNew[statCulSum[col]] = col;
            triMatRowNew[statCulSum[col]] = row;
            triMatElemNew[statCulSum[col]] = val;
            statCulSum[col] += 1;
            if (row != col)
            {
                triMatColNew[statCulSum[row]] = row;
                triMatRowNew[statCulSum[row]] = col;
                triMatElemNew[statCulSum[row]] = val;
                statCulSum[row] += 1;
            }
        }
        ASDP_FREE(sparse->triMatRow);
        ASDP_FREE(sparse->triMatCol);
        ASDP_FREE(sparse->triMatElem);
        sparse->triMatCol = triMatColNew;
        sparse->triMatRow = triMatRowNew;
        sparse->triMatElem = triMatElemNew;
        ASDP_FREE(statCulSum);
        ASDP_FREE(stat);
    }
    wSum->create = NULL;
    wSum->scal = NULL;
    wSum->norm = NULL;
    wSum->eig = NULL;
    wSum->getnnz = NULL;
    wSum->dump = NULL;
    wSum->getmatnz = NULL;
    wSum->add2buffer = NULL;
    wSum->view = NULL;
    wSum->iseye = NULL;
    wSum->isunitcol = NULL;
    wSum->mul_rk = NULL;
    wSum->mul_inner_rk_double = NULL;
    wSum->zeros = NULL;
    wSum->add_sdp_coeff = NULL;
    wSum->nrm1 = NULL;
    wSum->nrmInf = NULL;
}

static asdp_retcode AConeCreateUVt(sdp_coeff *sum, sdp_coeff **UVtPointer)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    sdp_coeff *UVt = NULL;
    ASDP_INIT(UVt, sdp_coeff, 1);
    ASDP_MEMCHECK(UVt);
    UVt->dataType = sum->dataType;
    UVt->destroy = sum->destroy;
    if (sum->dataType == SDP_COEFF_SPARSE)
    {
        sdp_coeff_sparse *sparse = (sdp_coeff_sparse *)sum->dataMat;
        sdp_coeff_sparse *sparseUVt;
        ASDP_INIT(sparseUVt, sdp_coeff_sparse, 1);
        sparseUVt->nTriMatElem = sparse->nTriMatElem;
        sparseUVt->nSDPCol = sparse->nSDPCol;
        ASDP_INIT(sparseUVt->triMatCol, int, sparseUVt->nTriMatElem);
        ASDP_INIT(sparseUVt->triMatRow, int, sparseUVt->nTriMatElem);
        ASDP_INIT(sparseUVt->triMatElem, double, sparseUVt->nTriMatElem);
        ASDP_MEMCPY(sparseUVt->triMatCol, sparse->triMatCol, int, sparseUVt->nTriMatElem);
        ASDP_MEMCPY(sparseUVt->triMatRow, sparse->triMatRow, int, sparseUVt->nTriMatElem);
        ASDP_MEMCPY(sparseUVt->triMatElem, sparse->triMatElem, double, sparseUVt->nTriMatElem);
        UVt->dataMat = (void *)sparseUVt;
    }
    else if (sum->dataType == SDP_COEFF_DENSE)
    {
        sdp_coeff_dense *dense = (sdp_coeff_dense *)sum->dataMat;
        sdp_coeff_dense *denseUVt;
        ASDP_INIT(denseUVt, sdp_coeff_dense, 1);
        denseUVt->nSDPCol = dense->nSDPCol;
        // Note: we store full matrix in UVt
        ASDP_INIT(denseUVt->dsMatElem, double, denseUVt->nSDPCol * denseUVt->nSDPCol);
        UVt->dataMat = (void *)denseUVt;
    }

    *UVtPointer = UVt;
exit_cleanup:
    return retcode;
}

extern void ASDPDetectSparsityOfSumSDP(asdp *ASolver)
{
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_cone *ACone = ASolver->ACones[iCone];
        AConeDenseDetectSparsity(&ACone->sdp_coeff_w_sum);
        ACone->addObjCoeffRand(ACone->coneData, ACone->sdp_obj_sum);
        AConeDenseDetectSparsity(&ACone->sdp_obj_sum);
        modifySlackVar(ACone->sdp_obj_sum, &ACone->sdp_slack_var);
    }
}

extern asdp_retcode ASDPSumSDPData(asdp *ASolver)
{
    // sum all sdp_coeff in ACone->sdp_coeff_w_sum
    // for detecting the row sparsity of SUM sdp_coeff in one cone
    // for choosing proper calculation in CG
    asdp_retcode retcode = ASDP_RETCODE_OK;
    double *weight;
    ASDP_INIT(weight, double, ASolver->nRows);
    srand(1);
    for (int i = 0; i < ASolver->nRows; ++i)
    {
        weight[i] = (double)rand() / RAND_MAX;
        weight[i] = weight[i] * 100;
#ifdef ASDP_SDPdata_DEBUG
        weight[i] = 1.0;
#endif
    }
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_cone *ACone = ASolver->ACones[iCone];
        ACone->sdpDataWSum(ACone->coneData, weight, ACone->sdp_coeff_w_sum);
        ACone->sdpDataWSum(ACone->coneData, weight, ACone->sdp_obj_sum);
    }
exit_cleanup:
    ASDP_FREE(weight);
    return retcode;
}

extern void ASDPNrm1Obj(asdp *ASolver)
{
    ASolver->cObjNrm1 = 0.0;
    if (ASolver->nLpCols > 0)
    {
        double temp = 0.0;
        asdpLPCone *lpCone = ASolver->lpCone;
        lpCone->coneObjNrm1(lpCone->coneData, &temp, ASolver->nLpCols);
        ASolver->cObjNrm1 += temp;
    }
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        double temp = 0.0;
        asdp_cone *ACone = ASolver->ACones[iCone];
        ACone->coneObjNrm1(ACone->coneData, &temp);
        ASolver->cObjNrm1 += temp;
    }
}

extern void ASDPNrm2Obj(asdp *ASolver)
{
    ASolver->cObjNrm2 = 0.0;
    if (ASolver->nLpCols > 0)
    {
        double temp = 0.0;
        asdpLPCone *lpCone = ASolver->lpCone;
        lpCone->coneObjNrm2Square(lpCone->coneData, &temp, ASolver->nLpCols);
        ASolver->cObjNrm2 += temp;
    }
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        double temp = 0.0;
        asdp_cone *ACone = ASolver->ACones[iCone];
        ACone->coneObjNrm2Square(ACone->coneData, &temp);
        ASolver->cObjNrm2 += temp;
    }
    ASolver->cObjNrm2 = pow(ASolver->cObjNrm2, 0.5);
}

extern void ASDPNrmInfObj(asdp *ASolver)
{
    ASolver->cObjNrmInf = 0.0;
    if (ASolver->nLpCols > 0)
    {
        asdpLPCone *lpCone = ASolver->lpCone;
        lpCone->coneObjNrmInf(lpCone->coneData, &ASolver->cObjNrmInf, ASolver->nLpCols);
    }
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        double temp = 0.0;
        asdp_cone *ACone = ASolver->ACones[iCone];
        ACone->coneObjNrmInf(ACone->coneData, &temp);
        ASolver->cObjNrmInf = ASDP_MAX(ASolver->cObjNrmInf, temp);
    }
}

extern void ASDPSet_w_coeff(asdp *ASolver)
{
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_cone *cone = ASolver->ACones[iCone];
        cone->addObjCoeff(cone->coneData, cone->sdp_coeff_w_sum);
    }
}

static void linSysProduct(asdp_cone *ACone, double *weight, asdp_rk_mat_dense *noUpdateVar, asdp_rk_mat_dense *updateVar, double *x, double *res)
{
    updateVar->matElem = x;
    updateVar->rank = noUpdateVar->rank;
    updateVar->nRows = noUpdateVar->nRows;
    ASDPUpdateConstrValCG(ACone, updateVar, noUpdateVar, weight, FLAG_INI);
    sdp_coeff *w_sum = ACone->sdp_coeff_w_sum;
    w_sum->zeros(w_sum->dataMat);
    ACone->sdpDataWSum(ACone->coneData, weight, w_sum);

    w_sum->mul_rk(w_sum->dataMat, noUpdateVar, res);
    int n = noUpdateVar->nRows * noUpdateVar->rank;
    double alpha = 1.0;
    int incx = 1;
    axpy(&n, &alpha, x, &incx, res, &incx);
}

extern void ADMMUpdateUVMvec(void *pData, double *x, double *res)
{
    admmCG *data = (admmCG *)pData;
    data->UpdateVarShell->matElem = x;
    linSysProduct(data->ACone, data->weight, data->noUpdateVar, data->UpdateVarShell, x, res);
}

extern void ASDPUpdateUV(asdp *ASolver, int iCone, int flagUV)
{
    asdp_cone *ACone = ASolver->ACones[iCone];
    double *b = ASolver->rowRHS;
    double *M1 = ASolver->M1temp;
    ASDP_MEMCPY(M1, b, double, ASolver->nRows);
    asdp_rk_mat_dense *M2 = ASolver->M2temp[iCone];

    // first write for V, hence for V just change pointer
    asdp_rk_mat_dense *noUpdateVar;
    asdp_rk_mat_dense *UpdateVar;
    if (flagUV == FLAG_U)
    {
        noUpdateVar = ASolver->V[iCone];
        UpdateVar = ASolver->U[iCone];
    }
    else if (flagUV == FLAG_V)
    {
        noUpdateVar = ASolver->U[iCone];
        UpdateVar = ASolver->V[iCone];
    }

    asdp_rk_mat_dense *S = ASolver->Vlag[iCone];
    double One = 1.0;
    double minusOne = -1.0;
    int incx = 1;
    // M1 = -b
    scal(&(ASolver->nRows), &(minusOne), M1, &incx);
    // M1 = sum_{j\neq \tilde{j}}A(X) - b
    axpy(&(ASolver->nRows), &(One), ASolver->constrValSum, &incx, M1, &incx);
    //    axpy(&(ASolver->nRows), &(minusOne), ASolver->constrVal[iCone], &incx, M1, &incx);
    ASolver->constrVal[iCone]->add(&minusOne, ASolver->constrVal[iCone]->constrVal, M1);
    // M1 = M1 * rho
    scal(&(ASolver->nRows), &(ASolver->rho), M1, &incx);
    axpy(&(ASolver->nRows), &(minusOne), ASolver->dualVar, &incx, M1, &incx);

    // Add C in the following function
    // and weighted sum of all sdp data in one cone
    ACone->sdp_obj_sum->zeros(ACone->sdp_obj_sum->dataMat);
    ACone->addObjCoeff(ACone->coneData, ACone->sdp_obj_sum);
    ACone->sdpDataWSum(ACone->coneData, M1, ACone->sdp_obj_sum);

    // (C-(lambda + weighted Sum A)) * V
    ACone->sdp_obj_sum->mul_rk(ACone->sdp_obj_sum->dataMat, noUpdateVar, M2->matElem);

    // M2 = (M2 -rho * V + S) or (M2 -rho * V - S)
    ASDPRkMatSub(M2, noUpdateVar, ASolver->rho, S, flagUV);

    // set bv in linear system
    double scalFactor = -1.0 / ASolver->rho;
    int n = M2->nRows * M2->rank;
    ASDP_ZERO(ASolver->bLinSys[iCone], double, n);
    axpy(&(n), &scalFactor, M2->matElem, &incx, ASolver->bLinSys[iCone], &incx);

    admmCG *MMat;
    ASDP_INIT(MMat, admmCG, 1);
    MMat->ACone = ACone;
    MMat->noUpdateVar = noUpdateVar;
    asdp_rk_mat_dense *shell;
    ASDP_INIT(shell, asdp_rk_mat_dense, 1);
    MMat->UpdateVarShell = shell;
    MMat->weight = M1;

    CGSetData(ASolver->CGLinsys[iCone], MMat, ADMMUpdateUVMvec);
    CGSolve(ASolver->CGLinsys[iCone], UpdateVar->matElem, ASolver->bLinSys[iCone], ASolver->dimacError[ASDP_DIMAC_ERROR_CONSTRVIO]);
    ASolver->cgTime += ASolver->CGLinsys[iCone]->cgDuration;
    ASolver->cgIter += ASolver->CGLinsys[iCone]->iter;

    ASDP_FREE(shell);
    ASDP_FREE(MMat);
}

extern void ASDPUpdateUVLP(asdp *ASolver, int iCol, int flagUV)
{
    double one = 1.0;
    double minusOne = -1.0;
    int incx = 1;
    asdpLPCone *lp_cone = ASolver->lpCone;
    double *b = ASolver->rowRHS;
    double *M1 = ASolver->M1temp;
    double zero = 0.0;
    scal(&ASolver->nRows, &zero, M1, &incx);
    axpy(&ASolver->nRows, &one, b, &incx, M1, &incx);
    double *noUpdateVar;
    double *UpdateVar;
    if (flagUV == FLAG_U)
    {
        noUpdateVar = &ASolver->vLp->matElem[iCol];
        UpdateVar = &ASolver->uLp->matElem[iCol];
    }
    else if (flagUV == FLAG_V)
    {
        noUpdateVar = &ASolver->uLp->matElem[iCol];
        UpdateVar = &ASolver->vLp->matElem[iCol];
    }

    double *s = &ASolver->vlagLp->matElem[iCol];
    // M1 = -b
    scal(&(ASolver->nRows), &(minusOne), M1, &incx);
    // M1 = sum_{j \neq \tilde{j}} A(X) - b
    axpy(&(ASolver->nRows), &(one), ASolver->constrValSum, &incx, M1, &incx);
    ASolver->constrValLP[iCol]->add(&minusOne, ASolver->constrValLP[iCol]->constrVal, M1);
    // M1 = M1 * rho
    scal(&(ASolver->nRows), &(ASolver->rho), M1, &incx);
    // M1 = M1 - lambd
    axpy(&(ASolver->nRows), &(minusOne), ASolver->dualVar, &incx, M1, &incx);

    double lpWObjSum = 0.0;
    lp_cone->objCoeffSum(lp_cone->coneData, &lpWObjSum, iCol);
    lp_cone->lpDataWSum(lp_cone->coneData, M1, &lpWObjSum, iCol);

    double M2 = lpWObjSum * noUpdateVar[0];
    if (flagUV == FLAG_U)
    {
        M2 = M2 - ASolver->rho * noUpdateVar[0] + s[0];
    }
    else if (flagUV == FLAG_V)
    {
        M2 = M2 - ASolver->rho * noUpdateVar[0] - s[0];
    }
    double blinSys = -1.0 * M2 / ASolver->rho;
    asdp_cone_lp *lpConeData = (asdp_cone_lp *)lp_cone->coneData;
    UpdateVar[0] = blinSys / (1 + lpConeData->nrm2Square[iCol] * noUpdateVar[0] * noUpdateVar[0]);
}

extern void ASDPUpdateS(asdp *ASolver, int iCone)
{
    asdp_cone *ACone = ASolver->ACones[iCone];
    asdp_rk_mat_dense *U = ASolver->U[iCone];
    asdp_rk_mat_dense *V = ASolver->V[iCone];
    asdp_rk_mat_dense *S = ASolver->Vlag[iCone];
    double rho = ASolver->rho;
    double minusRho = -rho;
    int rank = U->rank;
    int nRows = U->nRows;
    int n = rank * nRows;

    // S = S + rho * U
    int one = 1;
    axpy(&(n), &(rho), U->matElem, &(one), S->matElem, &(one));
    // S = S - rho * V
    axpy(&(n), &(minusRho), V->matElem, &(one), S->matElem, &(one));
}

extern void ASDPUpdateSLP(asdp *ASolver, int iCol)
{
    asdp_rk_mat_lp *s = ASolver->vlagLp;
    asdp_rk_mat_lp *u = ASolver->uLp;
    asdp_rk_mat_lp *v = ASolver->vLp;

    s->matElem[iCol] += ASolver->rho * u->matElem[iCol];
    s->matElem[iCol] -= ASolver->rho * v->matElem[iCol];
}

extern void ASDPUpdateDualVar(asdp *ASolver)
{
    double *dualVar = ASolver->dualVar;
    double rho = ASolver->rho;
    double minusRho = -rho;
    double *b = ASolver->rowRHS;
    double *constrSum = ASolver->constrValSum;
    int n = ASolver->nRows;

    // lambda = lambda + rho * b
    int incx = 1;
    axpy(&(n), &(rho), b, &(incx), dualVar, &(incx));
    // lambda = lambda - rho * constrVal
    axpy(&(n), &(minusRho), constrSum, &(incx), dualVar, &(incx));
}

static void MvecDense(void *MMat, double *xVec, double *yVec)
{
    sdp_coeff_dense *dense = (sdp_coeff_dense *)MMat;
    double alpha = -1.0;
    fds_symv_L(dense->nSDPCol, alpha, dense->dsMatElem, xVec, 0.0, yVec);
}

static void MvecSparse(void *MMat, double *xVec, double *yVec)
{
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *)MMat;
    int row = 0;
    int col = 0;
    double val;
    for (int i = 0; i < sparse->nTriMatElem; ++i)
    {
        row = sparse->triMatRow[i];
        col = sparse->triMatCol[i];
        val = sparse->triMatElem[i];
        yVec[row] -= val * xVec[col];
    }
}

extern asdp_retcode ASDPUpdateDimacError(asdp *ASolver)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    int nConstr = ASolver->nRows;
    ASDPInitConstrValAll(ASolver, ASolver->rLp, ASolver->rLp, ASolver->R, ASolver->R);
    ASDPInitConstrValSum(ASolver);
    double *constrVio = ASolver->constrVio;
    double one = 1.0;
    double minusOne = -1.0;
    axpbyAddition(&nConstr, &(one), ASolver->rowRHS, &(minusOne), ASolver->constrValSum, constrVio);
    int incx = 1;
    double diffNorm = 0;
    double rNorm = 0;
    double negOne = -1.0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone){
        int n = ASolver->R[iCone]->nRows * ASolver->R[iCone]->rank;
        axpy(&n, &negOne, ASolver->R[iCone]->matElem, &incx, ASolver->lag[iCone]->matElem, &incx);
        diffNorm += nrm2(&n, ASolver->lag[iCone]->matElem, &incx);
        rNorm += nrm2(&n, ASolver->R[iCone]->matElem, &incx);
    }
    ASolver->dimacError[ASDP_ERROR_DIFF_GAP] = ASolver->rho * diffNorm / (1 + rNorm);
    ASolver->dimacError[ASDP_DIMAC_ERROR_CONSTRVIO] = nrm2(&nConstr, constrVio, &incx) / (1 + ASolver->bRHSNrm1);
    if (ASolver->whichMethod == ADMMMethod)
    {
        double minUVnorm = INFINITY;
        double sumUVDis = 0.0;
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            asdp_rk_mat_dense *U = ASolver->U[iCone];
            asdp_rk_mat_dense *V = ASolver->V[iCone];

            int rank = U->rank;
            int nRows = U->nRows;
            int n = rank * nRows;
            double *UsubV = ASolver->UsubV[iCone];

            double Unorm = nrm2(&n, U->matElem, &incx);
            double Vnorm = nrm2(&n, V->matElem, &incx);
            minUVnorm = ASDP_MIN(minUVnorm, Unorm);
            minUVnorm = ASDP_MIN(minUVnorm, Vnorm);

            axpbyAddition(&n, &one, U->matElem, &minusOne, V->matElem, UsubV);
            sumUVDis += nrm2(&n, UsubV, &incx);
        }
        ASolver->dimacError[ASDP_DIMAC_ERROR_ASSYMMETRY] = sumUVDis / (1 + minUVnorm);
    }

    double gap = (ASolver->pObjVal - ASolver->dObjVal);
    ASolver->dimacError[ASDP_DIMAC_ERROR_PDGAP] = ASDP_ABS(gap) / (1 + ASDP_ABS(ASolver->pObjVal) + ASDP_ABS(ASolver->dObjVal));

    double *negLambd = ASolver->negLambd;
    ASDP_ZERO(negLambd, double, ASolver->nRows);
    axpy(&(ASolver->nRows), &minusOne, ASolver->dualVar, &incx, negLambd, &incx);

    ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] = 1.0;
    if (ASolver->dimacError[ASDP_DIMAC_ERROR_PDGAP] >= 1e-5 || ASolver->dimacError[ASDP_DIMAC_ERROR_CONSTRVIO] >= 1e-5)
    {
        goto exit_cleanup;
    }
#ifndef TERMINATION_SCHEME_OFFICIAL
    ;
#else
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
        if (Cone->sdp_slack_var->dataType == SDP_COEFF_DENSE)
        {
            ALanczosSetData(ASolver->LanczoSys[iCone], Cone->sdp_slack_var->dataMat, MvecDense);
        }
        else if (Cone->sdp_slack_var->dataType == SDP_COEFF_SPARSE)
        {
            ALanczosSetData(ASolver->LanczoSys[iCone], Cone->sdp_slack_var->dataMat, MvecSparse);
        }
        double dMaxStep[1] = {1.0};
        ALanczosSolve(ASolver->LanczoSys[iCone], ASolver->LanczosStart[iCone], dMaxStep);
        minEig = -1 * ASolver->LanczoSys[iCone]->dArray[1];
        double dualVio = ASDP_ABS(ASDP_MIN(minEig, 0));
        ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] += dualVio;
    }
    ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] /= (ASolver->cObjNrm1 + 1);
#endif
exit_cleanup:
    return retcode;
}

extern void ASDPNuclearNorm(asdp *ASolver)
{
    ASolver->traceSum = 0;
    int rank = 0;
    int nRows = 0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_rk_mat_dense *R = ASolver->R[iCone];
        rank = R->rank;
        nRows = R->nRows;
        for (int idx = 0; idx < nRows; ++idx)
        {
            ASolver->traceSum += dot(&rank, &R->matElem[idx], &nRows, &R->matElem[idx], &nRows);
        }
    }
}

extern asdp_retcode ASDPUpdateDimacErrorLP(asdp *ASolver)
{
    asdp_retcode retcode = ASDP_RETCODE_OK;
    int nConstr = ASolver->nRows;
    ASDPInitConstrValAllLP(ASolver, ASolver->rLp, ASolver->rLp, ASolver->R, ASolver->R);
    ASDPInitConstrValSumLP(ASolver);
    double *constrVio = ASolver->constrVio;
    double one = 1.0;
    double minusOne = -1.0;
    axpbyAddition(&nConstr, &one, ASolver->rowRHS, &minusOne, ASolver->constrValSum, constrVio);
    int incx = 1;
    ASolver->dimacError[ASDP_DIMAC_ERROR_CONSTRVIO] = nrm2(&nConstr, constrVio, &incx) / (1 + ASolver->bRHSNrm1);
    double diffNorm = 0;
    double rNorm = 0;
    double negOne = -1.0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone){
        int n = ASolver->R[iCone]->nRows * ASolver->R[iCone]->rank;
        axpy(&n, &negOne, ASolver->R[iCone]->matElem, &incx, ASolver->Vlag[iCone]->matElem, &incx);
        diffNorm += nrm2(&n, ASolver->Vlag[iCone]->matElem, &incx);
        rNorm += nrm2(&n, ASolver->R[iCone]->matElem, &incx);
    }
    ASolver->dimacError[ASDP_ERROR_DIFF_GAP] = ASolver->rho * diffNorm / (1 + rNorm);
    if (ASolver->nLpCols > 0){
        axpy(&ASolver->nLpCols, &negOne, ASolver->rLp->matElem, &incx, ASolver->vlagLp->matElem, &incx);
        diffNorm += nrm2(&ASolver->nLpCols, ASolver->vlagLp->matElem, &incx);
        rNorm += nrm2(&ASolver->nLpCols, ASolver->rLp->matElem, &incx);
    }
    ASolver->dimacError[ASDP_ERROR_DIFF_GAP] = ASolver->rho * diffNorm / rNorm;
    if (ASolver->whichMethod == ADMMMethod)
    {
        double minUVnorm = INFINITY;
        double sumUVDis = 0.0;

        double *uSubvLp = ASolver->uSubvLp;
        int idxMin = idamin(&(ASolver->uLp->nLPCols), ASolver->uLp->matElem, &incx);
        minUVnorm = ASDP_MIN(fabs(ASolver->uLp->matElem[idxMin]), minUVnorm);
        idxMin = idamin(&(ASolver->vLp->nLPCols), ASolver->vLp->matElem, &incx);
        minUVnorm = ASDP_MIN(fabs(ASolver->vLp->matElem[idxMin]), minUVnorm);

        axpbyAddition(&(ASolver->vLp->nLPCols), &one, ASolver->vLp->matElem, &minusOne, ASolver->uLp->matElem, uSubvLp);
        sumUVDis += nrm2(&(ASolver->vLp->nLPCols), uSubvLp, &incx);

        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            asdp_rk_mat_dense *U = ASolver->U[iCone];
            asdp_rk_mat_dense *V = ASolver->V[iCone];

            int rank = U->rank;
            int nRows = U->nRows;
            int n = rank * nRows;
            double *UsubV = ASolver->UsubV[iCone];

            double Unorm = nrm2(&n, U->matElem, &incx);
            double Vnorm = nrm2(&n, V->matElem, &incx);
            minUVnorm = ASDP_MIN(minUVnorm, Unorm);
            minUVnorm = ASDP_MIN(minUVnorm, Vnorm);

            axpbyAddition(&n, &one, U->matElem, &minusOne, V->matElem, UsubV);
            sumUVDis += nrm2(&n, UsubV, &incx);
        }
        ASolver->dimacError[ASDP_DIMAC_ERROR_ASSYMMETRY] = sumUVDis / (1 + minUVnorm);
    }
    double gap = (ASolver->pObjVal - ASolver->dObjVal);
    ASolver->dimacError[ASDP_DIMAC_ERROR_PDGAP] = ASDP_ABS(gap) / (1 + ASDP_ABS(ASolver->pObjVal) + ASDP_ABS(ASolver->dObjVal));

    double *negLambd = ASolver->negLambd;
    ASDP_ZERO(negLambd, double, ASolver->nRows);
    axpy(&(ASolver->nRows), &minusOne, ASolver->dualVar, &incx, negLambd, &incx);

    ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] = 1.0;
    if (ASolver->dimacError[ASDP_DIMAC_ERROR_PDGAP] >= 1e-5 || ASolver->dimacError[ASDP_DIMAC_ERROR_CONSTRVIO] >= 1e-5)
    {
        goto exit_cleanup;
    }

    ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] = 0.0;
    double lpObjWlpDataSum = 0.0;
    asdpLPCone *lp_cone = ASolver->lpCone;
    for (int iCol = 0; iCol < ASolver->nLpCols; ++iCol)
    {
        lpObjWlpDataSum = 0.0;
        lp_cone->objCoeffSum(lp_cone->coneData, &lpObjWlpDataSum, iCol);
        lp_cone->lpDataWSum(lp_cone->coneData, negLambd, &lpObjWlpDataSum, iCol);
        ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] += ASDP_ABS(ASDP_MIN(lpObjWlpDataSum, 0));
    }

    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_rk_mat_dense *U = ASolver->U[iCone];
        int nRows = U->nRows;

        asdp_cone *Cone = ASolver->ACones[iCone];
        Cone->sdp_slack_var->zeros(Cone->sdp_slack_var->dataMat);
        Cone->addObjCoeff(Cone->coneData, Cone->sdp_slack_var);
        Cone->sdpDataWSum(Cone->coneData, negLambd, Cone->sdp_slack_var);

        double minEig = 0.0;
        if (Cone->sdp_slack_var->dataType == SDP_COEFF_DENSE)
        {
            ALanczosSetData(ASolver->LanczoSys[iCone], Cone->sdp_slack_var->dataMat, MvecDense);
        }
        else if (Cone->sdp_slack_var->dataType == SDP_COEFF_SPARSE)
        {
            ALanczosSetData(ASolver->LanczoSys[iCone], Cone->sdp_slack_var->dataMat, MvecSparse);
        }
        double dMaxStep[1] = {1.0};
        ALanczosSolve(ASolver->LanczoSys[iCone], ASolver->LanczosStart[iCone], dMaxStep);
        minEig = -1 * ASolver->LanczoSys[iCone]->dArray[1];
        double dualVio = ASDP_ABS(ASDP_MIN(minEig, 0));
        ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] += dualVio;
    }
    ASolver->dimacError[ASDP_DIMAC_ERROR_DUALFEASIBLE] /= (ASolver->cObjNrm1 + 1);
exit_cleanup:
    return retcode;
}

extern void ASDPUpdateSDPVar(asdp *ASolver)
{
    int incx = 1;
    double minusOne = -1.0;
    double one = 1.0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        // update U
        ASDPUpdateUV(ASolver, iCone, FLAG_U);
        // ASolver->constrValSum - ASolver->constrVal[iCone]
        ASolver->constrVal[iCone]->add(&minusOne, ASolver->constrVal[iCone]->constrVal, ASolver->constrValSum);
        // update ASolver->constrVal[iCone]
        ASDPUpdateConstrVal(ASolver->ACones[iCone], ASolver->U[iCone], ASolver->V[iCone], ASolver->constrVal[iCone], FLAG_INI);
        // update ASolver->constrValSum
        ASolver->constrVal[iCone]->add(&one, ASolver->constrVal[iCone]->constrVal, ASolver->constrValSum);

        // update V
        ASDPUpdateUV(ASolver, iCone, FLAG_V);
        ASolver->constrVal[iCone]->add(&minusOne, ASolver->constrVal[iCone]->constrVal, ASolver->constrValSum);
        ASDPUpdateConstrVal(ASolver->ACones[iCone], ASolver->U[iCone], ASolver->V[iCone], ASolver->constrVal[iCone], FLAG_INI);
        ASolver->constrVal[iCone]->add(&one, ASolver->constrVal[iCone]->constrVal, ASolver->constrValSum);

        // update S
#ifdef NO_PENALTY_METHOD
        ASDPUpdateS(ASolver, iCone);
#endif
    }
}

extern void ASDPUpdateSDPLPVar(asdp *ASolver)
{
    int incx = 1;
    double minusOne = -1.0;
    double one = 1.0;
    ASDPUpdateSDPVar(ASolver);
    for (int iCol = 0; iCol < ASolver->nLpCols; ++iCol)
    {
        ASDPUpdateUVLP(ASolver, iCol, FLAG_U);
        ASolver->constrValLP[iCol]->add(&minusOne, ASolver->constrValLP[iCol]->constrVal, ASolver->constrValSum);
        ASDPUpdateConstrValLP(ASolver->lpCone, ASolver->uLp, ASolver->vLp, ASolver->constrValLP[iCol], iCol);
        ASolver->constrValLP[iCol]->add(&one, ASolver->constrValLP[iCol]->constrVal, ASolver->constrValSum);

        ASDPUpdateUVLP(ASolver, iCol, FLAG_V);
        ASolver->constrValLP[iCol]->add(&minusOne, ASolver->constrValLP[iCol]->constrVal, ASolver->constrValSum);
        ASDPUpdateConstrValLP(ASolver->lpCone, ASolver->uLp, ASolver->vLp, ASolver->constrValLP[iCol], iCol);
        ASolver->constrValLP[iCol]->add(&one, ASolver->constrValLP[iCol]->constrVal, ASolver->constrValSum);

#ifdef NO_PENALTY_METHOD
        ASDPUpdateSLP(ASolver, iCol);
#endif
    }
}

extern void BMSetGrad(asdp *ASolver, asdp_cone *ACone, asdp_rk_mat_dense *R, asdp_rk_mat_dense *Grad, int iCone)
{
    ASDP_ZERO(Grad->matElem, double, Grad->nRows * Grad->rank);
    double *M1 = ASolver->M1temp;
    int n = ASolver->nRows;
    int incx = 1;
    if (iCone == 0)
    {
        ASDP_ZERO(M1, double, n);
        // M1 = -lambda + rho * (constrVal - b)
        /*for (int iConstr = 0; iConstr < ASolver->nRows; ++iConstr){
            M1[iConstr] = -ASolver->dualVar[iConstr] + ASolver->rho * (ASolver->rowRHS[iConstr] - ASolver->constrValSum[iConstr]);
        }*/
        double minusOne = -1.0;
        axpy(&n, &minusOne, ASolver->dualVar, &incx, M1, &incx);
        double rho = ASolver->rho;
        double negRho = -1 * rho;
        axpy(&n, &negRho, ASolver->rowRHS, &incx, M1, &incx);
        axpy(&n, &rho, ASolver->constrValSum, &incx, M1, &incx);
    }

    ACone->sdp_obj_sum->zeros(ACone->sdp_obj_sum->dataMat);
    ACone->addObjCoeff(ACone->coneData, ACone->sdp_obj_sum);
    ACone->sdpDataWSum(ACone->coneData, M1, ACone->sdp_obj_sum);

    // Note: obj mul_rk no row sparse method
    ACone->sdp_obj_sum->mul_rk(ACone->sdp_obj_sum->dataMat, R, Grad->matElem);
    int len = Grad->nRows * Grad->rank;
    double TWO = 2.0;
    scal(&len, &TWO, Grad->matElem, &incx);
}

extern void BMSetGradLP(asdp *ASolver, asdpLPCone *lp_cone, asdp_rk_mat_lp *rLp, asdp_rk_mat_lp *Grad, int iCol)
{
    Grad->matElem[iCol] = 0.0;
    double *M1 = ASolver->M1temp;
    int n = ASolver->nRows;
    int incx = 1;

    if (iCol == 0)
    {
        if (n > 2000)
        {
            double zero = 0.0;
            scal(&n, &zero, M1, &incx);
        }
        else
        {
            ASDP_ZERO(M1, double, n);
        }

        // M1 = -labmd - rho * b + rho * (sum A(uv))
        double minusOne = -1.0;
        axpy(&n, &minusOne, ASolver->dualVar, &incx, M1, &incx);
        double rho = ASolver->rho;
        double negRho = -1 * rho;
        axpy(&n, &negRho, ASolver->rowRHS, &incx, M1, &incx);
        axpy(&n, &rho, ASolver->constrValSum, &incx, M1, &incx);
    }

    double lpObjWlpDataSum = 0.0;
    lp_cone->objCoeffSum(lp_cone->coneData, &lpObjWlpDataSum, iCol);
    lp_cone->lpDataWSum(lp_cone->coneData, M1, &lpObjWlpDataSum, iCol);
    Grad->matElem[iCol] = 2 * lpObjWlpDataSum * rLp->matElem[iCol];
}

extern void BMCalGrad(asdp *ASolver, asdp_rk_mat_lp *rLpDummy, asdp_rk_mat_lp *gradLpDummy, asdp_rk_mat_dense **R, asdp_rk_mat_dense **Grad, double *lagNormSquare)
{
    int incx = 1;
    double lagNorm2 = 0.0;
    lagNormSquare[0] = 0.0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        // S receive gradient info
        BMSetGrad(ASolver, ASolver->ACones[iCone], R[iCone], Grad[iCone], iCone);
        int n = Grad[iCone]->nRows * Grad[iCone]->rank;
        lagNorm2 = nrm2(&n, Grad[iCone]->matElem, &incx);
        lagNormSquare[0] += lagNorm2 * lagNorm2;
    }
}

extern void BMCalGradLP(asdp *ASolver, asdp_rk_mat_lp *rLp, asdp_rk_mat_lp *gradLp, asdp_rk_mat_dense **R, asdp_rk_mat_dense **Grad, double *lagNormSquare)
{
    int incx = 1;
    double lagNorm2 = 0.0;
    lagNormSquare[0] = 0.0;
    for (int iCol = 0; iCol < ASolver->nLpCols; ++iCol)
    {
        BMSetGradLP(ASolver, ASolver->lpCone, rLp, gradLp, iCol);
    }
    lagNorm2 = nrm2(&(gradLp->nLPCols), gradLp->matElem, &incx);
    lagNormSquare[0] += lagNorm2 * lagNorm2;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        // S receive gradient info
        BMSetGrad(ASolver, ASolver->ACones[iCone], R[iCone], Grad[iCone], iCone);
        int n = Grad[iCone]->nRows * Grad[iCone]->rank;
        lagNorm2 = nrm2(&n, Grad[iCone]->matElem, &incx);
        lagNormSquare[0] += lagNorm2 * lagNorm2;
    }
}

extern void LBFGSDirection(asdp *ASolver, lbfgs_node *head, asdp_rk_mat_lp *gradLpDummy, asdp_rk_mat_lp *dDummy, asdp_rk_mat_dense **Grad, asdp_rk_mat_dense **D, int innerIter)
{
    /*head is oldest info, update head node and move head pointer into next*/
    int incx = 1;
    int n = 0;
    double minusOne = -1.0;
    double *Dtemp;
    int idx = 0;
    if (innerIter == 0)
    {
        // D = -grad
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            n = Grad[iCone]->nRows * Grad[iCone]->rank;
            ASDP_ZERO(D[iCone]->matElem, double, n);
            axpy(&n, &minusOne, Grad[iCone]->matElem, &incx, D[iCone]->matElem, &incx);
        }
        return;
    }
    else
    {
        ASDP_INIT(Dtemp, double, head->allElem);
        idx = 0;
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            n = Grad[iCone]->nRows * Grad[iCone]->rank;
            ASDP_MEMCPY(&Dtemp[idx], Grad[iCone]->matElem, double, n);
            idx += n;
        }
    }

    int NodeNum = 0;
    if (innerIter <= ASolver->hisRecT + 1)
    {
        NodeNum = innerIter;
    }
    else
    {
        NodeNum = ASolver->hisRecT + 1;
    }
    lbfgs_node *node = head->prev;
    for (int k = 0; k < NodeNum; ++k)
    {
        // new to old, head -> previous is newst
        // calculate alpha
        double temp = 0.0;
        temp = dot(&(head->allElem), node->s, &incx, Dtemp, &incx);

        node->alpha = node->beta * temp;
        double alphaNeg = -1 * node->alpha;
        axpy(&(head->allElem), &alphaNeg, node->y, &incx, Dtemp, &incx);
        node = node->prev;
    }
    node = node->next;
    for (int k = 0; k < NodeNum; ++k)
    {
        double weight = node->alpha - node->beta * dot(&(head->allElem), node->y, &incx, Dtemp, &incx);
        axpy(&(head->allElem), &weight, node->s, &incx, Dtemp, &incx);
        node = node->next;
    }
    scal(&(head->allElem), &minusOne, Dtemp, &incx);
    idx = 0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        n = Grad[iCone]->nRows * Grad[iCone]->rank;
        ASDP_MEMCPY(D[iCone]->matElem, &Dtemp[idx], double, n);
        idx += n;
    }
    ASDP_FREE(Dtemp);
    return;
}

extern void LBFGSDirectionLP(asdp *ASolver, lbfgs_node *head, asdp_rk_mat_lp *gradLp, asdp_rk_mat_lp *d, asdp_rk_mat_dense **Grad, asdp_rk_mat_dense **D, int innerIter)
{
    /*head is oldest info, update head node and move head pointer into next*/
    int incx = 1;
    int n = 0;
    double minusOne = -1.0;
    double *Dtemp;
    int idx = 0;
    if (innerIter == 0)
    {
        // D = -grad
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            n = Grad[iCone]->nRows * Grad[iCone]->rank;
            ASDP_ZERO(D[iCone]->matElem, double, n);
            axpy(&n, &minusOne, Grad[iCone]->matElem, &incx, D[iCone]->matElem, &incx);
        }
        ASDP_ZERO(d->matElem, double, d->nLPCols);
        n = d->nLPCols;
        axpy(&n, &minusOne, gradLp->matElem, &incx, d->matElem, &incx);
        return;
    }
    else
    {
        // Dtemp = grad
        ASDP_INIT(Dtemp, double, head->allElem);
        idx = 0;
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            n = Grad[iCone]->nRows * Grad[iCone]->rank;
            ASDP_MEMCPY(&Dtemp[idx], Grad[iCone]->matElem, double, n);
            idx += n;
        }
        ASDP_MEMCPY(&Dtemp[idx], gradLp->matElem, double, gradLp->nLPCols);
    }

    int NodeNum = 0;
    if (innerIter <= ASolver->hisRecT + 1)
    {
        NodeNum = innerIter;
    }
    else
    {
        NodeNum = ASolver->hisRecT + 1;
    }
    lbfgs_node *node = head->prev;
    for (int k = 0; k < NodeNum; ++k)
    {
        double temp = 0.0;
        temp = dot(&(head->allElem), node->s, &incx, Dtemp, &incx);
        node->alpha = node->beta * temp;
        double alphaNeg = -1 * node->alpha;
        axpy(&(head->allElem), &alphaNeg, node->y, &incx, Dtemp, &incx);
        node = node->prev;
    }
    node = node->next;
    for (int k = 0; k < NodeNum; ++k)
    {
        double weight = node->alpha - node->beta * dot(&(head->allElem), node->y, &incx, Dtemp, &incx);
        axpy(&(head->allElem), &weight, node->s, &incx, Dtemp, &incx);
        node = node->next;
    }
    scal(&(head->allElem), &minusOne, Dtemp, &incx);
    idx = 0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        n = Grad[iCone]->nRows * Grad[iCone]->rank;
        ASDP_MEMCPY(D[iCone]->matElem, &Dtemp[idx], double, n);
        idx += n;
    }
    ASDP_MEMCPY(d->matElem, &Dtemp[idx], double, d->nLPCols);
    idx += d->nLPCols;
    ASDP_FREE(Dtemp);
    return;
}

extern void LBFGSDirectionUseGrad(asdp *ASolver, asdp_rk_mat_lp *dDummy, asdp_rk_mat_lp *gradLPDummy, asdp_rk_mat_dense **D, asdp_rk_mat_dense **Grad)
{
    double innerProduct = 0.0;
    int nElem = 0;
    int incx = 1;
    double minusOne = -1.0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        nElem = D[iCone]->nRows * D[iCone]->rank;
        innerProduct += dot(&nElem, D[iCone]->matElem, &incx, Grad[iCone]->matElem, &incx);
    }
    if (innerProduct >= 0)
    {
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            nElem = D[iCone]->nRows * D[iCone]->rank;
            ASDP_MEMCPY(D[iCone]->matElem, Grad[iCone]->matElem, double, nElem);
            scal(&nElem, &minusOne, D[iCone]->matElem, &incx);
        }
    }
}

extern void LBFGSDirectionUseGradLP(asdp *ASolver, asdp_rk_mat_lp *d, asdp_rk_mat_lp *gradLP, asdp_rk_mat_dense **D, asdp_rk_mat_dense **Grad)
{
    double innerProduct = 0.0;
    int nElem = 0;
    int incx = 1;
    double minusOne = -1.0;
    innerProduct += dot(&(gradLP->nLPCols), d->matElem, &incx, gradLP->matElem, &incx);
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        nElem = D[iCone]->nRows * D[iCone]->rank;
        innerProduct += dot(&nElem, D[iCone]->matElem, &incx, Grad[iCone]->matElem, &incx);
    }
    if (innerProduct >= 0)
    {
        ASDP_MEMCPY(d->matElem, gradLP->matElem, double, gradLP->nLPCols);
        scal(&(gradLP->nLPCols), &minusOne, d->matElem, &incx);
        for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
        {
            nElem = D[iCone]->nRows * D[iCone]->rank;
            ASDP_MEMCPY(D[iCone]->matElem, Grad[iCone]->matElem, double, nElem);
            scal(&nElem, &minusOne, D[iCone]->matElem, &incx);
        }
    }
}

extern void SetyAsNegGrad(asdp *ASolver, asdp_rk_mat_lp *gradLpDummy, asdp_rk_mat_dense **Grad)
{
    double minusOne = -1.0;
    int incx = 1;
    int n = 0;
    int idx = 0;
    lbfgs_node *head = ASolver->lbfgsHis;
    ASDP_ZERO(head->y, double, head->allElem);
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_rk_mat_dense *gradICone = Grad[iCone];
        n = gradICone->nRows * gradICone->rank;
        axpy(&n, &minusOne, gradICone->matElem, &incx, &head->y[idx], &incx);
        idx += n;
    }
}

extern void SetyAsNegGradLP(asdp *ASolver, asdp_rk_mat_lp *gradLp, asdp_rk_mat_dense **Grad)
{
    double minusOne = -1.0;
    int incx = 1;
    int n = 0;
    int idx = 0;
    lbfgs_node *head = ASolver->lbfgsHis;
    ASDP_ZERO(head->y, double, head->allElem);
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_rk_mat_dense *gradICone = Grad[iCone];
        n = gradICone->nRows * gradICone->rank;
        axpy(&n, &minusOne, gradICone->matElem, &incx, &head->y[idx], &incx);
        idx += n;
    }
    axpy(&gradLp->nLPCols, &minusOne, gradLp->matElem, &incx, &head->y[idx], &incx);
    // idx += gradLp->nLPCols;
}

extern void BMupdateVar(asdp *ASolver, asdp_rk_mat_lp *rLpDummy, asdp_rk_mat_lp *dDummy, asdp_rk_mat_dense **R, asdp_rk_mat_dense **D, double tau)
{
    int incx = 1;
    int nElem = 0;
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        nElem = R[iCone]->nRows * R[iCone]->rank;
        axpy(&nElem, &tau, D[iCone]->matElem, &incx, R[iCone]->matElem, &incx);
    }
}

extern void BMupdateVarLP(asdp *ASolver, asdp_rk_mat_lp *rLp, asdp_rk_mat_lp *d, asdp_rk_mat_dense **R, asdp_rk_mat_dense **D, double tau)
{
    int incx = 1;
    BMupdateVar(ASolver, rLp, d, R, D, tau);
    axpy(&(rLp->nLPCols), &tau, d->matElem, &incx, rLp->matElem, &incx);
}

extern void setlbfgsHisTwo(asdp *ASolver, asdp_rk_mat_lp *gradLpDummy, asdp_rk_mat_lp *dDummy, asdp_rk_mat_dense **Grad, asdp_rk_mat_dense **D, double tau)
{
    double minusOne = -1.0;
    double One = 1.0;
    int incx = 1;
    int n = 0;
    int idx = 0;
    lbfgs_node *head = ASolver->lbfgsHis;
    ASDP_ZERO(head->s, double, head->allElem);
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_rk_mat_dense *gradICone = Grad[iCone];
        n = gradICone->nRows * gradICone->rank;
        // head->s = tau * D
        axpy(&n, &tau, D[iCone]->matElem, &incx, &head->s[idx], &incx);
        // head->y += GradNew
        axpy(&n, &One, Grad[iCone]->matElem, &incx, &head->y[idx], &incx);
        idx += n;
    }
    head->beta = 1.0 / dot(&(head->allElem), head->y, &incx, head->s, &incx);
    ASolver->lbfgsHis = head->next;
}

extern void setlbfgsHisTwoLP(asdp *ASolver, asdp_rk_mat_lp *gradLp, asdp_rk_mat_lp *d, asdp_rk_mat_dense **Grad, asdp_rk_mat_dense **D, double tau)
{
    double minusOne = -1.0;
    double One = 1.0;
    int incx = 1;
    int n = 0;
    int idx = 0;
    lbfgs_node *head = ASolver->lbfgsHis;
    ASDP_ZERO(head->s, double, head->allElem);
    for (int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        asdp_rk_mat_dense *gradICone = Grad[iCone];
        n = gradICone->nRows * gradICone->rank;
        // head->s = tau * D
        axpy(&n, &tau, D[iCone]->matElem, &incx, &head->s[idx], &incx);
        // head->y += GradNew
        axpy(&n, &One, Grad[iCone]->matElem, &incx, &head->y[idx], &incx);
        idx += n;
    }
    axpy(&(d->nLPCols), &tau, d->matElem, &incx, &head->s[idx], &incx);
    axpy(&(d->nLPCols), &One, gradLp->matElem, &incx, &head->y[idx], &incx);
    // idx += d->nLPCols;
    head->beta = 1.0 / dot(&(head->allElem), head->y, &incx, head->s, &incx);
    //    if (fabs(head->beta) > 1e+30){
    //        if (head->beta > 0){
    //            head->beta = 10;
    //        }else{
    //            head->beta = -10;
    //        }
    //
    //    }
    ASolver->lbfgsHis = head->next;
}

extern void copyRtoV(asdp_rk_mat_lp *rLpDummy, asdp_rk_mat_lp *vlpDummy, asdp_rk_mat_dense **R, asdp_rk_mat_dense **V, int nCones)
{
    int n = 0;
    for (int iCone = 0; iCone < nCones; ++iCone)
    {
        n = R[iCone]->rank * R[iCone]->nRows;
        ASDP_MEMCPY(V[iCone]->matElem, R[iCone]->matElem, double, n);
    }
}

extern void copyRtoVLP(asdp_rk_mat_lp *rLp, asdp_rk_mat_lp *vlp, asdp_rk_mat_dense **R, asdp_rk_mat_dense **V, int nCones)
{
    int n = 0;
    ASDP_MEMCPY(vlp->matElem, rLp->matElem, double, vlp->nLPCols);
    for (int iCone = 0; iCone < nCones; ++iCone)
    {
        n = R[iCone]->rank * R[iCone]->nRows;
        ASDP_MEMCPY(V[iCone]->matElem, R[iCone]->matElem, double, n);
    }
}

extern int CheckAllRankMax(asdp *asolver, double aug_factor)
{
    int *rank_flags = malloc(asolver->nCones * sizeof(int));
    ASDP_ZERO(rank_flags, int, asolver->nCones);
    for (int iCone = 0; iCone < asolver->nCones; ++iCone)
    {
        int new_rank = ASDP_MIN(ceil(asolver->U[iCone]->rank * aug_factor), asolver->rank_max[iCone]);
        if (new_rank >= asolver->rank_max[iCone])
        {
            rank_flags[iCone] = 1;
        }
    }
    int sum_flag = 0;
    for (int i = 0; i < asolver->nCones; i++)
    {
        sum_flag += rank_flags[i];
    }
    free(rank_flags);
    return sum_flag == asolver->nCones;
}


#include <math.h>
#include "stdlib.h"
#include "def_lorads_solver.h"
#include "lorads.h"
#include "lorads_utils.h"
#include "lorads_sdp_conic.h"
#include "lorads_user_data.h"
#include "lorads_lp_conic.h"
#include "lorads_elements.h"
#include "lorads_vec_opts.h"
#include "lorads_cgs.h"
#include "lorads_alm.h"
#include "lorads_alg_common.h"
#include "lorads_admm.h"

extern void LORADSInitSolver(lorads_solver *ASolver, lorads_int nRows, lorads_int nCones, lorads_int *blkDims, lorads_int nLpCols)
{
    ASolver->nLpCols = nLpCols;
    // nRows: number of constraints
    // nCones: number of blks

    ASolver->nRows = nRows;   // constralorads_int number
    ASolver->nCones = nCones; // blk number

    LORADS_INIT(ASolver->rowRHS, double, nRows);
    LORADS_MEMCHECK(ASolver->rowRHS);

    // set dual variable
    LORADS_INIT(ASolver->var->dualVar, double, nRows);
    LORADS_ZERO(ASolver->var->dualVar, double, nRows);
    LORADS_MEMCHECK(ASolver->var->dualVar);

    LORADS_INIT(ASolver->var->bestDualVar, double, nRows);
    LORADS_ZERO(ASolver->var->bestDualVar, double, nRows);
    LORADS_MEMCHECK(ASolver->var->bestDualVar);

    // set Auxiliary variable
    LORADS_INIT(ASolver->var->constrVal, lorads_vec *, nCones);
    LORADS_MEMCHECK(ASolver->var->constrVal);
    for (lorads_int iCone = 0; iCone < nCones; ++iCone)
    {
        LORADS_INIT(ASolver->var->constrVal[iCone], lorads_vec, 1);
        LORADS_MEMCHECK(ASolver->var->constrVal[iCone]);
    }

    LORADS_INIT(ASolver->var->constrValSum, double, nRows);
    LORADS_MEMCHECK(ASolver->var->constrValSum);

    LORADS_INIT(ASolver->SDPCones, lorads_sdp_cone *, nCones); // allocate n cones
    LORADS_MEMCHECK(ASolver->SDPCones);

    for (lorads_int iCone = 0; iCone < nCones; ++iCone)
    {
        LORADS_INIT(ASolver->SDPCones[iCone], lorads_sdp_cone, 1);
    }
    double dim = (double)nLpCols;
    for (lorads_int iCone = 0; iCone < nCones; ++iCone)
    {
        dim += (double)blkDims[iCone];
    }

    LORADS_INIT(ASolver->dimacError, double, 5);
    LORADS_MEMCHECK(ASolver->dimacError);
    LORADS_ZERO(ASolver->dimacError, double, 4);

    if (ASolver->nLpCols > 0)
    {
        LORADS_INIT(ASolver->var->constrValLP, lorads_vec *, ASolver->nLpCols);
        LORADS_MEMCHECK(ASolver->var->constrValLP);
        for (lorads_int iLpCol = 0; iLpCol < ASolver->nLpCols; ++iLpCol)
        {
            LORADS_INIT(ASolver->var->constrValLP[iLpCol], lorads_vec, 1);
            LORADS_MEMCHECK(ASolver->var->constrValLP[iLpCol]);
        }
    }
    LORADS_INIT(ASolver->lpCone, lorads_lp_cone, 1); // allocate n cones
    LORADS_MEMCHECK(ASolver->lpCone);
}

extern void LORADSDestroySolver(lorads_solver *ASolver)
{
    // corresponding to `LORADSInitSolver`
    LORADS_FREE(ASolver->rowRHS);
    LORADS_FREE(ASolver->var->dualVar);
    LORADS_FREE(ASolver->var->bestDualVar);
    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        LORADS_FREE(ASolver->SDPCones[iCone]);
        if (ASolver->var->constrVal[iCone]->type == LORADS_DENSE_VEC){
            LORADS_FREE(ASolver->var->constrVal[iCone]);
        }
        else if (ASolver->var->constrVal[iCone]->type == LORADS_SPARSE_VEC)
        {
            LORADS_FREE(ASolver->var->constrVal[iCone]);
        }
    }
    LORADS_FREE(ASolver->var->constrVal);
    LORADS_FREE(ASolver->var->constrValSum);
    LORADS_FREE(ASolver->SDPCones);
    LORADS_FREE(ASolver->dimacError);
    if (ASolver->nLpCols > 0)
    {
        for  (lorads_int iLPCol = 0; iLPCol < ASolver->nLpCols; ++iLPCol)
        {
            if (ASolver->var->constrValLP[iLPCol]->type == LORADS_DENSE_VEC)
            {
                dense_vec *dense = (dense_vec *)ASolver->var->constrValLP[iLPCol]->data;
                LORADS_FREE(dense->val);
                LORADS_FREE(ASolver->var->constrValLP[iLPCol]);
            }
            else if (ASolver->var->constrValLP[iLPCol]->type == LORADS_SPARSE_VEC)
            {
                sparse_vec *sparse = (sparse_vec *)ASolver->var->constrValLP[iLPCol]->data;
                LORADS_FREE(sparse->nnzIdx);
                LORADS_FREE(sparse->val);
                LORADS_FREE(ASolver->var->constrValLP[iLPCol]);
            }
            LORADS_FREE(ASolver->var->constrVal);
        }
        LORADS_FREE(ASolver->var->constrValLP);
        LORADS_FREE(ASolver->lpCone);
    }
}

extern void LORADSSetDualObjective(lorads_solver *ASolver, double *dObj){
    LORADS_MEMCPY(ASolver->rowRHS, dObj, double, ASolver->nRows);
    return;
}

extern void LORADSInitConeData(lorads_solver *ASolver, user_data **SDPDatas,
                               double **coneMatElem, lorads_int **coneMatBeg, lorads_int **coneMatIdx,
                               lorads_int *BlkDims, lorads_int nConstrs, lorads_int nBlks,
                               lorads_int nLpCols, lorads_int *LpMatBeg, lorads_int *LpMatIdx, double *LpMatElem)
{
    if (nLpCols > 0){
        LORADSSetLpCone(ASolver->lpCone, ASolver->nRows, nLpCols, LpMatBeg, LpMatIdx, LpMatElem);
    }
    user_data *SDPData = NULL;
    for  (lorads_int iCone = 0; iCone < nBlks; ++iCone){
        LUserDataCreate(&SDPData); // initialize, set all equals to 0
        LUserDataSetConeData(SDPData, LORADS_CONETYPE_DENSE_SDP, nConstrs, BlkDims[iCone],
                             coneMatBeg[iCone], coneMatIdx[iCone], coneMatElem[iCone]);
        LORADSSetCone(ASolver, iCone, SDPData);
        SDPDatas[iCone] = SDPData;
        SDPData = NULL;
    }
}

extern void LORADSNrm1Obj(lorads_solver *ASolver)
{
    ASolver->cObjNrm1 = 0.0;
    if (ASolver->nLpCols > 0){
        double temp = 0.0;
        lorads_lp_cone *lpCone = ASolver->lpCone;
        lpCone->coneObjNrm1(lpCone->coneData, &temp, ASolver->nLpCols);
        ASolver->cObjNrm1 += temp;
    }
    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        double temp = 0.0;
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        ACone->coneObjNrm1(ACone->coneData, &temp);
        ASolver->cObjNrm1 += temp;
    }
}


extern void LORADSNrm2Obj(lorads_solver *ASolver)
{
    ASolver->cObjNrm2 = 0.0;
    if (ASolver->nLpCols > 0)
    {
        double temp = 0.0;
        lorads_lp_cone *lpCone = ASolver->lpCone;
        lpCone->coneObjNrm2Square(lpCone->coneData, &temp, ASolver->nLpCols);
        ASolver->cObjNrm2 += temp;
    }
    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        double temp = 0.0;
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        ACone->coneObjNrm2Square(ACone->coneData, &temp);
        ASolver->cObjNrm2 += temp;
    }
    ASolver->cObjNrm2 = pow(ASolver->cObjNrm2, 0.5);
}



extern void LORADSPreprocess(lorads_solver *ASolver, lorads_int *BlkDims)
{
    /* Start optimization */
    ASolver->dTimeBegin = LUtilGetTimeStamp();

    /* Process conic data */
    printf("Pre-solver starts \n");
    printf("  Processing the cones \n");
    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        AConeProcData(ASolver->SDPCones[iCone]);                     // set data
        AConePresolveData(ASolver->SDPCones[iCone], BlkDims[iCone]); // detect rank structure
    }
    printf("  End preprocess \n");

    // Detect `sum of sdp coeff` and `sum of sdp coeff and obj` strcuture
//    LORADSSumSDPData(ASolver);
//    LORADSDetectSparsityOfSumSDP(ASolver);


    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        if (ASolver->SDPCones[iCone]->type == LORADS_CONETYPE_DENSE_SDP){
            dense_vec *dense;
            LORADS_INIT(dense, dense_vec, 1);
            dense->nnz = ASolver->nRows;
            LORADS_INIT(dense->val, double, ASolver->nRows);
            ASolver->var->constrVal[iCone]->data = (void *)dense;
            ASolver->var->constrVal[iCone]->add = addDense;
            ASolver->var->constrVal[iCone]->zero = zeroDense;
            ASolver->var->constrVal[iCone]->type = LORADS_DENSE_VEC;
        }else if (ASolver->SDPCones[iCone]->type == LORADS_CONETYPE_SPARSE_SDP){
            sparse_vec *sparse;
            LORADS_INIT(sparse, sparse_vec , 1);
            lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)ASolver->SDPCones[iCone]->coneData;
            LORADS_INIT(sparse->val, double, cone->nRowElem);
            sparse->nnz = cone->nRowElem;
            LORADS_INIT(sparse->nnzIdx, lorads_int, cone->nRowElem);
            LORADS_MEMCPY(sparse->nnzIdx, cone->rowIdx, lorads_int, cone->nRowElem);
            ASolver->var->constrVal[iCone]->data = (void *)sparse;
            ASolver->var->constrVal[iCone]->add = addSparse;
            ASolver->var->constrVal[iCone]->zero = zeroSparse;
            ASolver->var->constrVal[iCone]->type = LORADS_SPARSE_VEC;
        }
    }

    if (ASolver->nLpCols > 0){
        lorads_lp_cone_data *lpCone = ASolver->lpCone->coneData;
        for  (lorads_int iCol = 0; iCol < ASolver->nLpCols; ++iCol){
            if (lpCone->lpCol[iCol]->dataType == LP_COEFF_DENSE){
                dense_vec *dense;
                LORADS_INIT(dense, dense_vec, 1);
                dense->nnz = ASolver->nRows;
                LORADS_INIT(dense->val, double, ASolver->nRows);
                ASolver->var->constrValLP[iCol]->data = (void *)dense;
                ASolver->var->constrValLP[iCol]->add = addDense;
                ASolver->var->constrValLP[iCol]->zero = zeroDense;
                ASolver->var->constrValLP[iCol]->type = LORADS_DENSE_VEC;
            }else if (lpCone->lpCol[iCol]->dataType == LP_COEFF_SPARSE){
                sparse_vec *sparse;
                LORADS_INIT(sparse, sparse_vec, 1);
                lp_coeff_sparse *lpCol = (lp_coeff_sparse *)lpCone->lpCol[iCol]->dataMat;
                LORADS_INIT(sparse->val, double, lpCol->nnz);
                sparse->nnz = lpCol->nnz;
                LORADS_INIT(sparse->nnzIdx, lorads_int, lpCol->nnz);
                LORADS_MEMCPY(sparse->nnzIdx, lpCol->rowPtr, lorads_int, lpCol->nnz);
                ASolver->var->constrValLP[iCol]->data = (void *)sparse;
                ASolver->var->constrValLP[iCol]->add = addSparse;
                ASolver->var->constrValLP[iCol]->zero = zeroSparse;
                ASolver->var->constrValLP[iCol]->type = LORADS_SPARSE_VEC;
            }else if (lpCone->lpCol[iCol]->dataType == LP_COEFF_ZERO){
                LORADS_ERROR_TRACE;
            }
        }
    }
}

extern void LORADSDestroyConeData(lorads_solver *ASolver)
{
    // destroy for function LORADSInitConeData_and_UVS
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        // SDPdata pointer is set as usrData here
        LORADS_FREE(ASolver->SDPCones[iCone]->usrData);
    }
    if (ASolver->nLpCols > 0)
    {
        ASolver->lpCone->destroyConeData(&ASolver->lpCone->coneData);
    }
    LORADS_FREE(ASolver->lpCone);
}

extern void destroyPreprocess(lorads_solver *ASolver)
{
    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        AConeDestroyPresolveData(ASolver->SDPCones[iCone]);
        ASolver->SDPCones[iCone]->coneDestroyData(&ASolver->SDPCones[iCone]->coneData);
    }
    LORADS_FREE(ASolver->constrVio);
}

extern void LORADSDetermineRank(lorads_solver *ASolver, lorads_int *blkDims, double timesRank)
{
    lorads_int nCones = ASolver->nCones;
    lorads_int *rankElem;
    LORADS_INIT(rankElem, lorads_int, nCones);
    LORADS_MEMCHECK(rankElem);
    LORADS_INIT(ASolver->rank_max, lorads_int, nCones);
    LORADS_MEMCHECK(ASolver->rank_max);
    for  (lorads_int iCone = 0; iCone < nCones; ++iCone)
    {
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        lorads_int nnzRows = 0;
        ACone->nnzStat(ACone->coneData, &nnzRows);
        if (timesRank <= 1e-6)
        {
            rankElem[iCone] = LORADS_MIN( (lorads_int)sqrt(2 * nnzRows) + 1, blkDims[iCone]);
        }
        else if (nnzRows / blkDims[iCone] >= 20 && blkDims[iCone] <= 400 && nCones <= 3){
            rankElem[iCone] = LORADS_MIN( (lorads_int)sqrt(2 * nnzRows) + 1, blkDims[iCone]);
        }
        else
        {
            rankElem[iCone] = LORADS_MIN(ceil(timesRank * log(blkDims[iCone])), LORADS_MIN( (lorads_int)sqrt(2 * nnzRows) + 1, blkDims[iCone]));
        }
        rankElem[iCone] = LORADS_MAX(1, rankElem[iCone]);

        ASolver->rank_max[iCone] = LORADS_MIN( (lorads_int)sqrt(2 * nnzRows) + 1, blkDims[iCone]);
    }
    ASolver->var->rankElem = rankElem;
}

extern void detectMaxCutProb(lorads_solver *ASolver, lorads_int *BlkDims, lorads_int *maxCut){
    maxCut[0] = -1;
    lorads_int stat = 0;
    lorads_sdp_cone *cone = ASolver->SDPCones[0];

    if (cone->type == LORADS_CONETYPE_DENSE_SDP)
    {
        lorads_cone_sdp_dense *denseCone = (lorads_cone_sdp_dense *)cone->coneData;
        for (lorads_int iRow = 0; iRow < denseCone->nRow; ++iRow)
        {
            sdp_coeff *sdpCoeff = denseCone->sdpRow[iRow];
            if (sdpCoeff->dataType == SDP_COEFF_SPARSE){
                sdp_coeff_sparse *sdpSparse = sdpCoeff->dataMat;
                if (sdpSparse->nTriMatElem == 1){
                    if (sdpSparse->triMatRow[0] == sdpSparse->triMatCol[0]){
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

extern void detectSparsitySDPCoeff(lorads_solver *ASolver)
{
    LORADS_INIT(ASolver->sparsitySDPCoeff, double, ASolver->nCones * ASolver->nRows);
    ASolver->nnzSDPCoeffSum = 0;
    ASolver->SDPCoeffSum = 0;
    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        ACone->nnzStatCoeff(ACone->coneData, &ASolver->sparsitySDPCoeff[iCone * ASolver->nRows], &ASolver->nnzSDPCoeffSum, &ASolver->SDPCoeffSum);
    }
    ASolver->overallSparse = (double)ASolver->nnzSDPCoeffSum / (double)ASolver->SDPCoeffSum / (double)ASolver->nRows;
}

extern void LORADS_RANDOM_rk_MAT(lorads_sdp_dense *U)
{
    lorads_int n = U->nRows * U->rank;
    LORADS_INIT(U->matElem, double, n);
    LORADS_MEMCHECK(U->matElem);
    for  (lorads_int i = 0; i < n; ++i)
    {
        U->matElem[i] = (double)rand() / RAND_MAX;
        U->matElem[i] -= (double)rand() / RAND_MAX;
    }
}

extern void LORADS_ONE_rk_MAT(lorads_sdp_dense *U){
    lorads_int n = U->nRows * U->rank;
    LORADS_INIT(U->matElem, double, n);
    LORADS_MEMCHECK(U->matElem);
    for  (lorads_int i = 0; i < n; ++i){
        U->matElem[i] = 1.0;
    }
}

extern void lpRandom(double *data, lorads_int n)
{
    for (lorads_int i = 0; i < n; ++i)
    {
        data[i] = (double)rand() / RAND_MAX;
        data[i] -= (double)rand() / RAND_MAX;
    }
}

extern void lpFix(double *data, lorads_int n)
{
    for (lorads_int i = 0; i < n; ++i)
    {
        if (i == 0)
        {
            data[i] = 1.0;
        }
        else
        {
        data[i] = 0.0;
        }
    }
}

extern void LORADSInitALMVars(lorads_solver *ASolver, lorads_int *rankElem, lorads_int *BlkDims, lorads_int nBlks, lorads_int nLpCols, lorads_int lbfgsHis)
{
    LORADS_INIT(ASolver->var->rLp, lorads_lp_dense, 1);
    LORADS_MEMCHECK(ASolver->var->rLp);
    ASolver->var->rLp->nCols = nLpCols;
    LORADS_INIT(ASolver->var->gradLp, lorads_lp_dense, 1);
    LORADS_MEMCHECK(ASolver->var->gradLp);
    ASolver->var->gradLp->nCols = nLpCols;

    srand(925);
    LORADS_INIT(ASolver->var->ADDSum, double, ASolver->nRows);
    LORADS_MEMCHECK(ASolver->var->ADDSum);
    LORADS_INIT(ASolver->var->ARDSum, double, ASolver->nRows);
    LORADS_MEMCHECK(ASolver->var->ARDSum);
    LORADS_INIT(ASolver->var->M1temp, double, ASolver->nRows);
    LORADS_MEMCHECK(ASolver->var->M1temp);
    // R is variable of ALMALM
    // R is also the average variable for ADMM split method
    LORADS_INIT(ASolver->var->R, lorads_sdp_dense *, nBlks);
    LORADS_MEMCHECK(ASolver->var->R);
    LORADS_INIT(ASolver->var->Grad, lorads_sdp_dense *, nBlks);
    LORADS_MEMCHECK(ASolver->var->Grad);
    lorads_int allElem = 0;
    if (nLpCols > 0){
        allElem += nLpCols;
    }
    for  (lorads_int iCone = 0; iCone < nBlks; ++iCone){
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        // set variable R
        lorads_sdp_dense *R;
        LORADS_INIT(R, lorads_sdp_dense, 1);
        LORADS_MEMCHECK(R);
        R->rank = rankElem[iCone];
        R->nRows = BlkDims[iCone];

#ifdef FIX_INI_POINT
        LORADS_ONE_rk_MAT(R);
#else
        LORADS_RANDOM_rk_MAT(R);
#endif
        ASolver->var->R[iCone] = R;
        allElem += R->nRows * R->rank;

        lorads_sdp_dense *Grad;
        LORADS_INIT(Grad, lorads_sdp_dense, 1);
        LORADS_MEMCHECK(Grad);
        Grad->rank = rankElem[iCone];
        Grad->nRows = BlkDims[iCone];
        lorads_int n = Grad->nRows * Grad->rank;
        LORADS_INIT(Grad->matElem, double, n);
        ASolver->var->Grad[iCone] = Grad;
    }
    if (nLpCols > 0)
    {
        LORADS_INIT(ASolver->var->rLp->matElem, double, nLpCols);
        LORADS_MEMCHECK(ASolver->var->rLp->matElem);
#ifdef FIX_INI_POINT
        lpFix(ASolver->var->rLp->matElem, nLpCols);
#else
        lpRandom(ASolver->var->rLp->matElem, nLpCols);
#endif
        LORADS_INIT(ASolver->var->gradLp->matElem, double, nLpCols);
        LORADS_MEMCHECK(ASolver->var->gradLp->matElem);
    }
    LORADS_INIT(ASolver->lbfgsHis, lbfgs_node, 1);
    LORADS_MEMCHECK(ASolver->lbfgsHis);
    // ALL  (hisRecT + 1) node
    lbfgs_node *head = ASolver->lbfgsHis;
    LORADS_INIT(head->s, double, allElem);
    LORADS_INIT(head->y, double, allElem);
    head->allElem = allElem;
    lbfgs_node *node = head;
    LORADS_INIT(ASolver->var->Dtemp, double, allElem);
    for  (lorads_int nodeNum = 0; nodeNum < lbfgsHis - 1; ++nodeNum)
    {
        lbfgs_node *nextNode;
        LORADS_INIT(nextNode, lbfgs_node, 1);
        LORADS_INIT(nextNode->s, double, allElem);
        LORADS_MEMCHECK(nextNode->s);
        LORADS_INIT(nextNode->y, double, allElem);
        LORADS_MEMCHECK(nextNode->y);
        nextNode->allElem = allElem;
        node->next = nextNode;
        nextNode->prev = node;
        node = node->next;
        if (nodeNum == lbfgsHis - 2)
        {
            node->next = head;
            head->prev = node;
        }
    }
    ASolver->lbfgsHis = head;
}

extern void LORADSDestroyALMVars(lorads_solver *ASolver)
{
    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        LORADS_FREE(ASolver->var->R[iCone]->matElem);
        LORADS_FREE(ASolver->var->R[iCone]);
        LORADS_FREE(ASolver->var->Grad[iCone]->matElem);
        LORADS_FREE(ASolver->var->Grad[iCone]);
    }
    lbfgs_node *node = ASolver->lbfgsHis;
    lbfgs_node *nextNode;
    for  (lorads_int nodeNum = 0; nodeNum < ASolver->hisRecT - 1; ++nodeNum)
    {
        LORADS_FREE(node->s);
        LORADS_FREE(node->y);
        nextNode = node->next;
        LORADS_FREE(node);
        node = nextNode;
    }
    LORADS_FREE(ASolver->var->Dtemp);
    LORADS_FREE(nextNode->s);
    LORADS_FREE(nextNode->y);
    LORADS_FREE(nextNode);
    LORADS_FREE(ASolver->var->R);
    LORADS_FREE(ASolver->var->Grad);
    if (ASolver->nLpCols > 0)
    {
        LORADS_FREE(ASolver->var->rLp->matElem);
    }
    LORADS_FREE(ASolver->var->rLp);
    LORADS_FREE(ASolver->var->gradLp);
    LORADS_FREE(ASolver->var->M1temp);
    LORADS_FREE(ASolver->var->ADDSum);
    LORADS_FREE(ASolver->var->ARDSum);
}

extern double normalRandom()
{
    double u1 = (double)rand() / RAND_MAX; // uniform [0, 1]
    double u2 = (double)rand() / RAND_MAX;

    // Box-Muller transform
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);

    return z;
}

extern void LORADSCGCreate(lorads_solver *ASolver)
{
    LORADS_INIT(ASolver->CGLinsys, lorads_cg_linsys *, ASolver->nCones);
    LORADS_MEMCHECK(ASolver->CGLinsys);
    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        LORADSCGSolverCreate(&ASolver->CGLinsys[iCone], ASolver->var->U[iCone]->nRows, ASolver->var->U[iCone]->rank, ASolver->nRows);
    }
}

extern void LORADSInitM1M2Temp(lorads_solver *ASolver)
{
    LORADS_INIT(ASolver->var->M2temp, lorads_sdp_dense *, ASolver->nCones);
    LORADS_MEMCHECK(ASolver->var->M2temp);
    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        LORADS_INIT(ASolver->var->M2temp[iCone], lorads_sdp_dense, 1);
        LORADS_MEMCHECK(ASolver->var->M2temp[iCone]);
        ASolver->var->M2temp[iCone]->rank = ASolver->var->U[iCone]->rank;
        ASolver->var->M2temp[iCone]->nRows = ASolver->var->U[iCone]->nRows;
        LORADS_INIT(ASolver->var->M2temp[iCone]->matElem, double, ASolver->var->M2temp[iCone]->rank * ASolver->var->M2temp[iCone]->nRows);
        LORADS_MEMCHECK(ASolver->var->M2temp[iCone]->matElem);
    }
}

extern void LORADSInitAuxiCri(lorads_solver *ASolver)
{
    lorads_int nConstr = ASolver->nRows;
    LORADS_INIT(ASolver->constrVio, double, nConstr);
    LORADS_MEMCHECK(ASolver->constrVio);
}


extern void LORADSInitADMMVars(lorads_solver *ASolver, lorads_int *rankElem, lorads_int *BlkDims, lorads_int nBlks, lorads_int nLpCols)
{
    LORADS_INIT(ASolver->var->uLp, lorads_lp_dense, 1);
    LORADS_MEMCHECK(ASolver->var->uLp);
    ASolver->var->uLp->nCols = nLpCols;
    LORADS_INIT(ASolver->var->vLp, lorads_lp_dense, 1);
    LORADS_MEMCHECK(ASolver->var->vLp);
    ASolver->var->vLp->nCols = nLpCols;
#ifdef DUAL_U_V
    LORADS_INIT(ASolver->var->sLp, lorads_lp_dense, 1);
    LORADS_MEMCHECK(ASolver->var->sLp);
    ASolver->var->sLp->nCols = nLpCols;
#endif
    if (nLpCols > 0)
    {
        LORADS_INIT(ASolver->var->uLp->matElem, double, nLpCols);
        LORADS_MEMCHECK(ASolver->var->uLp->matElem);
        LORADS_INIT(ASolver->var->vLp->matElem, double, nLpCols);
        LORADS_MEMCHECK(ASolver->var->vLp->matElem);
#ifdef FIX_INI_POINT
        lpFix(ASolver->var->uLp->matElem, nLpCols);
        lpFix(ASolver->var->vLp->matElem, nLpCols);
#else
        lpRandom(ASolver->var->uLp->matElem, nLpCols);
        lpRandom(ASolver->var->vLp->matElem, nLpCols);
#endif
#ifdef DUAL_U_V
        LORADS_INIT(ASolver->var->sLp->matElem, double, nLpCols);
        LORADS_MEMCHECK(ASolver->var->sLp->matElem);
        lpRandom(ASolver->var->sLp->matElem, nLpCols);
#endif

    }

    LORADS_INIT(ASolver->var->bLinSys, double *, nBlks);
    LORADS_MEMCHECK(ASolver->var->bLinSys);
    for  (lorads_int iCone = 0; iCone < nBlks; ++iCone)
    {
        LORADS_INIT(ASolver->var->bLinSys[iCone], double, BlkDims[iCone] * rankElem[iCone]);
        LORADS_MEMCHECK(ASolver->var->bLinSys[iCone]);
        LORADS_ZERO(ASolver->var->bLinSys[iCone], double, BlkDims[iCone] * rankElem[iCone]);
    }
    LORADS_INIT(ASolver->var->V, lorads_sdp_dense *, nBlks);
    LORADS_MEMCHECK(ASolver->var->V);
    LORADS_INIT(ASolver->var->U, lorads_sdp_dense *, nBlks);
    LORADS_MEMCHECK(ASolver->var->U);
#ifdef DUAL_U_V
    LORADS_INIT(ASolver->var->S, lorads_sdp_dense *, nBlks);
    LORADS_MEMCHECK(ASolver->var->S);
#endif

    for  (lorads_int iCone = 0; iCone < nBlks; ++iCone)
    {
        // set variable
        lorads_sdp_dense *U;
        LORADS_INIT(U, lorads_sdp_dense, 1);
        LORADS_MEMCHECK(U);
        U->rank = rankElem[iCone];
        U->nRows = BlkDims[iCone];

        lorads_sdp_dense *V;
        LORADS_INIT(V, lorads_sdp_dense, 1);
        LORADS_MEMCHECK(V);
        V->rank = rankElem[iCone];
        V->nRows = BlkDims[iCone];

#ifdef FIX_INI_POINT
        LORADS_ONE_rk_MAT(U);
        LORADS_ONE_rk_MAT(V);
//        LORADS_ONE_rk_MAT(S);
//        LORADS_ONE_rk_MAT(lag);
#else
        LORADS_RANDOM_rk_MAT(U);
        LORADS_RANDOM_rk_MAT(V);
#endif

        lorads_int incx = 1;
        ASolver->var->U[iCone] = U;
        ASolver->var->V[iCone] = V;
#ifdef DUAL_U_V
        lorads_sdp_dense *S;
        LORADS_INIT(S, lorads_sdp_dense, 1);
        LORADS_MEMCHECK(S);
        S->rank = rankElem[iCone];
        S->nRows = BlkDims[iCone];
        LORADS_INIT(S->matElem, double, S->rank * S->nRows);
        LORADS_MEMCHECK(S->matElem);
        LORADS_ZERO(S->matElem, double, S->rank * S->nRows);
        ASolver->var->S[iCone] = S;
#endif
    }

    LORADSCGCreate(ASolver);
    LORADSInitM1M2Temp(ASolver);
    LORADSInitAuxiCri(ASolver);
}

extern void LORADSDestroyADMMVars(lorads_solver *ASolver)
{
    // destroy for function LORADSInitConeData_and_UVS
    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        LORADS_FREE(ASolver->var->M2temp[iCone]->matElem);
        LORADS_FREE(ASolver->var->M2temp[iCone]);
        LORADS_FREE(ASolver->var->U[iCone]->matElem);
        LORADS_FREE(ASolver->var->U[iCone]);
        LORADS_FREE(ASolver->var->V[iCone]->matElem);
        LORADS_FREE(ASolver->var->V[iCone]);
#ifdef DUAL_U_V
        LORADS_FREE(ASolver->var->S[iCone]->matElem);
        LORADS_FREE(ASolver->var->S[iCone]);
#endif
        LORADS_FREE(ASolver->var->bLinSys[iCone]);
    }
    if (ASolver->nLpCols > 0)
    {
        LORADS_FREE(ASolver->var->uLp->matElem);
        LORADS_FREE(ASolver->var->vLp->matElem);
#ifdef DUAL_U_V
        LORADS_FREE(ASolver->var->sLp->matElem);
#endif
    }
    LORADS_FREE(ASolver->var->uLp);
    LORADS_FREE(ASolver->var->vLp);
#ifdef DUAL_U_V
    LORADS_FREE(ASolver->var->sLp);
#endif
    LORADS_FREE(ASolver->var->U);
    LORADS_FREE(ASolver->var->V);
#ifdef DUAL_U_V
    LORADS_FREE(ASolver->var->S);
#endif
    LORADS_FREE(ASolver->var->bLinSys);
    LORADS_FREE(ASolver->var->M2temp);
    LORADSCGDestroy(ASolver);
}

extern void LORADSInitFuncSet(lorads_func **pfunc, lorads_int nLpCols){
    lorads_func *func;
    LORADS_INIT(func, lorads_func, 1);
    if (nLpCols > 0){
        // function with lp cone
        func->ALMCalGrad = ALMCalGradLP;
        func->ALMCalq12p12 = ALMCalq12p12LP;
        func->InitConstrValAll = LORADSInitConstrValAllLP;
        func->InitConstrValSum = LORADSInitConstrValSumLP;
        func->LBFGSDirUseGrad = LBFGSDirectionUseGradLP;
        func->LBFGSDirection = LBFGSDirectionLP;
        func->admmUpdateVar = LORADSUpdateSDPLPVar;
        func->calObj_admm = LORADSCalObjUV_ADMM_LP;
        func->calObj_alm = LORADSCalObjRR_ALM_LP;
        func->copyRtoV = copyRtoVLP;
        func->setAsNegGrad = SetyAsNegGradLP;
        func->setlbfgsHisTwo = setlbfgsHisTwoLP;
        func->updateDimacsALM = LORADSUpdateDimacsErrorALMLP;
        func->updateDimacsADMM = LORADSUpdateDimacsErrorADMMLP;
        func->ALMupdateVar = ALMupdateVarLP;
    }else{
        // function with sdp cone only
        func->ALMCalGrad = ALMCalGrad;
        func->ALMCalq12p12 = ALMCalq12p12;
        func->InitConstrValAll = LORADSInitConstrValAll;
        func->InitConstrValSum = LORADSInitConstrValSum;
        func->LBFGSDirUseGrad = LBFGSDirectionUseGrad;
        func->LBFGSDirection = LBFGSDirection;
        func->admmUpdateVar = LORADSUpdateSDPVar;
        func->calObj_admm = LORADSCalObjUV_ADMM;
        func->calObj_alm = LORADSCalObjRR_ALM;
        func->copyRtoV = copyRtoV;
        func->setAsNegGrad = SetyAsNegGrad;
        func->setlbfgsHisTwo = setlbfgsHisTwo;
        func->updateDimacsALM = LORADSUpdateDimacsErrorALM;
        func->updateDimacsADMM = LORADSUpdateDimacsErrorADMM;
        func->ALMupdateVar = ALMupdateVar;
    }
    *pfunc = func;
}

extern lorads_int CheckAllRankMax(lorads_solver *asolver, double aug_factor)
{
    lorads_int *rank_flags = malloc(asolver->nCones * sizeof (lorads_int));
    LORADS_ZERO(rank_flags, lorads_int, asolver->nCones);
    for  (lorads_int iCone = 0; iCone < asolver->nCones; ++iCone){
        lorads_int new_rank = LORADS_MIN(ceil(asolver->var->U[iCone]->rank * aug_factor), asolver->rank_max[iCone]);
        if (new_rank >= asolver->rank_max[iCone]){
            rank_flags[iCone] = 1;
        }
    }
    lorads_int sum_flag = 0;
    for  (lorads_int i = 0; i < asolver->nCones; i++){
        sum_flag += rank_flags[i];
    }
    free(rank_flags);
    return sum_flag == asolver->nCones;
}

extern void lpRandomDiag(double *data, lorads_int n, lorads_int nRows, lorads_int nCols)
{
    if (n == 0){
        return;
    }
    lorads_int r = LORADS_MIN(nRows, nCols);
    LORADS_ZERO(data, double, n);
    for (lorads_int i = 0; i < r; ++i){
        data[i * nRows + i] = 1 / sqrt(r);
    }
}

extern void  LORADSCGReAllocate(lorads_solver *ASolver){
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
         LORADSCGSolverReCreate(&ASolver->CGLinsys[iCone], ASolver->var->U[iCone]->nRows, ASolver->var->U[iCone]->rank, ASolver->nRows);
    }
}

static void printinfo_aug(lorads_int iCone){
#ifdef LORADS_INT32
    printf("**Rank truncated to sqrt(2m) on SDP Cone No.%d.\n", iCone);
#endif
#ifdef UNIX_INT64
    printf("**Rank truncated to sqrt(2m) on SDP Cone No.%ld.\n", iCone);
#endif
#ifdef MAC_INT64
    printf("**Rank truncated to sqrt(2m) on SDP Cone No.%lld.\n", iCone);
#endif
}

extern lorads_int AUG_RANK(lorads_solver *ASolver, lorads_int *BlkDims, lorads_int nBlks, double aug_factor)
{
    // double rank

    lorads_int is_rank_max = CheckAllRankMax(ASolver, 1.0);
    // printf("AUG_RANK rank max? %d\n", is_rank_max);

    if (is_rank_max)
    {
        // printf("rank max!\n");
        return is_rank_max;
    }
    lorads_int allElem = 0;
    if (ASolver->nLpCols > 0){
        allElem += ASolver->nLpCols;
    }
    for  (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        lorads_int new_rank = LORADS_MIN(ceil(ASolver->var->U[iCone]->rank * aug_factor), ASolver->rank_max[iCone]);
        lorads_int aug_rank = new_rank - ASolver->var->U[iCone]->rank;
        lorads_int old_rank = ASolver->var->U[iCone]->rank;
        if (new_rank >= ASolver->rank_max[iCone]){
            printinfo_aug(iCone);
        }

        double *dataNewPtrU;
        LORADS_INIT(dataNewPtrU, double, new_rank * ASolver->var->U[iCone]->nRows);
        LORADS_MEMCPY(dataNewPtrU, ASolver->var->U[iCone]->matElem, double, old_rank * ASolver->var->U[iCone]->nRows);
        LORADS_FREE(ASolver->var->U[iCone]->matElem);
        LORADS_INIT(ASolver->var->U[iCone]->matElem, double, new_rank * ASolver->var->U[iCone]->nRows);
        LORADS_MEMCPY(ASolver->var->U[iCone]->matElem, dataNewPtrU, double, old_rank * ASolver->var->U[iCone]->nRows);
        LORADS_FREE(dataNewPtrU);
        lpRandomDiag(&ASolver->var->U[iCone]->matElem[ASolver->var->U[iCone]->nRows * old_rank], ASolver->var->U[iCone]->nRows * aug_rank, ASolver->var->U[iCone]->nRows, aug_rank);
        ASolver->var->U[iCone]->rank = new_rank;

        double *dataNewPtrV;
        LORADS_INIT(dataNewPtrV, double, new_rank * ASolver->var->V[iCone]->nRows);
        LORADS_MEMCPY(dataNewPtrV, ASolver->var->V[iCone]->matElem, double, old_rank * ASolver->var->V[iCone]->nRows);
        LORADS_FREE(ASolver->var->V[iCone]->matElem);
        LORADS_INIT(ASolver->var->V[iCone]->matElem, double, new_rank * ASolver->var->V[iCone]->nRows);
        LORADS_MEMCPY(ASolver->var->V[iCone]->matElem, dataNewPtrV, double, old_rank * ASolver->var->V[iCone]->nRows);
        LORADS_FREE(dataNewPtrV);
        lpRandomDiag(&ASolver->var->V[iCone]->matElem[ASolver->var->U[iCone]->nRows * old_rank], ASolver->var->V[iCone]->nRows * aug_rank, ASolver->var->V[iCone]->nRows, aug_rank);
        ASolver->var->V[iCone]->rank = new_rank;

        double *dataNewPtrR;
        LORADS_INIT(dataNewPtrR, double, new_rank * ASolver->var->R[iCone]->nRows);
        LORADS_MEMCPY(dataNewPtrR, ASolver->var->R[iCone]->matElem, double, old_rank * ASolver->var->R[iCone]->nRows);
        LORADS_FREE(ASolver->var->R[iCone]->matElem);
        LORADS_INIT(ASolver->var->R[iCone]->matElem, double, new_rank * ASolver->var->R[iCone]->nRows);
        LORADS_MEMCPY(ASolver->var->R[iCone]->matElem, dataNewPtrR, double, old_rank * ASolver->var->R[iCone]->nRows);
        LORADS_FREE(dataNewPtrR);
        lpRandomDiag(&ASolver->var->R[iCone]->matElem[ASolver->var->U[iCone]->nRows * old_rank], ASolver->var->U[iCone]->nRows * aug_rank, ASolver->var->U[iCone]->nRows, aug_rank);
        ASolver->var->R[iCone]->rank = new_rank;

        double *dataNewPtrGrad;
        LORADS_INIT(dataNewPtrGrad, double, new_rank * ASolver->var->Grad[iCone]->nRows);
        LORADS_MEMCPY(dataNewPtrGrad, ASolver->var->Grad[iCone]->matElem, double, old_rank * ASolver->var->Grad[iCone]->nRows);
        LORADS_FREE(ASolver->var->Grad[iCone]->matElem);
        LORADS_INIT(ASolver->var->Grad[iCone]->matElem, double, new_rank * ASolver->var->Grad[iCone]->nRows);
        LORADS_MEMCPY(ASolver->var->Grad[iCone]->matElem, dataNewPtrGrad, double, old_rank * ASolver->var->Grad[iCone]->nRows);
        LORADS_FREE(dataNewPtrGrad);
        lpRandomDiag(&ASolver->var->Grad[iCone]->matElem[ASolver->var->U[iCone]->nRows * old_rank], ASolver->var->U[iCone]->nRows * aug_rank, ASolver->var->U[iCone]->nRows, aug_rank);
        ASolver->var->Grad[iCone]->rank = new_rank;

        double *dataNewPtrM2;
        LORADS_INIT(dataNewPtrM2, double, new_rank * ASolver->var->M2temp[iCone]->nRows);
        LORADS_MEMCPY(dataNewPtrM2, ASolver->var->M2temp[iCone]->matElem, double, old_rank * ASolver->var->M2temp[iCone]->nRows);
        LORADS_FREE(ASolver->var->M2temp[iCone]->matElem);
        LORADS_INIT(ASolver->var->M2temp[iCone]->matElem, double, new_rank * ASolver->var->M2temp[iCone]->nRows);
        LORADS_MEMCPY(ASolver->var->M2temp[iCone]->matElem, dataNewPtrM2, double, old_rank * ASolver->var->M2temp[iCone]->nRows);
        LORADS_FREE(dataNewPtrM2);
        ASolver->var->M2temp[iCone]->rank = new_rank;

        allElem += ASolver->var->V[iCone]->rank * ASolver->var->V[iCone]->nRows;
        LORADS_FREE(ASolver->var->bLinSys[iCone]);
        LORADS_INIT(ASolver->var->bLinSys[iCone], double, new_rank * ASolver->var->U[iCone]->nRows);
        ASolver->var->rankElem[iCone] = new_rank;
    }
    LORADSCGReAllocate(ASolver);
    LORADS_FREE(ASolver->var->Dtemp);
    LORADS_INIT(ASolver->var->Dtemp, double, allElem);
    lbfgs_node *head = ASolver->lbfgsHis;
    head->allElem = allElem;
    LORADS_FREE(head->s);
    LORADS_INIT(head->s, double, allElem);
    LORADS_FREE(head->y);
    LORADS_INIT(head->y, double, allElem);
    lbfgs_node *node = head;
    for  (lorads_int nodeNum = 0; nodeNum < ASolver->hisRecT; ++nodeNum)
    {
        lbfgs_node *nextNode = node->next;
        nextNode->allElem = allElem;
        LORADS_FREE(nextNode->s);
        LORADS_INIT(nextNode->s, double, allElem);
        LORADS_FREE(nextNode->y);
        LORADS_INIT(nextNode->y, double, allElem);
        node = node->next;
    }

    return CheckAllRankMax(ASolver, aug_factor);
}

extern void printRes(double pObj, double dObj, double constrVio, double dualInfe, double pdgap, double constrVioInf, double dualInfeInf)
{
    printf("-----------------------------------------------------------------------\n");
    printf("Objective function Value are:\n");
    printf("\t 1.Primal Objective:            : %10.6e\n", pObj);
    printf("\t 2.Dual Objective:              : %10.6e\n", dObj);
    printf("Dimacs Error are:\n");
    printf("\t 1.Constraint Violation(1)      : %10.6e\n", constrVio);
    printf("\t 2.Dual Infeasibility(1)        : %10.6e\n", dualInfe);
    printf("\t 3.Primal Dual Gap              : %10.6e\n", pdgap);
    printf("\t 4.Primal Variable Semidefinite : %10.6e\n", 0.0);
    printf("\t 5.Constraint Violation(Inf)    : %10.6e\n", constrVioInf);
    printf("\t 6.Dual Infeasibility(Inf)      : %10.6e\n", dualInfeInf);
    printf("-----------------------------------------------------------------------\n");
}

extern void  LORADSEndProgram( lorads_solver *ASolver)
{
    printf("final rank: \n");
    for (int i = 0; i < ASolver->nCones; i++) {
        #ifdef LORADS_INT32
            printf("cone %d rank: %d", i, ASolver->var->U[i]->rank);
        #endif
        #ifdef UNIX_INT64
            printf("cone %d rank: %ld", i, ASolver->var->U[i]->rank);
        #endif
        #ifdef MAC_INT64
            printf("cone %d rank: %lld", i, ASolver->var->U[i]->rank);
        #endif
    }
    printf("\n");
    printf("-----------------------------------------------------------------------\n");
    if (ASolver->AStatus ==  LORADS_MAXITER)
    {
         printf("End Program due to reaching `the maximum number of iterations`:\n");
    }
    else if (ASolver->AStatus ==  LORADS_PRIMAL_DUAL_OPTIMAL)
    {
         printf("End Program due to reaching `Official terminate criteria`:\n");
    }
    else if (ASolver->AStatus ==  LORADS_PRIMAL_OPTIMAL)
    {
         printf("End Program due to reaching `final terminate criteria`:\n");
    }
    else if (ASolver->AStatus ==  LORADS_UNKNOWN)
    {
         printf("End Program but the status is unknown, please notify the authors\n");
    }
    else if (ASolver->AStatus == LORADS_TIME_LIMIT)
    {
        printf("End Program since time limit.\n");
    }
    printRes(ASolver->pObjVal, ASolver->dObjVal,
             ASolver->dimacError[ LORADS_DIMAC_ERROR_CONSTRVIO_L1],
             ASolver->dimacError[ LORADS_DIMAC_ERROR_DUALFEASIBLE_L1],
             ASolver->dimacError[ LORADS_DIMAC_ERROR_PDGAP],
             ASolver->dimacError[ LORADS_DIMAC_ERROR_CONSTRVIO_L1] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf),
             ASolver->dimacError[ LORADS_DIMAC_ERROR_DUALFEASIBLE_L1] * (1 + ASolver->cObjNrm1) / (1 + ASolver->cObjNrmInf));
}

extern void LORADS_ALMtoADMM(lorads_solver *ASolver, lorads_params *params, lorads_alm_state *alm_state, lorads_admm_state *admm_state)
{
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        lorads_sdp_dense *R = ASolver->var->R[iCone];
        lorads_sdp_dense *V = ASolver->var->V[iCone];
        LORADS_MEMCPY(V->matElem, R->matElem, double, V->nRows * V->rank);
    }

    // copy V -> U
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        lorads_sdp_dense *U = ASolver->var->U[iCone];
        lorads_sdp_dense *V = ASolver->var->V[iCone];
        LORADS_MEMCPY(U->matElem, V->matElem, double, V->nRows * V->rank);
    }
    if (ASolver->nLpCols > 0)
    {
        lorads_lp_dense *r = ASolver->var->rLp;
        lorads_lp_dense *u = ASolver->var->uLp;
        lorads_lp_dense *v = ASolver->var->vLp;
        LORADS_MEMCPY(u->matElem, r->matElem, double, v->nCols);
        LORADS_MEMCPY(v->matElem, r->matElem, double, v->nCols);
    }
    admm_state->l_1_dual_infeasibility = alm_state->l_1_dual_infeasibility;
    admm_state->l_1_primal_infeasibility = alm_state->l_1_primal_infeasibility;
    admm_state->l_2_dual_infeasibility = alm_state->l_2_dual_infeasibility;
    admm_state->l_inf_dual_infeasibility = alm_state->l_inf_dual_infeasibility;
    admm_state->l_inf_primal_infeasibility = alm_state->l_inf_primal_infeasibility;
    admm_state->l_2_primal_infeasibility = alm_state->l_2_primal_infeasibility;
    admm_state->primal_dual_gap = alm_state->primal_dual_gap;
    admm_state->rho = alm_state->rho * params->heuristicFactor;
    if (alm_state->rho > params->rhoMax){
        admm_state->rho = LORADS_MIN(sqrt(LORADS_MAX(params->rhoMax, alm_state->rho) / params->rhoMax) * params->rhoMax, alm_state->rho);
        params->rhoMax = admm_state->rho;
    }
}


extern void calculate_dual_infeasibility_solver(lorads_solver *ASolver){
    double *negLambd;
    LORADS_INIT(negLambd, double, ASolver->nRows);
    LORADS_ZERO(negLambd, double, ASolver->nRows);
    double minusOne = -1.0; lorads_int incx = 1;
    axpy(&ASolver->nRows, &minusOne, ASolver->var->dualVar, &incx, negLambd, &incx);
    double lp_obj_w_sum = 0.0;
    ASolver->dimacError[LORADS_DIMAC_ERROR_DUALFEASIBLE_L1] = 0.0;
    if (ASolver->nLpCols>0){
        lorads_lp_cone *lp_cone = ASolver->lpCone;
        for (lorads_int iCol = 0; iCol < ASolver->nLpCols; ++iCol){
            lp_obj_w_sum = 0.0;
            lp_cone->objCoeffSum(lp_cone->coneData, &lp_obj_w_sum, iCol);
            lp_cone->lpDataWSum(lp_cone->coneData, negLambd, &lp_obj_w_sum, iCol);
            ASolver->dimacError[LORADS_DIMAC_ERROR_DUALFEASIBLE_L1] += LORADS_ABS(LORADS_MIN(lp_obj_w_sum, 0));
        }
    }
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];

        ACone->sdp_slack_var->zeros(ACone->sdp_slack_var->dataMat);
        ACone->addObjCoeff(ACone->coneData, ACone->sdp_slack_var);
        ACone->sdpDataWSum(ACone->coneData, negLambd, ACone->sdp_slack_var);
//        ACone->sdp_slack_var->scaleData(ACone->sdp_slack_var->dataMat,  1e+6);
        dual_infeasible(ACone->sdp_slack_var->mv, ACone->sdp_slack_var->dataMat, &lp_obj_w_sum, ACone->sdp_slack_var->nSDPCol);
        ASolver->dimacError[LORADS_DIMAC_ERROR_DUALFEASIBLE_L1] += LORADS_ABS(LORADS_MIN(lp_obj_w_sum, 0));
    }
    ASolver->dimacError[LORADS_DIMAC_ERROR_DUALFEASIBLE_L1] /= ASolver->scaleObjHis;
    ASolver->dimacError[LORADS_DIMAC_ERROR_DUALFEASIBLE_L1] /= (ASolver->cObjNrm1 + 1);
    LORADS_FREE(negLambd);
}


extern void objScale_dualvar(lorads_solver *ASolver, double *scaleTemp, double *scaleHis){
    scaleHis[0] *= scaleTemp[0];
    if(ASolver->nLpCols > 0){
        lorads_lp_cone *lpCone = ASolver->lpCone;
        lpCone->scalObj(lpCone->coneData, scaleTemp[0]);
    }
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        ACone->objScale(ACone->coneData, scaleTemp[0]);
    }
    lorads_int incx = 1;
    scal(&ASolver->nRows, scaleTemp, ASolver->var->dualVar, &incx);
}

extern void cal_sdp_const(lorads_params *params, lorads_solver *ASolver, SDPConst *sdpConst){
    // Calculate obj nrm1, nrm2, nrmInf of obj
    LORADSNrm1Obj(ASolver);
    LORADSNrm2Obj(ASolver);
    LORADSNrmInfObj(ASolver);
    ASolver->bRHSNrm1 = nrm1(&ASolver->nRows, ASolver->rowRHS, &AIntConstantOne);
#ifdef UNDER_BLAS
    ASolver->bRHSNrmInf = fabs(ASolver->rowRHS[idamax_(&ASolver->nRows, ASolver->rowRHS, &AIntConstantOne)]);
#else
    ASolver->bRHSNrmInf = fabs(ASolver->rowRHS[idamax(&ASolver->nRows, ASolver->rowRHS, &AIntConstantOne)]);
#endif

    ASolver->bRHSNrm2 = nrm2(&ASolver->nRows, ASolver->rowRHS, &AIntConstantOne);
    sdpConst->l_1_norm_c = ASolver->cObjNrm1;
    sdpConst->l_inf_norm_c = ASolver->cObjNrmInf;
    sdpConst->l_2_norm_c = ASolver->cObjNrm2;
    sdpConst->l_1_norm_b = ASolver->bRHSNrm1;
    sdpConst->l_inf_norm_b = ASolver->bRHSNrmInf;
    sdpConst->l_2_norm_b = ASolver->bRHSNrm2;
}

extern double reopt(lorads_params *params, lorads_solver *ASolver, lorads_alm_state *alm_state_pointer, lorads_admm_state *admm_state_pointer, double *reopt_param, lorads_int *reopt_alm_iter, lorads_int *reopt_admm_iter, double timeSolveStart, int *admm_bad_iter_flag, int reopt_level) {
    lorads_int old_maxALMIter = params->maxALMIter;
    lorads_int old_maxADMMIter = params->maxADMMIter;
    double old_rhoMax = params->rhoMax;
    params->maxALMIter = reopt_alm_iter[0] - 1 + alm_state_pointer->outerIter;
    params->maxADMMIter = reopt_admm_iter[0];

    // scale the objective function and dual variable
    objScale_dualvar(ASolver, reopt_param, &ASolver->scaleObjHis);

    if (admm_state_pointer->rho <= params->rhoMax) {
        alm_state_pointer->rho = LORADS_MAX(admm_state_pointer->rho, alm_state_pointer->rho);
    }

    double start_time = LUtilGetTimeStamp();

    LORADS_ALMOptimize_reopt(params, ASolver, alm_state_pointer, true, sqrt(params->ALMRhoFactor), timeSolveStart);
    params->rhoMax = LORADS_MAX(
            sqrt(LORADS_MAX(admm_state_pointer->rho, alm_state_pointer->rho) / admm_state_pointer->rho) *
            admm_state_pointer->rho, params->rhoMax);
    LORADS_ALMtoADMM(ASolver, params, alm_state_pointer, admm_state_pointer);
    if ((*admm_bad_iter_flag == 0 || reopt_level < 2))
    {
        if (LORADSADMMOptimize_reopt(params, ASolver, admm_state_pointer, LORADS_MIN(admm_state_pointer->iter * 4, admm_state_pointer->iter + old_maxADMMIter), timeSolveStart) == RET_CODE_BAD_ITER){
            *admm_bad_iter_flag = 1;
        }
        else
        {
            *admm_bad_iter_flag = 0;
        }
    }
    // if (admm_state_pointer->primal_dual_gap >= params->phase1Tol)
    // {
    //     params->phase1Tol *= 0.1;
    //     params->phase1Tol = LORADS_MAX(params->phase1Tol, params->phase2Tol);
    // }

    double end_time = LUtilGetTimeStamp();
    params->maxALMIter = old_maxALMIter;
    params->maxADMMIter = old_maxADMMIter;
    params->rhoMax = old_rhoMax;
    return end_time - start_time;
}

extern void LORADSInitALMState(lorads_solver *ASolver, lorads_alm_state *alm_state_pointer, double rho, lorads_int outerIter, lorads_int innerIter){

    alm_state_pointer->dual_objective_value = 1e+30;
    alm_state_pointer->primal_objective_value = 1e+30;
    alm_state_pointer->l_1_dual_infeasibility = 1e+30;
    alm_state_pointer->l_1_primal_infeasibility = 1e+30;
    alm_state_pointer->l_inf_dual_infeasibility = 1e+30;
    alm_state_pointer->l_inf_primal_infeasibility = 1e+30;
    alm_state_pointer->rho = rho;
    alm_state_pointer->outerIter = outerIter;
    alm_state_pointer->innerIter = innerIter;
    alm_state_pointer->is_rank_updated = false;
    alm_state_pointer->tau = 0.0;
}

extern void LORADSInitADMMState(lorads_solver *ASolver, lorads_admm_state *admm_state_pointer, double rho){
    admm_state_pointer->dual_objective_value = 1e+30;
    admm_state_pointer->primal_objective_value = 1e+30;
    admm_state_pointer->primal_dual_gap = 1e+30;
    admm_state_pointer->l_1_dual_infeasibility = 1e+30;
    admm_state_pointer->l_1_primal_infeasibility = 1e+30;
    admm_state_pointer->l_inf_dual_infeasibility = 1e+30;
    admm_state_pointer->l_inf_primal_infeasibility = 1e+30;
    admm_state_pointer->l_2_dual_infeasibility = 1e+30;
    admm_state_pointer->l_2_primal_infeasibility = 1e+30;
    admm_state_pointer->rho = rho;
    admm_state_pointer->iter = 0;
}

extern void initial_solver_state(lorads_params *params, lorads_solver *ASolver, lorads_alm_state *alm_state_pointer, lorads_admm_state *admm_state_pointer, SDPConst *sdpConst){
    // Calculate obj nrm1, nrm2, nrmInf of obj
    cal_sdp_const(params, ASolver, sdpConst);
    // Initialize the ALM state
    lorads_int outerIter = 0;
    lorads_int innerIter = 0;
    double rho = 0.0;
    if (params->initRho == 0){
        lorads_int sumBlkDims = 0;
        for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
            sumBlkDims += ASolver->var->U[iCone]->nRows;
        }
        // rho = 1.0 / (sumBlkDims);
        rho = 1 / sqrt(sumBlkDims);
    }else{
        rho = params->initRho;
    }
    LORADSInitALMState(ASolver, alm_state_pointer, rho, outerIter, innerIter);
    // Initialize the ADMM state
    LORADSInitADMMState(ASolver, admm_state_pointer, rho);
    ASolver->scaleObjHis = 1;
    admm_state_pointer->nBlks = ASolver->nCones;
}


extern void printfProbInfo(lorads_solver *ASolver){
    printf("-----------------------------------------------------------------------\n");
    printf("Problem Information:\n");
#ifdef INT32
    printf("\t 1.Number of SDP Cones:         : %10d\n", ASolver->nCones);
    printf("\t 2.Number of LP Cones:          : %10d\n", ASolver->nLpCols);
    printf("\t 3.Number of Constraints:       : %10d\n", ASolver->nRows);
    printf("\t 4.sdp block dims:              : %3d,", ASolver->var->U[0]->nRows);
    for (lorads_int i = 1; i < ASolver->nCones; ++i){
        printf("%3d,", ASolver->var->U[i]->nRows);
    }
    printf("\n");
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        printf("iCone:%d\n", iCone);
        ASolver->SDPCones[iCone]->coneView(ASolver->SDPCones[iCone]->coneData);
    }
    printf("Initial rank:\n");
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        printf("iCone:%d, rank:%d\n", iCone, ASolver->var->U[iCone]->rank);
    }
#endif
#ifdef MAC_INT64
    printf("\t 1.Number of SDP Cones:         : %10lld\n", ASolver->nCones);
    printf("\t 2.Number of LP Cones:          : %10lld\n", ASolver->nLpCols);
    printf("\t 3.Number of Constraints:       : %10lld\n", ASolver->nRows);
    printf("\t 4.sdp block dims:              : %3lld,", ASolver->var->U[0]->nRows);
    for (lorads_int i = 1; i < ASolver->nCones; ++i){
        printf("%3lld,", ASolver->var->U[i]->nRows);
    }
    printf("\n");
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        printf("iCone:%lld\n", iCone);
        ASolver->SDPCones[iCone]->coneView(ASolver->SDPCones[iCone]->coneData);
    }
    printf("Initial rank:\n");
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        printf("iCone:%lld, rank:%lld\n", iCone, ASolver->var->U[iCone]->rank);
    }
#endif
#ifdef UNIX_INT64
    printf("\t 1.Number of SDP Cones:         : %10ld\n", ASolver->nCones);
    printf("\t 2.Number of LP Cones:          : %10ld\n", ASolver->nLpCols);
    printf("\t 3.Number of Constraints:       : %10ld\n", ASolver->nRows);
    printf("\t 4.sdp block dims:              : %3ld,", ASolver->var->U[0]->nRows);
    for (lorads_int i = 1; i < ASolver->nCones; ++i){
        printf("%3ld,", ASolver->var->U[i]->nRows);
    }
    printf("\n");
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        printf("iCone:%ld\n", iCone);
        ASolver->SDPCones[iCone]->coneView(ASolver->SDPCones[iCone]->coneData);
    }
    printf("Initial rank:\n");
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        printf("iCone:%ld, rank:%ld\n", iCone, ASolver->var->U[iCone]->rank);
    }
#endif
    printf("-----------------------------------------------------------------------\n");
}



#include <math.h>
#include "lorads_solver.h"
#include "lorads_vec_opts.h"
#include "lorads_dense_opts.h"
#include "lorads_admm.h"

extern void valRes(void *constrVal, vecType type, double **res)
{
    if (type == LORADS_DENSE_VEC){
        dense_vec *dense = (dense_vec *)constrVal;
        *res = dense->val;
    }
    else if (type == LORADS_SPARSE_VEC){
        sparse_vec *sparse = (sparse_vec *)constrVal;
        *res = sparse->val;
    }
}

extern void LORADSUVt(sdp_coeff *UVt_w_sum, lorads_sdp_dense *U, lorads_sdp_dense *V){
    /* sdp_coeff_w_sum is the sum of all sdp data in the cone,
     It has two data type
     - sparse: means most all sdp data coeff is sparse
     - dense:  there exist a sdp data coeff is dense

     calculate (UVt + VUt) / 2
     */
    if (UVt_w_sum->dataType == SDP_COEFF_SPARSE){
        sdp_coeff_sparse *sparse = (sdp_coeff_sparse *)UVt_w_sum->dataMat;
        // method1:
        lorads_int row = 0;
        lorads_int col = 0;
        lorads_int incx = 1;
        lorads_int incy = 1;
        for (lorads_int i = 0; i < sparse->nTriMatElem; ++i)
        {
            row = sparse->triMatRow[i];
            col = sparse->triMatCol[i];
            incx = U->nRows;
            incy = V->nRows;
            if (row != col){
                sparse->triMatElem[i] = 0.5 * dot(&(U->rank), &U->matElem[row], &incx, &V->matElem[col], &incy);
                sparse->triMatElem[i] += 0.5 *  dot(&(U->rank), &U->matElem[col], &incx, &V->matElem[row], &incy);
            }else{
                sparse->triMatElem[i] = dot(&(U->rank), &U->matElem[row], &incx, &V->matElem[col], &incy);
            }
        }
    }
    else if (UVt_w_sum->dataType == SDP_COEFF_DENSE){
        sdp_coeff_dense *dense = (sdp_coeff_dense *)UVt_w_sum->dataMat;
//        double *fullDataMat;
//        LORADS_INIT(fullDataMat, double, dense->nSDPCol * dense->nSDPCol);
        LORADS_ZERO(dense->fullMat, double, dense->nSDPCol * dense->nSDPCol);
//        printf("rank: %lld\n", U->rank);
        // alpha = 0.5, beta = 0.0;
        fds_syr2k(ACharConstantUploLow, 'N', U->nRows, U->rank, 0.5, U->matElem, V->matElem, 0.0, dense->fullMat);
        lorads_int idx = 0;
        lorads_int row = 0;
        for (lorads_int col = 0; col < dense->nSDPCol; ++col)
        {
            LORADS_MEMCPY(&dense->dsMatElem[idx], &dense->fullMat[dense->nSDPCol * col + row], double, dense->nSDPCol - col);
            row++;
            idx += (dense->nSDPCol - col);
        }
//        LORADS_FREE(fullDataMat);
    }
}


extern void LORADSInitConstrVal(lorads_sdp_cone *ACone, lorads_sdp_dense *U, lorads_sdp_dense *V, double *constrVal)
{
    LORADSUVt(ACone->sdp_coeff_w_sum, U, V);
    // Calculate Constraint Value for one cone
    ACone->coneAUV(ACone->coneData, U, V, constrVal, ACone->sdp_coeff_w_sum);
}

extern void LORADSInitConstrValAll(lorads_solver *ASolver, lorads_lp_dense *uLpDummy, lorads_lp_dense *vLpDummy, lorads_sdp_dense **U, lorads_sdp_dense **V){
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        double *constrVal;
        valRes(ASolver->var->constrVal[iCone]->data, ASolver->var->constrVal[iCone]->type, &constrVal);
        LORADSInitConstrVal(ASolver->SDPCones[iCone], U[iCone], V[iCone], constrVal);
    }
}

extern void LORADSInitConstrValAllLP(lorads_solver *ASolver, lorads_lp_dense *uLp, lorads_lp_dense *vLp, lorads_sdp_dense **U, lorads_sdp_dense **V)
{
    lorads_lp_cone *lpCone = ASolver->lpCone;
    LORADSInitConstrValAll(ASolver, uLp, vLp, U, V);
    for (lorads_int iCol = 0; iCol < lpCone->nCol; ++iCol){
        double *constrVal;
        valRes(ASolver->var->constrValLP[iCol]->data, ASolver->var->constrValLP[iCol]->type, &constrVal);
        lpCone->coneAUV(lpCone->coneData, uLp, vLp, constrVal, iCol);
    }
}

extern void LORADSInitConstrValObjVal(lorads_sdp_cone *ACone, lorads_sdp_dense *U, lorads_sdp_dense *V, double *constrVal, double *UVobjVal){
    LORADSUVt(ACone->sdp_obj_sum, U, V);
    // Calculate Constraint Value for one cone
    ACone->objAUV(ACone->coneData, U, V, UVobjVal, ACone->sdp_obj_sum);
//    printf("objVal: %f\n", UVobjVal[0]);
    ACone->coneAUV(ACone->coneData, U, V, constrVal, ACone->sdp_obj_sum);
}

extern void LORADSObjConstrValAll(lorads_solver *ASolver, lorads_sdp_dense **U, lorads_sdp_dense **V, double *UVobjVal){
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        ASolver->var->constrVal[iCone]->zero(ASolver->var->constrVal[iCone]->data);
        double *constrVal;
        valRes(ASolver->var->constrVal[iCone]->data, ASolver->var->constrVal[iCone]->type, &constrVal);
        LORADSInitConstrValObjVal(ASolver->SDPCones[iCone], U[iCone], V[iCone], constrVal, UVobjVal);
    }
}

extern void LORADSObjConstrValAllLP(lorads_solver *ASolver, lorads_lp_dense *uLp, lorads_lp_dense *vLp, lorads_sdp_dense **U, lorads_sdp_dense **V, double *UVobjVal)
{
    ASolver->lpCone->objAUV(ASolver->lpCone->coneData, uLp, vLp, UVobjVal);
    for (lorads_int iCol = 0; iCol < uLp->nCols; ++iCol)
    {
        double *constrVal;
        valRes(ASolver->var->constrValLP[iCol]->data, ASolver->var->constrValLP[iCol]->type, &constrVal);
        ASolver->lpCone->coneAUV(ASolver->lpCone->coneData, uLp, vLp, constrVal, iCol);
    }
    LORADSObjConstrValAll(ASolver, U, V, UVobjVal);
}

extern void LORADSUpdateConstrVal(lorads_sdp_cone *ACone, lorads_sdp_dense *U, lorads_sdp_dense *V, lorads_vec *constrVal)
{
    // constrVal = Acal(UVt)
    double *constrValRes;
    valRes(constrVal->data, constrVal->type, &constrValRes);
    LORADSInitConstrVal(ACone, U, V, constrValRes);
}

extern void LORADSInitConstrValSum(lorads_solver *ASolver)
{
    LORADS_ZERO(ASolver->var->constrValSum, double, ASolver->nRows);
    double alpha = 1.0;
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASolver->var->constrVal[iCone]->add(&alpha, ASolver->var->constrVal[iCone]->data, ASolver->var->constrValSum);
    }
}

extern void LORADSInitConstrValSumLP(lorads_solver *ASolver)
{
    LORADS_ZERO(ASolver->var->constrValSum, double, ASolver->nRows);
    double alpha = 1.0;
    for (lorads_int iCol = 0; iCol < ASolver->nLpCols; ++iCol)
    {
        ASolver->var->constrValLP[iCol]->add(&alpha, ASolver->var->constrValLP[iCol]->data, ASolver->var->constrValSum);
    }

    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        ASolver->var->constrVal[iCone]->add(&alpha, ASolver->var->constrVal[iCone]->data, ASolver->var->constrValSum);
    }
}


extern void copyRtoV(lorads_lp_dense *rLpDummy, lorads_lp_dense *vlpDummy, lorads_sdp_dense **R, lorads_sdp_dense **V, lorads_int nCones)
{
    lorads_int n = 0;
    for (lorads_int iCone = 0; iCone < nCones; ++iCone)
    {
        n = R[iCone]->rank * R[iCone]->nRows;
        LORADS_MEMCPY(V[iCone]->matElem, R[iCone]->matElem, double, n);
    }
}



extern void copyRtoVLP(lorads_lp_dense *rLp, lorads_lp_dense *vlp, lorads_sdp_dense **R, lorads_sdp_dense **V, lorads_int nCones)
{
    lorads_int n = 0;
    LORADS_MEMCPY(vlp->matElem, rLp->matElem, double, vlp->nCols);
    for (lorads_int iCone = 0; iCone < nCones; ++iCone)
    {
        n = R[iCone]->rank * R[iCone]->nRows;
        LORADS_MEMCPY(V[iCone]->matElem, R[iCone]->matElem, double, n);
    }
}





extern void LORADSUpdateSDPVar(lorads_solver *ASolver, double rho, double CG_tol, lorads_int CG_maxIter){
    double minusOne = -1.0;
    double one = 1.0;
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        // update U
#ifdef DUAL_U_V
        LORADSUpdateSDPVarOne_positive_S(ASolver, ASolver->var->U[iCone], ASolver->var->V[iCone], ASolver->var->S[iCone], iCone, rho, CG_tol, CG_maxIter);
#else
        LORADSUpdateSDPVarOne(ASolver, ASolver->var->U[iCone], ASolver->var->V[iCone], iCone, rho, CG_tol, CG_maxIter);
#endif
        // ASolver->constrValSum - ASolver->constrVal[iCone]
        ASolver->var->constrVal[iCone]->add(&minusOne, ASolver->var->constrVal[iCone]->data, ASolver->var->constrValSum);
        // update ASolver->constrVal[iCone]
        LORADSUpdateConstrVal(ASolver->SDPCones[iCone], ASolver->var->U[iCone], ASolver->var->V[iCone], ASolver->var->constrVal[iCone]);
        // update ASolver->constrValSum
        ASolver->var->constrVal[iCone]->add(&one, ASolver->var->constrVal[iCone]->data, ASolver->var->constrValSum);

        // update V
#ifdef DUAL_U_V
        LORADSUpdateSDPVarOne_negative_S(ASolver, ASolver->var->V[iCone], ASolver->var->U[iCone], ASolver->var->S[iCone], iCone, rho, CG_tol, CG_maxIter);
#else
        LORADSUpdateSDPVarOne(ASolver, ASolver->var->V[iCone], ASolver->var->U[iCone], iCone, rho, CG_tol, CG_maxIter);
#endif
        ASolver->var->constrVal[iCone]->add(&minusOne, ASolver->var->constrVal[iCone]->data, ASolver->var->constrValSum);
        LORADSUpdateConstrVal(ASolver->SDPCones[iCone], ASolver->var->U[iCone], ASolver->var->V[iCone], ASolver->var->constrVal[iCone]);
        ASolver->var->constrVal[iCone]->add(&one, ASolver->var->constrVal[iCone]->data, ASolver->var->constrValSum);
    }
}


extern void LORADSUpdateConstrValLP(lorads_lp_cone *lp_cone, lorads_lp_dense *uLp, lorads_lp_dense *vLp, lorads_vec *constrVal, lorads_int iCol)
{
    double *constrValRes;
    valRes(constrVal->data, constrVal->type, &constrValRes);
    lp_cone->coneAUV(lp_cone->coneData, uLp, vLp, constrValRes, iCol);
}

extern void LORADSUpdateSDPLPVar(lorads_solver *ASolver, double rho, double CG_tol, lorads_int CG_maxIter){
    double minusOne = -1.0;
    double one = 1.0;
    LORADSUpdateSDPVar(ASolver, rho, CG_tol, CG_maxIter);
    for (lorads_int iCol = 0; iCol < ASolver->nLpCols; ++iCol){
#ifdef DUAL_U_V
        LORADSUpdateLPVarOne_positive_S(ASolver, &ASolver->var->uLp->matElem[iCol], &ASolver->var->vLp->matElem[iCol],  iCol, rho, ASolver->var->sLp->matElem);
#else
        LORADSUpdateLPVarOne(ASolver, &ASolver->var->uLp->matElem[iCol], &ASolver->var->vLp->matElem[iCol],  iCol, rho);
#endif
        ASolver->var->constrValLP[iCol]->add(&minusOne, ASolver->var->constrValLP[iCol]->data, ASolver->var->constrValSum);
        LORADSUpdateConstrValLP(ASolver->lpCone, ASolver->var->uLp, ASolver->var->vLp, ASolver->var->constrValLP[iCol], iCol);
        ASolver->var->constrValLP[iCol]->add(&one, ASolver->var->constrValLP[iCol]->data, ASolver->var->constrValSum);

#ifdef DUAL_U_V
        LORADSUpdateLPVarOne_negative_S(ASolver, &ASolver->var->vLp->matElem[iCol], &ASolver->var->uLp->matElem[iCol],   iCol, rho, ASolver->var->sLp->matElem);
#else
        LORADSUpdateLPVarOne(ASolver, &ASolver->var->vLp->matElem[iCol], &ASolver->var->uLp->matElem[iCol],   iCol, rho);
#endif
        ASolver->var->constrValLP[iCol]->add(&minusOne, ASolver->var->constrValLP[iCol]->data, ASolver->var->constrValSum);
        LORADSUpdateConstrValLP(ASolver->lpCone, ASolver->var->uLp, ASolver->var->vLp, ASolver->var->constrValLP[iCol], iCol);
        ASolver->var->constrValLP[iCol]->add(&one, ASolver->var->constrValLP[iCol]->data, ASolver->var->constrValSum);
    }
}

extern void primalInfeasibility(lorads_solver *ASolver, lorads_sdp_dense **R, lorads_sdp_dense **R2, lorads_lp_dense *r, lorads_lp_dense *r2){
    LORADSInitConstrValAll(ASolver, r, r2, R, R2);
    LORADSInitConstrValSum(ASolver);
    double one = 1.0;
    double minusOne = -1.0;
    axpbyAddition(&ASolver->nRows, &(one), ASolver->rowRHS, &(minusOne), ASolver->var->constrValSum, ASolver->constrVio);
    lorads_int incx = 1;
    ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1] = nrm2(&ASolver->nRows, ASolver->constrVio, &incx) / (1 + ASolver->bRHSNrm1);
}

extern void primalInfeasibilityLP(lorads_solver *ASolver, lorads_sdp_dense **R, lorads_sdp_dense **R2, lorads_lp_dense *r, lorads_lp_dense *r2){
    LORADSInitConstrValAllLP(ASolver, r, r2, R, R2);
    LORADSInitConstrValSumLP(ASolver);
    double one = 1.0;
    double minusOne = -1.0;
    axpbyAddition(&ASolver->nRows, &(one), ASolver->rowRHS, &(minusOne), ASolver->var->constrValSum, ASolver->constrVio);
    lorads_int incx = 1;
    ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1] = nrm2(&ASolver->nRows, ASolver->constrVio, &incx) / (1 + ASolver->bRHSNrm1);
}

extern void LORADSUpdateDimacsErrorALM(lorads_solver *ASolver, lorads_sdp_dense **R, lorads_sdp_dense **R2, lorads_lp_dense *r, lorads_lp_dense *r2){
    primalInfeasibility(ASolver, R, R2, r, r2);
    double gap = (ASolver->pObjVal - ASolver->dObjVal);
    ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP] = LORADS_ABS(gap) / (1 + LORADS_ABS(ASolver->pObjVal) + LORADS_ABS(ASolver->dObjVal));
}

extern void LORADSUpdateDimacsErrorALMLP(lorads_solver *ASolver, lorads_sdp_dense **R, lorads_sdp_dense **R2, lorads_lp_dense *r, lorads_lp_dense *r2){
    primalInfeasibilityLP(ASolver, R, R2, r, r2);
    double gap = (ASolver->pObjVal - ASolver->dObjVal);
    ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP] = LORADS_ABS(gap) / (1 + LORADS_ABS(ASolver->pObjVal) + LORADS_ABS(ASolver->dObjVal));
}

extern void LORADSUpdateDimacsErrorADMM(lorads_solver *ASolver, lorads_sdp_dense **U, lorads_sdp_dense **V, lorads_lp_dense *u, lorads_lp_dense *v){
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        averageUV(U[iCone], V[iCone], ASolver->var->R[iCone]);
    }
    averageUVLP(u, v, ASolver->var->rLp);
    primalInfeasibility(ASolver, ASolver->var->R, ASolver->var->R, ASolver->var->rLp, ASolver->var->rLp);
    double gap = (ASolver->pObjVal - ASolver->dObjVal);
    ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP] = LORADS_ABS(gap) / (1 + LORADS_ABS(ASolver->pObjVal) + LORADS_ABS(ASolver->dObjVal));
}

extern void LORADSUpdateDimacsErrorADMMLP(lorads_solver *ASolver, lorads_sdp_dense **U, lorads_sdp_dense **V, lorads_lp_dense *u, lorads_lp_dense *v){
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        averageUV(U[iCone], V[iCone], ASolver->var->R[iCone]);
    }
    averageUVLP(u, v, ASolver->var->rLp);
    primalInfeasibilityLP(ASolver, ASolver->var->R, ASolver->var->R, ASolver->var->rLp, ASolver->var->rLp);
    double gap = (ASolver->pObjVal - ASolver->dObjVal);
    ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP] = LORADS_ABS(gap) / (1 + LORADS_ABS(ASolver->pObjVal) + LORADS_ABS(ASolver->dObjVal));
}

extern void LORADSNrmInfObj(lorads_solver *ASolver)
{
    ASolver->cObjNrmInf = 0.0;
    if (ASolver->nLpCols > 0)
    {
        lorads_lp_cone *lpCone = ASolver->lpCone;
        lpCone->coneObjNrmInf(lpCone->coneData, &ASolver->cObjNrmInf, ASolver->nLpCols);
    }
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        double temp = 0.0;
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        ACone->coneObjNrmInf(ACone->coneData, &temp);
        ASolver->cObjNrmInf = LORADS_MAX(ASolver->cObjNrmInf, temp);
    }
}

extern void LORADSUpdateDualVar(lorads_solver *ASolver, double rho)
{
    double *dualVar = ASolver->var->dualVar;
    double minusRho = -rho;
    double *b = ASolver->rowRHS;
    double *constrSum = ASolver->var->constrValSum;
    lorads_int n = ASolver->nRows;

    // lambda = lambda + rho * b
    lorads_int incx = 1;
    axpy(&(n), &(rho), b, &(incx), dualVar, &(incx));
    // lambda = lambda - rho * constrVal
    axpy(&(n), &(minusRho), constrSum, &(incx), dualVar, &(incx));
}

extern void LORADSCalDualObj(lorads_solver *ASolver)
{
    lorads_int n = ASolver->nRows;
    lorads_int one = 1;
    ASolver->dObjVal = dot(&n, ASolver->rowRHS, &one, ASolver->var->dualVar, &one);
    ASolver->dObjVal /= ASolver->scaleObjHis;
}


//extern void LORADSNuclearNorm(lorads_solver *ASolver)
//{
//    ASolver->traceSum = 0;
//    lorads_int rank = 0;
//    lorads_int nRows = 0;
//    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
//        lorads_sdp_dense *R = ASolver->var->R[iCone];
//        rank = R->rank;
//        nRows = R->nRows;
//        for(lorads_int idx = 0; idx < nRows; ++idx){
//            ASolver->traceSum += dot(&rank, &R->matElem[idx], &nRows, &R->matElem[idx], &nRows);
//        }
//    }
//}


extern lorads_int LORADSCheckSolverStatus(lorads_solver *ASolver){
    lorads_int retcode = LORADS_RETCODE_OK;
    if (ASolver->AStatus == LORADS_PRIMAL_OPTIMAL){
        retcode = LORADS_RETCODE_EXIT;
    }
    return retcode;
}
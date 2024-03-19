#ifdef HEADERPATH
#include "interface/asdp_conic_sdp.h"
#include "interface/def_asdp_user_data.h"
#include "interface/asdp_utils.h"
#include "src/asdp_sdpdata.h"
#include "src/sparse_opts.h"
#include "src/dense_opts.h"
#include "src/vec_opts.h"
#include "external/asdp_cs.h"
#include "external/def_asdp.h"
#include "src/def_asdp_sdpdata.h"
#include "src/asdp_debug.h"
#else
#include "asdp_conic_sdp.h"
#include "def_asdp_user_data.h"
#include "asdp_utils.h"
#include "asdp_sdpdata.h"
#include "sparse_opts.h"
#include "dense_opts.h"
#include "vec_opts.h"
#include "asdp_cs.h"
#include "def_asdp.h"
#include "def_asdp_sdpdata.h"
#include "asdp_debug.h"
#endif

#include <math.h>

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

#ifndef SMALL_DUAL_THRESHOLD
#define SMALL_DUAL_THRESHOLD (0)
#endif

#ifndef SPARSE_DUAL_THRESHOLD
#define SPARSE_DUAL_THRESHOLD (0.6)
#endif


#define SPARSE_EFFICIENCY (1.5)

extern asdp_retcode sdpSparseConeCreateImpl( void **pConeIn ) {
    asdp_cone_sdp_sparse **pCone = ( asdp_cone_sdp_sparse **)pConeIn;
    asdp_retcode retcode = ASDP_RETCODE_OK;
    ASDP_NULLCHECK(pCone);
    asdp_cone_sdp_sparse *cone = NULL;
    ASDP_INIT(cone, asdp_cone_sdp_sparse, 1);
    ASDP_MEMCHECK(cone);
    ASDP_ZERO(cone, asdp_cone_sdp_sparse, 1);
    *pCone = cone;
    
exit_cleanup:
    return retcode;
}

/** @brief Create a dense sdp cone
 *
 */
extern asdp_retcode sdpDenseConeCreateImpl( void **pConeIn ) {
    asdp_cone_sdp_dense **pCone = (asdp_cone_sdp_dense **)pConeIn;
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    ASDP_NULLCHECK(pCone);
    asdp_cone_sdp_dense *cone = NULL;
    
    ASDP_INIT(cone, asdp_cone_sdp_dense, 1);
    ASDP_MEMCHECK(cone);
    ASDP_ZERO(cone, asdp_cone_sdp_dense, 1);
    *pCone = cone;
    
exit_cleanup:
    return retcode;
}

extern void sdpDenseConeDestroyImpl( void **pConeIn ) {
    asdp_cone_sdp_dense **pCone = (asdp_cone_sdp_dense **)pConeIn;
    if ( !pCone ) {
        return;
    }
    
    sdpDenseConeClearImpl(*pCone);
    ASDP_FREE(*pCone);
    
    return;
}


extern void sdpDenseConeFeatureDetectImpl( void *coneIn, double *rowRHS,
                                           int coneIntFeatures[20], double coneDblFeatures[20] ) {
    /* When there is a single SDP cone. This routine detects if the SDP has the following structures
     
      1. (Almost) no primal interior point
      2. (Almost) no dual interior point
      3. Implied trace bound tr(X) == Z
      4. Implied dual upperbound and lower bound l <= y <= u
      5. No objective
      6. Dense cone
     
     Also the number of different cones will be records as the conic features (including objective)
     */
    
    /* We detect the case of no primal interior by seeking constraints that look like
       trace(a * a' * X) = b \approx 0.0
     */
    
    /* Get statistics */
    asdp_cone_sdp_dense *cone = (asdp_cone_sdp_dense *)coneIn;
    coneIntFeatures[INT_FEATURE_N_ZEORMATS] = cone->sdpConeStats[SDP_COEFF_ZERO];
    coneIntFeatures[INT_FEATURE_N_DSMATS] = cone->sdpConeStats[SDP_COEFF_DENSE];
    coneIntFeatures[INT_FEATURE_N_SPMATS] = cone->sdpConeStats[SDP_COEFF_SPARSE];
    coneIntFeatures[INT_FEATURE_N_SPR1MATS] = cone->sdpConeStats[SDP_COEFF_SPR1];
    coneIntFeatures[INT_FEATURE_N_DSR1MATS] = cone->sdpConeStats[SDP_COEFF_DSR1];
    
    int isNoPrimalInterior = 0;
    int isImpliedTraceX = 0;
    double dImpliedTraceX = 0.0;
    int iUnitCol = 0;
    int *iUnitColIdx = NULL;
    
    /* Detect if there is no primal interior point */
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        if ( sdpDataMatGetRank(cone->sdpRow[iRow]) == 1 &&
            fabs(rowRHS[iRow]) < 1e-03 * sdpDataMatNorm(cone->sdpRow[iRow], FRO_NORM) ) {
            isNoPrimalInterior = 1;
        }
    }
    
    if ( isNoPrimalInterior ) {
        coneIntFeatures[INT_FEATURE_I_NOPINTERIOR] = 1;
    }
    
    /* Detect implied trace bound
       We detect two cases of implied bound.
     
       The first case is the constraint trace(I * X) == a
       The second case is the diagonal constraint diag(X) = d */
    
    /* First detect is there is trace(I * X) == a */
    int isEyeMultiple = 0;
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        isEyeMultiple = sdpDataMatIsEye(cone->sdpRow[iRow], &dImpliedTraceX);
        if ( isEyeMultiple && rowRHS[iRow] / dImpliedTraceX > 0.0 ) {
            isImpliedTraceX = 1;
            dImpliedTraceX = rowRHS[iRow] / dImpliedTraceX;
            break;
        }
    }
    
    /* Then detect if there is diag(X) constraint. These constraints are of rank-one sparse coeffients */
    ASDP_INIT(iUnitColIdx, int, cone->nRow);
    dImpliedTraceX = 0.0;
    if ( iUnitColIdx && !isImpliedTraceX ) {
        for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
            if ( sdpDataMatIsUnitCol(cone->sdpRow[iRow], &iUnitCol) && !iUnitColIdx[iUnitCol] ) {
                iUnitColIdx[iUnitCol] = 1;
                dImpliedTraceX += rowRHS[iRow];
            }
        }
    }
    
    int nSumUnit = 0;
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        nSumUnit += iUnitColIdx[iRow];
    }
    
    if ( nSumUnit == cone->nCol ) {
        isImpliedTraceX = 1;
    }
    
    ASDP_FREE(iUnitColIdx);
    
    if ( isImpliedTraceX ) {
        coneIntFeatures[INT_FEATURE_I_IMPTRACE] = 1;
        coneDblFeatures[DBL_FEATURE_IMPTRACEX] = dImpliedTraceX;
    }
    
    /* Detect if there is no objective */
    if ( sdpDataMatGetType(cone->sdpObj) == SDP_COEFF_ZERO ) {
        coneIntFeatures[INT_FEATURE_I_NULLOBJ] = 1;
    }
    
    /* Detect dense cone */
    if ( cone->sdpConeStats[SDP_COEFF_DENSE] >= 0.7 * cone->nRow ) {
        coneIntFeatures[INT_FEATURE_I_VERYDENSE] = 1;
    }
    
    return;
}


extern asdp_retcode sdpDenseConePresolveImpl( void *coneIn ) {
    asdp_cone_sdp_dense *cone = (asdp_cone_sdp_dense *)coneIn;
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    double *preConeAuxi = NULL;
    /* cone->nCol * cone->nCol for a Full matrix */
    ASDP_INIT(preConeAuxi, double, cone->nCol);
    ASDP_MEMCHECK(preConeAuxi);
    
    /* no low rank*/
//    ASDP_CALL(sdpDataMatBuildUpEigs(cone->sdpObj, preConeAuxi));
    
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        /* Update conic statistics */
        cone->sdpConeStats[sdpDataMatGetType(cone->sdpRow[iRow])] -= 1;
//        ASDP_CALL(sdpDataMatBuildUpEigs(cone->sdpRow[iRow], preConeAuxi));
        cone->sdpConeStats[sdpDataMatGetType(cone->sdpRow[iRow])] += 1;
    }
    
exit_cleanup:
    
    ASDP_FREE(preConeAuxi);
    return retcode;
}

/** @brief Set data into a dense cone
 *
 */
extern asdp_retcode sdpDenseConeProcDataImpl( void *coneIn, int nRow, int nCol,
                                               int *coneMatBeg, int *coneMatIdx, double *coneMatElem ) {
    asdp_cone_sdp_dense *cone = (asdp_cone_sdp_dense *)coneIn;
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    /* Create the cone */
    cone->nRow = nRow;
    cone->nCol = nCol;
    
    ASDP_INIT(cone->sdpRow, sdp_coeff *, nRow);
    ASDP_MEMCHECK(cone->sdpRow);
    
    ASDP_ZERO(cone->sdpConeStats, int, 5);
    
    /* Primal objective */
    ASDP_CALL(sdpDataMatCreate(&cone->sdpObj));
    ASDP_CALL(sdpDataMatSetData(cone->sdpObj, nCol, coneMatBeg[1], coneMatIdx, coneMatElem));
    cone->sdpConeStats[sdpDataMatGetType(cone->sdpObj)] += 1;
    
    /* Constraint data */
    for ( int iRow = 0; iRow < nRow; ++iRow ) {
        ASDP_CALL(sdpDataMatCreate(&cone->sdpRow[iRow]));
        ASDP_CALL(sdpDataMatSetData(cone->sdpRow[iRow], nCol, coneMatBeg[iRow + 2] - coneMatBeg[iRow + 1],
                                     coneMatIdx + coneMatBeg[iRow + 1], coneMatElem + coneMatBeg[iRow + 1]));
        sdp_coeff_type coeffType = sdpDataMatGetType(cone->sdpRow[iRow]);
        cone->sdpConeStats[coeffType] += 1;
    }
        
exit_cleanup:
    return retcode;
}

extern asdp_retcode sdpSparseConeProcDataImpl( void *coneIn, int nRow, int nCol,
                                               int *coneMatBeg, int *coneMatIdx, double *coneMatElem ) {
    asdp_cone_sdp_sparse *cone = (asdp_cone_sdp_sparse *)coneIn;
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    /* Create the cone */
    cone->nRow = nRow;
    cone->nCol = nCol;
    
    ASDP_ZERO(cone->sdpConeStats, int, 5);
    
    /* Primal objective */
    ASDP_CALL(sdpDataMatCreate(&cone->sdpObj));
    ASDP_CALL(sdpDataMatSetData(cone->sdpObj, nCol, coneMatBeg[1], coneMatIdx, coneMatElem));
    cone->sdpConeStats[sdpDataMatGetType(cone->sdpObj)] += 1;
    
    /* Sparse constraint data */
    int nRowElem = 0;
    
    /* Count #Nonzeros */
    nRowElem = csp_nnz_cols(nRow, &coneMatBeg[1]);
    
    ASDP_INIT(cone->sdpRow, sdp_coeff *, nRowElem);
    ASDP_INIT(cone->rowIdx, int, nRowElem);
    ASDP_MEMCHECK(cone->sdpRow);
    ASDP_MEMCHECK(cone->rowIdx);
    
    int nRowElemTmp = 0;
    for ( int iRow = 0; iRow < nRow; ++iRow ) {
        if ( coneMatBeg[iRow + 2] - coneMatBeg[iRow + 1] > 0 ) {
            ASDP_CALL(sdpDataMatCreate(&cone->sdpRow[nRowElemTmp]));
            ASDP_CALL(sdpDataMatSetData(cone->sdpRow[nRowElemTmp], nCol, coneMatBeg[iRow + 2] - coneMatBeg[iRow + 1],
                                         coneMatIdx + coneMatBeg[iRow + 1], coneMatElem + coneMatBeg[iRow + 1]));
            sdp_coeff_type coeffType = sdpDataMatGetType(cone->sdpRow[nRowElemTmp]);
#ifdef ASDP_CONIC_DEBUG
            assert( coeffType != SDP_COEFF_ZERO );
#endif
            cone->sdpConeStats[coeffType] += 1;
            cone->rowIdx[nRowElemTmp] = iRow;
            nRowElemTmp += 1;
        } else {
            cone->sdpConeStats[SDP_COEFF_ZERO] += 1;
        }
    }
    
    cone->nRowElem = nRowElem;
    
#ifdef ASDP_CONIC_DEBUG
    assert( nRowElem == nRowElemTmp );
#endif
    
exit_cleanup:
    
    return retcode;
}

extern asdp_retcode sdpSparseConePresolveImpl( void *coneIn ) {
    asdp_cone_sdp_sparse *cone = (asdp_cone_sdp_sparse *)coneIn;
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    double *preConeAuxi = NULL;
    /* cone->nCol * cone->nCol for a Full matrix */
    ASDP_INIT(preConeAuxi, double, cone->nCol);
    ASDP_MEMCHECK(preConeAuxi);
    
//    ASDP_CALL(sdpDataMatBuildUpEigs(cone->sdpObj, preConeAuxi));
    for ( int iElem = 0; iElem < cone->nRowElem; ++iElem ) {
        cone->sdpConeStats[sdpDataMatGetType(cone->sdpRow[iElem])] -= 1;
//        ASDP_CALL(sdpDataMatBuildUpEigs(cone->sdpRow[iElem], preConeAuxi));
        cone->sdpConeStats[sdpDataMatGetType(cone->sdpRow[iElem])] += 1;
    }
exit_cleanup:
    
    ASDP_FREE(preConeAuxi);
    return retcode;
}
extern void sdpSparseConeFeatureDetectImpl( void *coneIn, double *rowRHS,
                                            int coneIntFeatures[20], double coneDblFeatures[20] ) {
    asdp_cone_sdp_dense *cone = (asdp_cone_sdp_dense *)coneIn;
    coneIntFeatures[INT_FEATURE_N_ZEORMATS] = cone->sdpConeStats[SDP_COEFF_ZERO];
    coneIntFeatures[INT_FEATURE_N_DSMATS] = cone->sdpConeStats[SDP_COEFF_DENSE];
    coneIntFeatures[INT_FEATURE_N_SPMATS] = cone->sdpConeStats[SDP_COEFF_SPARSE];
    coneIntFeatures[INT_FEATURE_N_SPR1MATS] = cone->sdpConeStats[SDP_COEFF_SPR1];
    coneIntFeatures[INT_FEATURE_N_DSR1MATS] = cone->sdpConeStats[SDP_COEFF_DSR1];
    
    return;
}


extern void sdpDenseConeViewImpl( void *coneIn ) {
    asdp_cone_sdp_dense *cone = (asdp_cone_sdp_dense *)coneIn;
    printf("- Dense SDP cone of %d rows. \n", cone->nRow);
    printf("- Objective: \n");
    sdpDataMatView(cone->sdpObj);
    
    printf("- Constraint: \n");
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        sdpDataMatView(cone->sdpRow[iRow]);
    }
    
    printf("- Conic statistics: Zero %d Sp %d Ds %d SpR1 %d DsR1 %d \n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE],
           cone->sdpConeStats[SDP_COEFF_SPR1], cone->sdpConeStats[SDP_COEFF_DSR1]);
    
    return;
}
extern void sdpDenseConeAUVImpl( void *coneIn, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *constrVal, sdp_coeff *UVt, int FLAG_UV ){
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        cone->sdpRow[iRow]->mul_inner_rk_double(cone->sdpRow[iRow]->dataMat, U, V, &constrVal[iRow], UVt->dataMat, UVt->dataType, FLAG_UV);
    }
}

extern void sdpDenseObjAUVImpl( void *coneIn, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *pObj, sdp_coeff *UVt, int FLAG){
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    double temp = 0.0;
    if (FLAG == FLAG_INI){
        cone->sdpObj->mul_inner_rk_double(cone->sdpObj->dataMat, U, V, &temp, UVt->dataMat, UVt->dataType, FLAG);
    }else{
        void *dummy;
        sdp_coeff_type dummyType;
        cone->sdpObj->mul_inner_rk_double(cone->sdpObj->dataMat, U, V, &temp, dummy, dummyType, FLAG);
    }
    
    pObj[0] += temp;
}

extern void sdpDenseConeObjNrm1(void *coneIn, double *objConeNrm1){
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrm1(cone->sdpObj->dataMat, &temp);
    objConeNrm1[0] += temp;
}

extern void sdpDenseConeObjNrm2Square(void *coneIn, double *objConeNrm2Square){
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrm2Square(cone->sdpObj->dataMat, &temp);
    objConeNrm2Square[0] += temp;
}

extern void sdpDenseConeObjNrmInf(void *coneIn, double *objConeNrmInf){
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrmInf(cone->sdpObj->dataMat, &temp);
    objConeNrmInf[0] = ASDP_MAX(temp, objConeNrmInf[0]);
}

extern void sdpDenseConeAddObjCoeff(void *coneIn, sdp_coeff *w_sum){
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    if (w_sum->dataType == SDP_COEFF_DENSE || w_sum->dataType == SDP_COEFF_SPARSE){
        cone->sdpObj->add_sdp_coeff(cone->sdpObj->dataMat, w_sum->dataMat, 1.0, w_sum->dataType);
    }else{
        ASDP_ERROR_TRACE;
        asdp_printf("copy obj sdp coeff to other but `dense, sparse` is developing!\n");
    }
}

extern void sdpDenseConeAddObjCoeffRand(void *coneIn, sdp_coeff *w_sum){
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    if (w_sum->dataType == SDP_COEFF_DENSE || w_sum->dataType == SDP_COEFF_SPARSE){
        cone->sdpObj->add_sdp_coeff(cone->sdpObj->dataMat, w_sum->dataMat, rand(), w_sum->dataType);
    }else{
        ASDP_ERROR_TRACE;
        asdp_printf("copy obj sdp coeff to other but `dense, sparse` is developing!\n");
    }
}


extern void sdpDenseConeDataScale(void *coneIn, double scaleFactorObj, double scaleFactorSDPData){
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    cone->sdpObj->scaleData(cone->sdpObj->dataMat, scaleFactorObj);
    for (int iRow = 0; iRow < cone->nRow; ++iRow){
        cone->sdpRow[iRow]->scaleData(cone->sdpRow[iRow]->dataMat, scaleFactorSDPData);
    }
}

extern void sdpDenseConeNnzStat(void *coneIn, int *stat){
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    for (int iRow = 0; iRow < cone->nRow; ++iRow){
        cone->sdpRow[iRow]->statNnz(stat);
    }
}

extern void sdpDenseConeNnzStatCoeff(void *coneIn, double *stat, int *nnzStat, int *eleStat){
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    int coeffNum = cone->sdpRow[0]->nSDPCol * (cone->sdpRow[0]->nSDPCol + 1) / 2;
    eleStat[0] += coeffNum;
    for (int iRow = 0; iRow < cone->nRow; ++iRow){
        int nnzCoeff = cone->sdpRow[iRow]->getnnz(cone->sdpRow[iRow]->dataMat);
        stat[iRow] = (double) nnzCoeff / (double) coeffNum;
        nnzStat[0] += nnzCoeff;
    }
}



extern void sdpDenseDataWeightSumImpl(void *coneIn, double *weight, sdp_coeff *sdp_data){
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    for (int iRow = 0; iRow < cone->nRow; ++iRow){
        // add data to `sdp_data`;
        cone->sdpRow[iRow]->add_sdp_coeff(cone->sdpRow[iRow]->dataMat, sdp_data->dataMat, weight[iRow], sdp_data->dataType);
    }
}



extern void sdpSparseConeViewImpl( void *coneIn ) {
    asdp_cone_sdp_sparse *cone = (asdp_cone_sdp_sparse *)coneIn;
    printf("- Sparse SDP cone of %d rows and %d nonzeros. \n", cone->nRow, cone->nRowElem);
    printf("- Objective: \n");
    sdpDataMatView(cone->sdpObj);
    
    printf("- Constraint: \n");
    for ( int iRow = 0; iRow < cone->nRowElem; ++iRow ) {
        printf("%d: ", cone->rowIdx[iRow]);
        sdpDataMatView(cone->sdpRow[iRow]);
    }
    
    printf("- Conic statistics: Zero %d Sp %d Ds %d SpR1 %d DsR1 %d \n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE],
           cone->sdpConeStats[SDP_COEFF_SPR1], cone->sdpConeStats[SDP_COEFF_DSR1]);

    return;
}
extern void sdpSparseConeScal( void *coneIn, double dScal ) {
    asdp_cone_sdp_sparse *cone = (asdp_cone_sdp_sparse *)coneIn;
    sdpDataMatScal(cone->sdpObj, dScal);
    return;
}

extern void sdpSparseConeAUVImpl(void *coneIn, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *constrVal, sdp_coeff *UVt, int FLAG_UV){
    asdp_cone_sdp_sparse *cone = (asdp_cone_sdp_sparse *)coneIn;
    for (int iRow = 0; iRow < cone->nRowElem; ++iRow){
        int idx = cone->rowIdx[iRow];
        cone->sdpRow[iRow]->mul_inner_rk_double(cone->sdpRow[iRow]->dataMat, U, V, &constrVal[iRow], UVt->dataMat, UVt->dataType, FLAG_UV);
    }
}

extern void sdpSparseObjAUVImpl(void *coneIn, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *pObj, sdp_coeff *UVt, int FLAG){
    asdp_cone_sdp_sparse *cone = (asdp_cone_sdp_sparse *)coneIn;
    double temp = 0.0;
    sdp_coeff_type dummy;
    cone->sdpObj->mul_inner_rk_double(cone->sdpObj->dataMat, U, V, &temp, UVt, dummy, FLAG_OBJ);
    pObj[0] += temp;
}

extern void sdpSparseConeObjNrm1(void *coneIn, double *objConeNrm1){
    asdp_cone_sdp_sparse *cone = (asdp_cone_sdp_sparse *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrm1(cone->sdpObj->dataMat, &temp);
    objConeNrm1[0] += temp;
}

extern void sdpSparseConeObjNrm2Square(void *coneIn, double *objConeNrm2Square){
    asdp_cone_sdp_sparse *cone = (asdp_cone_sdp_sparse *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrm2Square(cone->sdpObj->dataMat, &temp);
    objConeNrm2Square[0] += temp;
}

extern void sdpSparseConeObjNrmInf(void *coneIn, double *objConeNrmInf){
    asdp_cone_sdp_sparse *cone = (asdp_cone_sdp_sparse *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrmInf(cone->sdpObj->dataMat, &temp);
    objConeNrmInf[0] = ASDP_MAX(temp, objConeNrmInf[0]);
}

extern void sdpSparseConeAddObjCoeff(void *pConeIn, sdp_coeff *w_sum){
    asdp_cone_sdp_sparse *cone = ( asdp_cone_sdp_sparse *)pConeIn;
    cone->sdpObj->add_sdp_coeff(cone->sdpObj->dataMat, w_sum->dataMat, 1.0, w_sum->dataType);
}

extern void sdpSparseConeAddObjCoeffRand(void *pConeIn, sdp_coeff *w_sum){
    asdp_cone_sdp_sparse *cone = ( asdp_cone_sdp_sparse *)pConeIn;
    cone->sdpObj->add_sdp_coeff(cone->sdpObj->dataMat, w_sum->dataMat, rand(), w_sum->dataType);
}

extern void sdpSparseConeNnzStat(void *pConeIn, int *stat){
    asdp_cone_sdp_sparse *cone = ( asdp_cone_sdp_sparse *)pConeIn;
    for (int iRow = 0; iRow < cone->nRowElem; ++iRow){
        cone->sdpRow[iRow]->statNnz(stat);
    }
}

extern void sdpSparseConeNnzStatCoeff(void *pConeIn, double *stat, int *nnzStat, int *eleStat){
    asdp_cone_sdp_sparse *cone = ( asdp_cone_sdp_sparse *)pConeIn;
    int statSum = cone->sdpRow[0]->nSDPCol * (cone->sdpRow[0]->nSDPCol + 1) / 2;
    eleStat[0] += statSum;
    for (int iRow = 0; iRow < cone->nRowElem; ++iRow){
        int nnzCoeff = cone->sdpRow[iRow]->getnnz(cone->sdpRow[iRow]->dataMat);
        stat[cone->rowIdx[iRow]] = (double) nnzCoeff / (double) statSum;
        nnzStat[0] += nnzCoeff;
    }
}

extern void sdpSparseConeDataScale(void *pConeIn, double scaleFactorObj, double scaleFactorSDPData){
    asdp_cone_sdp_sparse *cone = ( asdp_cone_sdp_sparse *)pConeIn;
    cone->sdpObj->scaleData(cone->sdpObj->dataMat, scaleFactorObj);
    for (int iRow = 0; iRow < cone->nRow; ++iRow){
        cone->sdpRow[iRow]->scaleData(cone->sdpRow[iRow]->dataMat, scaleFactorSDPData);
    }
}


extern void sdpSparseDataWeightSumImpl(void *coneIn, double *weight, sdp_coeff *sdp_data){
    asdp_cone_sdp_sparse *cone = ( asdp_cone_sdp_sparse *)coneIn;
    for (int iRow = 0; iRow < cone->nRowElem; ++iRow){
        // add data to `sdp_data`;
        cone->sdpRow[iRow]->add_sdp_coeff(cone->sdpRow[iRow]->dataMat, sdp_data->dataMat, weight[cone->rowIdx[iRow]], sdp_data->dataType);
    }
}



extern void sdpDenseConeScal( void *coneIn, double dScal ) {
    asdp_cone_sdp_dense *cone = (asdp_cone_sdp_dense *)coneIn;
    sdpDataMatScal(cone->sdpObj, dScal);
    return;
}


extern void sdpSparseConeDestroyImpl( void **coneIn ) {
    asdp_cone_sdp_sparse **pCone = (asdp_cone_sdp_sparse **)coneIn;
    if ( !pCone ) {
        return;
    }
    
    sdpSparseConeClearImpl(*pCone);
    ASDP_FREE(*pCone);
    
    return;
}


extern void sdpDenseConeClearImpl( void *coneIn ) {
    asdp_cone_sdp_dense *cone = ( asdp_cone_sdp_dense *)coneIn;
    
    if ( !cone ) {
        return;
    }
    
    sdpDataMatDestroy(&cone->sdpObj);
    
    for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
        sdpDataMatDestroy(&cone->sdpRow[iRow]);
    }
    
    ASDP_FREE(cone->sdpRow);
    
    ASDP_ZERO(cone, asdp_cone_sdp_dense, 1);
    
    return;
}

extern void sdpSparseConeClearImpl( void *coneIn ) {
    asdp_cone_sdp_sparse *cone = (asdp_cone_sdp_sparse *)coneIn;
    
    if ( !cone ) {
        return;
    }
    
    ASDP_FREE(cone->rowIdx);
    sdpDataMatDestroy(&cone->sdpObj);
    
    for ( int iRow = 0; iRow < cone->nRowElem; ++iRow ) {
        sdpDataMatDestroy(&cone->sdpRow[iRow]);
    }
    
    ASDP_FREE(cone->sdpRow);
    
    ASDP_ZERO(cone, asdp_cone_sdp_sparse, 1);
    
    return;
}







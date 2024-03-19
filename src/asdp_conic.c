/** @file asdp\_conic.h
 *  @brief Implement ASDP general conic interface
 */

#ifdef HEADERPATH
#include "interface/asdp_utils.h"
#include "interface/def_asdp_user_data.h"
#include "interface/asdp_user_data.h"
#include "interface/def_asdp_conic.h"
#include "src/def_asdp_sdpdata.h"
#include "src/asdp_direct_linsys.h"
#include "src/asdp_conic.h"
#include "src/asdp_debug.h"
#else
#include "asdp_utils.h"
#include "def_asdp_user_data.h"
#include "asdp_user_data.h"
#include "def_asdp_conic.h"
#include "asdp_conic_sdp.h"
#include "def_asdp_sdpdata.h"
#include "asdp_direct_linsys.h"
#include "asdp_conic.h"
#include "asdp_debug.h"
#include "asdp_conic_lp.h"
#endif

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

extern asdp_retcode AConeCreate( asdp_cone **pACone ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    ASDP_NULLCHECK(pACone);
    asdp_cone *ACone;
    ASDP_INIT(ACone, asdp_cone, 1);
    ASDP_MEMCHECK(ACone);
    ASDP_ZERO(ACone, asdp_cone, 1);
    *pACone = ACone;
    
exit_cleanup:
    
    return retcode;
}

extern asdp_retcode LPSetConeFuncPointer( asdpLPCone *lpCone ){
    asdp_retcode retcode = ASDP_RETCODE_OK;
    lpCone->coneCreate = LPConeCreateImpl;
    lpCone->coneProcData = LPConeProcDataImpl;
    lpCone->conePresolveData = LPConePresolveImpl;
    lpCone->coneDestroyData = LPConeDestroyImpl;
    lpCone->coneView = LPConeViewImpl;
    lpCone->coneAUV = LPConeConeAUV;
    lpCone->objAUV = LPConeObjAUV;
    lpCone->coneObjNrm1 = LPConeObjNrm1;
    lpCone->coneObjNrm2Square = LPConeObjNrm2Square;
    lpCone->coneObjNrmInf = LPConeObjNrmInf;
    lpCone->lpDataWSum = LPConeDataWsum;
    lpCone->objCoeffSum = LPConeDataObjCoeffSum;
    lpCone->coneAUV2 = LPConeConeAUV2;
    return retcode;
}


extern asdp_retcode AConeSetData( asdp_cone *ACone, user_data *usrData) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    ACone->usrData = usrData;
    ACone->cone = AUserDataChooseCone(usrData);
    
     switch ( ACone->cone ) {
         case ASDP_CONETYPE_DENSE_SDP:
             ACone->coneCreate = sdpDenseConeCreateImpl;
             ACone->coneProcData = sdpDenseConeProcDataImpl;
             ACone->conePresolveData = sdpDenseConePresolveImpl;
             ACone->coneDestroyData = sdpDenseConeDestroyImpl;
             ACone->coneScal = sdpDenseConeScal;
             ACone->coneView = sdpDenseConeViewImpl;
             ACone->coneAUV = sdpDenseConeAUVImpl;
             ACone->objAUV = sdpDenseObjAUVImpl;
             ACone->sdpDataWSum = sdpDenseDataWeightSumImpl;
             ACone->getstat = sdpDenseConeFeatureDetectImpl;
             ACone->coneObjNrm1 = sdpDenseConeObjNrm1;
             ACone->coneObjNrm2Square = sdpDenseConeObjNrm2Square;
             ACone->coneObjNrmInf = sdpDenseConeObjNrmInf;
             ACone->addObjCoeff = sdpDenseConeAddObjCoeff;
             ACone->addObjCoeffRand = sdpDenseConeAddObjCoeffRand;
             ACone->nnzStat = sdpDenseConeNnzStat; // high level
             ACone->nnzStatCoeff = sdpDenseConeNnzStatCoeff; // more precision
             ACone->dataScale = sdpDenseConeDataScale;
             break;
         case ASDP_CONETYPE_SPARSE_SDP:
             ACone->coneCreate = sdpSparseConeCreateImpl;
             ACone->coneProcData = sdpSparseConeProcDataImpl;
             ACone->conePresolveData = sdpSparseConePresolveImpl;
             ACone->coneDestroyData = sdpSparseConeDestroyImpl;
             ACone->coneScal = sdpSparseConeScal;
             ACone->coneView = sdpSparseConeViewImpl;
             ACone->coneAUV = sdpSparseConeAUVImpl;
             ACone->objAUV = sdpSparseObjAUVImpl;
             ACone->sdpDataWSum = sdpSparseDataWeightSumImpl;
             ACone->getstat = sdpSparseConeFeatureDetectImpl;
             ACone->coneObjNrm1 = sdpSparseConeObjNrm1;
             ACone->coneObjNrm2Square = sdpSparseConeObjNrm2Square;
             ACone->coneObjNrmInf = sdpSparseConeObjNrmInf;
             ACone->addObjCoeff = sdpSparseConeAddObjCoeff;
             ACone->addObjCoeffRand = sdpSparseConeAddObjCoeffRand;
             ACone->nnzStat = sdpSparseConeNnzStat; // high level
             ACone->nnzStatCoeff = sdpSparseConeNnzStatCoeff; // more precision
             ACone->dataScale = sdpSparseConeDataScale;
             break;
         case ASDP_CONETYPE_SOCP:
             retcode = ASDP_RETCODE_FAILED;
             goto exit_cleanup;
         default:
             retcode = ASDP_RETCODE_FAILED;
             goto exit_cleanup;
     }
    
exit_cleanup:
    
    return retcode;
}

extern void AConeClear( asdp_cone *ACone ) {
    
    if ( !ACone ) {
        return;
    }
    ACone->sdp_slack_var->destroy(&ACone->sdp_slack_var->dataMat);
    ACone->sdp_coeff_w_sum->destroy(&ACone->sdp_coeff_w_sum->dataMat);
    ASDP_FREE(ACone->usrData);
    ASDP_FREE(ACone->sdp_slack_var);
    ASDP_FREE(ACone->sdp_coeff_w_sum);
    ACone->coneDestroyData(&ACone->coneData);
    
    ASDP_ZERO(ACone, asdp_cone, 1);
    return;
}

extern void AConeDestroy( asdp_cone **pACone ) {
    
    if ( !pACone ) {
        return;
    }
    ASDP_FREE(*pACone);
    return;
}

extern void ALpConeDestroy (asdpLPCone **pLPCone){
    if ( !pLPCone ) {
        return;
    }
    ASDP_FREE(*pLPCone);
    return;
}

extern void AConeView( asdp_cone *ACone ) {
     
    ACone->coneView(ACone->coneData);
    
    return;
}

extern asdp_retcode AConeProcData( asdp_cone *ACone ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    user_data *usrData = (user_data *) ACone->usrData;
    
    ASDP_CALL(ACone->coneCreate(&ACone->coneData)); // create a null coneData
    ASDP_CALL(ACone->coneProcData(ACone->coneData, usrData->nConicRow, usrData->nConicCol,
                                   usrData->coneMatBeg, usrData->coneMatIdx, usrData->coneMatElem));
    
exit_cleanup:
    return retcode;
}

extern void AConeDestroyProcData(asdp_cone *ACone){
    ACone->coneDestroyData(&ACone->coneData);
}

extern void destroyForAuxiDense(void **pA){
    sdp_coeff_dense *dense = (sdp_coeff_dense *)pA;
    ASDP_FREE(dense->dsMatElem);
    for (int i = 0; i < dense->nSDPCol; ++i){
        ASDP_FREE(dense->rowCol2NnzIdx[i]);
    }
    ASDP_FREE(dense->rowCol2NnzIdx);
    ASDP_FREE(pA);
}

extern void destroyForAuxiSparse(void **pA){
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *)pA;
    ASDP_FREE(sparse->triMatCol);
    ASDP_FREE(sparse->triMatRow);
    ASDP_FREE(sparse->triMatElem);
    for (int row = 0; row < sparse->nSDPCol; ++row){
        ASDP_FREE(sparse->rowCol2NnzIdx[row]);
    }
    ASDP_FREE(sparse->rowCol2NnzIdx);
    ASDP_FREE(pA);
}

extern asdp_retcode AConePresolveData( asdp_cone *ACone, int Dim) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    ASDP_CALL(ACone->conePresolveData(ACone->coneData));
    // set dense initially
    sdp_coeff *w_sum;
    ASDP_INIT(w_sum, sdp_coeff, 1);
    ASDP_MEMCHECK(w_sum);
    w_sum->dataType = SDP_COEFF_DENSE;
    w_sum->nSDPCol = Dim;
    sdp_coeff_dense *dataMat;
    ASDP_INIT(dataMat, sdp_coeff_dense, 1);
    dataMat->nSDPCol = w_sum->nSDPCol;
    ASDP_INIT(dataMat->dsMatElem, double, Dim * (Dim + 1)/2);
    ASDP_ZERO(dataMat->dsMatElem, double, Dim * (Dim + 1)/2);
    w_sum->dataMat = (void *)dataMat;
    w_sum->mul_rk = dataMatDenseMultiRkMat;
    w_sum->destroy = destroyForAuxiDense;
    w_sum->zeros = dataMatDenseZeros;
    ACone->sdp_coeff_w_sum = w_sum;
    
    sdp_coeff *sdp_obj_sum;
    ASDP_INIT(sdp_obj_sum, sdp_coeff, 1);
    ASDP_MEMCHECK(sdp_obj_sum);
    sdp_obj_sum->dataType = SDP_COEFF_DENSE;
    sdp_obj_sum->nSDPCol = Dim;
    sdp_coeff_dense *dataMatObj;
    ASDP_INIT(dataMatObj, sdp_coeff_dense, 1);
    dataMatObj->nSDPCol = sdp_obj_sum->nSDPCol;
    ASDP_INIT(dataMatObj->dsMatElem, double, Dim * (Dim + 1)/2);
    ASDP_ZERO(dataMatObj->dsMatElem, double, Dim * (Dim + 1)/2);
    sdp_obj_sum->dataMat = (void *)dataMatObj;
    sdp_obj_sum->mul_rk = dataMatDenseMultiRkMat;
    sdp_obj_sum->destroy = destroyForAuxiDense;
    sdp_obj_sum->zeros = dataMatDenseZeros;
    ACone->sdp_obj_sum = sdp_obj_sum;
    
    sdp_coeff *slackVarTemp;
    ASDP_INIT(slackVarTemp, sdp_coeff, 1);
    slackVarTemp->dataType = SDP_COEFF_DENSE;
    slackVarTemp->nSDPCol = Dim;
    sdp_coeff_dense *slackVarSdp_coeff;
    ASDP_INIT(slackVarSdp_coeff, sdp_coeff_dense, 1);
    slackVarSdp_coeff->nSDPCol = slackVarTemp->nSDPCol;
    ASDP_INIT(slackVarSdp_coeff->dsMatElem, double, Dim * (Dim + 1)/2);
    ASDP_ZERO(slackVarSdp_coeff->dsMatElem, double, Dim * (Dim + 1)/2);
    slackVarTemp->dataMat = slackVarSdp_coeff;
    slackVarTemp->mul_rk = dataMatDenseMultiRkMat;
    slackVarTemp->destroy = destroyForAuxiDense;
    slackVarTemp->zeros = dataMatDenseZeros;
    ACone->sdp_slack_var = slackVarTemp;

exit_cleanup:
    return retcode;
}

extern void AConeDestroyPresolveData(asdp_cone *ACone){
    ACone->sdp_coeff_w_sum->destroy(ACone->sdp_coeff_w_sum->dataMat);
    ASDP_FREE(ACone->sdp_coeff_w_sum);
    ACone->sdp_slack_var->destroy(ACone->sdp_slack_var->dataMat);
    ASDP_FREE(ACone->sdp_slack_var);
    ACone->sdp_obj_sum->destroy(ACone->sdp_obj_sum->dataMat);
    ASDP_FREE(ACone->sdp_obj_sum);
}

extern void AConeScalByConstant( asdp_cone *ACone, double dScal ) {
    
    ACone->coneScal(ACone->coneData, dScal);
    return;
}


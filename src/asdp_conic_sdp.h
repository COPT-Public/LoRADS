#ifndef ASDP_CONIC_SDP_H 
#define ASDP_CONIC_SDP_H

#ifdef HEADERPATH
#include "interface/def_asdp_conic.h"
#else
#include "def_asdp_conic.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Dense SDP cone */
extern asdp_retcode sdpDenseConeCreateImpl( void **pConeIn );
extern asdp_retcode sdpDenseConeProcDataImpl( void *cone, int nRow, int nCol,
                                               int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern asdp_retcode sdpDenseConePresolveImpl( void *cone );
extern void sdpDenseConeSetStartImpl( void *cone, double rResi );
extern double sdpDenseConeGetObjNorm( void *cone, int whichNorm );
extern double sdpDenseConeGetCoeffNorm( void *cone, int whichNorm );
extern void sdpDenseConeUpdateImpl( void *cone, double barHsdTau, double *rowDual );
extern asdp_retcode sdpDenseConeRatioTestImpl( void *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, int whichBuffer, double *maxStep );
extern void sdpDenseConeAddSymNnzImpl( void *cone, int iCol, int *schurMatCol );
extern void sdpDenseConeGetSymMapping( void *cone, int iCol, int *schurMatCol );
extern int sdpDenseConeGetDim( void *cone );
extern double sdpDenseConeGetObjNorm( void *cone, int whichNorm );
extern double sdpDenseConeGetCoeffNorm( void *cone, int whichNorm );
extern void sdpDenseConeScal( void *cone, double dScal );
extern asdp_retcode sdpDenseConeInteriorCheck( void *cone, double barHsdTau, double *rowDual, int *isInterior );
extern asdp_retcode sdpDenseConeInteriorCheckExpert( void *cone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef, int whichBuffer, int *isInterior );
extern void sdpDenseConeReduceResidual( void *cone, double resiReduction );
extern void sdpDenseConeSetPerturb( void *cone, double dDualPerturb );
extern asdp_retcode sdpDenseConeGetBarrier( void *cone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet );
extern asdp_retcode sdpDenseConeAddStepToBufferAndCheck( void *cone, double dStep, int whichBuffer, int *isInterior );
extern void sdpDenseConeClearImpl( void *cone );
extern void sdpDenseConeDestroyImpl( void **pCone );
extern void sdpDenseConeFeatureDetectImpl( void *cone, double *rowRHS, int coneIntFeatures[20], double coneDblFeatures[20] );
extern void sdpDenseConeViewImpl( void *cone );
extern void sdpDenseConeAUVImpl( void *coneIn, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *constrVal, sdp_coeff *UVt, int FLAG_UV );
extern void sdpDenseObjAUVImpl( void *coneIn, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *pObj, sdp_coeff *UVt, int FLAG);
extern void sdpDenseConeObjNrm1(void *cone, double *objConeNrm1);
extern void sdpDenseConeObjNrm2Square(void *coneIn, double *objConeNrm2Square);
extern void sdpDenseConeObjNrmInf(void *coneIn, double *objConeNrmInf);
extern void sdpDenseConeAddObjCoeff(void *cone, sdp_coeff *w_sum);
extern void sdpDenseConeAddObjCoeffRand(void *cone, sdp_coeff *w_sum);
extern void sdpDenseConeAddPreprocessRankOneConeDetect(void *coneIn, int *flagRankOne);
extern void sdpDenseConeDataScale(void *coneIn, double scaleFactorObj, double scaleFactorSDPData);
extern void sdpDenseConeNnzStat(void *coneIn, int *stat);
extern void sdpDenseConeNnzStatCoeff(void *coneIn, double *stat, int *nnzStat, int *eleStat);

/* Sparse SDP cone */
extern asdp_retcode sdpSparseConeCreateImpl( void **pCone );
extern asdp_retcode sdpSparseConeProcDataImpl( void *cone, int nRow, int nCol,
                                               int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern asdp_retcode sdpSparseConePresolveImpl( void *cone );
extern void sdpSparseConeSetStartImpl( void *cone, double rResi );
extern double sdpSparseConeGetObjNorm( void *cone, int whichNorm );
extern double sdpSparseConeGetCoeffNorm( void *cone, int whichNorm );
extern void sdpSparseConeScal( void *cone, double dScal );
extern void sdpSparseConeUpdateImpl( void *cone, double barHsdTau, double *rowDual );
extern asdp_retcode sdpSparseConeRatioTestImpl( void *cone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, int whichBuffer, double *maxStep );
extern int sdpSparseConeGetDim( void *cone );
extern asdp_retcode sdpSparseConeGetKKT( void *cone, void *kkt, int typeKKT );
extern asdp_retcode sdpSparseConeGetKKTByFixedStrategy( void *cone, void *kkt, int typeKKT, int kktStrategy );
extern int64_t sdpSparseConeGetSymNnzImpl( void *cone );
extern void sdpSparseConeAddSymNnzImpl( void *cone, int iCol, int *schurMatCol );
extern void sdpSparseConeGetSymMapping( void *cone, int iCol, int *schurMatCol );
extern asdp_retcode sdpSparseConeInteriorCheck( void *cone, double barHsdTau, double *rowDual, int *isInterior );
extern asdp_retcode sdpSparseConeInteriorCheckExpert( void *cone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef, int whichBuffer, int *isInterior );
extern void sdpSparseConeReduceResidual( void *cone, double resiReduction );
extern void sdpSparseConeSetPerturb( void *cone, double dDualPerturb );
extern asdp_retcode sdpSparseConeGetBarrier( void *cone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet );
extern asdp_retcode sdpSparseConeAddStepToBufferAndCheck( void *cone, double dStep, int whichBuffer, int *isInterior );
extern void sdpSparseConeClearImpl( void *cone );
extern void sdpSparseConeDestroyImpl( void **pCone );
extern void sdpSparseConeFeatureDetectImpl( void *cone, double *rowRHS, int coneIntFeatures[20], double coneDblFeatures[20] );
extern void sdpSparseConeViewImpl( void *cone );
extern void sdpSparseConeAUVImpl(void *coneIn, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *constrVal, sdp_coeff *UVt, int FLAG_UV);
extern void sdpSparseObjAUVImpl(void *coneIn, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *pObj, sdp_coeff *UVt, int FLAG);
extern void sdpSparseConeObjNrm1(void *cone, double *objConeNrm1);
extern void sdpSparseConeObjNrm2Square(void *coneIn, double *objConeNrm2Square);
extern void sdpSparseConeObjNrmInf(void *coneIn, double *objConeNrmInf);
extern void sdpSparseConeAddObjCoeff(void *cone, sdp_coeff *w_sum);
extern void sdpSparseConeAddObjCoeffRand(void *cone, sdp_coeff *w_sum);
extern void sdpSparseConeNnzStat(void *pConeIn, int *stat);
extern void sdpSparseConeNnzStatCoeff(void *pConeIn, double *stat, int *nnzStat, int *eleStat);
extern void sdpSparseConeFullDenseSDPData(void *pConeIn);
extern void sdpSparseConeAddPreprocessRankOneConeDetect(void *pConeIn, int *flagRankOne);
extern void sdpSparseConeDataScale(void *pConeIn, double scaleFactorObj, double scaleFactorSDPData);
extern void sdpSparseDataWeightSumImpl(void *cone, double *weight, sdp_coeff *sdp_data);
extern void sdpDenseDataWeightSumImpl(void *cone, double *weight, sdp_coeff *sdp_data);


#ifdef __cplusplus
}
#endif

#endif /* asdp_conic_sdp_h */

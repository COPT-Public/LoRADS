#ifndef ASDP_CONIC_LP_H
#define ASDP_CONIC_LP_H

#include <stdio.h>

#ifdef HEADERPATH
#include "interface/def_asdp_conic.h"
#else
#include "def_asdp_conic.h"
#endif



#ifdef __cplusplus
extern "C" {
#endif

extern asdp_retcode LPConeCreateImpl( void **pCone );
extern asdp_retcode LPConeProcDataImpl( void *cone, int nRow, int nCol, int *coneMatBeg, int *coneMatIdx, double *coneMatElem );
extern asdp_retcode LPConePresolveImpl( void *cone );
extern double LPConeGetObjNorm( asdp_cone_lp *cone, int whichNorm );
extern double LPConeGetCoeffNorm( asdp_cone_lp *cone, int whichNorm );
extern void LPConeClearImpl( void *cone );
extern void LPConeDestroyImpl( void **pCone );
extern void LPConeConeAUV(void *coneIn, asdp_rk_mat_lp *uLp, asdp_rk_mat_lp *vLp, double *AUV, int iCol);
extern void LPConeConeAUV2(void *coneIn, double *uvLp, double *AUVSum);
extern void LPConeObjAUV(void *coneIn, asdp_rk_mat_lp *uLp, asdp_rk_mat_lp *vLp, double *cUV);
extern void LPConeObjNrm1(void *coneIn, double *nrm1Val, int nLpCols);
extern void LPConeObjNrm2Square(void *coneIn, double *nrm2Val, int nLpCols);
extern void LPConeObjNrmInf(void *coneIn, double *nrmInf, int iCol);
extern void LPConeDataWsum(void *coneIn, double *weight, double *wSum, int iCol);
extern void LPConeDataObjCoeffSum(void *coneIn, double *res, int iCol);
extern void LPConeViewImpl( void *cone );
#ifdef __cplusplus
}
#endif


#endif /* ASDP_CONIC_LP_H */

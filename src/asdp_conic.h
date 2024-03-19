#ifndef ASDP_CONIC_H
#define ASDP_CONIC_H

#ifdef HEADERPATH
#include "interface/def_asdp_conic.h"
#include "interface/asdp_user_data.h"
#include "interface/def_asdp.h"
#else
#include "def_asdp_conic.h"
#include "asdp_user_data.h"
#include "def_asdp.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define KKT_TYPE_INFEASIBLE  (0)
#define KKT_TYPE_CORRECTOR   (1)
#define KKT_TYPE_HOMOGENEOUS (2)

#define ABS_NORM    (1)
#define FRO_NORM    (2)

#define BUFFER_DUALVAR   (0)
#define BUFFER_DUALCHECK (1)
#define BUFFER_DUALSTEP  (2)
extern asdp_retcode AConeCreate( asdp_cone **pACone );
extern asdp_retcode LPSetConeFuncPointer( asdpLPCone *lpCone );
extern asdp_retcode AConeSetData( asdp_cone *ACone, user_data *coneData );
extern asdp_retcode AConeProcData( asdp_cone *ACone );
extern void AConeDestroyProcData(asdp_cone *ACone);
extern void destroyUsrdata (user_data *usrData);
extern asdp_retcode AConePresolveData( asdp_cone *ACone, int Dim );
extern void AConeDestroyPresolveData(asdp_cone *ACone);
extern void destroyForAuxiSparse(void **pA);
extern void destroyForAuxiDense(void **pA);
extern void AConeClear( asdp_cone *ACone );
extern void AConeDestroy( asdp_cone **pACone );
extern void ALpConeDestroy (asdpLPCone **pLPCone);
extern void AConeDetectFeature( asdp_cone *ACone, double *rowRHS, int coneIntFeatures[20], double coneDblFeatures[20] );
extern void AConeView( asdp_cone *ACone );

extern void AConeSetStart( asdp_cone *ACone, double dConeStartVal );
extern void AConeUpdate( asdp_cone *ACone, double barHsdTau, double *rowDual );
extern asdp_retcode AConeRatioTest( asdp_cone *ACone, double barHsdTauStep, double *rowDualStep, double dAdaRatio, int whichBuffer, double *maxStep );
extern int64_t AConeGetSymNnz( asdp_cone *ACone );
extern void AConeAddSymNz( asdp_cone *ACone, int iCol, int *schurMatCol );
extern void AConeGetSymMapping( asdp_cone *ACone, int iCol, int *schurMatCol );
extern int AConeGetDim( asdp_cone *ACone );
extern double AConeGetCoeffNorm( asdp_cone *ACone, int whichNorm );
extern double AConeGetObjNorm( asdp_cone *ACone, int whichNorm );
extern asdp_retcode AConeBuildSchurComplement( asdp_cone *ACone, void *schurMat, int typeKKT );
extern asdp_retcode AConeBuildSchurComplementFixed( asdp_cone *ACone, void *schurMat, int typeKKT, int kktStrategy );
extern asdp_retcode AConeGetLogBarrier( asdp_cone *ACone, double barHsdTau, double *rowDual, int whichBuffer, double *logdet );
extern asdp_retcode AConeAddStepToBufferAndCheck( asdp_cone *ACone, double dStep, int whichBuffer, int *isInterior );
extern asdp_retcode AConeCheckIsInterior( asdp_cone *ACone, double barHsdTau, double *rowDual, int *isInterior );
extern asdp_retcode AConeCheckIsInteriorExpert( asdp_cone *ACone, double dCCoef, double dACoefScal, double *dACoef, double dEyeCoef, int whichBuffer, int *isInterior );
extern void AConeReduceResi( asdp_cone *ACone, double resiReduction );
extern void AConeSetPerturb( asdp_cone *ACone, double dPerturb );
extern int AConePFeasSolFound( asdp_cone *ACone, double barHsdTauStep, double *rowDualStep );
extern void AConeGetPrimal( asdp_cone *ACone, double dBarrierMu, double *dRowDual, double *dRowDualStep, double *dConePrimal, double *dConePrimal2 );
extern void AConeScalByConstant( asdp_cone *ACone, double dScal );

#ifdef __cplusplus
}
#endif

#endif /* asdp_conic_h */

#ifndef ASDP_ALGO_H
#define ASDP_ALGO_H

#ifdef HEADERPATH
#include "interface/asdp.h"
#else
#include "asdp.h"
#include "asdp_user_data.h"
#include "def_asdp.h"
#endif

typedef struct
{
    asdp_cone *ACone;
    asdp_rk_mat_dense *noUpdateVar;
    asdp_rk_mat_dense *UpdateVarShell;
    double *weight;
} admmCG;

#ifdef __cplusplus
extern "C"
{
#endif

#define ASDP_DIMAC_ERROR_CONSTRVIO (0)
#define ASDP_DIMAC_ERROR_ASSYMMETRY (1)
#define ASDP_DIMAC_ERROR_DUALFEASIBLE (2)
#define ASDP_DIMAC_ERROR_PDGAP (3)
#define ASDP_ERROR_DIFF_GAP (4)
#define PI (3.1415926)
#define pi01 (0.31415926)

    extern asdp_retcode ASDPSolverCreate(asdp **pHSolver);
    extern void ASDPDestroy(asdp **pHSolver);
    extern asdp_retcode ASDPSetLpCone(asdpLPCone *lpCone, int nRows, int nLpCols, int *LpMatBeg, int *LpMatIdx, double *LpMatElem);
    extern asdp_retcode ASDPInitConeData(asdp *ASolver, user_data **SDPDatas, double **coneMatElem, int **coneMatBeg, int **coneMatIdx, int *BlkDims, int nConstrs, int nBlks, int nLpCols, int *LpMatBeg, int *LpMatIdx, double *LpMatElem);
    extern double normalRandom();
    extern void ASDPDestroyConeData(asdp *ASolver);
    extern void lpRandom(double *data, int n);
    extern asdp_retcode ASDPInitBMVars(asdp *ASolver, int *BlkDims, int nBlks, int nLpCols);
    extern void ASDPDestroyBMVars(asdp *ASolver);
    extern asdp_retcode ASDPInitADMMVars(asdp *ASolver, int *BlkDims, int nBlks, int nLpCols);
    extern void ASDPDestroyADMMVars(asdp *ASolver);
    extern asdp_retcode ASDPInitSolver(asdp *ASolver, int nRows, int nCones, int *blkDims, int nLpCols, double rho, double rhoMax, int maxIter, int strategy);
    extern asdp_retcode ASDPDetermineRank(asdp *ASolver, int *blkDims, double timesRank, int rankSpecify);
    extern void ASDPDestroySolver(asdp *ASolver);
    extern void ASDPDestroyRankElem(asdp *ASolver);
    extern asdp_retcode ASDPCGCreate(asdp *ASolver);
    extern void ASDPCGReAllocate(asdp *ASolver);
    extern asdp_retcode ASDPLanczosCreate(asdp *ASolver);
    extern asdp_retcode ASDPDirectCreate(asdp *ASolver);
    extern void ASDPCGDetroy(asdp *ASolver);
    extern void ASDPLanczosDetroy(asdp *ASolver);
    extern asdp_retcode ASDPSetCone(asdp *ASolver, int iCone, void *userCone);
    extern asdp_retcode ASDPSetVarUVS(asdp *ASolver, int iBlk, void *U, int FALG_UV);
    extern asdp_retcode ASDP_ONE_rk_MAT(asdp_rk_mat_dense *U);
    extern asdp_retcode ASDP_HALF_rk_MAT(asdp_rk_mat_dense *U);
    extern asdp_retcode ASDP_RANDOM_rk_MAT(asdp_rk_mat_dense *U);
    extern void ASDPSetDualObjective(asdp *ASolver, double *dObj);
    extern void ASDP_linSysSetPPrint(asdp_cone *ACone, int iBlk);
    extern void ASDPRkMatSub(asdp_rk_mat_dense *A, asdp_rk_mat_dense *B, double alpha, asdp_rk_mat_dense *S, int flag_UV);
    extern void ASDPInitConstrVal(asdp_cone *ACone, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *constrVal, int flag_UV);
    extern void ASDPInitConstrValAll(asdp *ASolver, asdp_rk_mat_lp *uLpDummy, asdp_rk_mat_lp *vLpDummy, asdp_rk_mat_dense **U, asdp_rk_mat_dense **V);
    extern void ASDPInitConstrValAllLP(asdp *ASolver, asdp_rk_mat_lp *uLp, asdp_rk_mat_lp *vLp, asdp_rk_mat_dense **U, asdp_rk_mat_dense **V);
    extern void ASDPUpdateConstrVal(asdp_cone *ACone, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, constrValStruct *constrVal, int FLAG_UV);
    extern void ASDPUpdateConstrValCG(asdp_cone *ACone, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *weight, int FLAG_UV);
    extern void ASDPUpdateConstrValLP(asdpLPCone *lp_cone, asdp_rk_mat_lp *uLp, asdp_rk_mat_lp *vLp, constrValStruct *constrVal, int iCol);
    extern void ASDPInitConstrValSum(asdp *ASolver);
    extern void ASDPInitConstrValSumLP(asdp *ASolver);
    extern void ASDPConstrValSumBMtemp(asdp *ASolver, double *q12);
    extern void ASDPDetectSparsity(asdp *ASolver);
    extern asdp_retcode ASDPInitM1M2Temp(asdp *ASolver);
    extern asdp_retcode ASDPInitAuxiCri(asdp *ASolver);
    extern asdp_retcode ASDPInitConeAUorV(asdp *ASolver);
    extern void AConeDestroyAUorV(asdp *ASolver);
    extern void ASDPCalObj(asdp *ASolver, int FLAG_BM_USEV);
    extern void ASDPCalObjLP(asdp *ASolver, int FLAG_BM_USEV);
    extern void BMCalq12p12(asdp *ASolver, asdp_rk_mat_lp *rLpDummy, asdp_rk_mat_lp *dDummy, asdp_rk_mat_dense **R, asdp_rk_mat_dense **D, double *q1, double *q2, double *p12);
    extern void BMCalq12p12LP(asdp *ASolver, asdp_rk_mat_lp *rLp, asdp_rk_mat_lp *d, asdp_rk_mat_dense **R, asdp_rk_mat_dense **D, double *q1, double *q2, double *p12);
    extern int BMLineSearch(double rho, int n, double *lambd, double p1, double p2, double *q0, double *q1, double *q2, double *tau);
    extern void ASDPCalDualObj(asdp *ASolver);
    extern asdp_retcode ASDPSumSDPData(asdp *ASolver);
    extern void ASDPDetectSparsityOfSumSDP(asdp *ASolver);
    extern void ASDPAddPreprocessRankOneConeDetect(asdp *ASolver);
    extern void ASDPDetectSparsityOfSumObjAndSDPData(asdp *ASolver);
    extern void ASDPModifySDPDataType(asdp *ASolver);
    extern void ASDPNrm1Obj(asdp *ASolver);
    extern void ASDPNrm2Obj(asdp *ASolver);
    extern void ASDPNrmInfObj(asdp *ASolver);
    extern void ASDPSet_w_coeff(asdp *ASolver);
    extern asdp_retcode ASDPUpdateSdpDataType(asdp *ASolver);
    extern void ASDPLinSetB(double *b, asdp_rk_mat_dense *M2, double rho);
    extern void ASDPUpdateUV(asdp *ASolver, int iBlk, int flagUV);
    extern void ASDPUpdateUVLP(asdp *ASolver, int iCol, int flagUV);
    extern void valRes(void *constrVal, constrType type, double **res);
    extern void addDense(double *alpha, void *constrVal, double *vec);
    extern void addSparse(double *alpha, void *constrVal, double *vec);
    extern void zeroDense(void *constrVal);
    extern void zeroSparse(void *constrVal);
    extern void ASDPUpdateS(asdp *ASolver, int iBlk);
    extern void ASDPUpdateSLP(asdp *ASolver, int iCol);
    extern void ASDPUpdateDualVar(asdp *ASolver);
    extern asdp_retcode ASDPUpdateDimacError(asdp *ASolver);
    extern void ASDPNuclearNorm(asdp *ASolver);
    extern asdp_retcode ASDPUpdateDimacErrorLP(asdp *ASolver);
    extern void ASDPUpdateSDPLPVar(asdp *ASolver);
    extern void ASDPUpdateSDPVar(asdp *ASolver);
    extern void BMSetGrad(asdp *ASolver, asdp_cone *ACone, asdp_rk_mat_dense *R, asdp_rk_mat_dense *Grad, int iCone);
    extern void BMSetGradLP(asdp *ASolver, asdpLPCone *lp_cone, asdp_rk_mat_lp *rLp, asdp_rk_mat_lp *Grad, int iCol);
    extern void BMCalGrad(asdp *ASolver, asdp_rk_mat_lp *rLpDummy, asdp_rk_mat_lp *gradLpDummy, asdp_rk_mat_dense **R, asdp_rk_mat_dense **Grad, double *lagNormSquare);
    extern void BMCalGradLP(asdp *ASolver, asdp_rk_mat_lp *rLp, asdp_rk_mat_lp *gradLp, asdp_rk_mat_dense **R, asdp_rk_mat_dense **Grad, double *lagNormSquare);
    extern void LBFGSDirection(asdp *ASolver, lbfgs_node *head, asdp_rk_mat_lp *gradLpDummy, asdp_rk_mat_lp *dDummy, asdp_rk_mat_dense **Grad, asdp_rk_mat_dense **D, int innerIter);
    extern void LBFGSDirectionLP(asdp *ASolver, lbfgs_node *head, asdp_rk_mat_lp *gradLp, asdp_rk_mat_lp *d, asdp_rk_mat_dense **Grad, asdp_rk_mat_dense **D, int innerIter);
    extern void LBFGSDirectionUseGrad(asdp *ASolver, asdp_rk_mat_lp *dDummy, asdp_rk_mat_lp *gradLPDummy, asdp_rk_mat_dense **D, asdp_rk_mat_dense **Grad);
    extern void LBFGSDirectionUseGradLP(asdp *ASolver, asdp_rk_mat_lp *d, asdp_rk_mat_lp *gradLP, asdp_rk_mat_dense **D, asdp_rk_mat_dense **Grad);
    extern void SetyAsNegGrad(asdp *ASolver, asdp_rk_mat_lp *gradLpDummy, asdp_rk_mat_dense **Grad);
    extern void SetyAsNegGradLP(asdp *ASolver, asdp_rk_mat_lp *gradLp, asdp_rk_mat_dense **Grad);
    extern void BMupdateVar(asdp *ASolver, asdp_rk_mat_lp *rLpDummy, asdp_rk_mat_lp *dDummy, asdp_rk_mat_dense **R, asdp_rk_mat_dense **D, double tau);
    extern void BMupdateVarLP(asdp *ASolver, asdp_rk_mat_lp *rLp, asdp_rk_mat_lp *d, asdp_rk_mat_dense **R, asdp_rk_mat_dense **D, double tau);
    extern void setlbfgsHisTwo(asdp *ASolver, asdp_rk_mat_lp *gradLpDummy, asdp_rk_mat_lp *dDummy, asdp_rk_mat_dense **Grad, asdp_rk_mat_dense **D, double tau);
    extern void setlbfgsHisTwoLP(asdp *ASolver, asdp_rk_mat_lp *gradLp, asdp_rk_mat_lp *d, asdp_rk_mat_dense **Grad, asdp_rk_mat_dense **D, double tau);

    extern int AUG_RANK(asdp *ASolver, int *BlkDims, int nBlks, double aug_factor);

    extern void copyRtoV(asdp_rk_mat_lp *rLpDummy, asdp_rk_mat_lp *vlpDummy, asdp_rk_mat_dense **R, asdp_rk_mat_dense **V, int nCones);

    extern void copyRtoVLP(asdp_rk_mat_lp *rLp, asdp_rk_mat_lp *vlp, asdp_rk_mat_dense **R, asdp_rk_mat_dense **V, int nCones);

    extern void dualGradStandardSDP(asdp *ASolver, double *lambda, double *grad, asdp_rk_mat_dense **Udummy, asdp_rk_mat_dense **Vdummy);

    extern void dualGradStandardSDPLP(asdp *ASolver, double *lambda, double *grad, asdp_rk_mat_dense **Udummy, asdp_rk_mat_dense **Vdummy);

    extern void dualGradNoncvxSDP(asdp *ASolver, double *lambda, double *grad, asdp_rk_mat_dense **Udummy, asdp_rk_mat_dense **V);

    extern void dualGradNoncvxSDPLP(asdp *ASolver, double *lambda, double *grad, asdp_rk_mat_dense **Udummy, asdp_rk_mat_dense **V);
    extern int CheckAllRankMax(asdp *asolver, double aug_factor);
#ifdef __cplusplus
}
#endif

#endif /* asdp_algo_h */

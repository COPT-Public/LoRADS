#ifndef LORADS_SDP_CONIC
#define LORADS_SDP_CONIC


#include "def_lorads_sdp_conic.h"

extern void LORADSSetCone(lorads_solver *ASolver, lorads_int iCone, void *userCone);
extern void AConeProcData(lorads_sdp_cone *ACone );
extern void AConePresolveData( lorads_sdp_cone *ACone, lorads_int Dim);
extern void LORADSSumSDPData(lorads_solver *ASolver);
extern void destroyForAuxiSparse(void **pA);
extern void AConeDenseDetectSparsity(sdp_coeff **sdp_coeff_w_sum_pointer);
extern void LORADSDetectSparsityOfSumSDP(lorads_solver *ASolver);
extern void AConeDestroyPresolveData(lorads_sdp_cone *ACone);
extern void sdpDenseConeClearImpl( void *cone );
extern void sdpDenseConeDestroyImpl( void **pCone );
extern void sdpDenseConeFeatureDetectImpl( void *cone, double *rowRHS, lorads_int coneIntFeatures[20], double coneDblFeatures[20] );
extern void sdpDenseConeViewImpl( void *cone );
extern void sdpDenseConeAUVImpl( void *coneIn, lorads_sdp_dense *U, lorads_sdp_dense *V, double *constrVal, sdp_coeff *UVt);
extern void sdpDenseObjAUVImpl( void *coneIn, lorads_sdp_dense *U, lorads_sdp_dense *V, double *pObj, sdp_coeff *UVt);
extern void sdpDenseConeObjNrm1(void *cone, double *objConeNrm1);
extern void sdpDenseConeObjNrm2Square(void *coneIn, double *objConeNrm2Square);
extern void sdpDenseConeObjNrmInf(void *coneIn, double *objConeNrmInf);
extern void sdpDenseConeAddObjCoeff(void *cone, sdp_coeff *w_sum);
extern void sdpDenseConeAddObjCoeffRand(void *cone, sdp_coeff *w_sum);
extern void sdpDenseConeDataScale(void *coneIn, double scaleFactorSDPData);
extern void sdpDenseConeNnzStat(void *coneIn, lorads_int *stat);
extern void sdpDenseConeNnzStatCoeff(void *coneIn, double *stat, lorads_int *nnzStat, lorads_int *eleStat);

/* Declaration of ARPACK functions */
extern void dsaupd_(int *ido, char *bmat,  lorads_int *n, char *which,  int *nev, double *tol, double *resid,
                     int *ncv, double *v,  lorads_int *ldv,  int *iparam,  int *ipntr, double *workd,
                    double *workl,  int *lworkl,  int *info);

extern void dseupd_(int *rvec, char *HowMny,  int *select, double *d, double *z,  lorads_int *ldz, double *sigma,
                    char *bmat,  lorads_int *n, char *which,  int *nev, double *tol, double *resid,
                     int *ncv, double *v,  lorads_int *ldv,  int *iparam,  int *ipntr, double *workd,
                    double *workl,  int *lworkl,  int *info);

int dual_infeasible(void (*matvec) (void *M, double *x, double *y, lorads_int n), void *M, double *res, lorads_int n);

#endif
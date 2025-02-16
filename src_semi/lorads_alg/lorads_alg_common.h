#ifndef LORADS_ALG_COMMON_H
#define LORADS_ALG_COMMON_H

#include "def_lorads_elements.h"
#include "lorads_solver.h"


extern void LORADSInitConstrValAll(lorads_solver *ASolver, lorads_lp_dense *uLpDummy, lorads_lp_dense *vLpDummy, lorads_sdp_dense **U, lorads_sdp_dense **V);
extern void LORADSInitConstrValAllLP(lorads_solver *ASolver, lorads_lp_dense *uLp, lorads_lp_dense *vLp, lorads_sdp_dense **U, lorads_sdp_dense **V);
extern void LORADSInitConstrValSum(lorads_solver *ASolver);
extern void LORADSInitConstrValSumLP(lorads_solver *ASolver);
extern void LORADSUpdateSDPVar(lorads_solver *ASolver, double rho, double CG_tol, lorads_int CG_maxIter);
extern void LORADSUpdateSDPLPVar(lorads_solver *ASolver, double rho, double CG_tol, lorads_int CG_maxIter);
extern void LORADSUVt(sdp_coeff *UVt_w_sum, lorads_sdp_dense *U, lorads_sdp_dense *V);
extern void copyRtoV(lorads_lp_dense *rLpDummy, lorads_lp_dense *vlpDummy, lorads_sdp_dense **R, lorads_sdp_dense **V, lorads_int nCones);
extern void copyRtoVLP(lorads_lp_dense *rLp, lorads_lp_dense *vlp, lorads_sdp_dense **R, lorads_sdp_dense **V, lorads_int nCones);
extern void LORADSUpdateDimacsErrorALM(lorads_solver *ASolver, lorads_sdp_dense **R, lorads_sdp_dense **R2, lorads_lp_dense *r, lorads_lp_dense *r2);
extern void LORADSUpdateDimacsErrorALMLP(lorads_solver *ASolver, lorads_sdp_dense **R, lorads_sdp_dense **R2, lorads_lp_dense *r, lorads_lp_dense *r2);
extern void LORADSUpdateDimacsErrorADMM(lorads_solver *ASolver, lorads_sdp_dense **U, lorads_sdp_dense **V, lorads_lp_dense *u, lorads_lp_dense *v);
extern void LORADSUpdateDimacsErrorADMMLP(lorads_solver *ASolver, lorads_sdp_dense **U, lorads_sdp_dense **V, lorads_lp_dense *u, lorads_lp_dense *v);
extern void LORADSNrmInfObj(lorads_solver *ASolver);
extern void LORADSUpdateDualVar(lorads_solver *ASolver, double rho);
extern void LORADSCalDualObj(lorads_solver *ASolver);
extern void LORADSNuclearNorm(lorads_solver *ASolver);
extern void LORADSCheckDimacErrALMCriteria(lorads_solver *ASolver);
extern void LORADSObjConstrValAll(lorads_solver *ASolver, lorads_sdp_dense **U, lorads_sdp_dense **V, double *UVobjVal);
extern void LORADSObjConstrValAllLP(lorads_solver *ASolver, lorads_lp_dense *uLp, lorads_lp_dense *vLp, lorads_sdp_dense **U, lorads_sdp_dense **V, double *UVobjVal);
extern lorads_int LORADSCheckSolverStatus(lorads_solver *ASolver);


#endif
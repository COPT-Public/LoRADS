#ifndef LORADS_ADMM_H
#define LORADS_ADMM_H


#include "lorads_solver.h"
#include "lorads.h"

extern lorads_int LORADSADMMOptimize(lorads_params *params, lorads_solver *ASolver, lorads_admm_state *admm_iter_state, lorads_int iter_celling, double timeSolveStart);
extern lorads_int LORADSADMMOptimize_reopt(lorads_params *params, lorads_solver *ASolver, lorads_admm_state *admm_iter_state, lorads_int iter_celling, double timeSolveStart);
extern void LORADSCalObjUV_ADMM(lorads_solver *ASolver);
extern void LORADSCalObjUV_ADMM_LP(lorads_solver *ASolver);
extern void averageUV(lorads_sdp_dense *U, lorads_sdp_dense *V, lorads_sdp_dense *UVavg);
extern void averageUVLP(lorads_lp_dense *ulp, lorads_lp_dense *vlp, lorads_lp_dense *uvavg);
extern void LORADSUpdateSDPVarOne(lorads_solver *ASolver, lorads_sdp_dense *updateVar, lorads_sdp_dense *noUpdateVar, lorads_int iCone, double rho, double CG_tol, lorads_int CG_maxIter);
extern void LORADSUpdateSDPVarOne_positive_S(lorads_solver *ASolver, lorads_sdp_dense *updateVar, lorads_sdp_dense *noUpdateVar, lorads_sdp_dense *S, lorads_int iCone, double rho, double CG_tol, lorads_int CG_maxIter);
extern void LORADSUpdateSDPVarOne_negative_S(lorads_solver *ASolver, lorads_sdp_dense *updateVar, lorads_sdp_dense *noUpdateVar, lorads_sdp_dense *S, lorads_int iCone, double rho, double CG_tol, lorads_int CG_maxIter);


extern void LORADSUpdateLPVarOne(lorads_solver *ASolver,  double *UpdateVar, double *noUpdateVar, lorads_int iCol, double rho);
extern void LORADSUpdateLPVarOne_positive_S(lorads_solver *ASolver,  double *UpdateVar, double *noUpdateVar, lorads_int iCol, double rho, double *sLp);
extern void LORADSUpdateLPVarOne_negative_S(lorads_solver *ASolver,  double *UpdateVar, double *noUpdateVar, lorads_int iCol, double rho, double *sLp);
#endif
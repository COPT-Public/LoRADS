#ifndef LORADS_ALM_H
#define LORADS_ALM_H


#include "def_lorads_elements.h"
#include "lorads_solver.h"

extern void ALMCalGradLP(lorads_solver *ASolver, lorads_lp_dense *rLp, lorads_lp_dense *gradLp, lorads_sdp_dense **R, lorads_sdp_dense **Grad, double *lagNormSquare, double rho);
extern void ALMSetGradLP(lorads_solver *ASolver, lorads_lp_cone *lp_cone, lorads_lp_dense *rLp, lorads_lp_dense *Grad, lorads_int iCol, double rho);
extern void ALMCalGrad(lorads_solver *ASolver, lorads_lp_dense *rLpDummy, lorads_lp_dense *gradLpDummy, lorads_sdp_dense **R, lorads_sdp_dense **Grad, double *lagNormSquare, double rho);
extern void ALMSetGrad(lorads_solver *ASolver, lorads_sdp_cone *ACone, lorads_sdp_dense *R, lorads_sdp_dense *Grad, lorads_int iCone, double rho);

double LORADSnthroot(double base, lorads_int n);
extern lorads_int LORADScubic_equation(double a, double b, double c, double d, double *res);
extern lorads_int ALMLineSearch(double rho, lorads_int n, double *lambd, double p1, double p2, double *q0, double *q1, double *q2, double *tau);
extern void LBFGSDirection(lorads_params *params, lorads_solver *ASolver, lbfgs_node *head, lorads_lp_dense *gradLpDummy, lorads_lp_dense *dDummy, lorads_sdp_dense **Grad, lorads_sdp_dense **D, lorads_int innerIter);
extern void LBFGSDirectionLP(lorads_params *params, lorads_solver *ASolver, lbfgs_node *head, lorads_lp_dense *gradLp, lorads_lp_dense *d, lorads_sdp_dense **Grad, lorads_sdp_dense **D, lorads_int innerIter);
extern void LBFGSDirectionUseGrad(lorads_solver *ASolver, lorads_lp_dense *dDummy, lorads_lp_dense *gradLPDummy, lorads_sdp_dense **D, lorads_sdp_dense **Grad);
extern void LBFGSDirectionUseGradLP(lorads_solver *ASolver, lorads_lp_dense *d, lorads_lp_dense *gradLP, lorads_sdp_dense **D, lorads_sdp_dense **Grad);
extern void ALMCalq12p12LP(lorads_solver *ASolver, lorads_lp_dense *rLp, lorads_lp_dense *d, lorads_sdp_dense **R, lorads_sdp_dense **D, double *q1, double *q2, double *p12);
extern void SetyAsNegGrad(lorads_solver *ASolver, lorads_lp_dense *gradLpDummy, lorads_sdp_dense **Grad);
extern void SetyAsNegGradLP(lorads_solver *ASolver, lorads_lp_dense *gradLp, lorads_sdp_dense **Grad);
extern void setlbfgsHisTwo(lorads_solver *ASolver, lorads_lp_dense *gradLpDummy, lorads_lp_dense *dDummy, lorads_sdp_dense **Grad, lorads_sdp_dense **D, double tau);
extern void setlbfgsHisTwoLP(lorads_solver *ASolver, lorads_lp_dense *gradLp, lorads_lp_dense *d, lorads_sdp_dense **Grad, lorads_sdp_dense **D, double tau);
extern void ALMupdateVar(lorads_solver *ASolver, lorads_lp_dense *rLpDummy, lorads_lp_dense *dDummy, lorads_sdp_dense **R, lorads_sdp_dense **D, double tau);
extern void ALMupdateVarLP(lorads_solver *ASolver, lorads_lp_dense *rLp, lorads_lp_dense *d, lorads_sdp_dense **R, lorads_sdp_dense **D, double tau);
extern void ALMCalq12p12(lorads_solver *ASolver, lorads_lp_dense *rLpDummy, lorads_lp_dense *dDummy, lorads_sdp_dense **R, lorads_sdp_dense **D, double *q1, double *q2, double *p12);
extern void ALMCalq12p12LP(lorads_solver *ASolver, lorads_lp_dense *rLp, lorads_lp_dense *d, lorads_sdp_dense **R, lorads_sdp_dense **D, double *q1, double *q2, double *p12);
extern void LORADSCalObjRR_ALM(lorads_solver *ASolver);
extern void LORADSCalObjRR_ALM_LP(lorads_solver *ASolver);
extern lorads_int LORADS_ALMOptimize_reopt(lorads_params *params, lorads_solver *ASolver, lorads_alm_state *alm_iter_state, bool early_stop, double rho_update_factor, double timeSolveStart);
extern lorads_int LORADS_ALMOptimize(lorads_params *params, lorads_solver *ASolver, lorads_alm_state *alm_iter_state, double rho_update_factor, double timeSolveStart);
#endif
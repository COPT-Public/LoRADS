#ifndef LORADS_SOLVER
#define LORADS_SOLVER

#include "def_lorads_solver.h"
#include "lorads.h"
#include "lorads_user_data.h"
extern void LORADSInitSolver(lorads_solver *ASolver, lorads_int nRows, lorads_int nCones, lorads_int *blkDims, lorads_int nLpCols);
extern void LORADSDestroySolver(lorads_solver *ASolver);
extern void LORADSSetDualObjective(lorads_solver *ASolver, double *dObj);
extern void LORADSInitConeData(lorads_solver *ASolver, user_data **SDPDatas,
                               double **coneMatElem, lorads_int **coneMatBeg, lorads_int **coneMatIdx,
                               lorads_int *BlkDims, lorads_int nConstrs, lorads_int nBlks,
                               lorads_int nLpCols, lorads_int *LpMatBeg, lorads_int *LpMatIdx, double *LpMatElem);
extern void LORADSPreprocess(lorads_solver *ASolver, lorads_int *BlkDims);
extern void LORADSDestroyConeData(lorads_solver *ASolver);
extern void destroyPreprocess(lorads_solver *ASolver);
extern void LORADSDetermineRank(lorads_solver *ASolver, lorads_int *blkDims, double timesRank);
extern void detectSparsitySDPCoeff(lorads_solver *ASolver);
extern void LORADSInitALMVars(lorads_solver *ASolver, lorads_int *rankElem, lorads_int *BlkDims, lorads_int nBlks, lorads_int nLpCols, lorads_int lbfgsHis);
extern void LORADSDestroyALMVars(lorads_solver *ASolver);
extern void LORADSInitADMMVars(lorads_solver *ASolver, lorads_int *rankElem, lorads_int *BlkDims, lorads_int nBlks, lorads_int nLpCols);
extern void LORADSDestroyADMMVars(lorads_solver *ASolver);
extern void LORADSInitFuncSet(lorads_func **pfunc, lorads_int nLpCols);
extern lorads_int CheckAllRankMax(lorads_solver *asolver, double aug_factor);
extern lorads_int AUG_RANK(lorads_solver *ASolver, lorads_int *BlkDims, lorads_int nBlks, double aug_factor);
extern void  LORADSEndProgram( lorads_solver *ASolver);
extern void printRes(double pObj, double dObj, double constrVio, double dualInfe, double pdgap, double constrVioInf, double dualInfeInf);
extern void LORADS_ALMtoADMM(lorads_solver *ASolver, lorads_params *params, lorads_alm_state *alm_state, lorads_admm_state *admm_state);
extern void calculate_dual_infeasibility_solver(lorads_solver *ASolver);
extern void objScale_dualvar(lorads_solver *ASolver, double *scaleTemp, double *scaleHis);
extern void cal_sdp_const(lorads_params *params, lorads_solver *ASolver, SDPConst *sdpConst);
extern void initial_solver_state(lorads_params *params, lorads_solver *ASolver, lorads_alm_state *alm_state_pointer, lorads_admm_state *admm_state_pointer, SDPConst *sdpConst);
extern double reopt(lorads_params *params, lorads_solver *ASolver, lorads_alm_state *alm_state_pointer, lorads_admm_state *admm_state_pointer, double *reopt_param, lorads_int *reopt_alm_iter, lorads_int *reopt_admm_iter, double timeSolveStart, int *admm_bad_iter_flag, int reopt_level);
extern void printfProbInfo(lorads_solver *ASolver);
#endif
#ifndef DEF_LORADS_SOLVER
#define DEF_LORADS_SOLVER

#include <stdbool.h>
#include "def_lorads_lp_conic.h"
#include "def_lorads_sdp_conic.h"
#include "def_lorads_elements.h"
#include "def_lorads_cgs.h"
#include "def_lorads_lbfgs.h"
#include "lorads.h"

typedef struct{
    /* Variables SDPCone */
    lorads_sdp_dense **U;    // admm variable, and lbfgs descent direction D
    lorads_sdp_dense **V;    // admm variable only
#ifdef DUAL_U_V
    lorads_sdp_dense **S;    // admm dual variable
#endif
    lorads_sdp_dense **R;    // average variable for storage and ALM variable
    lorads_sdp_dense **Grad; // grad of R

    /* Variable LPCone */
    lorads_lp_dense *rLp;
    lorads_lp_dense *uLp;
    lorads_lp_dense *vLp;
#ifdef DUAL_U_V
    lorads_lp_dense *sLp;
#endif
    lorads_lp_dense *gradLp;

    /* ALM lbfgs and ADMM Variables */
    double *dualVar;

    /* Auxiliary variable */
    lorads_vec **constrVal; // constraint violation [iCone][iConstr]
    double *constrValSum;        // constraint violation [iConstr]
    double *ARDSum;              // q1 in line search of ALM
    double *ADDSum;              // q2 in line search of ALM
    double **bLinSys;            // for solving linear system, bInLinSys[iCone]
    lorads_int *rankElem;               // all rank
    double *M1temp;              // M1 for solving linear system
    double *bestDualVar;
    lorads_sdp_dense **M2temp; // M2 for solving linear system
    double *Dtemp;
    lorads_vec **constrValLP;
}lorads_variable;

typedef struct{
    /* User data */
    lorads_int nRows;      // constraint number
    double *rowRHS; // b of Ax = b

    /* Cones */
    lorads_int nCones; // sdp cones block number
    lorads_sdp_cone **SDPCones;
    lorads_cg_linsys **CGLinsys;

    // variable
    lorads_variable *var;

    /* Auxiliary variable LPCone (SDPCone but rank is 1, dim is 1) */
    lorads_lp_cone *lpCone;
    lorads_int nLpCols;

    /* ALM lbfgs Variables */
    lorads_int hisRecT;
    lbfgs_node *lbfgsHis; // record difference of primal variable history and gradient history, lbfgsHis[iCone]

    /* Monitor */
    lorads_int nIterCount;
    double cgTime;
    lorads_int cgIter;
    lorads_int checkSolTimes;
//    double traceSum;

    /* Convergence criterion */
    double pObjVal;
    double dObjVal;
    double pInfeas;
    double dInfeas;
    double cObjNrm1;
    double cObjNrm2;
    double cObjNrmInf;
    double bRHSNrm1;
    double bRHSNrmInf;
    double bRHSNrm2;
    double *dimacError;
    double *constrVio;

    /* Starting time */
    double dTimeBegin;

    // scale
    double cScaleFactor;
    double bScaleFactor;

    // check exit alm
    lorads_int *rank_max;
    double *sparsitySDPCoeff;
    double overallSparse;
    lorads_int nnzSDPCoeffSum;
    lorads_int SDPCoeffSum;
    lorads_status AStatus;

    double scaleObjHis;
} lorads_solver;


typedef struct {
    void (*InitConstrValAll)    (lorads_solver *, lorads_lp_dense *, lorads_lp_dense *, lorads_sdp_dense **, lorads_sdp_dense **);
    void (*InitConstrValSum)    (lorads_solver *);
    void (*ALMCalGrad)          (lorads_solver *, lorads_lp_dense *, lorads_lp_dense *, lorads_sdp_dense **, lorads_sdp_dense **, double *, double);
    void (*LBFGSDirection)      (lorads_params *, lorads_solver *, lbfgs_node *, lorads_lp_dense *, lorads_lp_dense *, lorads_sdp_dense **, lorads_sdp_dense **, lorads_int);
    void (*LBFGSDirUseGrad)     (lorads_solver *, lorads_lp_dense *, lorads_lp_dense *, lorads_sdp_dense **, lorads_sdp_dense **);
    void (*copyRtoV)            (lorads_lp_dense *rlp, lorads_lp_dense *vlp, lorads_sdp_dense **, lorads_sdp_dense **, lorads_int);
    void (*ALMCalq12p12)        (lorads_solver *, lorads_lp_dense *, lorads_lp_dense *, lorads_sdp_dense **, lorads_sdp_dense **, double *, double *, double *);
    void (*setAsNegGrad)        (lorads_solver *, lorads_lp_dense *, lorads_sdp_dense **);
    void (*ALMupdateVar)        (lorads_solver *, lorads_lp_dense *, lorads_lp_dense *, lorads_sdp_dense **, lorads_sdp_dense **, double);
    void (*setlbfgsHisTwo)      (lorads_solver *, lorads_lp_dense *, lorads_lp_dense *, lorads_sdp_dense **, lorads_sdp_dense **, double);
    void (*updateDimacsALM)     (lorads_solver *, lorads_sdp_dense **, lorads_sdp_dense **, lorads_lp_dense *, lorads_lp_dense *);
    void (*updateDimacsADMM)    (lorads_solver *, lorads_sdp_dense **, lorads_sdp_dense **, lorads_lp_dense *, lorads_lp_dense *);
    void (*calObj_admm)         (lorads_solver *);
    void (*calObj_alm)          (lorads_solver *);
    void (*updateDimac)         (lorads_solver *);

    void (*admmUpdateVar)       (lorads_solver *, double, double, lorads_int);
}lorads_func;


typedef struct{
    bool is_rank_updated;
    lorads_int outerIter;
    lorads_int innerIter;
    double rho;
    double l_inf_primal_infeasibility;
    double l_1_primal_infeasibility;
    double l_2_primal_infeasibility;
    double primal_dual_gap;
    double primal_objective_value;
    double dual_objective_value;
    double l_inf_dual_infeasibility;
    double l_1_dual_infeasibility;
    double l_2_dual_infeasibility;
    double tau;
}lorads_alm_state;

typedef struct{
    lorads_int iter;
    lorads_int nBlks;
    lorads_int cg_iter;
    double rho;
    double l_1_dual_infeasibility;
    double l_inf_dual_infeasibility;
    double l_1_primal_infeasibility;
    double l_inf_primal_infeasibility;
    double l_2_primal_infeasibility;
    double l_2_dual_infeasibility;
    double primal_objective_value;
    double dual_objective_value;
    double primal_dual_gap;
}lorads_admm_state;


typedef struct{
    double l_1_norm_c;
    double l_2_norm_c;
    double l_inf_norm_c;
    double l_1_norm_b;
    double l_2_norm_b;
    double l_inf_norm_b;
}SDPConst;

#endif
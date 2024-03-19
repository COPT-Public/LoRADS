#ifndef DEF_ASDP_FUNC_SET_H
#define DEF_ASDP_FUNC_SET_H


#include "asdp_lbfgs.h"
#include "def_asdp.h"
#include "def_asdp_rk_mat.h"


typedef struct {
    void (*InitConstrValAll)    (asdp *, asdp_rk_mat_lp *, asdp_rk_mat_lp *, asdp_rk_mat_dense **, asdp_rk_mat_dense **);
    void (*InitConstrValSum)    (asdp *);
    void (*BMCalGrad)           (asdp *, asdp_rk_mat_lp *, asdp_rk_mat_lp *, asdp_rk_mat_dense **, asdp_rk_mat_dense **, double *);
    void (*LBFGSDirection)      (asdp *, lbfgs_node *, asdp_rk_mat_lp *, asdp_rk_mat_lp *, asdp_rk_mat_dense **, asdp_rk_mat_dense **, int);
    void (*LBFGSDirUseGrad)     (asdp *, asdp_rk_mat_lp *, asdp_rk_mat_lp *, asdp_rk_mat_dense **, asdp_rk_mat_dense **);
    void (*copyRtoV)            (asdp_rk_mat_lp *rlp, asdp_rk_mat_lp *vlp, asdp_rk_mat_dense **, asdp_rk_mat_dense **, int);
    void (*BMCalq12p12)         (asdp *, asdp_rk_mat_lp *, asdp_rk_mat_lp *, asdp_rk_mat_dense **, asdp_rk_mat_dense **, double *, double *, double *);
    void (*setAsNegGrad)        (asdp *, asdp_rk_mat_lp *, asdp_rk_mat_dense **);
    void (*BMupdateVar)         (asdp *, asdp_rk_mat_lp *, asdp_rk_mat_lp *, asdp_rk_mat_dense **, asdp_rk_mat_dense **, double);
    void (*setlbfgsHisTwo)      (asdp *, asdp_rk_mat_lp *, asdp_rk_mat_lp *, asdp_rk_mat_dense **, asdp_rk_mat_dense **, double);
    void (*calObj)              (asdp *, int);
    asdp_retcode (*updateDimac) (asdp *);
//    void (*dualGradStand)       (asdp *, double *, double *, asdp_rk_mat_dense **, asdp_rk_mat_dense **);
//    void (*dualGradNcvx)        (asdp *, double *, double *, asdp_rk_mat_dense **, asdp_rk_mat_dense **);
    
    void (*admmUpdateVar)       (asdp *);
}asdp_func;

#ifdef __cplusplus
extern "C" {
#endif
extern void ASDPInitFuncSet(asdp_func **pfunc, int nLpCols);
#ifdef __cplusplus
}
#endif
#endif

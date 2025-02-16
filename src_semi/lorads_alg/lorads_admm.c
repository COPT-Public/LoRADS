#include <math.h>
#include "lorads.h"
#include "lorads_solver.h"
#include "lorads_alg_common.h"
#include "lorads_vec_opts.h"
#include "lorads_cgs.h"

extern void ADMMPrintLog(lorads_admm_state *admm_iter_state, double time){
#ifdef INT32
    printf("ADMM Iter:%d pObj:%5.5e dObj:%5.5e pInfea(1):%5.5e pInfea(Inf):%5.5e pdGap:%5.5e rho:%3.2f cgIter:%d Time:%3.2f\n",
            admm_iter_state->iter, admm_iter_state->primal_objective_value,
            admm_iter_state->dual_objective_value, admm_iter_state->l_1_primal_infeasibility,
            admm_iter_state->l_inf_primal_infeasibility, admm_iter_state->primal_dual_gap,
            admm_iter_state->rho, (int) ((double)admm_iter_state->cg_iter / (double)admm_iter_state->nBlks), time);
#endif
#ifdef MAC_INT64
    printf("ADMM Iter:%lld pObj:%5.5e dObj:%5.5e pInfea(1):%5.5e pInfea(Inf):%5.5e pdGap:%5.5e rho:%3.2f cgIter:%d Time:%3.2f\n",
           admm_iter_state->iter, admm_iter_state->primal_objective_value,
           admm_iter_state->dual_objective_value, admm_iter_state->l_1_primal_infeasibility,
           admm_iter_state->l_inf_primal_infeasibility, admm_iter_state->primal_dual_gap,
           admm_iter_state->rho, (int) ((double)admm_iter_state->cg_iter / (double)admm_iter_state->nBlks), time);
#endif
#ifdef UNIX_INT64
    printf("ADMM Iter:%ld pObj:%5.5e dObj:%5.5e pInfea(1):%5.5e pInfea(Inf):%5.5e pdGap:%5.5e rho:%3.2f cgIter:%d Time:%3.2f\n",
           admm_iter_state->iter, admm_iter_state->primal_objective_value,
           admm_iter_state->dual_objective_value, admm_iter_state->l_1_primal_infeasibility,
           admm_iter_state->l_inf_primal_infeasibility, admm_iter_state->primal_dual_gap,
           admm_iter_state->rho, (int) ((double)admm_iter_state->cg_iter / (double)admm_iter_state->nBlks), time);
#endif
}


extern lorads_int LORADSADMMOptimize(lorads_params *params, lorads_solver *ASolver, lorads_admm_state *admm_iter_state,
                                           lorads_int iter_celling, double timeSolveStart){
    if (admm_iter_state->primal_dual_gap <= params->phase2Tol && admm_iter_state->l_1_primal_infeasibility <= params->phase2Tol){
        return RET_CODE_OK;
    }
    double endCG_tol = LORADS_MIN(admm_iter_state->l_1_primal_infeasibility * 1e-4, 1e-8);
    lorads_int CGMaxIter = 800;
    admm_iter_state->rho = LORADS_MIN(admm_iter_state->rho, params->rhoMax);
    lorads_func *aFunc;
    LORADSInitFuncSet(&aFunc, ASolver->nLpCols);
    ASolver->cgIter = 0;
    lorads_int retcode = RET_CODE_OK;
    // Initialize info
    // Initialize Constraint function value <A, UV^T>
    LORADSInitConstrValAll(ASolver, ASolver->var->uLp, ASolver->var->vLp, ASolver->var->U, ASolver->var->V);
    LORADSInitConstrValSum(ASolver);
    // Initialize Objective function value
    aFunc->calObj_admm(ASolver);
    LORADSCalDualObj(ASolver);
    aFunc->updateDimacsADMM(ASolver, ASolver->var->U, ASolver->var->V, ASolver->var->uLp, ASolver->var->vLp);
    admm_iter_state->primal_objective_value = ASolver->pObjVal;
    admm_iter_state->dual_objective_value = ASolver->dObjVal;
    admm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
    admm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
    admm_iter_state->l_inf_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf);
    admm_iter_state->l_2_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrm2);
    ASolver->cgTime = 0.0;

    double cur_rho_max = params->rhoMax;

    double old_l_inf_primal_infeasibility_mean = 1e30;

    double primal_infeasibility_buffer[10];
    LORADS_ZERO(primal_infeasibility_buffer, double, 10);
    lorads_int buffer_size = 10;
    int bad_pd = 0;
    int count = 0;
    double origTime = LUtilGetTimeStamp();
    while (admm_iter_state->iter <= params->maxADMMIter || admm_iter_state->primal_dual_gap >= params->phase2Tol ||
           admm_iter_state->l_1_primal_infeasibility >= params->phase2Tol){
        if (admm_iter_state->iter >= iter_celling){
            break;
        }
        endCG_tol = LORADS_MIN(admm_iter_state->l_1_primal_infeasibility * 1e-2, 1e-8);
        aFunc->admmUpdateVar(ASolver, admm_iter_state->rho, endCG_tol, CGMaxIter);
        admm_iter_state->cg_iter = ASolver->cgIter;
        aFunc->calObj_admm(ASolver);
        LORADSCalDualObj(ASolver);
        aFunc->updateDimacsADMM(ASolver, ASolver->var->U, ASolver->var->V, ASolver->var->uLp, ASolver->var->vLp);
        admm_iter_state->primal_objective_value = ASolver->pObjVal;
        admm_iter_state->dual_objective_value = ASolver->dObjVal;
        admm_iter_state->l_inf_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf);
        admm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
        if (admm_iter_state->l_inf_primal_infeasibility >= 1e10 || admm_iter_state->primal_dual_gap >= 1-1e-8){
            printf("Numerical Error!");
            return RET_CODE_NUM_ERR;
        }
        if (admm_iter_state->primal_dual_gap <= params->phase2Tol * 5)
        {
            bad_pd -= 5;
            bad_pd = LORADS_MAX(0, bad_pd);
        }
        else if (admm_iter_state->primal_dual_gap <= params->phase2Tol)
        {
            bad_pd -= 10;
            bad_pd = LORADS_MAX(0, bad_pd);
        }

        if (admm_iter_state->primal_dual_gap >= params->phase1Tol * 1e2)
        {
            bad_pd += 2;
        }

        if (bad_pd >= 800){
            return RET_CODE_OK;
        }

        primal_infeasibility_buffer[count % 10] = admm_iter_state->l_inf_primal_infeasibility;
        if (admm_iter_state->l_inf_primal_infeasibility <= params->phase2Tol){
            aFunc->updateDimacsADMM(ASolver, ASolver->var->U, ASolver->var->V, ASolver->var->uLp, ASolver->var->vLp);
            admm_iter_state->primal_objective_value = ASolver->pObjVal;
            admm_iter_state->dual_objective_value = ASolver->dObjVal;
            admm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
            admm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
            ADMMPrintLog(admm_iter_state, LUtilGetTimeStamp() - origTime);
            return RET_CODE_OK;
        }
        LORADSUpdateDualVar(ASolver, admm_iter_state->rho);
        if ((admm_iter_state->iter + 1) % params->rhoFreq == 0){
            admm_iter_state->rho *= params->rhoFactor;
            if (admm_iter_state->rho >= cur_rho_max){
                admm_iter_state->rho = cur_rho_max;
                if ((admm_iter_state->iter + 1) % (params->rhoFreq * 100) == 0){
                    double l_inf_primal_infeasibility_mean =
                            nrm1(&buffer_size, primal_infeasibility_buffer, &AIntConstantOne) / (double) buffer_size;
                    if (l_inf_primal_infeasibility_mean / old_l_inf_primal_infeasibility_mean >= 0.65){
                        admm_iter_state->rho *= pow(params->rhoFactor, round(log(params->rhoFreq * 100) / log(params->rhoFreq)));
                        cur_rho_max = admm_iter_state->rho;
                    }
                    old_l_inf_primal_infeasibility_mean = l_inf_primal_infeasibility_mean;
                }
            }
            if (admm_iter_state->rho >= params->rhoCellingADMM){
                admm_iter_state->rho = params->rhoCellingADMM;
            }
        }
        if (admm_iter_state->iter % 50 == 0){
            aFunc->updateDimacsADMM(ASolver, ASolver->var->U, ASolver->var->V, ASolver->var->uLp, ASolver->var->vLp);
            admm_iter_state->primal_objective_value = ASolver->pObjVal;
            admm_iter_state->dual_objective_value = ASolver->dObjVal;
            admm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
            admm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
            ADMMPrintLog(admm_iter_state, LUtilGetTimeStamp() - origTime);
            if (LUtilGetTimeStamp() - timeSolveStart >= params->timeSecLimit){
                return RET_CODE_TIME_OUT;
            }
        }
        if (admm_iter_state->primal_dual_gap <= params->phase2Tol*1e-3 && admm_iter_state->l_1_primal_infeasibility <= params->phase2Tol*1e-3){
            printf("Early Stop When DIMACS Errors Are Well-Satisfied");
            return RET_CODE_OK;
        }
        admm_iter_state->iter++;
    }
    return RET_CODE_OK;
}


extern lorads_int LORADSADMMOptimize_reopt(lorads_params *params, lorads_solver *ASolver, lorads_admm_state *admm_iter_state,
                                     lorads_int iter_celling, double timeSolveStart){
    if (admm_iter_state->primal_dual_gap <= params->phase2Tol && admm_iter_state->l_1_primal_infeasibility <= params->phase2Tol){
        return RET_CODE_OK;
    }
    double endCG_tol = LORADS_MIN(admm_iter_state->l_1_primal_infeasibility * 1e-4, 1e-8);
    lorads_int CGMaxIter = 800;
    admm_iter_state->rho = LORADS_MIN(admm_iter_state->rho, params->rhoMax);
    lorads_func *aFunc;
    LORADSInitFuncSet(&aFunc, ASolver->nLpCols);
    ASolver->cgIter = 0;
    lorads_int retcode = RET_CODE_OK;
    // Initialize info
    // Initialize Constraint function value <A, UV^T>
    LORADSInitConstrValAll(ASolver, ASolver->var->uLp, ASolver->var->vLp, ASolver->var->U, ASolver->var->V);
    LORADSInitConstrValSum(ASolver);
    // Initialize Objective function value
    aFunc->calObj_admm(ASolver);
    LORADSCalDualObj(ASolver);
    aFunc->updateDimacsADMM(ASolver, ASolver->var->U, ASolver->var->V, ASolver->var->uLp, ASolver->var->vLp);
    admm_iter_state->primal_objective_value = ASolver->pObjVal;
    admm_iter_state->dual_objective_value = ASolver->dObjVal;
    admm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
    admm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
    admm_iter_state->l_inf_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf);
    admm_iter_state->l_2_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrm2);
    printf("enter admm reopt \n");
    ADMMPrintLog(admm_iter_state, 0);


    ASolver->cgTime = 0.0;

    double cur_rho_max = params->rhoMax;

    double old_l_inf_primal_infeasibility_mean = 1e30;

    double primal_infeasibility_buffer[10];
    LORADS_ZERO(primal_infeasibility_buffer, double, 10);
    lorads_int buffer_size = 10;
    int count = 0;
    int bad_pd = 0;
    double origTime = LUtilGetTimeStamp();
    while (admm_iter_state->iter <= params->maxADMMIter || admm_iter_state->primal_dual_gap >= params->phase2Tol ||
            admm_iter_state->l_1_primal_infeasibility >= params->phase2Tol){
        if (admm_iter_state->iter >= iter_celling){
            ADMMPrintLog(admm_iter_state, 0);
#ifdef INI32
            printf("exit admm since maxiter greater than iter_celling:%d\n", iter_celling);
#endif
#ifdef MAC_INT64
            printf("exit admm since maxiter greater than iter_celling:%lld\n", iter_celling);
#endif
#ifdef UNIX_INT64
            printf("exit admm since maxiter greater than iter_celling:%ld\n", iter_celling);
#endif
            
            break;
        }
        endCG_tol = LORADS_MIN(admm_iter_state->l_1_primal_infeasibility * 1e-4, 1e-8);
        aFunc->admmUpdateVar(ASolver, admm_iter_state->rho, endCG_tol, CGMaxIter);
        admm_iter_state->cg_iter = ASolver->cgIter;
        aFunc->calObj_admm(ASolver);
        LORADSCalDualObj(ASolver);
        aFunc->updateDimacsADMM(ASolver, ASolver->var->U, ASolver->var->V, ASolver->var->uLp, ASolver->var->vLp);
        admm_iter_state->primal_objective_value = ASolver->pObjVal;
        admm_iter_state->dual_objective_value = ASolver->dObjVal;
        admm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
        admm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
        admm_iter_state->l_inf_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrmInf);
        admm_iter_state->l_2_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1] * (1 + ASolver->bRHSNrm1) / (1 + ASolver->bRHSNrm2);
        if (admm_iter_state->l_inf_primal_infeasibility >= 1e10 || admm_iter_state->primal_dual_gap >= 1-1e-8){
            ADMMPrintLog(admm_iter_state, LUtilGetTimeStamp() - origTime);
            printf("Numerical Error!\n");
            return RET_CODE_NUM_ERR;
        }

        if (admm_iter_state->primal_dual_gap <= params->phase2Tol * 5)
        {
            bad_pd -= 5;
            bad_pd = LORADS_MAX(0, bad_pd);
        }
        else if (admm_iter_state->primal_dual_gap <= params->phase2Tol)
        {
            bad_pd -= 10;
            bad_pd = LORADS_MAX(0, bad_pd);
        }

        if (admm_iter_state->primal_dual_gap >= params->phase1Tol * 1e2)
        {
            bad_pd += 2;
        }

        if (bad_pd >= 200){
            printf("------");
            return RET_CODE_OK;
        }

        primal_infeasibility_buffer[count % 10] = admm_iter_state->l_inf_primal_infeasibility;
        if (admm_iter_state->l_1_primal_infeasibility <= params->phase2Tol){
            aFunc->updateDimacsADMM(ASolver, ASolver->var->U, ASolver->var->V, ASolver->var->uLp, ASolver->var->vLp);
            admm_iter_state->primal_objective_value = ASolver->pObjVal;
            admm_iter_state->dual_objective_value = ASolver->dObjVal;
            admm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
            admm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
            if (admm_iter_state->primal_dual_gap <= params->phase2Tol){
                ADMMPrintLog(admm_iter_state, LUtilGetTimeStamp() - origTime);
                return RET_CODE_OK;
            }
        }
        LORADSUpdateDualVar(ASolver, admm_iter_state->rho);
        if (admm_iter_state->iter % params->rhoFreq == 0){
            admm_iter_state->rho *= params->rhoFactor;
            if (admm_iter_state->rho >= cur_rho_max){
                admm_iter_state->rho = cur_rho_max;
                if (admm_iter_state->iter % (params->rhoFreq * 100) == 0){
                    double l_inf_primal_infeasibility_mean =
                            nrm1(&buffer_size, primal_infeasibility_buffer, &AIntConstantOne) / buffer_size;
                    if (l_inf_primal_infeasibility_mean / old_l_inf_primal_infeasibility_mean >= 0.65){
                        admm_iter_state->rho *= pow(params->rhoFactor, round(log(params->rhoFreq * 100) / log(params->rhoFreq)));
                        cur_rho_max = admm_iter_state->rho;
                    }
                    old_l_inf_primal_infeasibility_mean = l_inf_primal_infeasibility_mean;
                }
            }
            if (admm_iter_state->rho >= params->rhoCellingADMM){
                admm_iter_state->rho = params->rhoCellingADMM;
            }
        }
        if (admm_iter_state->iter % 50 == 0){
            aFunc->updateDimacsADMM(ASolver, ASolver->var->U, ASolver->var->V, ASolver->var->uLp, ASolver->var->vLp);
            admm_iter_state->primal_objective_value = ASolver->pObjVal;
            admm_iter_state->dual_objective_value = ASolver->dObjVal;
            admm_iter_state->primal_dual_gap = ASolver->dimacError[LORADS_DIMAC_ERROR_PDGAP];
            admm_iter_state->l_1_primal_infeasibility = ASolver->dimacError[LORADS_DIMAC_ERROR_CONSTRVIO_L1];
            ADMMPrintLog(admm_iter_state, LUtilGetTimeStamp() - origTime);
            if (LUtilGetTimeStamp() - timeSolveStart >= params->timeSecLimit){
                return RET_CODE_TIME_OUT;
            }
        }
        if (admm_iter_state->primal_dual_gap <= params->phase2Tol*1e-3 && admm_iter_state->l_1_primal_infeasibility <= params->phase2Tol*1e-3){
            printf("Early Stop When DIMACS Errors Are Well-Satisfied");
            return RET_CODE_OK;
        }
        admm_iter_state->iter++;
    }
    ADMMPrintLog(admm_iter_state, LUtilGetTimeStamp() - origTime);
    return RET_CODE_OK;
}


extern void averageUV(lorads_sdp_dense *U, lorads_sdp_dense *V, lorads_sdp_dense *UVavg){
    lorads_int n = U->rank * U->nRows;
    for (lorads_int i = 0; i < n; ++i){
        UVavg->matElem[i] = (U->matElem[i] + V->matElem[i]) / 2;
    }
}

extern void averageUVLP(lorads_lp_dense *ulp, lorads_lp_dense *vlp, lorads_lp_dense *uvavg){
    lorads_int n = ulp->nCols;
    for (lorads_int i = 0; i < n; ++i){
        uvavg->matElem[i] = (ulp->matElem[i] + vlp->matElem[i]) / 2;
    }
}


extern void LORADSCalObjUV_ADMM(lorads_solver *ASolver){
    ASolver->pObjVal = 0.0;
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        lorads_sdp_dense *R = ASolver->var->R[iCone];
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        lorads_sdp_dense *U = ASolver->var->U[iCone];
        lorads_sdp_dense *V = ASolver->var->V[iCone];
        averageUV(U, V, R);
        LORADSUVt(ACone->sdp_obj_sum, R, R);
        ACone->objAUV(ACone->coneData, R, R, &ASolver->pObjVal, ACone->sdp_obj_sum);
    }
    ASolver->pObjVal /= ASolver->scaleObjHis;
}

extern void LORADSCalObjUV_ADMM_LP(lorads_solver *ASolver){
    ASolver->pObjVal = 0.0;
    averageUVLP(ASolver->var->uLp, ASolver->var->vLp, ASolver->var->rLp);
    ASolver->lpCone->objAUV(ASolver->lpCone->coneData, ASolver->var->rLp, ASolver->var->rLp, &ASolver->pObjVal);
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone){
        lorads_sdp_dense *R = ASolver->var->R[iCone];
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        lorads_sdp_dense *U = ASolver->var->U[iCone];
        lorads_sdp_dense *V = ASolver->var->V[iCone];
        averageUV(U, V, R);
        LORADSUVt(ACone->sdp_obj_sum, R, R);
        ACone->objAUV(ACone->coneData, R, R, &ASolver->pObjVal, ACone->sdp_obj_sum);
//        printf("iCone:%d, %.20f\n", iCone, ASolver->pObjVal);
    }
    ASolver->pObjVal /= ASolver->scaleObjHis;
}

extern void LORADSUpdateConstrValCG(lorads_sdp_cone *ACone, lorads_sdp_dense *U, lorads_sdp_dense *V, double *weight)
{
    LORADSUVt(ACone->sdp_coeff_w_sum, U, V);
    ACone->coneAUV(ACone->coneData, U, V, weight, ACone->sdp_coeff_w_sum);
    if (ACone->type == LORADS_CONETYPE_SPARSE_SDP){
        // put result into right position
        double temp = 0.0;
        int idx = 0;
        lorads_cone_sdp_sparse *sparse = (lorads_cone_sdp_sparse *)ACone->coneData;
        for (lorads_int i = sparse->nRowElem - 1; i >= 0; i--)
        {
            temp = weight[i];
            weight[i] = 0.0;
            idx = sparse->rowIdx[i];
            weight[idx] = temp;
        }
    }
}


static void linSysProduct(lorads_sdp_cone *ACone, double *weight, lorads_sdp_dense *noUpdateVar, lorads_sdp_dense *updateVar, double *x, double *res)
{
    updateVar->matElem = x;
    updateVar->rank = noUpdateVar->rank;
    updateVar->nRows = noUpdateVar->nRows;
    LORADSUpdateConstrValCG(ACone, updateVar, noUpdateVar, weight);
    sdp_coeff *w_sum = ACone->sdp_coeff_w_sum;
    w_sum->zeros(w_sum->dataMat);
    ACone->sdpDataWSum(ACone->coneData, weight, w_sum);

    w_sum->mul_rk(w_sum->dataMat, noUpdateVar, res);
    lorads_int n = noUpdateVar->nRows * noUpdateVar->rank;
    double alpha = 1.0;
    lorads_int incx = 1;
    axpy(&n, &alpha, x, &incx, res, &incx);
}

extern void LORADSRkMatSub(lorads_sdp_dense *A,  lorads_sdp_dense *B, double alpha){
    // A = A - alpha * B
    double negAlpha = -1 * alpha;
    lorads_int n = A->nRows * A->rank;
    lorads_int incx = 1;
    axpy(&n, &negAlpha, B->matElem, &incx, A->matElem, &incx);
}

extern void LORADSRkMatSub_positive_S(lorads_sdp_dense *A,  lorads_sdp_dense *B, double alpha, lorads_sdp_dense *S){
    // A = A - alpha * B + S
    double negAlpha = -1 * alpha;
    lorads_int n = A->nRows * A->rank;
    lorads_int incx = 1;
    double one = 1.0;
    axpy(&n, &negAlpha, B->matElem, &incx, A->matElem, &incx);
    axpy(&n, &one, S->matElem, &incx, A->matElem, &incx);
}

extern void LORADSRkMatSub_negative_S(lorads_sdp_dense *A,  lorads_sdp_dense *B, double alpha, lorads_sdp_dense *S){
    // A = A - alpha * B - S
    double negAlpha = -1 * alpha;
    lorads_int n = A->nRows * A->rank;
    lorads_int incx = 1;
    double negative_one = -1.0;
    axpy(&n, &negAlpha, B->matElem, &incx, A->matElem, &incx);
    axpy(&n, &negative_one, S->matElem, &incx, A->matElem, &incx);
}

extern void ADMMUpdateUVMvec(void *pData, double *x, double *res)
{
    admmCG *data = (admmCG *)pData;
    data->UpdateVarShell->matElem = x;
    linSysProduct(data->ACone, data->weight, data->noUpdateVar, data->UpdateVarShell, x, res);
}

extern void LORADSUpdateSDPVarOne(lorads_solver *ASolver, lorads_sdp_dense *updateVar, lorads_sdp_dense *noUpdateVar, lorads_int iCone, double rho, double CG_tol, lorads_int CG_maxIter){
    lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
    double *b = ASolver->rowRHS;
    double *M1 = ASolver->var->M1temp;
    LORADS_MEMCPY(M1, b, double, ASolver->nRows);
    lorads_sdp_dense *M2 = ASolver->var->M2temp[iCone];
    double One = 1.0;
    double minusOne = -1.0;
    lorads_int incx = 1;
    // M1 = -b
    scal(&(ASolver->nRows), &(minusOne), M1, &incx);
    // M1 = sum_{j\neq \tilde{j}}A(X) - b
    axpy(&(ASolver->nRows), &(One), ASolver->var->constrValSum, &incx, M1, &incx);
    // axpy(&(ASolver->nRows), &(minusOne), ASolver->constrVal[iCone], &incx, M1, &incx);
    ASolver->var->constrVal[iCone]->add(&minusOne, ASolver->var->constrVal[iCone]->data, M1);
    // M1 = M1 * rho
    scal(&(ASolver->nRows), &(rho), M1, &incx);
    axpy(&(ASolver->nRows), &(minusOne), ASolver->var->dualVar, &incx, M1, &incx);

    // Add C in the following function
    // and weighted sum of all sdp data in one cone
    ACone->sdp_obj_sum->zeros(ACone->sdp_obj_sum->dataMat);
    ACone->addObjCoeff(ACone->coneData, ACone->sdp_obj_sum);
    ACone->sdpDataWSum(ACone->coneData, M1, ACone->sdp_obj_sum);

    // (C-(lambda + weighted Sum A)) * V
    ACone->sdp_obj_sum->mul_rk(ACone->sdp_obj_sum->dataMat, noUpdateVar, M2->matElem);

    // M2 = (M2 -rho * V + S) or (M2 -rho * V - S)
    LORADSRkMatSub(M2, noUpdateVar, rho);

    // set bv in linear system
    double scalFactor = -1.0 / rho;
    lorads_int n = M2->nRows * M2->rank;
    LORADS_ZERO(ASolver->var->bLinSys[iCone], double, n);
    axpy(&(n), &scalFactor, M2->matElem, &incx, ASolver->var->bLinSys[iCone], &incx);

    admmCG *MMat;
    LORADS_INIT(MMat, admmCG, 1);
    MMat->ACone = ACone;
    MMat->noUpdateVar = noUpdateVar;
    lorads_sdp_dense *shell;
    LORADS_INIT(shell, lorads_sdp_dense, 1);
    MMat->UpdateVarShell = shell;
    MMat->weight = M1;

    CGSetData(ASolver->CGLinsys[iCone], MMat, ADMMUpdateUVMvec);
    CGSolve(ASolver->CGLinsys[iCone], updateVar->matElem, ASolver->var->bLinSys[iCone], CG_tol, CG_maxIter);
    ASolver->cgTime += ASolver->CGLinsys[iCone]->cgDuration;
    ASolver->cgIter += ASolver->CGLinsys[iCone]->iter;
    LORADS_FREE(shell);
    LORADS_FREE(MMat);
}


extern void LORADSUpdateSDPVarOne_positive_S(lorads_solver *ASolver, lorads_sdp_dense *updateVar, lorads_sdp_dense *noUpdateVar, lorads_sdp_dense *S, lorads_int iCone, double rho, double CG_tol, lorads_int CG_maxIter){
    lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
    double *b = ASolver->rowRHS;
    double *M1 = ASolver->var->M1temp;
    LORADS_MEMCPY(M1, b, double, ASolver->nRows);
    lorads_sdp_dense *M2 = ASolver->var->M2temp[iCone];
    double One = 1.0;
    double minusOne = -1.0;
    lorads_int incx = 1;
    // M1 = -b
    scal(&(ASolver->nRows), &(minusOne), M1, &incx);
    // M1 = sum_{j\neq \tilde{j}}A(X) - b
    axpy(&(ASolver->nRows), &(One), ASolver->var->constrValSum, &incx, M1, &incx);
    // axpy(&(ASolver->nRows), &(minusOne), ASolver->constrVal[iCone], &incx, M1, &incx);
    ASolver->var->constrVal[iCone]->add(&minusOne, ASolver->var->constrVal[iCone]->data, M1);
    // M1 = M1 * rho
    scal(&(ASolver->nRows), &(rho), M1, &incx);
    axpy(&(ASolver->nRows), &(minusOne), ASolver->var->dualVar, &incx, M1, &incx);

    // Add C in the following function
    // and weighted sum of all sdp data in one cone
    ACone->sdp_obj_sum->zeros(ACone->sdp_obj_sum->dataMat);
    ACone->addObjCoeff(ACone->coneData, ACone->sdp_obj_sum);
    ACone->sdpDataWSum(ACone->coneData, M1, ACone->sdp_obj_sum);

    // (C-(lambda + weighted Sum A)) * V
    ACone->sdp_obj_sum->mul_rk(ACone->sdp_obj_sum->dataMat, noUpdateVar, M2->matElem);

    // M2 = (M2 -rho * V + S) or (M2 -rho * V - S)
    LORADSRkMatSub_positive_S(M2, noUpdateVar, rho, S);

    // set bv in linear system
    double scalFactor = -1.0 / rho;
    lorads_int n = M2->nRows * M2->rank;
    LORADS_ZERO(ASolver->var->bLinSys[iCone], double, n);
    axpy(&(n), &scalFactor, M2->matElem, &incx, ASolver->var->bLinSys[iCone], &incx);

    admmCG *MMat;
    LORADS_INIT(MMat, admmCG, 1);
    MMat->ACone = ACone;
    MMat->noUpdateVar = noUpdateVar;
    lorads_sdp_dense *shell;
    LORADS_INIT(shell, lorads_sdp_dense, 1);
    MMat->UpdateVarShell = shell;
    MMat->weight = M1;

    CGSetData(ASolver->CGLinsys[iCone], MMat, ADMMUpdateUVMvec);
    CGSolve(ASolver->CGLinsys[iCone], updateVar->matElem, ASolver->var->bLinSys[iCone], CG_tol, CG_maxIter);
    ASolver->cgTime += ASolver->CGLinsys[iCone]->cgDuration;
    ASolver->cgIter += ASolver->CGLinsys[iCone]->iter;
    LORADS_FREE(shell);
    LORADS_FREE(MMat);
}



extern void LORADSUpdateSDPVarOne_negative_S(lorads_solver *ASolver, lorads_sdp_dense *updateVar, lorads_sdp_dense *noUpdateVar, lorads_sdp_dense *S, lorads_int iCone, double rho, double CG_tol, lorads_int CG_maxIter){
    lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
    double *b = ASolver->rowRHS;
    double *M1 = ASolver->var->M1temp;
    LORADS_MEMCPY(M1, b, double, ASolver->nRows);
    lorads_sdp_dense *M2 = ASolver->var->M2temp[iCone];
    double One = 1.0;
    double minusOne = -1.0;
    lorads_int incx = 1;
    // M1 = -b
    scal(&(ASolver->nRows), &(minusOne), M1, &incx);
    // M1 = sum_{j\neq \tilde{j}}A(X) - b
    axpy(&(ASolver->nRows), &(One), ASolver->var->constrValSum, &incx, M1, &incx);
    // axpy(&(ASolver->nRows), &(minusOne), ASolver->constrVal[iCone], &incx, M1, &incx);
    ASolver->var->constrVal[iCone]->add(&minusOne, ASolver->var->constrVal[iCone]->data, M1);
    // M1 = M1 * rho
    scal(&(ASolver->nRows), &(rho), M1, &incx);
    axpy(&(ASolver->nRows), &(minusOne), ASolver->var->dualVar, &incx, M1, &incx);

    // Add C in the following function
    // and weighted sum of all sdp data in one cone
    ACone->sdp_obj_sum->zeros(ACone->sdp_obj_sum->dataMat);
    ACone->addObjCoeff(ACone->coneData, ACone->sdp_obj_sum);
    ACone->sdpDataWSum(ACone->coneData, M1, ACone->sdp_obj_sum);

    // (C-(lambda + weighted Sum A)) * V
    ACone->sdp_obj_sum->mul_rk(ACone->sdp_obj_sum->dataMat, noUpdateVar, M2->matElem);

    // M2 = (M2 -rho * V + S) or (M2 -rho * V - S)
    LORADSRkMatSub_negative_S(M2, noUpdateVar, rho, S);

    // set bv in linear system
    double scalFactor = -1.0 / rho;
    lorads_int n = M2->nRows * M2->rank;
    LORADS_ZERO(ASolver->var->bLinSys[iCone], double, n);
    axpy(&(n), &scalFactor, M2->matElem, &incx, ASolver->var->bLinSys[iCone], &incx);

    admmCG *MMat;
    LORADS_INIT(MMat, admmCG, 1);
    MMat->ACone = ACone;
    MMat->noUpdateVar = noUpdateVar;
    lorads_sdp_dense *shell;
    LORADS_INIT(shell, lorads_sdp_dense, 1);
    MMat->UpdateVarShell = shell;
    MMat->weight = M1;

    CGSetData(ASolver->CGLinsys[iCone], MMat, ADMMUpdateUVMvec);
    CGSolve(ASolver->CGLinsys[iCone], updateVar->matElem, ASolver->var->bLinSys[iCone], CG_tol, CG_maxIter);
    ASolver->cgTime += ASolver->CGLinsys[iCone]->cgDuration;
    ASolver->cgIter += ASolver->CGLinsys[iCone]->iter;
    LORADS_FREE(shell);
    LORADS_FREE(MMat);
}



extern void LORADSUpdateLPVarOne(lorads_solver *ASolver,  double *UpdateVar, double *noUpdateVar, lorads_int iCol, double rho)
{
    double one = 1.0;
    double minusOne = -1.0;
    lorads_int incx = 1;
    lorads_lp_cone *lp_cone = ASolver->lpCone;
    double *b = ASolver->rowRHS;
    double *M1 = ASolver->var->M1temp;
    double zero = 0.0;
    scal(&ASolver->nRows, &zero, M1, &incx);
    axpy(&ASolver->nRows, &one, b, &incx, M1, &incx);

//    double *s = &ASolver->var->vlagLp->matElem[iCol];
    // M1 = -b
    scal(&(ASolver->nRows), &(minusOne), M1, &incx);
    // M1 = sum_{j \neq \tilde{j}} A(X) - b
    axpy(&(ASolver->nRows), &(one), ASolver->var->constrValSum, &incx, M1, &incx);
    ASolver->var->constrValLP[iCol]->add(&minusOne, ASolver->var->constrValLP[iCol]->data, M1);
    // M1 = M1 * rho
    scal(&(ASolver->nRows), &(rho), M1, &incx);
    // M1 = M1 - lambd
    axpy(&(ASolver->nRows), &(minusOne), ASolver->var->dualVar, &incx, M1, &incx);

    double lpWObjSum = 0.0;
    lp_cone->objCoeffSum(lp_cone->coneData, &lpWObjSum, iCol);
    lp_cone->lpDataWSum(lp_cone->coneData, M1, &lpWObjSum, iCol);

    double M2 = lpWObjSum * noUpdateVar[0];
    // flag = 1 for U, -1 for V
    M2 = M2 - rho * noUpdateVar[0];
    double blinSys = -1.0 * M2 / rho;
    lorads_lp_cone_data *lpConeData = (lorads_lp_cone_data *)lp_cone->coneData;
    UpdateVar[0] = blinSys / (1 + lpConeData->nrm2Square[iCol] * noUpdateVar[0] * noUpdateVar[0]);
}


extern void LORADSUpdateLPVarOne_positive_S(lorads_solver *ASolver,  double *UpdateVar, double *noUpdateVar, lorads_int iCol, double rho, double *sLp)
{
    double one = 1.0;
    double minusOne = -1.0;
    lorads_int incx = 1;
    lorads_lp_cone *lp_cone = ASolver->lpCone;
    double *b = ASolver->rowRHS;
    double *M1 = ASolver->var->M1temp;
    double zero = 0.0;
    scal(&ASolver->nRows, &zero, M1, &incx);
    axpy(&ASolver->nRows, &one, b, &incx, M1, &incx);

//    double *s = &ASolver->var->vlagLp->matElem[iCol];
    // M1 = -b
    scal(&(ASolver->nRows), &(minusOne), M1, &incx);
    // M1 = sum_{j \neq \tilde{j}} A(X) - b
    axpy(&(ASolver->nRows), &(one), ASolver->var->constrValSum, &incx, M1, &incx);
    ASolver->var->constrValLP[iCol]->add(&minusOne, ASolver->var->constrValLP[iCol]->data, M1);
    // M1 = M1 * rho
    scal(&(ASolver->nRows), &(rho), M1, &incx);
    // M1 = M1 - lambd
    axpy(&(ASolver->nRows), &(minusOne), ASolver->var->dualVar, &incx, M1, &incx);

    double lpWObjSum = 0.0;
    lp_cone->objCoeffSum(lp_cone->coneData, &lpWObjSum, iCol);
    lp_cone->lpDataWSum(lp_cone->coneData, M1, &lpWObjSum, iCol);

    double M2 = lpWObjSum * noUpdateVar[0];
    // flag = 1 for U, -1 for V
    M2 = M2 - rho * noUpdateVar[0] + sLp[iCol];
    double blinSys = -1.0 * M2 / rho;
    lorads_lp_cone_data *lpConeData = (lorads_lp_cone_data *)lp_cone->coneData;
    UpdateVar[0] = blinSys / (1 + lpConeData->nrm2Square[iCol] * noUpdateVar[0] * noUpdateVar[0]);
}


extern void LORADSUpdateLPVarOne_negative_S(lorads_solver *ASolver,  double *UpdateVar, double *noUpdateVar, lorads_int iCol, double rho, double *sLp)
{
    double one = 1.0;
    double minusOne = -1.0;
    lorads_int incx = 1;
    lorads_lp_cone *lp_cone = ASolver->lpCone;
    double *b = ASolver->rowRHS;
    double *M1 = ASolver->var->M1temp;
    double zero = 0.0;
    scal(&ASolver->nRows, &zero, M1, &incx);
    axpy(&ASolver->nRows, &one, b, &incx, M1, &incx);

//    double *s = &ASolver->var->vlagLp->matElem[iCol];
    // M1 = -b
    scal(&(ASolver->nRows), &(minusOne), M1, &incx);
    // M1 = sum_{j \neq \tilde{j}} A(X) - b
    axpy(&(ASolver->nRows), &(one), ASolver->var->constrValSum, &incx, M1, &incx);
    ASolver->var->constrValLP[iCol]->add(&minusOne, ASolver->var->constrValLP[iCol]->data, M1);
    // M1 = M1 * rho
    scal(&(ASolver->nRows), &(rho), M1, &incx);
    // M1 = M1 - lambd
    axpy(&(ASolver->nRows), &(minusOne), ASolver->var->dualVar, &incx, M1, &incx);

    double lpWObjSum = 0.0;
    lp_cone->objCoeffSum(lp_cone->coneData, &lpWObjSum, iCol);
    lp_cone->lpDataWSum(lp_cone->coneData, M1, &lpWObjSum, iCol);

    double M2 = lpWObjSum * noUpdateVar[0];
    // flag = 1 for U, -1 for V
    M2 = M2 - rho * noUpdateVar[0] - sLp[iCol];
    double blinSys = -1.0 * M2 / rho;
    lorads_lp_cone_data *lpConeData = (lorads_lp_cone_data *)lp_cone->coneData;
    UpdateVar[0] = blinSys / (1 + lpConeData->nrm2Square[iCol] * noUpdateVar[0] * noUpdateVar[0]);
}
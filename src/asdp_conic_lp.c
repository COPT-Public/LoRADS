#ifdef HEADERPATH
#include "interface/def_asdp.h"
#include "interface/asdp_conic_sdp.h"
#include "interface/def_asdp_user_data.h"
#include "interface/asdp_utils.h"
#include "interface/def_asdp_schur.h"
#include "src/asdp_sdpdata.h"
#include "src/sparse_opts.h"
#include "src/dense_opts.h"
#include "src/vec_opts.h"
#include "src/asdp_linsolver.h"
#include "external/asdp_cs.h"
#else
#include "def_asdp.h"
#include "asdp_conic_sdp.h"
#include "def_asdp_user_data.h"
#include "asdp_utils.h"
#include "asdp_sdpdata.h"
#include "sparse_opts.h"
#include "dense_opts.h"
#include "vec_opts.h"
#include "asdp_cs.h"
#include "asdp_lpdata.h"
#endif

#include <math.h>


extern asdp_retcode LPConeCreateImpl( void **pConeIn ) {
    asdp_cone_lp **pCone = (asdp_cone_lp **)pConeIn;
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    ASDP_NULLCHECK(pCone);
    asdp_cone_lp *cone = NULL;
    
    ASDP_INIT(cone, asdp_cone_lp, 1);
    ASDP_MEMCHECK(cone);
    ASDP_ZERO(cone, asdp_cone_lp, 1);
    *pCone = cone;
    
exit_cleanup:
    return retcode;
}

extern asdp_retcode LPConeProcDataImpl( void *coneIn, int nRow, int nCol, int *coneMatBeg, int *coneMatIdx, double *coneMatElem ) {
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    cone->nRow = nRow;
    cone->nCol = nCol;
    
    ASDP_INIT(cone->objMatElem, double, nCol);
    ASDP_MEMCHECK(cone->objMatElem);
    ASDP_ZERO(cone->objMatElem, double, nCol);
    int nObjNnz = coneMatBeg[1];
    
    ASDP_INIT(cone->rowMatBeg, int, nRow + 1);
    ASDP_MEMCHECK(cone->rowMatBeg);
    ASDP_INIT(cone->rowMatIdx, int, coneMatBeg[nRow + 1] - nObjNnz);
    ASDP_MEMCHECK(cone->rowMatIdx);
    ASDP_INIT(cone->rowMatElem, double, coneMatBeg[nRow + 1] - nObjNnz);
    ASDP_MEMCHECK(cone->rowMatElem);
    
    /* Copy LP coefficient */
    for ( int iElem = 0; iElem < coneMatBeg[1]; ++iElem ) {
        cone->objMatElem[coneMatIdx[iElem]] = coneMatElem[iElem];
    }
    
    /* Copy LP constraint matrix */
    ASDP_MEMCPY(cone->rowMatBeg, coneMatBeg + 1, int, nRow + 1);
    ASDP_MEMCPY(cone->rowMatIdx, coneMatIdx + nObjNnz, int, coneMatBeg[nRow + 1] - nObjNnz);
    ASDP_MEMCPY(cone->rowMatElem, coneMatElem + nObjNnz, double, coneMatBeg[nRow + 1] - nObjNnz);
    
    /* Remove the shift of objective */
    for ( int iCol = 0; iCol < nRow + 1; ++iCol ) {
        cone->rowMatBeg[iCol] -= nObjNnz;
    }
    
exit_cleanup:
    return retcode;
}

extern asdp_retcode LPConePresolveImpl( void *coneIn ) {
    asdp_retcode retcode = ASDP_RETCODE_OK;
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    
    ASDP_INIT(cone->lpCol, lp_coeff *, cone->nCol);
    int nnz = 0;
    ASDP_INIT(cone->nrm2Square, double, cone->nCol);
    double temp = 0.0;
    int incx = 1;
    
    
    int *nnzStat;
    ASDP_INIT(nnzStat, int, cone->nCol);
    ASDP_ZERO(nnzStat, int, cone->nCol);
    double **lpCoeffTemp;
    ASDP_INIT(lpCoeffTemp, double *, cone->nCol);
    int **lpCoeffIdxTemp;
    ASDP_INIT(lpCoeffIdxTemp, int *, cone->nCol);
    
    for (int iRow = 0; iRow < cone->nRow; ++iRow){
        int nnzRow = cone->rowMatBeg[iRow + 1] - cone->rowMatBeg[iRow];
        int *nnzRowIdx = &cone->rowMatIdx[cone->rowMatBeg[iRow]];
        for (int i = 0; i < nnzRow; ++i){
            nnzStat[nnzRowIdx[i]] += 1;
        }
    }
    
    for (int iCol = 0; iCol < cone->nCol; ++iCol){
        ASDP_INIT(lpCoeffTemp[iCol], double, nnzStat[iCol]);
        ASDP_ZERO(lpCoeffTemp[iCol], double, nnzStat[iCol]);
        ASDP_INIT(lpCoeffIdxTemp[iCol], int, nnzStat[iCol]);
        ASDP_ZERO(lpCoeffIdxTemp[iCol], int, nnzStat[iCol]);
    }
    ASDP_ZERO(nnzStat, int, cone->nCol);
    for (int iRow = 0; iRow < cone->nRow; ++iRow){
        int nnzRow = cone->rowMatBeg[iRow + 1] - cone->rowMatBeg[iRow];
        int *nnzRowIdx = &cone->rowMatIdx[cone->rowMatBeg[iRow]];
        double *elem = &cone->rowMatElem[cone->rowMatBeg[iRow]];
        for (int i = 0; i < nnzRow; ++i){
            int iCol = nnzRowIdx[i];
            lpCoeffIdxTemp[iCol][nnzStat[iCol]] = iRow;
            lpCoeffTemp[iCol][nnzStat[iCol]] = elem[i];
            nnzStat[iCol] += 1;
        }
    }
    
    
    for (int iCol = 0; iCol < cone->nCol; ++iCol){
        ASDP_INIT(cone->lpCol[iCol], lp_coeff, 1);
        // choose datatype
        nnz = nnzStat[iCol];
        if ( (double)nnz / (double)cone->nRow < 0.25 && nnz != 0){
            LPDataMatIChooseType(cone->lpCol[iCol], LP_COEFF_SPARSE);
        }else if (nnz == 0){
            LPDataMatIChooseType(cone->lpCol[iCol], LP_COEFF_ZERO);
        }else{
            LPDataMatIChooseType(cone->lpCol[iCol], LP_COEFF_DENSE);
        }
        cone->lpCol[iCol]->create(&cone->lpCol[iCol]->dataMat, cone->nRow, nnz, lpCoeffIdxTemp[iCol], lpCoeffTemp[iCol]);
        
        temp = nrm2(&nnzStat[iCol], lpCoeffTemp[iCol], &incx);
        cone->nrm2Square[iCol] = temp * temp;
    }
    // no need the three element any more
    ASDP_FREE(cone->rowMatBeg);
    ASDP_FREE(cone->rowMatIdx);
    ASDP_FREE(cone->rowMatElem);
    // free auxi var
    for (int iCol = 0; iCol < cone->nCol; ++iCol){
        ASDP_FREE(lpCoeffTemp[iCol]);
        ASDP_FREE(lpCoeffIdxTemp[iCol]);
    }
    ASDP_FREE(lpCoeffTemp);
    ASDP_FREE(lpCoeffIdxTemp);
    ASDP_FREE(nnzStat);
exit_cleanup:
    return retcode;
}

extern double LPConeGetObjNorm( asdp_cone_lp *cone, int whichNorm ) {
    
    double dNorm = 0.0;
    if ( whichNorm == ABS_NORM ) {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            dNorm += fabs(cone->objMatElem[iCol]);
        }
    } else {
        for ( int iCol = 0; iCol < cone->nCol; ++iCol ) {
            dNorm += cone->objMatElem[iCol] * cone->objMatElem[iCol];
        }
        
        dNorm = sqrt(dNorm);
    }
    
    return dNorm;
}

extern double LPConeGetCoeffNorm( asdp_cone_lp *cone, int whichNorm ) {
    
    if ( whichNorm == ABS_NORM ) {
        return csp_sum_abs(cone->nRow, cone->rowMatBeg, cone->rowMatIdx, cone->rowMatElem);
    } else {
        return csp_fro_norm(cone->nRow, cone->rowMatBeg, cone->rowMatIdx, cone->rowMatElem);
    }
    
    assert( 0 );
}

extern void LPConeScal( void *coneIn, double dScal ) {
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    rscl(&cone->nCol, &dScal, cone->objMatElem, &AIntConstantOne);
    return;
}

extern void LPConeConeAUV(void *coneIn, asdp_rk_mat_lp *uLp, asdp_rk_mat_lp *vLp, double *AUV, int iCol){
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    double uv = uLp->matElem[iCol] * vLp->matElem[iCol];
    cone->lpCol[iCol]->mul_inner_rk_double(cone->lpCol[iCol]->dataMat, &uv, AUV);
}

extern void LPConeConeAUV2(void *coneIn, double *uvLp, double *AUVSum){
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    int nCol = cone->nCol;
    int nRow = cone->nRow;
    ASDP_ZERO(AUVSum, double, nRow);
    double *AUVSumTemp;
    ASDP_INIT(AUVSumTemp, double, nRow);
    double alpha = 1.0;
    int incx = 1;
    for (int iCol = 0; iCol < nCol; ++iCol){
        cone->lpCol[iCol]->mul_inner_rk_double(cone->lpCol[iCol]->dataMat, &uvLp[iCol], AUVSumTemp);
        axpy(&nRow, &alpha, AUVSumTemp, &incx, AUVSum, &incx);
    }
    ASDP_FREE(AUVSumTemp);
}

extern void LPConeObjAUV(void *coneIn, asdp_rk_mat_lp *uLp, asdp_rk_mat_lp *vLp, double *cUV){
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    double *uv;
    ASDP_INIT(uv, double, vLp->nLPCols);
    ASDP_MEMCPY(uv, vLp->matElem, double, vLp->nLPCols);
    
    vvscl(&(vLp->nLPCols), uLp->matElem, uv);
    int incx = 1;
    cUV[0] += dot(&(vLp->nLPCols), cone->objMatElem, &incx, uv, &incx);
    
    ASDP_FREE(uv);
}

extern void LPConeObjNrm1(void *coneIn, double *nrm1Val, int nLpCols){
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    int incx = 1;
    nrm1Val[0] = nrm1(&nLpCols, cone->objMatElem, &incx);
}

extern void LPConeObjNrm2Square(void *coneIn, double *nrm2Val, int nLpCols){
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    int incx = 1;
    nrm2Val[0] = nrm1(&nLpCols, cone->objMatElem, &incx);
    nrm2Val[0] = nrm2Val[0] * nrm2Val[0];
}

extern void LPConeObjNrmInf(void *coneIn, double *nrmInf, int nLpCols){
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    int incx = 1;
    int idx = idamax(&nLpCols, cone->objMatElem, &incx);
    nrmInf[0] = ASDP_MAX(fabs(cone->objMatElem[idx]), nrmInf[0]);
}

extern void LPConeDataWsum(void *coneIn, double *weight, double *wSum, int iCol){
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    cone->lpCol[iCol]->weight_sum(cone->lpCol[iCol]->dataMat, weight, wSum);
}

extern void LPConeDataObjCoeffSum(void *coneIn, double *res, int iCol){
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    res[0] += cone->objMatElem[iCol];
}

extern void LPConeClearImpl( void *coneIn ) {
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    if ( !cone ) {
        return;
    }
    
    ASDP_FREE(cone->objMatElem);
    for (int iCol = 0; iCol < cone->nCol; ++iCol){
        cone->lpCol[iCol]->destroy(&cone->lpCol[iCol]->dataMat);
        ASDP_FREE(cone->lpCol[iCol]);
    }
    
    ASDP_FREE(cone->lpCol);
    ASDP_FREE(cone->nrm2Square);
    ASDP_ZERO(cone, asdp_cone_lp, 1);
    
    return;
}

extern void LPConeDestroyImpl( void **pConeIn ) {
    asdp_cone_lp **pCone = (asdp_cone_lp **)pConeIn;
    if ( !pCone ) {
        return;
    }
    
    LPConeClearImpl(*pCone);
    ASDP_FREE(*pCone);
    
    return;
}

extern void LPConeViewImpl( void *coneIn ) {
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    asdp_printf("LP Cone of %d variables and %d constraints \n", cone->nRow, cone->nCol);
    return;
}

extern void LPConeGetStatsImpl( void *coneIn, double *rowRHS, int coneIntFeatures[20], double coneDblFeatures[20] ) {
    asdp_cone_lp *cone = (asdp_cone_lp *)coneIn;
    /* Detect two special structures for LPs:
     
     1. Implied dual equality (free primal variable)
     2. Implied dual box constraint
    
     */
    
    if ( cone->nCol % 2 != 0 || cone->nCol < 100 ) {
        return;
    }
    
    /* Use two auxiliary arrays from the LP datas structure */
    double *dDualBoundUpperTmp = NULL;
    double *dDualBoundLowerTmp = NULL;
    
    ASDP_INIT(dDualBoundLowerTmp, double, cone->nRow);
    ASDP_INIT(dDualBoundUpperTmp, double, cone->nRow);
    ASDP_ZERO(dDualBoundLowerTmp, double, cone->nRow);
    ASDP_ZERO(dDualBoundUpperTmp, double, cone->nRow);
    
    if ( !dDualBoundLowerTmp || !dDualBoundUpperTmp ) {
        ASDP_FREE(dDualBoundLowerTmp);
        ASDP_FREE(dDualBoundUpperTmp);
        return;
    }
    
    int nHalfCols = (int) cone->nCol / 2;
    
    int isImpliedDual = 1;
    int isImpliedUpper = 0;
    int isImpliedLower = 0;
    
    double dMaxDualBoundUpper = 1.0;
    double dMinDualBoundLower = -1.0;
    
    
    /* First detect if l <= y <= u */
    for ( int iCol = 0; iCol < cone->nRow; ++iCol ) {
        for ( int iElem = cone->rowMatBeg[iCol]; iElem < cone->rowMatBeg[iCol + 1]; ++iElem ) {
            
            if ( cone->rowMatBeg[iCol + 1] - cone->rowMatBeg[iCol] > 2 ) {
                isImpliedDual = 0;
                break;
            }
            
            if ( cone->rowMatElem[iElem] > 0.0 ) {
                if ( dDualBoundUpperTmp[iCol] ) {
                    isImpliedDual = 0;
                    break;
                }
                isImpliedUpper = 1;
                double dUpperBound = cone->objMatElem[cone->rowMatIdx[iElem]] / cone->rowMatElem[iElem];
                dDualBoundUpperTmp[iCol] = ASDP_MAX(dDualBoundUpperTmp[iCol], dUpperBound);
            } else {
                if ( dDualBoundLowerTmp[iCol] ) {
                    isImpliedDual = 0;
                    break;
                }
                isImpliedLower = 1;
                double dLowerBound = cone->objMatElem[cone->rowMatIdx[iElem]] / cone->rowMatElem[iElem];
                dDualBoundLowerTmp[iCol] = ASDP_MIN(dDualBoundLowerTmp[iCol], dLowerBound);
            }
        }
        
        if ( !isImpliedDual ) {
            break;
        }
    }
    
    if ( isImpliedDual ) {
        
        coneIntFeatures[INT_FEATURE_I_IMPYBOUND] = 1;
        
        if ( isImpliedUpper ) {
            for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
                dMaxDualBoundUpper = ASDP_MAX(dMaxDualBoundUpper, dDualBoundUpperTmp[iRow]);
            }
            
            if ( dMaxDualBoundUpper <= 0.0 ) {
                dMaxDualBoundUpper = 1.0;
            }
            
            coneDblFeatures[DBL_FEATURE_IMPYBOUNDUP] = dMaxDualBoundUpper;
        }
        
        if ( isImpliedLower ) {
            for ( int iRow = 0; iRow < cone->nRow; ++iRow ) {
                dMinDualBoundLower = ASDP_MIN(dMinDualBoundLower, dDualBoundLowerTmp[iRow]);
            }
            
            if ( dMinDualBoundLower >= 0.0 ) {
                dMinDualBoundLower = -1.0;
            }
            
            coneDblFeatures[DBL_FEATURE_IMPYBOUNDLOW] = dMinDualBoundLower;
        }
    }
    
    ASDP_FREE(dDualBoundLowerTmp);
    ASDP_FREE(dDualBoundUpperTmp);
    
    for ( int iCol = 0; iCol < nHalfCols; ++iCol ) {
        if ( cone->objMatElem[iCol] + cone->objMatElem[iCol + nHalfCols] != 0.0 ) {
            return;
        }
    }
    
    for ( int iCol = 0; iCol < cone->nRow; ++iCol ) {
        int colNnz = cone->rowMatBeg[iCol + 1] - cone->rowMatBeg[iCol];
        int nHalfNnz = (int) colNnz / 2;
        if ( colNnz % 2 != 0 ) {
            return;
        }
        
        for ( int iRow = 0; iRow < nHalfNnz; ++iRow ) {
            if ( cone->rowMatElem[cone->rowMatBeg[iCol] + iRow] + \
                cone->rowMatElem[cone->rowMatBeg[iCol] + iRow + nHalfNnz] != 0.0 ) {
                return;
            }
        }
    }
    
    coneIntFeatures[INT_FEATURE_I_NODINTERIOR] = 1;
    
    return;
}

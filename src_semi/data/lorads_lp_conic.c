


#include "lorads.h"
#include "def_lorads_lp_conic.h"
#include "lorads_vec_opts.h"
#include "lorads_sdp_data.h"
#include "lorads_lp_data.h"

extern void lp_cone_create(lorads_lp_cone_data **coneData){
    lorads_lp_cone_data *conePtr;
    LORADS_INIT(conePtr, lorads_lp_cone_data, 1);
    LORADS_MEMCHECK(conePtr);
    *coneData = conePtr;
}

extern void lp_cone_proc(lorads_lp_cone_data *cone, lorads_int nRow, lorads_int nCol, lorads_int *lpMatBeg, lorads_int *lpMatIdx, double *lpMatElem){
    cone->nRow = nRow;
    cone->nCol = nCol;

    LORADS_INIT(cone->objMatElem, double, nCol);
    LORADS_MEMCHECK(cone->objMatElem);
    LORADS_ZERO(cone->objMatElem, double, nCol);
    lorads_int nObjNnz = lpMatBeg[1];

    LORADS_INIT(cone->rowMatBeg, lorads_int, nRow + 1);
    LORADS_MEMCHECK(cone->rowMatBeg);
    LORADS_INIT(cone->rowMatIdx, lorads_int, lpMatBeg[nRow + 1] - nObjNnz);
    LORADS_MEMCHECK(cone->rowMatIdx);
    LORADS_INIT(cone->rowMatElem, double, lpMatBeg[nRow + 1] - nObjNnz);
    LORADS_MEMCHECK(cone->rowMatElem);

    /* Copy LP coefficient */
    for ( lorads_int iElem = 0; iElem < lpMatBeg[1]; ++iElem ) {
        cone->objMatElem[lpMatIdx[iElem]] = lpMatElem[iElem];
    }

    /* Copy LP constraint matrix */
    LORADS_MEMCPY(cone->rowMatBeg, lpMatBeg + 1, lorads_int, nRow + 1);
    LORADS_MEMCPY(cone->rowMatIdx, lpMatIdx + nObjNnz, lorads_int, lpMatBeg[nRow + 1] - nObjNnz);
    LORADS_MEMCPY(cone->rowMatElem, lpMatElem + nObjNnz, double, lpMatBeg[nRow + 1] - nObjNnz);

    /* Remove the shift of objective */
    for ( lorads_int iCol = 0; iCol < nRow + 1; ++iCol ) {
        cone->rowMatBeg[iCol] -= nObjNnz;
    }
}

extern void lp_cone_presolve(lorads_lp_cone_data *cone){
    LORADS_INIT(cone->lpCol, lp_coeff *, cone->nCol);
    lorads_int nnz = 0;
    LORADS_INIT(cone->nrm2Square, double, cone->nCol);
    double temp = 0.0;
    lorads_int incx = 1;


    lorads_int *nnzStat;
    LORADS_INIT(nnzStat, lorads_int, cone->nCol);
    LORADS_ZERO(nnzStat, lorads_int, cone->nCol);
    double **lpCoeffTemp;
    LORADS_INIT(lpCoeffTemp, double *, cone->nCol);
    lorads_int **lpCoeffIdxTemp;
    LORADS_INIT(lpCoeffIdxTemp, lorads_int *, cone->nCol);

    for (int iRow = 0; iRow < cone->nRow; ++iRow){
        lorads_int nnzRow = cone->rowMatBeg[iRow + 1] - cone->rowMatBeg[iRow];
        lorads_int *nnzRowIdx = &cone->rowMatIdx[cone->rowMatBeg[iRow]];
        for (int i = 0; i < nnzRow; ++i){
            nnzStat[nnzRowIdx[i]] += 1;
        }
    }

    for (int iCol = 0; iCol < cone->nCol; ++iCol){
        LORADS_INIT(lpCoeffTemp[iCol], double, nnzStat[iCol]);
        LORADS_ZERO(lpCoeffTemp[iCol], double, nnzStat[iCol]);
        LORADS_INIT(lpCoeffIdxTemp[iCol], lorads_int, nnzStat[iCol]);
        LORADS_ZERO(lpCoeffIdxTemp[iCol], lorads_int, nnzStat[iCol]);
    }
    LORADS_ZERO(nnzStat, lorads_int, cone->nCol);
    for (int iRow = 0; iRow < cone->nRow; ++iRow){
        lorads_int nnzRow = cone->rowMatBeg[iRow + 1] - cone->rowMatBeg[iRow];
        lorads_int *nnzRowIdx = &cone->rowMatIdx[cone->rowMatBeg[iRow]];
        double *elem = &cone->rowMatElem[cone->rowMatBeg[iRow]];
        for (int i = 0; i < nnzRow; ++i){
            lorads_int iCol = nnzRowIdx[i];
            lpCoeffIdxTemp[iCol][nnzStat[iCol]] = iRow;
            lpCoeffTemp[iCol][nnzStat[iCol]] = elem[i];
            nnzStat[iCol] += 1;
        }
    }


    for (int iCol = 0; iCol < cone->nCol; ++iCol){
        LORADS_INIT(cone->lpCol[iCol], lp_coeff, 1);
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
    LORADS_FREE(cone->rowMatBeg);
    LORADS_FREE(cone->rowMatIdx);
    LORADS_FREE(cone->rowMatElem);
    // free auxi var
    for (int iCol = 0; iCol < cone->nCol; ++iCol){
        LORADS_FREE(lpCoeffTemp[iCol]);
        LORADS_FREE(lpCoeffIdxTemp[iCol]);
    }
    LORADS_FREE(lpCoeffTemp);
    LORADS_FREE(lpCoeffIdxTemp);
    LORADS_FREE(nnzStat);
}

extern void lp_cone_obj_nrm1(lorads_lp_cone_data *coneData, double *nrm1Val, lorads_int nLpCols){
    lorads_int incx = 1;
    nrm1Val[0] = nrm1(&nLpCols, coneData->objMatElem, &incx);
}

extern void lp_cone_obj_nrm2Square(lorads_lp_cone_data *coneData, double *nrm2Val, lorads_int nLpCols){
    lorads_int incx = 1;
    double temp;
    temp = nrm1(&nLpCols, coneData->objMatElem, &incx);
    nrm2Val[0] = temp * temp;
}

extern void lp_cone_obj_nrmInf(lorads_lp_cone_data *coneData, double *nrmInf, lorads_int nLpCols){
    lorads_int incx = 1;
#ifdef UNDER_BLAS
    lorads_int idx = idamax_(&nLpCols, coneData->objMatElem, &incx);
#else
    lorads_int idx = idamax(&nLpCols, coneData->objMatElem, &incx);
#endif
    nrmInf[0] = LORADS_MAX(fabs(coneData->objMatElem[idx]), nrmInf[0]);
}

extern void destroy_lp_cone_data(lorads_lp_cone_data **coneData){
    lorads_lp_cone_data *data = *coneData;
    LORADS_FREE(data->objMatElem);
    LORADS_FREE(data->rowMatBeg);
    LORADS_FREE(data->rowMatIdx);
    LORADS_FREE(data->rowMatElem);
    LORADS_FREE(data->nrm2Square);
    for (lorads_int i = 0; i < data->nCol; ++i){
        data->lpCol[i]->destroy(&data->lpCol[i]->dataMat);
        LORADS_FREE(data->lpCol[i]);
    }
    LORADS_FREE(data->lpCol);
}

extern void lp_cone_view(lorads_lp_cone_data *cone ) {
#ifdef lorads_int32
    printf("LP Cone of %d variables and %d constraints \n", cone->nRow, cone->nCol);
#endif
#ifdef UNIX_INT64
    printf("LP Cone of %ld variables and %ld constraints \n", cone->nRow, cone->nCol);
#endif
#ifdef MAC_INT64
    printf("LP Cone of %lld variables and %lld constraints \n", cone->nRow, cone->nCol);
#endif
    return;
}

extern void lp_cone_AUV(lorads_lp_cone_data *coneData , lorads_lp_dense *uLp, lorads_lp_dense *vLp, double *AUV, lorads_int iCol){
    double uv = uLp->matElem[iCol] * vLp->matElem[iCol];
    coneData->lpCol[iCol]->mul_inner_rk_double(coneData->lpCol[iCol]->dataMat, &uv, AUV);
}

extern void lp_cone_AUV2(lorads_lp_cone_data *cone, double *uvLp, double *AUVSum){
    lorads_int nCol = cone->nCol;
    lorads_int nRow = cone->nRow;
    LORADS_ZERO(AUVSum, double, nRow);
//    double *AUVSumTemp;
//    LORADS_INIT(AUVSumTemp, double, nRow);
//    double alpha = 1.0;
//    lorads_int incx = 1;
    for (int iCol = 0; iCol < nCol; ++iCol){
        cone->lpCol[iCol]->mul_inner_rk_double(cone->lpCol[iCol]->dataMat, &uvLp[iCol], AUVSum);
//        cone->lpCol[iCol]->mul_inner_rk_double(cone->lpCol[iCol]->dataMat, &uvLp[iCol], AUVSumTemp);
//        axpy(&nRow, &alpha, AUVSumTemp, &incx, AUVSum, &incx);
    }
//    LORADS_FREE(AUVSumTemp);
}

extern void lp_cone_objAUV(lorads_lp_cone_data *cone, lorads_lp_dense *uLp, lorads_lp_dense *vLp, double *cUV){
    double *uv;
    LORADS_INIT(uv, double, vLp->nCols);
    LORADS_MEMCPY(uv, vLp->matElem, double, vLp->nCols);

    vvscl(&(vLp->nCols), uLp->matElem, uv);
    lorads_int incx = 1;
    cUV[0] += dot(&(vLp->nCols), cone->objMatElem, &incx, uv, &incx);
    LORADS_FREE(uv);
}

extern void lp_cone_scalObj(lorads_lp_cone_data *cone, double alpha){
    lorads_int nCol = cone->nCol;
    lorads_int incx = 1;
    double *objMatElem = cone->objMatElem;
    scal(&nCol, &alpha, objMatElem, &incx);
}

extern void lp_cone_Wsum(lorads_lp_cone_data *cone, double *weight, double *wSum, lorads_int iCol){
    cone->lpCol[iCol]->weight_sum(cone->lpCol[iCol]->dataMat, weight, wSum);
}

extern void lp_cone_ObjCoeffSum(lorads_lp_cone_data *cone, double *res, lorads_int iCol){
    res[0] += cone->objMatElem[iCol];
}

extern void LORADSSetLpCone(lorads_lp_cone *lp_cone, lorads_int nRows,
                            lorads_int nLpCols, lorads_int *lpMatBeg,
                            lorads_int *lpMatIdx, double *LpMatElem){
    lp_cone->nCol = nLpCols;
    lp_cone->coneObjNrm1 = lp_cone_obj_nrm1;
    lp_cone->coneObjNrm2Square = lp_cone_obj_nrm2Square;
    lp_cone->coneView = lp_cone_view;
    lp_cone->coneAUV = lp_cone_AUV;
    lp_cone->coneAUV2 = lp_cone_AUV2;
    lp_cone->objAUV = lp_cone_objAUV;
    lp_cone->coneObjNrmInf = lp_cone_obj_nrmInf;
    lp_cone->lpDataWSum = lp_cone_Wsum;
    lp_cone->objCoeffSum = lp_cone_ObjCoeffSum;
    lp_cone->destroyConeData = destroy_lp_cone_data;
    lp_cone->scalObj = lp_cone_scalObj;

    lp_cone_create(&lp_cone->coneData);
    lp_cone_proc(lp_cone->coneData, nRows, nLpCols, lpMatBeg, lpMatIdx, LpMatElem);
    lp_cone_presolve(lp_cone->coneData);
}

#include <assert.h>
#include "lorads_lp_data.h"
#include "def_lorads_lp_data.h"
#include "lorads_utils.h"
#include "lorads_vec_opts.h"


static void LPdataMatCreateZeroImpl(void **pA, lorads_int nRows, lorads_int nnz, lorads_int *dataMatIdx, double *dataMatElem){

    assert(nnz == 0);

    if ( !pA ) {
        LORADS_ERROR_TRACE;
    }

    lp_coeff_zero *zero;
    LORADS_INIT(zero, lp_coeff_zero, 1);
    LORADS_MEMCHECK(zero);

    zero->nRows = nRows;
    *pA = (void *) zero;
}


static void LPdataMatDestroyZeroImpl(void **pA){
    return;
}


static void LPdataMatZeroMulInnerRkDouble(void *pA, double *puv, double *constrVal){
    lp_coeff_zero *zero = (lp_coeff_zero *)pA;
    LORADS_ZERO(constrVal, double, zero->nRows);
    return;
}

static void LPdataMatZeroWeightSum(void *pA, double *weight, double *res){
    return;
}

static void LPdataMatZeroScaleData(void *pA, double scaleFactor){
    return;
}

static void LPdataMatCreateSparseImpl(void **pA, lorads_int nRows, lorads_int nnz, lorads_int *dataMatIdx, double *dataMatElem){
    if ( !pA ) {
        LORADS_ERROR_TRACE;
    }
    lp_coeff_sparse *sparse;
    LORADS_INIT(sparse, lp_coeff_sparse, 1);
    LORADS_MEMCHECK(sparse);

    sparse->nRows = nRows;
    sparse->nnz = nnz;
    LORADS_INIT(sparse->rowPtr, lorads_int, nnz);
    LORADS_MEMCHECK(sparse->rowPtr);
    LORADS_INIT(sparse->val, double, nnz);
    LORADS_MEMCHECK(sparse->val);

    LORADS_MEMCPY(sparse->rowPtr, dataMatIdx, lorads_int, nnz);
    LORADS_MEMCPY(sparse->val, dataMatElem, double, nnz);

    *pA = (void *) sparse;
}


static void LPdataMatDestroySparseImpl(void **pA){
    lp_coeff_sparse *sparse = (lp_coeff_sparse *) *pA;
    LORADS_FREE(sparse->rowPtr);
    LORADS_FREE(sparse->val);
    LORADS_FREE(sparse);
    return;
}


static void LPdataMatSparseMulInnerRkDouble(void *pA, double *puv, double *constrVal){
    lp_coeff_sparse *sparse = (lp_coeff_sparse *)pA;
//    double uv = puv[0];
//    for (lorads_int i = 0; i < sparse->nnz; ++i){
//        constrVal[i] = sparse->val[i] * uv;
//    }
    if (sparse->nnz < 5){
        for (lorads_int i = 0; i < sparse->nnz; ++i){
            constrVal[i] = puv[0] * sparse->val[i];
        }
    }else{
        LORADS_MEMCPY(constrVal, sparse->val, double, sparse->nnz);
        lorads_int incx = 1;
        scal(&sparse->nnz, puv, constrVal, &incx);
    }
    return;
}

static void LPdataMatSparseWeightSum(void *pA, double *weight, double *res){
    lp_coeff_sparse *sparse = (lp_coeff_sparse *)pA;
    for (lorads_int i = 0; i < sparse->nnz; ++i){
        res[0] += sparse->val[i] * weight[sparse->rowPtr[i]];
    }
    return;
}


static void LPdataMatSparseScaleData(void *pA, double scaleFactor){
    lp_coeff_sparse *sparse = (lp_coeff_sparse *)pA;
    vvscl(&sparse->nnz, &scaleFactor, sparse->val);
    return;
}


static void LPdataMatCreateDenseImpl(void **pA, lorads_int nRows, lorads_int nnz, lorads_int *dataMatIdx, double *dataMatElem){
    if ( !pA ) {
        LORADS_ERROR_TRACE;
    }

    lp_coeff_dense *dense;
    LORADS_INIT(dense, lp_coeff_dense, 1);
    LORADS_MEMCHECK(dense);
    LORADS_INIT(dense->val, double, nRows);
    LORADS_MEMCHECK(dense->val);

    for (lorads_int i = 0; i < dense->nRows; ++i){
        dense->val[dataMatIdx[i]] = dataMatElem[i];
    }
    *pA = (void *) dense;
}

static void LPdataMatDestroyDenseImpl(void **pA){
    lp_coeff_dense *dense = (lp_coeff_dense *) *pA;
    LORADS_FREE(dense->val);
    LORADS_ZERO(dense, lp_coeff_dense, 1);
    return;
}


static void LPdataMatDenseMulInnerRkDouble(void *pA, double *puv, double *constrVal){
    lp_coeff_dense *dense = (lp_coeff_dense *)pA;
    LORADS_MEMCPY(constrVal, dense->val, double, dense->nRows);
    double uv = puv[0];
    lorads_int incx = 1;
    scal(&(dense->nRows), &uv, constrVal, &incx);
    return;
}


static void LPdataMatDenseWeightSum(void *pA, double *weight, double *res){
    lp_coeff_dense *dense = (lp_coeff_dense *)pA;
    lorads_int incx = 1;
    res[0] += dot(&(dense->nRows), dense->val, &incx, weight, &incx);
    return;
}

static void LPdataMatDenseScaleData(void *pA, double scaleFactor){
    lp_coeff_dense *dense = (lp_coeff_dense *)pA;
    lorads_int incx = 1;
    vvscl(&dense->nRows, &scaleFactor, dense->val);
    return;
}


extern void LPDataMatIChooseType(lp_coeff *lpCoeff, lp_coeff_type dataType){
    lpCoeff->dataType = dataType;
    switch (dataType) {
        case LP_COEFF_ZERO:
            lpCoeff->create = LPdataMatCreateZeroImpl;
            lpCoeff->destroy = LPdataMatDestroyZeroImpl;
            lpCoeff->mul_inner_rk_double = LPdataMatZeroMulInnerRkDouble;
            lpCoeff->weight_sum = LPdataMatZeroWeightSum;
//            lpCoeff->scaleData = LPdataMatZeroScaleData;
            break;
        case LP_COEFF_DENSE:
            lpCoeff->create = LPdataMatCreateDenseImpl;
            lpCoeff->destroy = LPdataMatDestroyDenseImpl;
            lpCoeff->mul_inner_rk_double = LPdataMatDenseMulInnerRkDouble;
            lpCoeff->weight_sum = LPdataMatDenseWeightSum;
//            lpCoeff->scaleData = LPdataMatDenseScaleData;
            break;
        case LP_COEFF_SPARSE:
            lpCoeff->create = LPdataMatCreateSparseImpl;
            lpCoeff->destroy = LPdataMatDestroySparseImpl;
            lpCoeff->mul_inner_rk_double = LPdataMatSparseMulInnerRkDouble;
            lpCoeff->weight_sum = LPdataMatSparseWeightSum;
            break;
        default:
            assert(0);
            break;
    }
    return;
}

/* @file asdp_lpdata.h
   @brief Implement the lp coefficient data operations
 */


#ifdef HEADERPATH
#include "interface/asdp_conic.h"
#include "interface/asdp_utils.h"
#include "interface/def_asdp_sdpdata.h"
#include "src/asdp_sdpdata.h"
#include "src/vec_opts.h"
#include "src/sparse_opts.h"
#include "src/dense_opts.h"
#include "src/r1_opts.h"
#include "external/asdp_cs.h"
#else
#include "asdp_conic.h"
#include "asdp_utils.h"
#include "def_asdp_sdpdata.h"
#include "def_asdp_lpdata.h"
#include "asdp_sdpdata.h"
#include "vec_opts.h"
#include "sparse_opts.h"
#include "dense_opts.h"
#include "r1_opts.h"
#include "asdp_cs.h"
#include "def_asdp_rk_mat.h"
#include "asdp_debug.h"
#endif

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

#include <math.h>


static asdp_retcode LPdataMatCreateZeroImpl(void **pA, int nRows, int nnz, int *dataMatIdx, double *dataMatElem){
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    assert(nnz == 0);
    
    if ( !pA ) {
        retcode = ASDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    lp_coeff_zero *zero;
    ASDP_INIT(zero, lp_coeff_zero, 1);
    ASDP_MEMCHECK(zero);
    
    zero->nRows = nRows;
    *pA = (void *) zero;
    
exit_cleanup:
    return retcode;
}

static void LPdataMatDestroyZeroImpl(void **pA){
    return;
}


static void LPdataMatZeroMulInnerRkDouble(void *pA, double *puv, double *constrVal){
    lp_coeff_zero *zero = (lp_coeff_zero *)pA;
    ASDP_ZERO(constrVal, double, zero->nRows);
    return;
}

static void LPdataMatZeroWeightSum(void *pA, double *weight, double *res){
    return;
}

static void LPdataMatZeroScaleData(void *pA, double scaleFactor){
    return;
}

static asdp_retcode LPdataMatCreateSparseImpl(void **pA, int nRows, int nnz, int *dataMatIdx, double *dataMatElem){
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    if ( !pA ) {
        retcode = ASDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    lp_coeff_sparse *sparse;
    ASDP_INIT(sparse, lp_coeff_sparse, 1);
    ASDP_MEMCHECK(sparse);
    
    sparse->nRows = nRows;
    sparse->nnz = nnz;
    ASDP_INIT(sparse->rowPtr, int, nnz);
    ASDP_MEMCHECK(sparse->rowPtr);
    ASDP_INIT(sparse->val, double, nnz);
    ASDP_MEMCHECK(sparse->val);
    
    ASDP_MEMCPY(sparse->rowPtr, dataMatIdx, int, nnz);
    ASDP_MEMCPY(sparse->val, dataMatElem, double, nnz);
    
    *pA = (void *) sparse;
    
exit_cleanup:
    return retcode;
}


static void LPdataMatDestroySparseImpl(void **pA){
    lp_coeff_sparse *sparse = (lp_coeff_sparse *) *pA;
    ASDP_FREE(sparse->rowPtr);
    ASDP_FREE(sparse->val);
    ASDP_FREE(sparse);
    return;
}


static void LPdataMatSparseMulInnerRkDouble(void *pA, double *puv, double *constrVal){
    lp_coeff_sparse *sparse = (lp_coeff_sparse *)pA;
//    double uv = puv[0];
//    for (int i = 0; i < sparse->nnz; ++i){
//        constrVal[i] = sparse->val[i] * uv;
//    }
    if (sparse->nnz < 5){
        for (int i = 0; i < sparse->nnz; ++i){
            constrVal[i] = puv[0] * sparse->val[i];
        }
    }else{
        ASDP_MEMCPY(constrVal, sparse->val, double, sparse->nnz);
        int incx = 1;
        scal(&sparse->nnz, puv, constrVal, &incx);
    }
    return;
}

static void LPdataMatSparseWeightSum(void *pA, double *weight, double *res){
    lp_coeff_sparse *sparse = (lp_coeff_sparse *)pA;
    for (int i = 0; i < sparse->nnz; ++i){
        res[0] += sparse->val[i] * weight[sparse->rowPtr[i]];
    }
    return;
}


static void LPdataMatSparseScaleData(void *pA, double scaleFactor){
    lp_coeff_sparse *sparse = (lp_coeff_sparse *)pA;
    vvscl(&sparse->nnz, &scaleFactor, sparse->val);
    return;
}


static asdp_retcode LPdataMatCreateDenseImpl(void **pA, int nRows, int nnz, int *dataMatIdx, double *dataMatElem){
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    if ( !pA ) {
        retcode = ASDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    lp_coeff_dense *dense;
    ASDP_INIT(dense, lp_coeff_dense, 1);
    ASDP_MEMCHECK(dense);
    ASDP_INIT(dense->val, double, nRows);
    ASDP_MEMCHECK(dense->val);
    
    for (int i = 0; i < dense->nRows; ++i){
        dense->val[dataMatIdx[i]] = dataMatElem[i];
    }
    
    *pA = (void *) dense;
    
exit_cleanup:
    return retcode;
}

static void LPdataMatDestroyDenseImpl(void **pA){
    lp_coeff_dense *dense = (lp_coeff_dense *) *pA;
    ASDP_FREE(dense->val);
    ASDP_ZERO(dense, lp_coeff_dense, 1);
    return;
}


static void LPdataMatDenseMulInnerRkDouble(void *pA, double *puv, double *constrVal){
    lp_coeff_dense *dense = (lp_coeff_dense *)pA;
    ASDP_MEMCPY(constrVal, dense->val, double, dense->nRows);
    double uv = puv[0];
    int incx = 1;
    scal(&(dense->nRows), &uv, constrVal, &incx);
    return;
}


static void LPdataMatDenseWeightSum(void *pA, double *weight, double *res){
    lp_coeff_dense *dense = (lp_coeff_dense *)pA;
    int incx = 1;
    res[0] += dot(&(dense->nRows), dense->val, &incx, weight, &incx);
    return;
}

static void LPdataMatDenseScaleData(void *pA, double scaleFactor){
    lp_coeff_dense *dense = (lp_coeff_dense *)pA;
    int incx = 1;
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
            lpCoeff->scaleData = LPdataMatZeroScaleData;
            break;
        case LP_COEFF_DENSE:
            lpCoeff->create = LPdataMatCreateDenseImpl;
            lpCoeff->destroy = LPdataMatDestroyDenseImpl;
            lpCoeff->mul_inner_rk_double = LPdataMatDenseMulInnerRkDouble;
            lpCoeff->weight_sum = LPdataMatDenseWeightSum;
            lpCoeff->scaleData = LPdataMatDenseScaleData;
            break;
        case LP_COEFF_SPARSE:
            lpCoeff->create = LPdataMatCreateSparseImpl;
            lpCoeff->destroy = LPdataMatDestroySparseImpl;
            lpCoeff->mul_inner_rk_double = LPdataMatSparseMulInnerRkDouble;
            lpCoeff->weight_sum = LPdataMatSparseWeightSum;
            lpCoeff->scaleData = LPdataMatSparseScaleData;
            break;
        default:
            assert(0);
            break;
    }
    return;
}

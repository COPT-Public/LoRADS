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
static asdp_retcode dataMatCreateZeroImpl( void **pA, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    (void) dataMatNnz;
    (void) dataMatIdx;
    (void) dataMatElem;
    
    if ( !pA ) {
        retcode = ASDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff_zero *zero;
    ASDP_INIT(zero, sdp_coeff_zero, 1);
    ASDP_MEMCHECK(zero);
    
    zero->nSDPCol = nSDPCol;
    *pA = (void *) zero;
    
exit_cleanup:
    return retcode;
}

static void dataMatViewZeroImpl( void *A );
static void dataMatViewSparseImpl( void *A );
static void dataMatViewDenseImpl( void *A );
static void dataMatViewRankOneSparseImpl( void *A );
static void dataMatViewRankOneDenseImpl( void *A );

static asdp_retcode dataMatCreateSparseImpl( void **pA, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    if ( !pA ) {
        retcode = ASDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff_sparse *sparse;
    ASDP_INIT(sparse, sdp_coeff_sparse, 1);
    ASDP_MEMCHECK(sparse);
    
    sparse->nSDPCol = nSDPCol;
    sparse->nTriMatElem = dataMatNnz;
    
    ASDP_INIT(sparse->triMatRow, int, dataMatNnz);
    ASDP_INIT(sparse->triMatCol, int, dataMatNnz);
    ASDP_INIT(sparse->triMatElem, double, dataMatNnz);
    
    ASDP_MEMCHECK(sparse->triMatRow);
    ASDP_MEMCHECK(sparse->triMatCol);
    ASDP_MEMCHECK(sparse->triMatElem);
    
    // Note tsp_decompress work when ascending sort
    if ( !AUtilCheckIfAscending(dataMatNnz, dataMatIdx) ) {
        AUtilAscendSortDblByInt(dataMatElem, dataMatIdx, 0, dataMatNnz - 1);
    }
    
    tsp_decompress(sparse->nSDPCol, sparse->nTriMatElem, dataMatIdx, dataMatElem,
                   sparse->triMatRow, sparse->triMatCol, sparse->triMatElem);
    
    *pA = (void *) sparse;
    
exit_cleanup:
    return retcode;
}

static asdp_retcode dataMatCreateDenseImpl( void **pA, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    if ( !pA ) {
        retcode = ASDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff_dense *dense;
    ASDP_INIT(dense, sdp_coeff_dense, 1);
    ASDP_MEMCHECK(dense);
    
    dense->nSDPCol = nSDPCol;
    ASDP_INIT(dense->dsMatElem, double, PACK_NNZ(nSDPCol));
    ASDP_MEMCHECK(dense->dsMatElem);
    
    pds_decompress(dataMatNnz, dataMatIdx, dataMatElem, dense->dsMatElem);

    *pA = (void *) dense;
    
exit_cleanup:
    return retcode;
}


static void dataMatScalZeroImpl( void *A, double alpha ) {
    
    return;
}

/** @brief A = alpha \* A for sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[in] alpha Scale parameter
 */
static void dataMatScalSparseImpl( void *A, double alpha ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    int incx = 1;
    
    scal(&spA->nTriMatElem, &alpha, spA->triMatElem, &incx);
    
    return;
}

/** @brief A = alpha \* A for dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[in] alpha Scale parameter
 */
static void dataMatScalDenseImpl( void *A, double alpha ) {
    
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    int nElem = PACK_NNZ(dsA->nSDPCol);
    int incx = 1;
    
    scal(&nElem, &alpha, dsA->dsMatElem, &incx);
    
    return;
}


static double dataMatNormZeroImpl( void *A, int type ) {
    
    return 0.0;
}

/** @brief Calculate norm of sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[in] type Type of norm
 *  @return Norm of A
 */
static double dataMatNormSparseImpl( void *A, int type ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    double nrmA = 0.0;
    
    if ( type == FRO_NORM ) {
        for ( int i = 0; i < spA->nTriMatElem; ++i ) {
            nrmA += ( spA->triMatCol[i] == spA->triMatRow[i] ) ?
                      spA->triMatElem[i] * spA->triMatElem[i] :
                      spA->triMatElem[i] * spA->triMatElem[i] * 2.0;
        }
        nrmA = sqrt(nrmA);
    } else if ( type == ABS_NORM ) {
        for ( int i = 0; i < spA->nTriMatElem; ++i ) {
            nrmA += ( spA->triMatCol[i] == spA->triMatRow[i] ) ?
                      fabs(spA->triMatElem[i]) :
                      fabs(spA->triMatElem[i]) * 2.0;
        }
    }
    
    return nrmA;
}

static void dataMatSparseNrm1(void *A, double *res){
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    double nrmA = 0.0;
    for ( int i = 0; i < spA->nTriMatElem; ++i ) {
        nrmA += ( spA->triMatCol[i] == spA->triMatRow[i] ) ?
                  fabs(spA->triMatElem[i]) :
                  fabs(spA->triMatElem[i]) * 2.0;
    }
    res[0] = nrmA;
}

static void dataMatSparseNrm2Square(void *A, double *res){
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    double nrmA = 0.0;
    for ( int i = 0; i < spA->nTriMatElem; ++i ) {
        nrmA += ( spA->triMatCol[i] == spA->triMatRow[i] ) ?
                  pow(spA->triMatElem[i], 2.0) :
                pow(spA->triMatElem[i], 2.0) * 2.0;
    }
    res[0] = nrmA;
}

static void dataMatSparseNrmInf(void *A, double *res){
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    double nrmA = 0.0;
    for ( int i = 0; i < spA->nTriMatElem; ++i ) {
        nrmA = ASDP_MAX(ASDP_ABS(spA->triMatElem[i]), nrmA);
    }
    res[0] = nrmA;
}

extern void dataMatSparseZeros(void *A){
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    ASDP_ZERO(spA->triMatElem, double, spA->nTriMatElem);
}

static void dataMatSparseStatNnz(int *nnzStat){
    nnzStat[0]++;
    return;
}

static void dataMatSparseAddPreprocessRankOneConeDetect(void *A, int *flagRankOne){
    return;
}

static void dataMatSparseScale(void *A, double scaleFactor){
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *)A;
    vvscl(&sparse->nTriMatElem, &scaleFactor, sparse->triMatElem);
}


extern void dataMatSparsefullDenseSDP(void *A){
    return;
}


/** @brief Calculate norm of dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[in] type Type of norm
 *  @return Norm of A
 */
static double dataMatNormDenseImpl( void *A, int type ) {
    
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    double nrmA = 0.0;
    int nCol = dsA->nSDPCol;
    int incx = 1;
    int colLen; ///< Exclude diagonal
    double colNrm; ///< Exclude diagonal
    double *p = dsA->dsMatElem;
    
    if ( type == FRO_NORM ) {
        for ( int i = 0; i < nCol; ++i ) {
            nrmA += p[0] * p[0];
            colLen = nCol - i - 1;
            colNrm = nrm2(&colLen, p + 1, &incx);
            nrmA += colNrm * colNrm * 2.0;
            p = p + colLen + 1;
        }
        nrmA = sqrt(nrmA);
    } else if ( type == ABS_NORM ) {
        for ( int i = 0; i < nCol; ++i ) {
            nrmA += fabs(p[0]);
            colLen = nCol - i - 1;
            nrmA += nrm1(&colLen, p + 1, &incx) * 2;
            p = p + colLen + 1;
        }
    }
    
    return nrmA;
}

static void dataMatDenseNrm1(void *A, double *res){
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    double nrmA = 0.0;
    int nCol = dsA->nSDPCol;
    int incx = 1;
    int colLen; ///< Exclude diagonal
    double colNrm; ///< Exclude diagonal
    double *p = dsA->dsMatElem;
    for ( int i = 0; i < nCol; ++i ) {
        nrmA += fabs(p[0]);
        colLen = nCol - i - 1;
        nrmA += nrm1(&colLen, p + 1, &incx) * 2;
        p = p + colLen + 1;
    }
    res[0] = nrmA;
}

static void dataMatDenseNrm2Square(void *A, double *res){
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    double nrmA = 0.0;
    int nCol = dsA->nSDPCol;
    int incx = 1;
    int colLen; ///< Exclude diagonal
    double colNrm; ///< Exclude diagonal
    double *p = dsA->dsMatElem;
    for ( int i = 0; i < nCol; ++i ) {
        nrmA += pow(p[0], 2);
        colLen = nCol - i - 1;
        nrmA += pow(nrm2(&colLen, p + 1, &incx), 2) * 2;
        p = p + colLen + 1;
    }
    res[0] = nrmA;
}


static void dataMatDenseNrmInf(void *A, double *res){
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    double nrmA = 0.0;
    int nCol = dsA->nSDPCol;
    int incx = 1;
    int colLen; ///< Exclude diagonal
    double colNrm; ///< Exclude diagonal
    double *p = dsA->dsMatElem;
    res[0] = 0;
    for (int i = 0; i < (nCol+1)*nCol /2; ++i ){
        res[0] = ASDP_MAX(res[0], ASDP_ABS(p[i]));
    }
    
    return;
}

extern void dataMatDenseZeros(void *A){
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    ASDP_ZERO(dsA->dsMatElem, double, dsA->nSDPCol * (dsA->nSDPCol + 1)/2);
}

static void dataMatDenseAddPreprocessRankOneConeDetect(void *A, int *flagRankOne){
    return;
}

static void dataMatDenseScale(void *A, double scaleFactor){
    // A = A * scaleFactor
    sdp_coeff_dense *dense = (sdp_coeff_dense *)A;
    int n = dense->nSDPCol * (dense->nSDPCol + 1)/2;
    vvscl(&n, &scaleFactor, dense->dsMatElem);
}

static void dataMatDenseStatNnz(int *nnzStat){
    nnzStat[0]++;
    return;
}


static asdp_retcode dataMatBuildUpEigZeroImpl( void *A, int *rank, double *auxiFullMat, double **eVals, double **eVecs ) {
    
    *rank = 0;
    
    return ASDP_RETCODE_OK;
}

static asdp_retcode dataMatBuildUpEigSparseImpl( void *A, int *rank, double *auxiFullMat, double **eVals, double **eVecs ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    if ( !eVals || !eVecs ) {
        retcode = ASDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    double sgn = 0.0;
    ASDP_ZERO(auxiFullMat, double, sparse->nSDPCol);
    int isRankOne = tsp_r1_extract(sparse->nSDPCol, sparse->nTriMatElem, sparse->triMatRow,
                                   sparse->triMatCol, sparse->triMatElem, &sgn, auxiFullMat);
    
    if ( !isRankOne ) {
        goto exit_cleanup;
    }
    
    *rank = 1;
    ASDP_INIT(*eVals, double, 1);
    ASDP_INIT(*eVecs, double, sparse->nSDPCol);
    ASDP_MEMCHECK(*eVecs);
    (*eVals)[0] = sgn;
    ASDP_MEMCPY(*eVecs, auxiFullMat, double, sparse->nSDPCol);
    
exit_cleanup:
    return retcode;
}

static asdp_retcode dataMatBuildUpEigDenseImpl( void *A, int *rank, double *auxiFullMat, double **eVals, double **eVecs ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    if ( !eVals || !eVecs ) {
        retcode = ASDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    
    double sgn = 0.0;
    int isRankOne = pds_r1_extract(dense->nSDPCol, dense->dsMatElem, &sgn, auxiFullMat);
    
    if ( !isRankOne ) {
        goto exit_cleanup;
    }
    
    *rank = 1;
    ASDP_INIT(*eVals, double, 1);
    ASDP_INIT(*eVecs, double, dense->nSDPCol);
    ASDP_MEMCHECK(*eVecs);
    (*eVals)[0] = sgn;
    ASDP_MEMCPY(*eVecs, auxiFullMat, double, dense->nSDPCol);
    
exit_cleanup:
    return retcode;
}



static int dataMatGetNnzZeroImpl( void *A ) {
    
    return 0;
}

/** @brief Calculate number of nonzero elements in sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @return Number of nonzero elements in A
 */
static int dataMatGetNnzSparseImpl( void *A ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    
    return spA->nTriMatElem;
}

/** @brief Calculate number of nonzero elements in dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @return Number of nonzero elements in A
 */
static int dataMatGetNnzDenseImpl( void *A ) {
    
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    
    return PACK_NNZ(dsA->nSDPCol);
}


static void dataMatDumpZeroImpl( void *A, double *v ) {
    
    sdp_coeff_zero *zeroA = (sdp_coeff_zero *) A;
    ASDP_ZERO(v, double, zeroA->nSDPCol * zeroA->nSDPCol);
    
    return;
}

/** @brief Construct full matrix from sparse matrix A
 *  @param[in] A sdpSparseData pointer
 *  @param[out] v Store full matrix
 */
static void dataMatDumpSparseImpl( void *A, double *v ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    ASDP_ZERO(v, double, spA->nSDPCol * spA->nSDPCol);
    int nCol = spA->nSDPCol;
    
    for ( int i = 0; i < spA->nTriMatElem; ++i ) {
        v[spA->triMatCol[i] * nCol + spA->triMatRow[i]] = spA->triMatElem[i];
        if ( spA->triMatCol[i] != spA->triMatRow[i] ) {
            v[spA->triMatRow[i] * nCol + spA->triMatCol[i]] = spA->triMatElem[i];
        }
    }
    
    return;
}

/** @brief Construct full matrix from dense matrix A
 *  @param[in] A sdpDenseData pointer
 *  @param[out] v Store full matrix
 */
static void dataMatDumpDenseImpl( void *A, double *v ) {
    
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    int nCol = dsA->nSDPCol;
    int dsIdx; ///< Index for (col, row) = (i, j) in dsMatElem
    
    for ( int i = 0, j; i < nCol; ++i) { // Column index
        dsIdx = i * (nCol * 2 - i + 1) / 2;
        v[i * nCol + i] = dsA->dsMatElem[dsIdx];
        for ( j = i + 1; j < nCol; ++j ) { // Row index
            dsIdx = i * (nCol * 2 - i - 1) / 2 + j;
            v[i * nCol + j] = dsA->dsMatElem[dsIdx];
            v[j * nCol + i] = v[i * nCol + j];
        }
    }
    
    return;
}


static void dataMatGetZeroSparsityImpl( void *A, int *spout ) {
    return;
}

static void dataMatGetSparseSparsityImpl( void *A, int *spout ) {
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    
    for ( int k = 0; k < spA->nTriMatElem; ++k ) {
        spout[PACK_IDX(spA->nSDPCol, spA->triMatRow[k], spA->triMatCol[k])] = 1;
    }

    return;
}

static void dataMatGetDenseSparsityImpl( void *A, int *spout ) {
    
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    
    for ( int i = 0; i < PACK_NNZ(dsA->nSDPCol); ++i ) {
        spout[i] = 1;
    }
    
    /* This method should not be invoked */
    assert( 0 );
    return;
}


static void dataMatAddZeroToBufferImpl( void *A, double a, int *spmap, double *buffer ) {
    /* Let buffer <- buffer + a * A */
    return;
}

static void dataMatAddSparseToBufferImpl( void *A, double a, int *spmat, double *buffer ) {
    /* If spmat exists, then the buffer uses sparse data structure,
       Otherwise (if spmat is NULL), buffer is a dense array
     */
    
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    
    if ( !spmat ) {
        /* Simply add it */
        for ( int i = 0; i < spA->nTriMatElem; ++i ) {
            FULL_ENTRY(buffer, spA->nSDPCol, spA->triMatRow[i], spA->triMatCol[i]) += \
            a * spA->triMatElem[i];
        }
    } else {
        /* Add according to the mapping */
        for ( int i = 0; i < spA->nTriMatElem; ++i ) {
            buffer[spmat[PACK_IDX(spA->nSDPCol, spA->triMatRow[i], spA->triMatCol[i])]] += \
            a * spA->triMatElem[i];
        }
    }
    
    return;
}

static void dataMatAddDenseToBufferImpl( void *A, double a, int *spmat, double *buffer ) {
    /* Whenever a dense matrix exists, the dual matrix must be dense */
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    
    if ( !spmat ) {
        /* We are adding a dense packed matrix to a dense full matrix */
        double *dPackElem = dsA->dsMatElem;
        double *dFullElem = buffer;
        for ( int i = 0; i < dsA->nSDPCol; ++i ) {
            for ( int j = 0; j < dsA->nSDPCol - i; ++j ) {
                dFullElem[j] += a * dPackElem[j];
            }
            dFullElem += dsA->nSDPCol + 1;
            dPackElem += (dsA->nSDPCol - i);
        }
        
    } else {
        /* This case should never happen */
        assert( 0 );
    }
    
    return;
}



static void dataMatClearZeroImpl( void *A ) {
    
    if ( !A ) {
        return;
    }
    
    ASDP_ZERO(A, sdp_coeff_zero, 1);
    
    return;
}

static void dataMatClearSparseImpl( void *A) {
    
    if ( !A ) {
        return;
    }
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    ASDP_FREE(sparse->triMatRow);
    ASDP_FREE(sparse->triMatCol);
    ASDP_FREE(sparse->triMatElem);

    ASDP_ZERO(A, sdp_coeff_sparse, 1);
    
    return;
}

static void dataMatClearDenseImpl( void *A ) {
    
    if ( !A ) {
        return;
    }
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    
    ASDP_FREE(dense->dsMatElem);
    ASDP_ZERO(A, sdp_coeff_dense, 1);
    
    return;
}

static void dataMatClearRankOneSparseImpl( void *A ) {
    
    if ( !A ) {
        return;
    }
    
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    
    ASDP_FREE(spr1->spR1MatElem);
    ASDP_FREE(spr1->spR1MatIdx);
    ASDP_FREE(spr1->spR1MatFactor);
    
    if (spr1->nnzVal != NULL){
        ASDP_FREE(spr1->nnzVal);
    }
    
    ASDP_ZERO(A, sdp_coeff_spr1, 1);
    
    return;
}

static void dataMatDestroyZeroImpl( void **pA ) {
    
    if ( !pA ) {
        return;
    }
    
    dataMatClearZeroImpl(*pA);
    sdp_coeff *p = (sdp_coeff *)*pA;
    ASDP_FREE(*pA);
    
    return;
}

static void dataMatDestroySparseImpl( void **pA ) {
    
    if ( !pA ) {
        return;
    }
    sdp_coeff *p = (sdp_coeff *)*pA;
    dataMatClearSparseImpl(*pA);
    ASDP_FREE(*pA);
    
    return;
}

extern void dataMatDestroyDenseImpl( void **pA ) {
    
    if ( !pA ) {
        return;
    }
    
    dataMatClearDenseImpl(*pA);
    sdp_coeff *p = (sdp_coeff *)*pA;
    ASDP_FREE(*pA);
    
    return;
}

static void dataMatViewZeroImpl( void *A ) {
    
    printf("Zero matrix of size %d L1 = [%5.3e] L2 = [%5.3e] \n",
           ((sdp_coeff_zero *) A)->nSDPCol, 0.0, 0.0);
    
    return;
}

static void dataMatViewSparseImpl( void *A ) {
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    printf("Sparse matrix of size %d and %d nnzs L1 = [%5.3e] L2 = [%5.3e]. \n",
           sparse->nSDPCol, sparse->nTriMatElem,
           dataMatNormSparseImpl(sparse, ABS_NORM),
           dataMatNormSparseImpl(A, FRO_NORM));
    return;
}



static void dataMatViewDenseImpl( void *A ) {
    
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    printf("Dense matrix of size %d L1 = [%5.3e] L2 = [%5.3e]. \n",
           dense->nSDPCol, dataMatNormDenseImpl(dense, ABS_NORM),
           dataMatNormDenseImpl(dense, FRO_NORM));
    
    return;
}

static int dataMatZeroIsEye( void *A, double *dEyeMultiple ) {
    
    return 0;
}

static int dataMatSparseIsEye( void *A, double *dEyeMultiple ) {
    
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    if ( sparse->nTriMatElem != sparse->nSDPCol ) {
        return 0;
    }
    
    double dDiagElem = sparse->triMatElem[0];
    
    for ( int iElem = 0; iElem < sparse->nTriMatElem; ++iElem ) {
        if ( sparse->triMatCol[iElem] != sparse->triMatRow[iElem] ) {
            return 0;
        }
        
        if ( sparse->triMatElem[iElem] != dDiagElem ) {
            return 0;
        }
    }
    
    if ( dEyeMultiple ) {
        *dEyeMultiple = dDiagElem;
    }
    
    return 1;
}

static int dataMatDenseIsEye( void *A, double *dEyeMultiple ) {
    
    return 0;
}

static int dataMatZeroIsUnitCol( void *A, int *iUnitCol ) {
    return 0;
}

static asdp_retcode dataMatZeroMultiRkMat(void *A, asdp_rk_mat_dense *X, double *AX){
    // AX = A * X = 0
    asdp_retcode retcode = ASDP_RETCODE_OK;
exit_cleanup:
    return retcode;
}

static asdp_retcode dataMatZeroMultiRkMatInnerProRkMat(void *A, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *res, void *UVt, sdp_coeff_type UVtType, int dummy2){
    // res = <AU, V> = 0
    res[0] = 0.0; // no change
    return ASDP_RETCODE_OK;
}

static void dataMatZeroAddDenseSDPCoeff(void *A, void *B, double weight, sdp_coeff_type Btype){
    return;
}

static void dataMatZeroNrm1(void *A, double *res){
    res[0] = 0.0;
    return;
}

static void dataMatZeroNrm2Square(void *A, double *res){
    res[0] = 0.0;
    return;
}

static void dataMatZeroNrmInf (void *A, double *res){
    res[0] = 0.0;
    return;
}

static void dataMatZeroStatNnz (int *nnzStat){
    return;
}

static void dataMatZeroAddPreprocessRankOneConeDetect(void *A, int *flagRankOne){
    return;
}

static void dataMatZeroScale(void *A, double scaleFactor){
    return;
}


static void dataMatZerofullDenseSDP(void *A){
    return;
}


extern asdp_retcode dataMatSparseMultiRkMat(void *A, asdp_rk_mat_dense *X, double *AX){
    asdp_retcode retcode = ASDP_RETCODE_OK;
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    
    ASDP_ZERO(AX, double, X->rank * X->nRows);
    double a = 1.0;
    tsp_ASemiLowbyMul(sparse->nSDPCol, sparse->nTriMatElem, a, sparse->triMatRow, sparse->triMatCol,
                      sparse->triMatElem, X->matElem, X->rank, AX);
    return retcode;
}



static void tsp_ASpSemiLowby(int n, int nnz, double a, int *idxInv, int *Ai, int *Aj, double *Ax, double *x, double *y){
    if (a == 0.0){
        return;
    }
    for (int i = 0; i < nnz; ++i){
        int row = Ai[i];
        int col = Aj[i];
        // find nnz index order
        int nnzIdx = idxInv[row];
        if (nnzIdx >= n){
            ASDP_ERROR_TRACE;
            asdp_printf("There exists a bug, the index should smaller than %d", n);
        }
        y[nnzIdx] += a * x[col] * Ax[i];
        if (row != col){
            nnzIdx = idxInv[col];
            y[nnzIdx] += a * x[row] * Ax[i];
        }
    }
}

static void tsp_AUV(int n, int nnzA, int *Ai, int *Aj, double *Ax, int nnzUVt, int *UVti, int *UVtj, double *UVtx, int **nnzIdx, double *res){
    // A is sparse, UVtx
    int incx = 1;
    int rowA = 0; int colA = 0;
    int rowUVt = 0; int colUVt = 0;
    double temp;
    if (nnzA == nnzUVt){
        res[0] += dot(&nnzA, Ax, &incx, UVtx, &incx);
        for (int i = 0; i < nnzA; ++i){
            rowA = Ai[i]; colA = Aj[i];
            rowUVt = UVti[i]; colUVt = UVtj[i];
            if (rowA == colA && rowA == rowUVt && rowUVt == colUVt){
                res[0] -= 0.5 * Ax[i] * UVtx[i];
            }
        }
    }else if (nnzA != nnzUVt){
        assert(nnzUVt > nnzA);
        int i = 0; int j = 0; int idx = 0;
        for (i = 0; i < nnzA; ++i){
            rowA = Ai[i];
            colA = Aj[i];
            idx = nnzIdx[rowA][colA];
            temp = Ax[i] * UVtx[idx];
            res[0] += temp;
            if (rowA == colA){
                res[0] -= 0.5 * temp;
            }
        }
    }
}

static void tsp_ADenseUV(int n, int nnzA, int *Ai, int *Aj, double *Ax, int nUVt, double *UVtx, int **nnzIdx, double *res){
    int rowA = 0; int colA = 0; int idx = 0;
    double temp;
    for (int i = 0; i < nnzA; ++i){
        rowA = Ai[i]; colA = Aj[i];
        idx = nnzIdx[rowA][colA];
        temp = Ax[i] * UVtx[idx];
        res[0] += temp;
        if (rowA == colA){
            res[0] -= 0.5 * temp;
        }
    }
}


extern void tsp_AUVObj(int n, int nnz, int *Ai, int *Aj, double *Ax, int rank, int nRows, double *Ux, double *Vx, double *res){
    /*
     Triplet matrix multiplication:
        res = <AU, V>
     where A is a sparse matrix in triplet format
     Ai is the row index of each non-zero element
     Aj is the column index of each non-zero element
     Ax is the value of each non-zero element
     n is the dimension of the matrix
     */
    int incxNew = 1;
    int row; int col;
    double *UVt;
    double temp = 0.0;
    ASDP_INIT(UVt, double, nnz);
    int nSquare = pow(n, 2);
    double sparseRatio = (double) nnz / (double) nSquare;
    if ( sparseRatio < 0.08){
        int incx = nRows;
        for (int i = 0; i < nnz; ++i){
            row = Ai[i]; col = Aj[i];
            UVt[i] = dot(&rank, &Ux[row], &incx, &Vx[col], &incx);
            UVt[i] += dot(&rank, &Ux[col], &incx, &Vx[row], &incx);
        }
    }else{
        double *UVtFull;
        ASDP_INIT(UVtFull, double, nSquare);
        fds_syr2k('L', 'N', nRows, rank, 1.0, Ux, Vx, 0.0, UVtFull);
        for (int i = 0; i < nnz; ++i){
            row = Ai[i]; col = Aj[i];
            UVt[i] = UVtFull[col * n + row];
        }
        ASDP_FREE(UVtFull);
    }
    temp += dot(&nnz, Ax, &incxNew, UVt, &incxNew);
    for (int i = 0; i < nnz; ++i){
        if (Ai[i] == Aj[i]){
            temp -= 0.5 * Ax[i] * UVt[i];
        }
    }
    res[0] += temp;
    ASDP_FREE(UVt);
}


static asdp_retcode dataMatSparseMultiRkMatInnerProRkMat(void *A, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *res, void *UVtIn, sdp_coeff_type UVtType, int FLAG_UV){
    /*
      res = <AU, V>
      UVt is (UVt + VUt) in fact
     */
    asdp_retcode retcode = ASDP_RETCODE_OK;
    res[0] = 0.0;

    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    if(FLAG_UV == FLAG_INI || FLAG_UV == FLAG_UVt){
        // Note: A is lower triangular matrix
        if (UVtType == SDP_COEFF_SPARSE){
            sdp_coeff_sparse *UVt = (sdp_coeff_sparse *)UVtIn;
            tsp_AUV(sparse->nSDPCol, sparse->nTriMatElem, sparse->triMatRow, sparse->triMatCol, sparse->triMatElem, UVt->nTriMatElem, UVt->triMatRow, UVt->triMatCol, UVt->triMatElem, UVt->rowCol2NnzIdx, res);
        }else{
            sdp_coeff_dense *UVt = (sdp_coeff_dense *)UVtIn;
             tsp_ADenseUV(sparse->nSDPCol, sparse->nTriMatElem, sparse->triMatRow, sparse->triMatCol, sparse->triMatElem, UVt->nSDPCol, UVt->dsMatElem, UVt->rowCol2NnzIdx, res);
        }
    }else if (FLAG_UV == FLAG_OBJ){
        tsp_AUVObj(sparse->nSDPCol, sparse->nTriMatElem, sparse->triMatRow, sparse->triMatCol, sparse->triMatElem, U->rank, U->nRows, U->matElem, V->matElem, res);
    }
exit_cleanup:
    return  retcode;
}


static void dataMatSparseAddDenseSDPCoeff(void *A, void *B, double weight){
    // B = B + A * weight
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    sdp_coeff_dense *dense = (sdp_coeff_dense *) B;
    int n = dense->nSDPCol;
    int row, col;
    if (dense->rowCol2NnzIdx == NULL){
        for (int i = 0; i < sparse->nTriMatElem; ++i){
            int row = sparse->triMatRow[i];
            int col = sparse->triMatCol[i];
            if (row >= col){
                dense->dsMatElem[(n * col -col*(col+1)/2 + row)] += weight * sparse->triMatElem[i];
            }
        }
    }else{
        int idx; int row; int col;
        for (int i = 0; i < sparse->nTriMatElem; ++i){
            row = sparse->triMatRow[i];
            col = sparse->triMatCol[i];
            assert(row >= col);
            idx = dense->rowCol2NnzIdx[row][col];
            assert(idx != -1);
            dense->dsMatElem[idx] += weight * sparse->triMatElem[i];
        }
    }

}



static void dataMatSparseAddSparseSDPCoeff(void *A, void *B, double weight){
    // B = B + weight * A
    sdp_coeff_sparse *sparseA = (sdp_coeff_sparse *) A;
    sdp_coeff_sparse *sparseB = (sdp_coeff_sparse *) B;
    int nnzA = sparseA->nTriMatElem;
    int nnzB = sparseB->nTriMatElem;
    assert(nnzB >= nnzA);
    
    int rowA, colA; int j;
    for (int i = 0; i < nnzA; ++i){
        rowA = sparseA->triMatRow[i];
        colA = sparseA->triMatCol[i];
        j = sparseB->rowCol2NnzIdx[rowA][colA];
        sparseB->triMatElem[j] += weight * sparseA->triMatElem[i];
    }
}

static void dataMatSparseAddSDPCoeff(void *A, void *B, double weight, sdp_coeff_type B_type){
    if (B_type == SDP_COEFF_DENSE){
        dataMatSparseAddDenseSDPCoeff(A, B, weight);
    }else if (B_type == SDP_COEFF_SPARSE){
        dataMatSparseAddSparseSDPCoeff(A, B, weight);
    }else{
        ASDP_ERROR_TRACE;
    }
}



static int dataMatDenseIsUnitCol( void *A, int *iUnitCol ) {
    return 0;
}

extern asdp_retcode dataMatDenseMultiRkMat(void *A, asdp_rk_mat_dense *X, double *AX){
    /*
     AX = A * X, A is dense meanse result is dense
     */
    asdp_retcode retcode = ASDP_RETCODE_OK;
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    double alpha = 1.0;
    double beta = 0.0;
    // construct full matrix
    double *fullMat;
    ASDP_INIT(fullMat, double, dense->nSDPCol * dense->nSDPCol);
    int idx = 0;
    for (int col = 0; col < dense->nSDPCol; ++col){
        for (int row = col; row < dense->nSDPCol; ++row){
            fullMat[dense->nSDPCol * col + row] = dense->dsMatElem[idx];
            fullMat[dense->nSDPCol * row + col] = dense->dsMatElem[idx];
            idx++;
        }
    }
    char side = 'L'; // C:= alpha * A * B + beta * C;
    char uplo = 'L';
    int m = dense->nSDPCol;
    int n = X->rank;
#ifdef UNDER_BLAS
    dsymm_(&side, &uplo, &m, &n, &alpha, fullMat, &m, X->matElem, &m, &beta, AX, &m );
#else
    dsymm(&side, &uplo, &m, &n, &alpha, fullMat, &m, X->matElem, &m, &beta, AX, &m );
#endif
    ASDP_FREE(fullMat);
exit_cleanup:
    return retcode;
}

static asdp_retcode dataMatDenseMultiRkMatInnerProRkMat(void *A, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *res, void *UVtIn, sdp_coeff_type UVtType, int FLAG_UV){
    /*
     res = <AU, V>
     dense A, dense U, dense V
     */
    res[0] = 0.0;
    asdp_retcode retcode = ASDP_RETCODE_OK;
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    
    double alpha = 1.0;
    double beta = 0.0;
    int n = 0;
    int incx = 1;
    if(FLAG_UV == FLAG_INI || FLAG_UV == FLAG_UVt){
        assert(UVtType == SDP_COEFF_DENSE);
        sdp_coeff_dense *UVt = (sdp_coeff_dense *)UVtIn;
        n = UVt->nSDPCol * (UVt->nSDPCol + 1) / 2;
        res[0] += dot(&n, dense->dsMatElem, &incx, UVt->dsMatElem, &incx);
        int idx = 0; int colNum = UVt->nSDPCol;
        for (int i = 0; i < UVt->nSDPCol; ++i){
            res[0] -= 0.5 * dense->dsMatElem[idx] * UVt->dsMatElem[idx];
            idx += colNum;
            colNum -= 1;
        }
    }else if (FLAG_UV == FLAG_OBJ){
        double *UVt;
        ASDP_INIT(UVt, double, dense->nSDPCol * dense->nSDPCol);
        fds_syr2k(ACharConstantUploLow, 'N', U->nRows, U->rank, 1.0, U->matElem, V->matElem, 1.0, UVt);
        double *UVtlow;
        n = dense->nSDPCol * (dense->nSDPCol + 1)/2;
        ASDP_INIT(UVtlow, double, n);
        int idx = 0; int row = 0;
        for (int col = 0; col < dense->nSDPCol; ++col){
            ASDP_MEMCPY(&UVtlow[idx], &UVt[dense->nSDPCol * col + row], double, dense->nSDPCol - col);
            row++;
            idx += (dense->nSDPCol - col);
        }
        res[0] += dot(&n, dense->dsMatElem, &incx, UVtlow, &incx);
        idx = 0; int colNum = dense->nSDPCol;
        for (int i = 0; i < dense->nSDPCol; ++i){
            res[0] -= 0.5 * dense->dsMatElem[idx] * UVt[idx];
            idx += colNum;
            colNum -= 1;
        }
        ASDP_FREE(UVt);
        ASDP_FREE(UVtlow);
    }

exit_cleanup:
    return retcode;
}

static void dataMatDenseAddDenseSDPCoeff(void *A, void *B, double weight){
    // B += A * weight
    sdp_coeff_dense *denseA = (sdp_coeff_dense *) A;
    sdp_coeff_dense *denseB = (sdp_coeff_dense *) B;
    int n = denseA->nSDPCol * (denseA->nSDPCol + 1) / 2;
    int incx = 1;
    axpy(&n, &weight, denseA->dsMatElem, &incx, denseB->dsMatElem, &incx);
        
}

static void dataMatDenseAddSDPCoeff(void *A, void *B, double weight, sdp_coeff_type B_type){
    if (B_type == SDP_COEFF_DENSE){
        dataMatDenseAddDenseSDPCoeff(A, B, weight);
    }else{
        ASDP_ERROR_TRACE;
    }
}

extern asdp_retcode dataMatDenseAddSparseSDPCoeff(void *A, void *B, double weight){
    // B += A * weight
    asdp_retcode retcode = ASDP_RETCODE_FAILED;
    asdp_printf("The function is not used");
exit_cleanup:
    return retcode;
        
}

static int dataMatRankOneSparseIsUnitCol( void *A, int *iUnitCol ) {
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    if ( spr1->nSpR1FactorElem == 1 && spr1->spR1MatElem[0] == 1.0 && spr1->spR1FactorSign == 1.0 ) {
        if ( iUnitCol ) {
            *iUnitCol = spr1->spR1MatIdx[0];
        }
        return 1;
    }
    return 0;
}


static asdp_retcode dataMatRankOneSparseCopyTo(sdp_coeff *dst, sdp_coeff *src, sdp_coeff_type srcType){
    asdp_retcode retcode = ASDP_RETCODE_OK;
    if (srcType == SDP_COEFF_DENSE){
        sdp_coeff_dense *temp;
        ASDP_INIT(temp, sdp_coeff_dense, 1);
        sdp_coeff_spr1 *srcData = (sdp_coeff_spr1 *)src->dataMat;
        int n = srcData->nSDPCol;
        ASDP_INIT(temp->dsMatElem, double, (n+1) * n /2);
        ASDP_MEMCHECK(temp->dsMatElem);
        for(int i = 0; i < srcData->nSpR1FactorElem; ++i){
            int row = srcData->spR1MatIdx[i];
            for (int j = 0; j <= i; ++j){
                int col = srcData->spR1MatIdx[j];
                temp->dsMatElem[(n * col - col * (col + 1)/2 + row)] = srcData->spR1FactorSign * srcData->spR1MatElem[i] * srcData->spR1MatElem[j];
            }
        }
        dst->dataMat = (void *)temp;
    }else{
        asdp_printf("Copy rank one sparse sdp_coeff to other type sdp_coeff but dense is developing!");
    }
    
exit_cleanup:
    return retcode;
}

static void dataMatRankOneSparseAddDenseSDPCoeff(void *A, void *B, double weight){
    // B += weight * A
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    sdp_coeff_dense *dense = (sdp_coeff_dense *) B;
    int n = spr1->nSDPCol;
    int row, col;
    if (dense->rowCol2NnzIdx == NULL){
        for (int i = 0; i<spr1->nSpR1FactorElem; ++i){
            for (int j = 0; j <= i; ++j){
                row = spr1->spR1MatIdx[i];
                col = spr1->spR1MatIdx[j];
                dense->dsMatElem[(n * col -col*(col+1)/2 + row)] += weight * spr1->spR1FactorSign * spr1->spR1MatElem[i] * spr1->spR1MatElem[j];
            }
        }
    }else{
        int idx = 0; int count = 0;
        for (int j = 0; j < spr1->nSpR1FactorElem; ++j){
            for (int i = j; i < spr1->nSpR1FactorElem; ++i){
                row = spr1->spR1MatIdx[i];
                col = spr1->spR1MatIdx[j];
                idx = dense->rowCol2NnzIdx[row][col];
                dense->dsMatElem[idx] += weight * spr1->nnzVal[count];
                count++;
            }
        }
    }

}


static void dataMatRankOneSparseAddSparseSDPCoeff(void *A, void *B, double weight){
    // B += weight * A
    sdp_coeff_spr1 *spr1 = (sdp_coeff_spr1 *) A;
    sdp_coeff_sparse *sparseB= (sdp_coeff_sparse *) B;
    int row, col, idx; int count = 0;
    for (int j = 0; j < spr1->nSpR1FactorElem; ++j){
        for (int i = j; i < spr1->nSpR1FactorElem; ++i){
            row = spr1->spR1MatIdx[i];
            col = spr1->spR1MatIdx[j];
            idx = sparseB->rowCol2NnzIdx[row][col];
            sparseB->triMatElem[idx] += weight * spr1->nnzVal[count];
            count++;
        }
    }
}

extern asdp_retcode dataMatRankOneDenseAddSparseSDPCoeff(void *A, void *B, double weight){
    // not used
    asdp_printf("It is impossible add rank one dense matrix to sparse matrix");
    ASDP_ERROR_TRACE;
    // B += weight * A
    asdp_retcode retcode = ASDP_RETCODE_OK;
    sdp_coeff_dsr1 *dr1 = (sdp_coeff_dsr1 *) A;
    sdp_coeff_sparse *sparseB = (sdp_coeff_sparse *) B;
    if (sparseB->nTriMatElem < 3){
        int flag = 0;
        for (int row = 0; row < dr1->nSDPCol; ++row){
            for (int col = 0; col <= row; ++col){
                flag = 0;
                for (int iElem = 0; iElem < sparseB->nTriMatElem; ++iElem){
                    int rowB = sparseB->triMatRow[iElem];
                    int colB = sparseB->triMatCol[iElem];
                    if (row == rowB && col == colB){
                        sparseB->triMatElem[iElem] += weight * dr1->r1FactorSign * dr1->r1MatFactor[row] * dr1->r1MatFactor[col];
                        flag = 1;
                    }
                }
                if (flag == 0){
                    asdp_printf("No find proper position in sparseB at function dataMatRankOneDenseAddSparseSDPCoeff");
                    retcode = ASDP_RETCODE_FAILED;
                    goto exit_cleanup;
                }
            }
        }
    }else{
        int iA = 0;
        int jA = 0;
        int j = 0;
        while ( iA >= jA && j < sparseB->nTriMatElem && iA < dr1->nSDPCol){
            int rowA = iA;
            int colA = jA;
            int rowB = sparseB->triMatRow[j];
            int colB = sparseB->triMatCol[j];
            if (colA < colB||(colA == colB && rowA < rowB)){
                // Move the A pointer
                iA++;
                if (iA == dr1->nSDPCol){
                    jA++;
                    iA = jA;
                }
            }else if (colA == colB && rowA == rowB){
                sparseB->triMatElem[j] += weight * dr1->r1FactorSign * dr1->r1MatFactor[rowA] * dr1->r1MatFactor[colA];
                iA++;
                if (iA == dr1->nSDPCol){
                    jA++;
                    iA = jA;
                }
                j++;
            }else{
                j++;
            }
        }
    }
exit_cleanup:
    return retcode;
}



static void sdpDataMatIChooseType( sdp_coeff *sdpCoeff, sdp_coeff_type dataType ) {
    
    sdpCoeff->dataType = dataType;
     
    switch ( dataType ) {
        case SDP_COEFF_ZERO:
            sdpCoeff->create = dataMatCreateZeroImpl;
            sdpCoeff->scal = dataMatScalZeroImpl;
            sdpCoeff->norm = dataMatNormZeroImpl;
            sdpCoeff->eig = dataMatBuildUpEigZeroImpl;
            sdpCoeff->getnnz = dataMatGetNnzZeroImpl;
            sdpCoeff->dump = dataMatDumpZeroImpl;
            sdpCoeff->getmatnz = dataMatGetZeroSparsityImpl;
            sdpCoeff->add2buffer = dataMatAddZeroToBufferImpl;
            sdpCoeff->destroy = dataMatDestroyZeroImpl;
            sdpCoeff->view = dataMatViewZeroImpl;
            sdpCoeff->iseye = dataMatZeroIsEye;
            sdpCoeff->isunitcol = dataMatZeroIsUnitCol;
            sdpCoeff->mul_rk = dataMatZeroMultiRkMat;
            sdpCoeff->mul_inner_rk_double = dataMatZeroMultiRkMatInnerProRkMat;
            sdpCoeff->add_sdp_coeff = dataMatZeroAddDenseSDPCoeff;
            sdpCoeff->nrm1 = dataMatZeroNrm1;
            sdpCoeff->nrm2Square = dataMatZeroNrm2Square;
            sdpCoeff->nrmInf = dataMatZeroNrmInf;
            sdpCoeff->zeros = NULL;
            sdpCoeff->statNnz = dataMatZeroStatNnz;
            sdpCoeff->addPreprocessRankOneConeDetect = dataMatZeroAddPreprocessRankOneConeDetect;
            sdpCoeff->scaleData = dataMatZeroScale;
            break;
        case SDP_COEFF_SPARSE:
            sdpCoeff->create = dataMatCreateSparseImpl;
            sdpCoeff->scal = dataMatScalSparseImpl;
            sdpCoeff->norm = dataMatNormSparseImpl;
            sdpCoeff->eig = dataMatBuildUpEigSparseImpl;
            sdpCoeff->getnnz = dataMatGetNnzSparseImpl;
            sdpCoeff->dump = dataMatDumpSparseImpl;
            sdpCoeff->getmatnz = dataMatGetSparseSparsityImpl;
            sdpCoeff->add2buffer = dataMatAddSparseToBufferImpl;
            sdpCoeff->destroy = dataMatDestroySparseImpl;
            sdpCoeff->view = dataMatViewSparseImpl;
            sdpCoeff->mul_rk = dataMatSparseMultiRkMat;
            sdpCoeff->mul_inner_rk_double = dataMatSparseMultiRkMatInnerProRkMat;
            sdpCoeff->add_sdp_coeff = dataMatSparseAddSDPCoeff;
            sdpCoeff->nrm1 = dataMatSparseNrm1;
            sdpCoeff->nrm2Square = dataMatSparseNrm2Square;
            sdpCoeff->nrmInf = dataMatSparseNrmInf;
            sdpCoeff->zeros = dataMatSparseZeros;
            sdpCoeff->statNnz = dataMatSparseStatNnz;
            sdpCoeff->addPreprocessRankOneConeDetect = dataMatSparseAddPreprocessRankOneConeDetect;
            sdpCoeff->scaleData = dataMatSparseScale;
            break;
        case SDP_COEFF_DENSE:
            sdpCoeff->create = dataMatCreateDenseImpl;
            sdpCoeff->scal = dataMatScalDenseImpl;
            sdpCoeff->norm = dataMatNormDenseImpl;
            sdpCoeff->eig = dataMatBuildUpEigDenseImpl;
            sdpCoeff->getnnz = dataMatGetNnzDenseImpl;
            sdpCoeff->dump = dataMatDumpDenseImpl;
            sdpCoeff->getmatnz = dataMatGetDenseSparsityImpl;
            sdpCoeff->add2buffer = dataMatAddDenseToBufferImpl;
            sdpCoeff->destroy = dataMatDestroyDenseImpl;
            sdpCoeff->view = dataMatViewDenseImpl;
            sdpCoeff->iseye = dataMatDenseIsEye;
            sdpCoeff->isunitcol = dataMatDenseIsUnitCol;
            sdpCoeff->mul_rk = dataMatDenseMultiRkMat;
            sdpCoeff->mul_inner_rk_double = dataMatDenseMultiRkMatInnerProRkMat;
            sdpCoeff->add_sdp_coeff = dataMatDenseAddSDPCoeff;
            sdpCoeff->nrm1 = dataMatDenseNrm1;
            sdpCoeff->nrm2Square = dataMatDenseNrm2Square;
            sdpCoeff->nrmInf = dataMatDenseNrmInf;
            sdpCoeff->zeros = dataMatDenseZeros;
            sdpCoeff->statNnz = dataMatDenseStatNnz;
            sdpCoeff->addPreprocessRankOneConeDetect = dataMatDenseAddPreprocessRankOneConeDetect;
            sdpCoeff->scaleData = dataMatDenseScale;
            break;
        case SDP_COEFF_SPR1:{
            printf("Not used!\n");
            break;
        }
        case SDP_COEFF_DSR1:{
            printf("Not used!\n");
            break;
        }
           
        default:
            assert( 0 );
            break;
    }
    
    return;
}

/* External methods for the SDP data */
extern asdp_retcode sdpDataMatCreate( sdp_coeff **psdpCoeff ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    if ( !psdpCoeff ) {
        retcode = ASDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    sdp_coeff *sdpCoeff;
    ASDP_INIT(sdpCoeff, sdp_coeff, 1);
    ASDP_MEMCHECK(sdpCoeff);
    
    /* -1 means there is no eigen-decomposition available */
    sdpCoeff->eigRank = -1;
    *psdpCoeff = sdpCoeff;
    
exit_cleanup:
    return retcode;
}

/** @brief Set SDP data matrix. This routine selects data type and assiciate structure with their operations
 * 
 *
 */
extern asdp_retcode sdpDataMatSetData( sdp_coeff *sdpCoeff, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem ) {
    // sdpCoeff  : data mat struct
    // nSDPCol   : data mat dimension
    // dataMatNnz: data mat non-zeros number
    asdp_retcode retcode = ASDP_RETCODE_OK;
    sdpCoeff->nSDPCol = nSDPCol;
    
    /* Choose data matrix type */
    int nPack = PACK_NNZ(nSDPCol);
    
    /* At this stage, only sparse, zero and dense matrices are classified */
    if ( dataMatNnz == 0 ) {
        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_ZERO);
    } else if ( dataMatNnz > 0.3 * nPack ) {
        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_DENSE);
    } else {
        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_SPARSE);
    }
    
    /* Create data */
    ASDP_CALL(sdpCoeff->create(&sdpCoeff->dataMat, nSDPCol,
                                dataMatNnz, dataMatIdx, dataMatElem));
    
exit_cleanup:
    
    return retcode;
}

extern int sdpDataMatGetRank( sdp_coeff *sdpCoeff ) {
    
    if ( sdpCoeff->dataType == SDP_COEFF_ZERO ) {
        return 0;
    } else if ( sdpCoeff->dataType == SDP_COEFF_DSR1 || sdpCoeff->dataType == SDP_COEFF_SPR1 ) {
        return 1;
    } else if ( sdpCoeff->eigRank != -1 ) {
        return sdpCoeff->eigRank;
    }
    
    return sdpCoeff->nSDPCol;
}

extern void sdpDataMatScal( sdp_coeff *sdpCoeff, double scal ) {
    
    sdpCoeff->scal(sdpCoeff->dataMat, scal);
    
    return;
}

extern double sdpDataMatNorm( sdp_coeff *sdpCoeff, int type ) {
    
    return sdpCoeff->norm(sdpCoeff->dataMat, type);
}

/* Used in presolving and low-rank (rank one) structure detection */
//extern asdp_retcode sdpDataMatBuildUpEigs( sdp_coeff *sdpCoeff, double *dAuxFullMatrix ) {
//
//    asdp_retcode retcode = ASDP_RETCODE_OK;
//
//    ASDP_ZERO(dAuxFullMatrix, double, sdpCoeff->nSDPCol);
//    ASDP_CALL(sdpCoeff->eig(sdpCoeff->dataMat, &sdpCoeff->eigRank, dAuxFullMatrix,
//                             &sdpCoeff->eigVals, &sdpCoeff->eigVecs));
//
//    if ( sdpCoeff->eigRank != 1 ) {
//        goto exit_cleanup;
//    }
//
//    /* The matrix is rank-one */
//    /* Count nnz in the decomposition */
//    int nRankOneNz = 0;
//    for ( int iRow = 0; iRow < sdpCoeff->nSDPCol; ++iRow ) {
//        if ( fabs(sdpCoeff->eigVecs[iRow]) > 1e-10 ) {
//            nRankOneNz += 1;
//        }
//    }
//
//    /* We no longer need it */
//    sdpCoeff->destroy(&sdpCoeff->dataMat);
//
//    int useDense = 0;
//    if ( nRankOneNz > 0.5 * sdpCoeff->nSDPCol ) {
//        useDense = 1;
//    }
//
//    if ( useDense ) {
//
//        sdp_coeff_dsr1 *dsr1;
//
//        ASDP_INIT(dsr1, sdp_coeff_dsr1, 1);
//        ASDP_MEMCHECK(dsr1);
//
//        dsr1->nSDPCol = sdpCoeff->nSDPCol;
//        dsr1->r1FactorSign = sdpCoeff->eigVals[0];
//
//        ASDP_INIT(dsr1->r1MatFactor, double, sdpCoeff->nSDPCol);
//        ASDP_MEMCHECK(dsr1->r1MatFactor);
//        ASDP_MEMCPY(dsr1->r1MatFactor, sdpCoeff->eigVecs, double, sdpCoeff->nSDPCol);
//
//        dataMatNormalizeRankOneDenseImpl(dsr1);
//        sdpCoeff->dataMat = dsr1;
//        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_DSR1);
//
//    } else {
//
//        sdp_coeff_spr1 *spr1;
//
//        ASDP_INIT(spr1, sdp_coeff_spr1, 1);
//        ASDP_MEMCHECK(spr1);
//
//        spr1->nSDPCol = sdpCoeff->nSDPCol;
//        spr1->spR1FactorSign = sdpCoeff->eigVals[0];
//        spr1->nSpR1FactorElem = nRankOneNz;
//
//        ASDP_INIT(spr1->spR1MatIdx, int, nRankOneNz);
//        ASDP_INIT(spr1->spR1MatElem, double, nRankOneNz);
//        ASDP_INIT(spr1->spR1MatFactor, double, sdpCoeff->nSDPCol);
//
//        int iNz = 0;
//        for ( int iRow = 0; iRow < sdpCoeff->nSDPCol; ++iRow ) {
//            if ( fabs(sdpCoeff->eigVecs[iRow] ) > 1e-10 ) {
//                spr1->spR1MatIdx[iNz] = iRow;
//                spr1->spR1MatElem[iNz] = sdpCoeff->eigVecs[iRow];
//                spr1->spR1MatFactor[iRow] = sdpCoeff->eigVecs[iRow];
//                iNz += 1;
//            }
//        }
//
//        assert( iNz == nRankOneNz );
//        dataMatNormalizeRankOneSparseImpl(spr1);
//        sdpCoeff->dataMat = spr1;
//        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_SPR1);
//    }
//
//    /* Reset rank */
//    sdpCoeff->eigRank = -1;
//    ASDP_FREE(sdpCoeff->eigVals);
//    ASDP_FREE(sdpCoeff->eigVecs);
//
//exit_cleanup:
//    return retcode;
//}

extern inline int sdpDataMatGetNnz( sdp_coeff *sdpCoeff ) {
    
    return sdpCoeff->getnnz(sdpCoeff->dataMat);
}

extern inline void sdpDataMatDump( sdp_coeff *sdpCoeff, double *dFullMatrix ) {
    
    return sdpCoeff->dump(sdpCoeff->dataMat, dFullMatrix);
}

extern void sdpDataMatGetMatNz( sdp_coeff *sdpCoeff, int *iMatSpsPattern ) {
    
    sdpCoeff->getmatnz(sdpCoeff->dataMat, iMatSpsPattern);
    
    return;
}

extern void sdpDataMatAddToBuffer( sdp_coeff *sdpCoeff, double dElem, int *iMatSpsPattern, double *dBuffer ) {
    
    if ( dElem ) {
        sdpCoeff->add2buffer(sdpCoeff->dataMat, dElem, iMatSpsPattern, dBuffer);
    }
    
    return;
}

extern inline sdp_coeff_type sdpDataMatGetType( sdp_coeff *sdpCoeff ) {
    
    return sdpCoeff->dataType;
}

extern inline void sdpDataMatClear( sdp_coeff *sdpCoeff ) {
    
    if ( !sdpCoeff ) {
        return;
    }
    
    sdpCoeff->destroy(&sdpCoeff->dataMat);
    
    if ( sdpCoeff->eigRank != -1 ) {
        ASDP_FREE(sdpCoeff->eigVals);
        ASDP_FREE(sdpCoeff->eigVecs);
    }
    
    ASDP_ZERO(sdpCoeff, sdp_coeff, 1);
    
    return;
}

extern inline void sdpDataMatDestroy( sdp_coeff **psdpCoeff ) {
    
    if ( !psdpCoeff ) {
        return;
    }
    
    sdpDataMatClear(*psdpCoeff);
    ASDP_FREE(*psdpCoeff);
    
    return;
}

extern inline void sdpDataMatView( sdp_coeff *sdpCoeff ) {
    
    sdpCoeff->view(sdpCoeff->dataMat);
    
    return;
}

extern int sdpDataMatIsEye( sdp_coeff *sdpCoeff, double *dEyeMultiple ) {
    
    return sdpCoeff->iseye(sdpCoeff->dataMat, dEyeMultiple);
}

extern int sdpDataMatIsUnitCol( sdp_coeff *sdpCoeff, int *iUnitCol ) {
    
    return sdpCoeff->isunitcol(sdpCoeff->dataMat, iUnitCol) ;
}

static asdp_retcode copyDenseToDense(sdp_coeff *dst, sdp_coeff *src){
    asdp_retcode retcode = ASDP_RETCODE_OK;
    sdp_coeff_dense *temp;
    ASDP_INIT(temp, sdp_coeff_dense, 1);
    ASDP_MEMCHECK(temp);
    sdp_coeff_dense *srcData = (sdp_coeff_dense *)src->dataMat;
    temp->nSDPCol = srcData->nSDPCol;
    int n = srcData->nSDPCol;
    ASDP_INIT(temp->dsMatElem, double, (n+1)*n/2);
    ASDP_MEMCPY(temp->dsMatElem, srcData->dsMatElem, double, (n+1)*n/2);
    dst->dataMat = (void *)temp;
exit_cleanup:
    return retcode;
}

static asdp_retcode copyZeroToDense(sdp_coeff *dst, sdp_coeff *src){
    asdp_retcode retcode = ASDP_RETCODE_OK;
    sdp_coeff_dense *temp;
    ASDP_INIT(temp, sdp_coeff_dense, 1);
    temp->nSDPCol = src->nSDPCol;
    int n = src->nSDPCol;
    ASDP_INIT(temp->dsMatElem, double, (n+1)*n/2);
    ASDP_MEMCHECK(temp->dsMatElem);
    ASDP_ZERO(temp->dsMatElem, double, (n+1)*n/2);
    dst->dataMat = (void *)temp;
exit_cleanup:
    return retcode;
}

static asdp_retcode copySparseToDense(sdp_coeff *dst, sdp_coeff *src){
    asdp_retcode retcode = ASDP_RETCODE_OK;
    sdp_coeff_dense *temp;
    ASDP_INIT(temp, sdp_coeff_dense, 1);
    sdp_coeff_sparse *srcData = (sdp_coeff_sparse *)src->dataMat;
    temp->nSDPCol = srcData->nSDPCol;
    int n = srcData->nSDPCol;
    ASDP_INIT(temp->dsMatElem, double, (n+1)*n/2);
    ASDP_MEMCHECK(temp->dsMatElem);
    ASDP_ZERO(temp->dsMatElem, double, (n+1)*n/2);
    for (int i = 0; i<srcData->nTriMatElem; ++i){
        int row = srcData->triMatRow[i];
        int col = srcData->triMatCol[i];
        temp->dsMatElem[(n * col - col * (col + 1) / 2 + row)] = srcData->triMatElem[i];
    }
    dst->dataMat = (void *)temp;
exit_cleanup:
    return retcode;
}



#include "lorads_utils.h"
#include "def_lorads_sdp_data.h"
#include "lorads_sdp_data.h"
#include "def_lorads_elements.h"
#include "lorads_vec_opts.h"
#include "lorads_sparse_opts.h"
#include "lorads_dense_opts.h"


extern void statNnz(double *vec, lorads_int n){
    lorads_int nnz = 0;
    for (lorads_int i = 0; i < n; ++i){
        if (fabs(vec[i]) > 1e-10){
            nnz += 1;
        }
    }
#ifdef INT32
    printf("nnz = %d\n", nnz);
#endif
#ifdef UNIX_INT64
    printf("nnz = %ld\n", nnz);
#endif
#ifdef MAC_INT64
    printf("nnz = %lld\n", nnz);
#endif
}


extern lorads_int hash_function(lorads_int row, lorads_int col, lorads_int size) {
    return (row + col) % size;
}

lorads_int find_index(Dict *dict, lorads_int row, lorads_int col) {
    lorads_int hash_index = hash_function(row, col, dict->size);
    DictNode *node = dict->table[hash_index];
    while (node != NULL) {
        if (node->element.row == row && node->element.col == col) {
            return node->element.index;
        }
        node = node->next;
    }
    return -1;
}

extern void sdpDataMatClear( sdp_coeff *sdpCoeff ) {
    if ( !sdpCoeff ) {
        return;
    }
    sdpCoeff->destroy(&sdpCoeff->dataMat);
    LORADS_ZERO(sdpCoeff, sdp_coeff, 1);
}

extern void sdpDataMatDestroy( sdp_coeff **psdpCoeff ) {
    if ( !psdpCoeff ) {
        return;
    }
    sdpDataMatClear(*psdpCoeff);
    LORADS_FREE(*psdpCoeff);
}

static void dataMatCreateZeroImpl( void **pA, lorads_int nSDPCol, lorads_int dataMatNnz, lorads_int *dataMatIdx, double *dataMatElem ) {
    (void) dataMatNnz;
    (void) dataMatIdx;
    (void) dataMatElem;
    if ( !pA ) {
        LORADS_ERROR_TRACE;
    }
    sdp_coeff_zero *zero;
    LORADS_INIT(zero, sdp_coeff_zero, 1);
    LORADS_MEMCHECK(zero);
    if (zero){
        zero->nSDPCol = nSDPCol;
    }else{
        LORADS_ERROR_TRACE;
    }

    *pA = (void *) zero;
}

static void dataMatViewZeroImpl( void *A );
static void dataMatViewSparseImpl( void *A );
static void dataMatViewDenseImpl( void *A );

static void dataMatCreateSparseImpl( void **pA, lorads_int nSDPCol, lorads_int dataMatNnz, lorads_int *dataMatIdx, double *dataMatElem ) {
    if ( !pA ) {
        goto exit_cleanup;
    }

    sdp_coeff_sparse *sparse;
    LORADS_INIT(sparse, sdp_coeff_sparse, 1);
    LORADS_MEMCHECK(sparse);
    sparse->nSDPCol = nSDPCol;
    sparse->nTriMatElem = dataMatNnz;

    LORADS_INIT(sparse->triMatRow, lorads_int, dataMatNnz);
    LORADS_INIT(sparse->triMatCol, lorads_int, dataMatNnz);
    LORADS_INIT(sparse->triMatElem, double, dataMatNnz);

    LORADS_MEMCHECK(sparse->triMatRow);
    LORADS_MEMCHECK(sparse->triMatCol);
    LORADS_MEMCHECK(sparse->triMatElem);

    // Note tsp_decompress work when ascending sort
    if ( !LUtilCheckIfAscending(dataMatNnz, dataMatIdx) ) {
        LUtilAscendSortDblByInt(dataMatElem, dataMatIdx, 0, dataMatNnz - 1);
    }

    tsp_decompress(sparse->nSDPCol, sparse->nTriMatElem, dataMatIdx, dataMatElem,
                   sparse->triMatRow, sparse->triMatCol, sparse->triMatElem);

    *pA = (void *) sparse;
    exit_cleanup:
        return;
}

static void dataMatCreateDenseImpl( void **pA, lorads_int nSDPCol, lorads_int dataMatNnz, lorads_int *dataMatIdx, double *dataMatElem ) {
    if ( !pA ) {
        LORADS_ERROR_TRACE;
    }

    sdp_coeff_dense *dense;
    LORADS_INIT(dense, sdp_coeff_dense, 1);
    LORADS_MEMCHECK(dense);

    dense->nSDPCol = nSDPCol;
    LORADS_INIT(dense->dsMatElem, double, PACK_NNZ(nSDPCol));
    LORADS_MEMCHECK(dense->dsMatElem);

    pds_decompress(dataMatNnz, dataMatIdx, dataMatElem, dense->dsMatElem);
    *pA = (void *) dense;
}


static void dataMatScalZeroImpl( void *A, double alpha ) {
    return;
}


static void dataMatScalSparseImpl( void *A, double alpha ) {
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    lorads_int incx = 1;
    scal(&spA->nTriMatElem, &alpha, spA->triMatElem, &incx);
    return;
}

static void dataMatSparseNrm1(void *A, double *res){
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    double nrmA = 0.0;
    for (lorads_int i = 0; i < spA->nTriMatElem; ++i){
        nrmA += 2 * fabs(spA->triMatElem[i]);
        lorads_int row = spA->triMatRow[i];
        lorads_int col = spA->triMatCol[i];
        if (row == col){
            nrmA -= fabs(spA->triMatElem[i]);
        }
    }
    res[0] = nrmA;
}

static void dataMatSparseNrm2Square(void *A, double *res){
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    double nrmA = 0.0;
    for (lorads_int i = 0; i < spA->nTriMatElem; ++i){
        nrmA += 2 * pow(spA->triMatElem[i], 2);
        lorads_int row = spA->triMatRow[i];
        lorads_int col = spA->triMatCol[i];
        if (row == col){
            nrmA -= pow(spA->triMatElem[i], 2);
        }
    }
    res[0] = nrmA;
}

static void dataMatSparseNrmInf(void *A, double *res){
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    double nrmA = 0.0;
    for ( lorads_int i = 0; i < spA->nTriMatElem; ++i ) {
        nrmA =  LORADS_MAX( LORADS_ABS(spA->triMatElem[i]), nrmA);
    }
    res[0] = nrmA;
}

extern void dataMatSparseZeros(void *A){
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    LORADS_ZERO(spA->triMatElem, double, spA->nTriMatElem);
}

static void dataMatSparseStatNnz(lorads_int *nnzStat){
    nnzStat[0]++;
    return;
}

extern void dataMatSparseScale(void *A, double scaleFactor){
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *)A;
    lorads_int incx = 1;
    scal(&sparse->nTriMatElem, &scaleFactor, sparse->triMatElem, &incx);
}

static void dataMatSparseCollectNnzPos(void *A, lorads_int *nnzRow, lorads_int *nnzCol, lorads_int *nnzIdx){
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    LORADS_MEMCPY(&nnzRow[nnzIdx[0]], spA->triMatRow, lorads_int, spA->nTriMatElem);
    LORADS_MEMCPY(&nnzCol[nnzIdx[0]], spA->triMatCol, lorads_int, spA->nTriMatElem);
    nnzIdx[0] += spA->nTriMatElem;
    return;
}

static void dataMatSparseReConstructIndex(void *A, Dict *dict){
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    lorads_int nnz = spA->nTriMatElem;
    lorads_int row, col;
    lorads_int idx = 0;
    LORADS_INIT(spA->nnzIdx2ResIdx, lorads_int, spA->nTriMatElem);
    for (lorads_int i = 0; i < nnz; ++i){
        row = spA->triMatRow[i];
        col = spA->triMatCol[i];
        idx = find_index(dict, row, col);
        if (idx == -1){
            LORADS_ERROR_TRACE;
        }
        spA->nnzIdx2ResIdx[i] = idx;
    }
    return;
}

static void dataMatDenseNrm1(void *A, double *res){
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    double nrmA = 0.0;
    lorads_int nElem = dsA->nSDPCol * (dsA->nSDPCol + 1) / 2;
    lorads_int colLen; ///< Exclude diagonal
    double colNrm; ///< Exclude diagonal
    double *p = dsA->dsMatElem;
    lorads_int nCol = dsA->nSDPCol;
    lorads_int incx = 1;
    for ( lorads_int i = 0; i < nCol; ++i ) {
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
    lorads_int nCol = dsA->nSDPCol;
    lorads_int incx = 1;
    lorads_int colLen; ///< Exclude diagonal
    double colNrm; ///< Exclude diagonal
    double *p = dsA->dsMatElem;
    for ( lorads_int i = 0; i < nCol; ++i ) {
        nrmA += pow(p[0], 2);
        colLen = nCol - i - 1;
        nrmA += pow(nrm2(&colLen, p + 1, &incx), 2) * 2;
        p = p + colLen + 1;
    }
    res[0] = nrmA;
}


static void dataMatDenseNrmInf(void *A, double *res){
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    lorads_int nElem = dsA->nSDPCol * (dsA->nSDPCol + 1) / 2;
    double *p = dsA->dsMatElem;
    res[0] = 0;
    for (lorads_int i = 0; i < nElem; ++i ){
        res[0] =  LORADS_MAX(res[0],  LORADS_ABS(p[i]));
    }
    return;
}

extern void dataMatDenseZeros(void *A){
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    LORADS_ZERO(dsA->dsMatElem, double, dsA->nSDPCol * (dsA->nSDPCol + 1) / 2);
}


extern void dataMatDenseScale(void *A, double scaleFactor){
    // A = A * scaleFactor
    sdp_coeff_dense *dense = (sdp_coeff_dense *)A;
    lorads_int n = dense->nSDPCol * (dense->nSDPCol + 1) / 2;
    scal(&n, &scaleFactor, dense->dsMatElem, &AIntConstantOne);
}

static void dataMatDenseStatNnz(lorads_int *nnzStat){
    nnzStat[0]++;
    return;
}

static void dataMatDenseReConstructIndex(void *A, Dict *dict){
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    return;
}

static lorads_int dataMatGetNnzZeroImpl( void *A ) {
    return 0;
}

static lorads_int dataMatGetNnzSparseImpl( void *A ) {
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    return spA->nTriMatElem;
}

static lorads_int dataMatGetNnzDenseImpl( void *A ) {
    sdp_coeff_dense *dsA = (sdp_coeff_dense *) A;
    return dsA->nSDPCol * (dsA->nSDPCol + 1) / 2;
}


static void dataMatGetZeroSparsityImpl( void *A, lorads_int *spout ) {
    return;
}

static void dataMatGetSparseSparsityImpl( void *A, lorads_int *spout ) {
    sdp_coeff_sparse *spA = (sdp_coeff_sparse *) A;
    for ( lorads_int k = 0; k < spA->nTriMatElem; ++k ) {
        spout[FULL_IDX(spA->nSDPCol, spA->triMatRow[k], spA->triMatCol[k])] = 1;
    }

    return;
}

static void dataMatGetDenseSparsityImpl( void *A, lorads_int *spout ) {
    /* This method should not be invoked */
    assert( 0 );
}


static void dataMatClearZeroImpl( void *A ) {
    if ( !A ) {
        return;
    }
    LORADS_ZERO(A, sdp_coeff_zero, 1);
    return;
}

static void dataMatClearSparseImpl( void *A) {
    if ( !A ) {
        return;
    }
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    LORADS_FREE(sparse->triMatRow);
    LORADS_FREE(sparse->triMatCol);
    LORADS_FREE(sparse->triMatElem);
    if (sparse->nnzIdx2ResIdx){
        LORADS_FREE(sparse->nnzIdx2ResIdx);
    }
    LORADS_ZERO(A, sdp_coeff_sparse, 1);
    return;
}

static void dataMatClearDenseImpl( void *A ) {
    if ( !A ) {
        return;
    }
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    LORADS_FREE(dense->dsMatElem);
    LORADS_ZERO(A, sdp_coeff_dense, 1);
    return;
}

static void dataMatDestroyZeroImpl( void **pA ) {
    if ( !pA ) {
        return;
    }
    dataMatClearZeroImpl(*pA);
    sdp_coeff *p = (sdp_coeff *)*pA;
    LORADS_FREE(*pA);
    return;
}

static void dataMatDestroySparseImpl( void **pA ) {
    if ( !pA ) {
        return;
    }
    sdp_coeff *p = (sdp_coeff *)*pA;
    dataMatClearSparseImpl(*pA);
    LORADS_FREE(*pA);
    return;
}

extern void dataMatDestroyDenseImpl( void **pA ) {
    if ( !pA ) {
        return;
    }
    dataMatClearDenseImpl(*pA);
    sdp_coeff *p = (sdp_coeff *)*pA;
    LORADS_FREE(*pA);
    return;
}

static void dataMatViewZeroImpl( void *A ) {
#ifdef INT32
    printf("Zero matrix of size %d L1 = [%5.3e] L2 = [%5.3e] \n",
           ((sdp_coeff_zero *) A)->nSDPCol, 0.0, 0.0);
#endif
#ifdef UNIX_INT64
    printf("Zero matrix of size %ld L1 = [%5.3e] L2 = [%5.3e] \n",
           ((sdp_coeff_zero *) A)->nSDPCol, 0.0, 0.0);
#endif
#ifdef MAC_INT64
    printf("Zero matrix of size %lld L1 = [%5.3e] L2 = [%5.3e] \n",
           ((sdp_coeff_zero *) A)->nSDPCol, 0.0, 0.0);
#endif
    return;
}

static void dataMatViewSparseImpl( void *A ) {
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
#ifdef INT32
    printf("Sparse matrix of size %d and %d nnzs. \n",
           sparse->nSDPCol, sparse->nTriMatElem);
#endif
#ifdef MAC_INT64
    printf("Sparse matrix of size %lld and %lld nnzs. \n",
           sparse->nSDPCol, sparse->nTriMatElem);
#endif

#ifdef UNIX_INT64
    printf("Sparse matrix of size %ld and %ld nnzs. \n",
           sparse->nSDPCol, sparse->nTriMatElem);
#endif
    return;
}



static void dataMatViewDenseImpl( void *A ) {
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
#ifdef INT32
    printf("Dense matrix of size %d. \n", dense->nSDPCol);
#endif
#ifdef UNIX_INT64
    printf("Dense matrix of size %ld. \n", dense->nSDPCol);
#endif
#ifdef MAC_INT64
    printf("Dense matrix of size %lld. \n", dense->nSDPCol);
#endif
    return;
}

static void dataMatZeroMultiRkMat(void *A, lorads_sdp_dense *X, double *AX){
    // AX = A * X = 0
    return;
}

static void dataMatZeroMultiRkMatInnerProRkMat(void *A, lorads_sdp_dense *U, lorads_sdp_dense *V, double *res, void *UVt, sdp_coeff_type UVtType){
    // res = <AU, V> = 0
    res[0] = 0.0; // no change
    return;
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

static void dataMatZeroNrmInf(void *A, double *res){
    res[0] = 0.0;
    return;
}

static void dataMatZeroStatNnz (lorads_int *nnzStat){
    return;
}

static void dataMatZeroScale(void *A, double scaleFactor){
    return;
}

static void dataMatZeroCollectNnzPos(void *A, lorads_int *nnzRow, lorads_int *nnzCol, lorads_int *nnzIdx){
    return;
}

static void dataMatReConstructIndex(void *A, Dict *dict){
    return;
}


extern void dataMatSparseMultiRkMat(void *A, lorads_sdp_dense *X, double *AX){
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    LORADS_ZERO(AX, double, X->rank * X->nRows);
    lorads_int row, col = 0;
    for (int i = 0; i < sparse->nTriMatElem; ++i){
        row = sparse->triMatRow[i];
        col = sparse->triMatCol[i];
        axpy(&X->rank, &sparse->triMatElem[i], &X->matElem[col], &sparse->nSDPCol, &AX[row], &sparse->nSDPCol);
        if (row != col){
            axpy(&X->rank, &sparse->triMatElem[i], &X->matElem[row], &sparse->nSDPCol, &AX[col], &sparse->nSDPCol);
        }
    }
    return;
}

extern void dataMatSparseMV(void *A, double *x, double *y, lorads_int n){
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    assert(n == sparse->nSDPCol);
    LORADS_ZERO(y, double, n);
    lorads_int row, col = 0;
    lorads_int one = 1;
    for (int i = 0; i < sparse->nTriMatElem; ++i){
        row = sparse->triMatRow[i];
        col = sparse->triMatCol[i];
        axpy(&one, &sparse->triMatElem[i], &x[col], &sparse->nSDPCol, &y[row], &sparse->nSDPCol);
        if (row != col){
            axpy(&one, &sparse->triMatElem[i], &x[row], &sparse->nSDPCol, &y[col], &sparse->nSDPCol);
        }
    }
    return;
}


static void sparseAUV(lorads_int n, lorads_int nnzA, lorads_int *Ai, lorads_int *Aj, double *Ax,
                      lorads_int nnzUVt, lorads_int *UVti, lorads_int *UVtj, double *UVtx, lorads_int *nnzIdx, double *res){
    lorads_int incx = 1;
    lorads_int rowA, colA, rowUVt, colUVt = 0;
    if (nnzA == nnzUVt){
        res[0] += 2 * dot(&nnzA, Ax, &incx, UVtx, &incx);
        for (lorads_int i = 0; i < nnzA; ++i){
            rowA = Ai[i]; colA = Aj[i];
            rowUVt = UVti[i]; colUVt = UVtj[i];
            if (rowA == colA && rowA == rowUVt && rowUVt == colUVt){
                res[0] -= Ax[i] * UVtx[i];
            }
        }
    }else{
        double temp;
        assert(nnzUVt >  nnzA);
        lorads_int i, j, idx = 0;
        for (i = 0; i < nnzA; ++i){
            rowA = Ai[i];
            colA = Aj[i];
            idx = nnzIdx[i];
            temp = 2 * Ax[i] * UVtx[idx];
            res[0] += temp;
            if (rowA == colA){
                res[0] -= 0.5 * temp;
            }
        }
    }
}

static void denseAUV(lorads_int nnzA, const lorads_int *Ai, const lorads_int *Aj, double *Ax,
                     const double *UVtx, lorads_int *nnzIdx, double *res){
    lorads_int rowA, colA, idx = 0;
    double temp = 0.0;
    for (lorads_int i = 0; i < nnzA; ++i){
        rowA = Ai[i]; colA = Aj[i];
        idx = nnzIdx[i];
        temp = 2 * Ax[i] * UVtx[idx];
        res[0] += temp;
        if (rowA == colA){
            res[0] -= 0.5 * temp;
        }
    }
}

static void dataMatSparseMultiRkMatInnerProRkMat(void *A, lorads_sdp_dense *U, lorads_sdp_dense *V, double *res, void *UVtIn, sdp_coeff_type UVtType){
    /*
    res = <AU, V>
    UVt is (UVt + VUt) / 2 in fact
     */
    res[0] = 0.0;

    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;

    if (UVtType == SDP_COEFF_SPARSE){
        sdp_coeff_sparse *UVt = (sdp_coeff_sparse *)UVtIn;
        sparseAUV(sparse->nSDPCol, sparse->nTriMatElem, sparse->triMatRow, sparse->triMatCol, sparse->triMatElem,
                  UVt->nTriMatElem, UVt->triMatRow, UVt->triMatCol, UVt->triMatElem, sparse->nnzIdx2ResIdx, res);
    }else{
        sdp_coeff_dense *UVt = (sdp_coeff_dense *)UVtIn;
        denseAUV(sparse->nTriMatElem, sparse->triMatRow, sparse->triMatCol, sparse->triMatElem,
                 UVt->dsMatElem, sparse->nnzIdx2ResIdx, res);
    }
}

static void dataMatSparseAddDenseSDPCoeff(void *A, void *B, double weight){
    // B = B + A * weight
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *) A;
    sdp_coeff_dense *dense = (sdp_coeff_dense *) B;
    lorads_int n = dense->nSDPCol;
    if (sparse->nnzIdx2ResIdx == NULL){
        for (lorads_int i = 0; i < sparse->nTriMatElem; ++i){
            lorads_int row = sparse->triMatRow[i];
            lorads_int col = sparse->triMatCol[i];
            if (row >= col){
                dense->dsMatElem[(n * col -col*(col+1)/2 + row)] += weight * sparse->triMatElem[i];
            }
        }
    }else{
        lorads_int idx; lorads_int row; lorads_int col;
        for (lorads_int i = 0; i < sparse->nTriMatElem; ++i){
            row = sparse->triMatRow[i];
            col = sparse->triMatCol[i];
            assert(row >= col);
//            idx = dense->rowCol2NnzIdx[row][col];
            idx = sparse->nnzIdx2ResIdx[i];
            assert(idx != -1);
            dense->dsMatElem[idx] += weight * sparse->triMatElem[i];
        }
    }
}



static void dataMatSparseAddSparseSDPCoeff(void *A, void *B, double weight){
    // B = B + weight * A
    sdp_coeff_sparse *sparseA = (sdp_coeff_sparse *) A;
    sdp_coeff_sparse *sparseB = (sdp_coeff_sparse *) B;
    lorads_int nnzA = sparseA->nTriMatElem;
    lorads_int nnzB = sparseB->nTriMatElem;
    assert(nnzB >= nnzA);

    lorads_int rowA, colA; lorads_int idx;
    for (lorads_int i = 0; i < nnzA; ++i){
        rowA = sparseA->triMatRow[i];
        colA = sparseA->triMatCol[i];
//        idx = sparseB->rowCol2NnzIdx[rowA][colA];
        idx = sparseA->nnzIdx2ResIdx[i];
        sparseB->triMatElem[idx] += weight * sparseA->triMatElem[i];
    }
}

static void dataMatSparseAddSDPCoeff(void *A, void *B, double weight, sdp_coeff_type B_type){
    if (B_type == SDP_COEFF_DENSE){
        dataMatSparseAddDenseSDPCoeff(A, B, weight);
    }else if (B_type == SDP_COEFF_SPARSE){
        dataMatSparseAddSparseSDPCoeff(A, B, weight);
    }else{
        LORADS_ERROR_TRACE;
    }
}

extern void dataMatDenseMultiRkMat(void *A, lorads_sdp_dense *X, double *AX){
    /*
     AX = A * X, A is dense meanse result is dense
     */
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    double alpha = 1.0;
    double beta = 0.0;
    LORADS_ZERO(dense->fullMat, double, dense->nSDPCol * dense->nSDPCol);
    lorads_int idx = 0;
    for (lorads_int col = 0; col < dense->nSDPCol; ++col){
        for (lorads_int row = col; row < dense->nSDPCol; ++row){
            dense->fullMat[dense->nSDPCol * col + row] = dense->dsMatElem[idx];
            dense->fullMat[dense->nSDPCol * row + col] = dense->dsMatElem[idx];
            idx++;
        }
    }
    char side = 'L'; // C:= alpha * A * B + beta * C;
    char uplo = 'L';
    lorads_int m = dense->nSDPCol;
    lorads_int n = X->rank;
#ifdef UNDER_BLAS
    dsymm_(&side, &uplo, &m, &n, &alpha, dense->fullMat, &m, X->matElem, &m, &beta, AX, &m );
#else
    dsymm(&side, &uplo, &m, &n, &alpha, dense->fullMat, &m, X->matElem, &m, &beta, AX, &m );
#endif
}

extern void dataMatDenseMV(void *A, double *x, double *y, lorads_int n){
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    double alpha = 1.0;
    double beta = 0.0;
    LORADS_ZERO(dense->fullMat, double, dense->nSDPCol * dense->nSDPCol);
    lorads_int idx = 0;
    for (lorads_int col = 0; col < dense->nSDPCol; ++col){
        for (lorads_int row = col; row < dense->nSDPCol; ++row){
            dense->fullMat[dense->nSDPCol * col + row] = dense->dsMatElem[idx];
            dense->fullMat[dense->nSDPCol * row + col] = dense->dsMatElem[idx];
            idx++;
        }
    }
    lorads_int m = dense->nSDPCol;
    assert(n == m);
    char side = 'L'; // C:= alpha * A * B + beta * C;
    char uplo = 'L';
    lorads_int k = 1;
#ifdef UNDER_BLAS
    dsymm_(&side, &uplo, &m, &k, &alpha, dense->fullMat, &m, x, &m, &beta, y, &m );
#else
    dsymm(&side, &uplo, &m, &k, &alpha, dense->fullMat, &m, x, &m, &beta, y, &m );
#endif
}

static void dataMatDenseMultiRkMatInnerProRkMat(void *A, lorads_sdp_dense *U, lorads_sdp_dense *V, double *res, void *UVtIn, sdp_coeff_type UVtType){
    /* res = <AU, V>
     dense A, dense U, dense V
     */
    res[0] = 0.0;
    sdp_coeff_dense *dense = (sdp_coeff_dense *) A;
    lorads_int n = 0;
    lorads_int incx = 1;
    sdp_coeff_dense *UVt = (sdp_coeff_dense *)UVtIn;
    n = UVt->nSDPCol * (UVt->nSDPCol + 1) / 2;
    res[0] += 2 * dot(&n, dense->dsMatElem, &incx, UVt->dsMatElem, &incx);
    lorads_int idx = 0; lorads_int colNum = UVt->nSDPCol;
    for (lorads_int i = 0; i < UVt->nSDPCol; ++i){
        res[0] -= dense->dsMatElem[idx] * UVt->dsMatElem[idx];
        idx += colNum;
        colNum -= 1;
    }
}

static void dataMatDenseAddDenseSDPCoeff(void *A, void *B, double weight){
    // B += A * weight
    sdp_coeff_dense *denseA = (sdp_coeff_dense *) A;
    sdp_coeff_dense *denseB = (sdp_coeff_dense *) B;
    lorads_int n = denseA->nSDPCol * (denseA->nSDPCol + 1) / 2;
    lorads_int incx = 1;
    axpy(&n, &weight, denseA->dsMatElem, &incx, denseB->dsMatElem, &incx);
}

static void dataMatDenseAddSDPCoeff(void *A, void *B, double weight, sdp_coeff_type B_type){
    if (B_type == SDP_COEFF_DENSE){
        dataMatDenseAddDenseSDPCoeff(A, B, weight);
    }else{
        LORADS_ERROR_TRACE;
    }
}

static void sdpDataMatIChooseType(sdp_coeff *sdpCoeff, sdp_coeff_type dataType ) {
    sdpCoeff->dataType = dataType;
    switch (dataType) {
        case SDP_COEFF_ZERO:
            sdpCoeff->create = dataMatCreateZeroImpl;
            sdpCoeff->getnnz = dataMatGetNnzZeroImpl;
            sdpCoeff->getmatnz = dataMatGetZeroSparsityImpl;
            sdpCoeff->destroy = dataMatDestroyZeroImpl;
            sdpCoeff->view = dataMatViewZeroImpl;
            sdpCoeff->mul_rk = dataMatZeroMultiRkMat;
            sdpCoeff->mv = NULL;
            sdpCoeff->mul_inner_rk_double = dataMatZeroMultiRkMatInnerProRkMat;
            sdpCoeff->add_sdp_coeff = dataMatZeroAddDenseSDPCoeff;
            sdpCoeff->nrm1 = dataMatZeroNrm1;
            sdpCoeff->nrm2Square = dataMatZeroNrm2Square;
            sdpCoeff->nrmInf = dataMatZeroNrmInf;
            sdpCoeff->zeros = NULL;
            sdpCoeff->statNnz = dataMatZeroStatNnz;
            sdpCoeff->scaleData = dataMatZeroScale;
            sdpCoeff->collectNnzPos = dataMatZeroCollectNnzPos;
            sdpCoeff->reConstructIndex = dataMatReConstructIndex;
            break;
        case SDP_COEFF_SPARSE:
            sdpCoeff->create = dataMatCreateSparseImpl;
            sdpCoeff->getnnz = dataMatGetNnzSparseImpl;
            sdpCoeff->getmatnz = dataMatGetSparseSparsityImpl;
            sdpCoeff->destroy = dataMatDestroySparseImpl;
            sdpCoeff->view = dataMatViewSparseImpl;
            sdpCoeff->mul_rk = dataMatSparseMultiRkMat;
            sdpCoeff->mv = dataMatSparseMV; // for dual infeasibility
            sdpCoeff->mul_inner_rk_double = dataMatSparseMultiRkMatInnerProRkMat;
            sdpCoeff->add_sdp_coeff = dataMatSparseAddSDPCoeff;
            sdpCoeff->nrm1 = dataMatSparseNrm1;
            sdpCoeff->nrm2Square = dataMatSparseNrm2Square;
            sdpCoeff->nrmInf = dataMatSparseNrmInf;
            sdpCoeff->zeros = dataMatSparseZeros;
            sdpCoeff->statNnz = dataMatSparseStatNnz;
            sdpCoeff->scaleData = dataMatSparseScale;
            sdpCoeff->collectNnzPos = dataMatSparseCollectNnzPos;
            sdpCoeff->reConstructIndex = dataMatSparseReConstructIndex;
            break;
        case SDP_COEFF_DENSE:
            sdpCoeff->create = dataMatCreateDenseImpl;
            sdpCoeff->getnnz = dataMatGetNnzDenseImpl;
            sdpCoeff->getmatnz = dataMatGetDenseSparsityImpl;
            sdpCoeff->destroy = dataMatDestroyDenseImpl;
            sdpCoeff->view = dataMatViewDenseImpl;
            sdpCoeff->mul_rk = dataMatDenseMultiRkMat;
            sdpCoeff->mv = dataMatDenseMV; // for dual infeasibility
            sdpCoeff->mul_inner_rk_double = dataMatDenseMultiRkMatInnerProRkMat;
            sdpCoeff->add_sdp_coeff = dataMatDenseAddSDPCoeff;
            sdpCoeff->nrm1 = dataMatDenseNrm1;
            sdpCoeff->nrm2Square = dataMatDenseNrm2Square;
            sdpCoeff->nrmInf = dataMatDenseNrmInf;
            sdpCoeff->zeros = dataMatDenseZeros;
            sdpCoeff->statNnz = dataMatDenseStatNnz;
            sdpCoeff->scaleData = dataMatDenseScale;
            sdpCoeff->collectNnzPos = NULL; // not used
            sdpCoeff->reConstructIndex = dataMatDenseReConstructIndex; // not used
            break;
        default:
            assert(0);
            break;
    }
}

/* External methods for the SDP data */
extern void sdpDataMatCreate( sdp_coeff **psdpCoeff ) {
    if ( !psdpCoeff ) {
        LORADS_ERROR_TRACE;
    }
    sdp_coeff *sdpCoeff;
    LORADS_INIT(sdpCoeff, sdp_coeff, 1);
    LORADS_MEMCHECK(sdpCoeff);
    *psdpCoeff = sdpCoeff;
}

extern void sdpDataMatSetData( sdp_coeff *sdpCoeff, lorads_int nSDPCol, lorads_int dataMatNnz, lorads_int *dataMatIdx, double *dataMatElem ) {
    // sdpCoeff  : data mat struct
    // nSDPCol   : data mat dimension
    // dataMatNnz: data mat non-zeros number
    sdpCoeff->nSDPCol = nSDPCol;

    /* At this stage, only sparse, zero and dense matrices are classified */
    if ( dataMatNnz == 0 ) {
        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_ZERO);
    } else if ( (double) dataMatNnz > 0.1 * (double)(nSDPCol * (nSDPCol + 1) / 2) ) {
        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_DENSE);
    } else {
        sdpDataMatIChooseType(sdpCoeff, SDP_COEFF_SPARSE);
    }
    /* Create data */
    sdpCoeff->create(&sdpCoeff->dataMat, nSDPCol,
                               dataMatNnz, dataMatIdx, dataMatElem);
}

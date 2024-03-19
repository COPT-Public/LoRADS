#ifdef HEADERPATH
#include "src/sparse_opts.h"
#include "src/vec_opts.h"
#include "src/asdp_debug.h"
#include "interface/asdp_utils.h"
#else
#include "sparse_opts.h"
#include "vec_opts.h"
#include "asdp_utils.h"
#include "asdp_sdpdata.h"
#include "asdp_debug.h"
#include "dense_opts.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <math.h>

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

/* Compressed column operations */
extern void csp_Axpby( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y ) {
    // CSC matrix multiplication
    // y = a * A * x + y
    // where A is a sparse matrix in CSC format
    // Ap is the pointer to the start of each column
    // Ai is the row index of each non-zero element
    // Ax is the value of each non-zero element
    // n is the dimension of the matrix
    if ( a == 0.0 ) {
        return;
    }
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            y[Ai[j]] += a * x[i] * Ax[j];
        }
    }
    
    return;
}


extern void csp_ATxpby( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y ) {
    // CSC matrix multiplication
    // y = a * A^T * x + y
    // where A is a sparse matrix in CSC format
    // Ap is the pointer to the start of each column
    // Ai is the row index of each non-zero element
    // Ax is the value of each non-zero element
    // n is the dimension of the matrix
    if ( a == 0.0 ) {
        return;
    }
    
    double aTy = 0.0;
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i], aTy = 0.0; j < Ap[i + 1]; ++j ) {
            aTy += x[Ai[j]] * Ax[j];
        }
        y[i] += a * aTy;
    }
    
    return;
}

extern double csp_sum_abs( int n, int *Ap, int *Ai, double *Ax ) {
    // Sum of absolute values of a sparse matrix
    // where A is a sparse matrix in CSC format
    // Ap is the pointer to the start of each column
    // Ai is the row index of each non-zero element
    // Ax is the value of each non-zero element
    // n is the dimension of the matrix
    double sabs = 0.0;
    for ( int i = 0, j; i < Ap[n]; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++i ) {
            sabs += ( Ai[i] == i ) ? 0.5 * fabs(Ax[i]) : fabs(Ax[i]);
        }
    }
    
    return sabs;
}

extern double csp_fro_norm( int n, int *Ap, int *Ai, double *Ax ) {
    // Frobenius norm of a sparse matrix
    // where A is a sparse matrix in CSC format
    // Ap is the pointer to the start of each column
    // Ai is the row index of each non-zero element
    // Ax is the value of each non-zero element
    // n is the dimension of the matrix
    double nrm = 0.0;
    for ( int i = 0, j; i < Ap[n]; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++i ) {
            nrm += ( Ai[i] == i ) ? 0.5 * Ax[i] * Ax[i]: Ax[i] * Ax[i];
        }
    }
    
    return sqrt(nrm);
}

/** @brief Column sparse matrix aApB
 *  This function is called in the following contexts:
 *
 *  1. Updating dual variable S + a \* dS
 *
 */
extern void csp_aApB( int n, int nnz, double a, int *Al, double *Ax, double *Bx ) {
    
    if ( Al ) {
        /* In this case B is sparse and Al indicates where each Ax is located in Bx */
        for ( int i = 0; i < nnz; ++i ) {
            Bx[i] += a * Ax[Al[i]];
        }
        
    } else {
        /* In this case B is dense and axpy is sufficient*/
        
    }
    
    return;
}

/** @brief Get the number of `nonzero columns` in an csc matrix
 *
 */
extern int csp_nnz_cols ( int n, int *Ap ) {
    
    int nzcols = 0;
    
    for ( int i = 0; i < n; ++i ) {
        nzcols += ( Ap[i + 1] - Ap[i] > 0 );
    }
    
    return nzcols;
}

extern void csp_dump( int n, int *Ap, int *Ai, double *Ax, double *v ) {
    // csc matrix to dense matrix
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            v[i * n + Ai[j]] = v[Ai[j] * n + i] = Ax[j];
        }
    }
    
    return;
}

/* Decompress a column to triplet matrix */
extern void tsp_decompress( int n, int nnz, int *Ci, double *Cx, int *Ai, int *Aj, double *Ax ) {
    
    int j = 0, idthresh = n;
    
    for ( int k = 0; k < nnz; ++k ) {
        while ( Ci[k] >= idthresh ) {
            j += 1;
            idthresh += n - j;
        }
        Ai[k] = Ci[k] - idthresh + n;
        Aj[k] = j;
        Ax[k] = Cx[k];
    }
    
    return;
}

/** @brief Check if the matrix is rank-one
 *
 * a is an n by one auxiliary vector and on exit, the returned integer values tells
 * if the matrix is rank-one and a would contain the rank-one element
 *
 * In a word A = sgn \* a \* a'
 *
 */
extern int tsp_r1_extract( int n, int nnz, int *Ai, int *Aj, double *Ax, double *sgn, double *a ) {
    // rank 1 decomposition:  A = sgn \* a \* a'
    assert( nnz > 0 );
    
    /* Get the first nonzero */
    int i = Ai[0];
    int j = Aj[0];
    double v = Ax[0];
    
    if ( i != j ) {
        return 0;
    }
    
    if ( nnz == 1 ) {
        *sgn = Ax[0];
        a[i] = 1.0;
        return 1;
    }
    
    double s = ( v > 0 ) ? 1.0 : -1.0;
    v = sqrt(fabs(v));
    
    /* Assume that the matrix is rank-one */
    int k = 0;
    int anz = 0;
    for ( k = 0; k < nnz; ++k ) {
        if ( Aj[k] > i ) {
            break;
        }
        a[Ai[k]] = Ax[k] / v;
        anz += 1;
    }
    
    /* The number of nnzs in the submatrix must match */
    if ( nnz != (int) (anz * ( anz + 1 ) / 2) ) {
        return 0;
    }
    
    if ( k == n ) {
        return 0;
    }
    
    /* Now a contains the rank-one components */
    double eps = 0.0;
    
    if ( s == 1.0 ) {
        /* If condition is broken into two cases */
        for ( int k = 0; k < nnz; ++k ) {
            eps += fabs(Ax[k] - a[Ai[k]] * a[Aj[k]]);
        }
    } else {
        for ( int k = 0; k < nnz; ++k ) {
            eps += fabs(Ax[k] + a[Ai[k]] * a[Aj[k]]);
        }
    }
    
    if ( eps > 1e-10 ) {
        return 0;
    }
    
    *sgn = s;
    
    return 1;
}


/* Lower triplet operations */
/** @brief Scale a triplet matrix
 *
 */
extern void tsp_scal( double a, int nnz, double *Ax ) {
    // scale a triplet matrix
    int incx = 1;
    scal(&nnz, &a, Ax, &incx);
    
    return;
}

/** @brief Compute one norm of a triplet matrix
 *
 */
extern double tsp_sum_abs( int nnz, int *Ai, int *Aj, double *Ax ) {
    // sum of absolute values
    double nrm = 0.0;
    for ( int i = 0; i < nnz; ++i ) {
        nrm += ( Ai[i] == Aj[i] ) ? 0.5 * fabs(Ax[i]) : fabs(Ax[i]);
    }
    
    return 2.0 * nrm;
}

extern double tsp_fro_norm( int nnz, int *Ai, int *Aj, double *Ax ) {
    // compute the Frobenius norm of a triplet matrix
    double nrm = 0.0;
    for ( int i = 0; i < nnz; ++i ) {
        nrm += ( Ai[i] == Aj[i] ) ? 0.5 * Ax[i] * Ax[i]: Ax[i] * Ax[i];
    }
    
    return sqrt(2.0 * nrm);
}

extern void tsp_dump( int n, int nnz, int *Ai, int *Aj, double *Ax, double *v ) {
    // triplet matrix to dense matrix
    int i, j;
    for ( int k = 0; k < nnz; ++k ) {
        i = Ai[k]; j = Aj[k];
        v[n * i + j] = v[n * j + i] = Ax[k];
    }
    
    return;
}

extern double tsp_quadform( int n, int nnz, int *Ai, int *Aj, double *Ax, double *v ) {
    // quadratic form: v' * A * v
    double quadform = 0.0, tmp = 0.0;
    
    for ( int k = 0; k < nnz; ++k ) {
        tmp = Ax[k] * v[Ai[k]] * v[Aj[k]];
        if ( Ai[k] == Aj[k] ) {
            // diagonal element
            quadform += 0.5 * tmp;
        } else {
            // non-diagonal element
            quadform += tmp;
        }
    }
    
    return 2.0 * quadform;
}

extern void tsp_Axpby( int n, int nnz, double a, int *Ai, int *Aj, double *Ax, double *x, double *y ) {
    /*
      Triplet matrix multiplication
      y = a * A * x + y
      where A is a sparse matrix in triplet format
      Ai is the row index of each non-zero element
      Aj is the column index of each non-zero element
      Ax is the value of each non-zero element
      n is the dimension of the matrix
     */
    if ( a == 0.0 ) {
        return;
    }
    
    for ( int i = 0; i < nnz; ++i ) {
        y[Ai[i]] += a * x[Aj[i]] * Ax[i];
    }
    
    return;
}

extern void tsp_ASemiLowby( int n, int nnz, double a, int *Ai, int *Aj, double *Ax, double *x, double *y ){
    /*
     Triplet matrix vector multiplication
     y = a * A * x + y;
     where A is a sparse matrix in triplet format
     Ai is the row index of each non-zero element
     Aj is the column index of each non-zero element
     Ax is the value of each non-zero element
     n is the dimension of the matrix
     Note: only store lower semi triangular part
     */
    if ( a == 0.0 ) {
        return;
    }
    for ( int i = 0; i < nnz; ++i){
        int row = Ai[i];
        int col = Aj[i];
        y[row] += a * x[col] * Ax[i];
        if (row != col){
            y[col] += a * x[row] * Ax[i];
        }
    }
}

extern void tsp_ASemiLowbyMul(int n, int nnz, double a, int *Ai, int *Aj, double *Ax, double *x, int k, double *y){
    /*
     Triplet matrix matrix multiplication
     */
    if (a == 0.0){
        return;
    }
    int row = 0;
    int col = 0;
    for (int i = 0; i < nnz; ++i){
        row = Ai[i];
        col = Aj[i];
        axpy(&k, &Ax[i], &x[col], &n, &y[row], &n);
        if (row != col){
            axpy(&k, &Ax[i], &x[row], &n, &y[col], &n);
        }
    }
}




extern void ATspx(int m, int n, double *A, int nnz, int *idx, double *x, double *res){
    /*
      res = A' * x
      A is dense matrix (m x n)
      x is sparse with index set idx
     */
    double zero = 0.0;
    int incx = 1;
    int incy = 1;
    double alpha = 1.0;
    int index = 1;
    ASDP_ZERO(res, double, n);
    
    if (n < 20){
        for (int i = 0; i < nnz; ++i){
            index = idx[i];
            for (int row = 0; row < n; ++row){
                // row is transpose matrix's row
                res[row] += A[row * m + index] * x[i];
            }
        }
    }else{
        double *ATemp;
        int num = nnz * n;
        ASDP_INIT(ATemp, double, num);
        ASDP_ZERO(ATemp, double, num);
        for (int i = 0; i < nnz; ++i){
            index = idx[i];
            axpy(&n, &alpha, &A[index], &m, &ATemp[i * n], &incy);
        }
        double beta = 0.0;
#ifdef UNDER_BLAS
        dgemv_(&ACharConstantNoTrans, &n, &nnz, &alpha, ATemp, &n, x, &incx, &beta, res, &incy);
#else
        dgemv(&ACharConstantNoTrans, &n, &nnz, &alpha, ATemp, &n, x, &incx, &beta, res, &incy);
#endif
        ASDP_FREE(ATemp);
    }
}


static void coo2csr(int n, int nnz, int *cooRowInd, int *cooColInd, double *cooVal,
             int **cscColPtrRES, int **cscRowIndRES, double **cscValRES) {
    // NOTE: csc and csr is the same in symmetry matrix
    // Initialize cscColPtr
    int *cscColPtr = NULL;
    ASDP_INIT(cscColPtr, int, n+1);

    // nonzero in each column
    int count = 0;
    for (int i = 0; i < nnz; i++) {
        count += 2;
        cscColPtr[cooColInd[i] + 1]++;
        cscColPtr[cooRowInd[i] + 1]++;
        if (cooRowInd[i] == cooColInd[i]){
            cscColPtr[cooRowInd[i] + 1]--;
            count -= 1;
        }
    }
    
    int *cscRowInd;
    double *cscVal;
    
    ASDP_INIT(cscRowInd, int, count);
    ASDP_INIT(cscVal, double, count);
    
    // calculate cscColPtr
    for (int i = 0; i < n; i++) {
        cscColPtr[i + 1] += cscColPtr[i];
    }

    // fill in cscRowInd and cscVal
    for (int i = 0; i < nnz; i++) {
        int row = cooRowInd[i];
        int col = cooColInd[i];
        
        if (row != col){
            int dest2 = cscColPtr[row];
            cscRowInd[dest2] = cooColInd[i];
            cscVal[dest2] = cooVal[i];
            cscColPtr[row]++;
        }
        
        int dest1 = cscColPtr[col];
        

        cscRowInd[dest1] = cooRowInd[i];
        cscVal[dest1] = cooVal[i];
    
        cscColPtr[col]++;
        
    }

    // recover cscColPtr
    for (int i = n; i > 0; i--) {
        cscColPtr[i] = cscColPtr[i - 1];
    }
    cscColPtr[0] = 0;
    
    *cscRowIndRES = cscRowInd;
    *cscColPtrRES = cscColPtr;
    *cscValRES = cscVal;
}


/* mkl function
 extern void sparse_Min_Eig(sdp_coeff *slackVar, double *minEig){
    sdp_coeff_sparse *slackVarMat = (sdp_coeff_sparse *)slackVar->dataMat;
    int n = slackVarMat->nSDPCol;
    int nnz = slackVarMat->nTriMatElem;
    int *csrColPtr;
    int *csrRowInd;
    double *csrVal;
    coo2csr(slackVarMat->nSDPCol, nnz, slackVarMat->triMatRow, slackVarMat->triMatCol, slackVarMat->triMatElem, &csrColPtr, &csrRowInd, &csrVal);
    
    char which = 'S';

    int  pm[128];
    mkl_sparse_ee_init(pm);
    
    sparse_status_t retcode;
    sparse_matrix_t A = NULL;
    retcode = mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ZERO, n, n, csrColPtr, csrColPtr + 1, csrRowInd, csrVal);
    
    if (retcode != 0){
        printf("csr create failed!\n");
    }
    
    struct matrix_descr descrA;
    descrA.type = SPARSE_MATRIX_TYPE_SYMMETRIC;
//    descrA.type = SPARSE_MATRIX_TYPE_GENERAL;
    descrA.mode = SPARSE_FILL_MODE_LOWER;
//    descrA.mode = SPARSE_FILL_MODE_UPPER;
    descrA.diag = SPARSE_DIAG_NON_UNIT;
    
    int k = 0;
    // how many eig you need
    int k0 = 1;
    // receive minimum eig
    double *E;
    ASDP_INIT(E, double, k0);
    // receive eig vector
    double *X;
    ASDP_INIT(X, double, k0 * n);
    // receive residual
    double *res;
    ASDP_INIT(res, double, k0);
    
    
    retcode = mkl_sparse_d_ev(&which, pm, A,  descrA, k0, &k, E, X, res);
    
    if (retcode != 0){
        printf("error code:%d\n", retcode);
        printf("mkl_sparse_d_ev failed!\n");
    }
    
    *minEig = E[0];
    
    mkl_sparse_destroy(A);
    free(E);
    free(X);
    free(res);
    ASDP_FREE(csrColPtr);
    ASDP_FREE(csrRowInd);
    ASDP_FREE(csrVal);
}
*/

//// not used
//
//#ifdef HEADERPATH
//#include "linalg/def_ asdp_lanczos.h"
//#include "linalg/ asdp_lanczos.h"
//#include "linalg/vec_opts.h"
//#include "linalg/dense_opts.h"
//#include "interface/asdp_utils.h"
//#else
//#include "def_asdp.h"
//#include "asdp.h"
//#include <stdio.h>
//#include "asdp_utils.h"
//#include "asdp_user_data.h"
//#include "def_asdp_user_data.h"
//#include "asdp_conic.h"
//#include "asdp_algo.h"
//#include "asdp_sdpdata.h"
//#include "asdp_debug.h"
//#include "vec_opts.h"
//#include "dense_opts.h"
//#include "sparse_opts.h"
//#include "asdp_direct_linsys.h"
//#include "asdp_lbfgs.h"
//#include "asdp_chol.h"
//#include "asdp_rank_reduce.h"
//#endif
//
//#include <math.h>
//
//extern void dgeqp3(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
//extern void dgemm(char *transA, char *transB, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
//
//extern void transpose(double *A, int m, int n, double **AtPointer){
//    // transpose A
//    if (*AtPointer == A){
//        double *newAt;
//        ASDP_INIT(newAt, double, m * n);
//        for (int col = 0; col < n; ++col){
//            for (int row = 0; row < m; ++row){
//                newAt[row * n + col] = A[col * m + row];
//            }
//        }
//        ASDP_FREE(A);
//        *AtPointer = newAt;
//    }else{
//        double *At = *AtPointer;
//        for (int col = 0; col < n; ++col){
//            for (int row = 0; row < m; ++row){
//                At[row * n + col] = A[col * m + row];
//            }
//        }
//    }
//}
//
//extern void LQ_decomposition(double *Umat, int nRows, int rank,  double **L, int **P){
//    assert(nRows >= rank);
//    double *UmatTran;
//    ASDP_INIT(UmatTran, double, nRows * rank);
//    transpose(Umat, nRows, rank, &UmatTran);
//    int rowNum = rank;
//    int colNum = nRows;
//    // matrix is (m x n), Q is (m x m) matrix and R is (m x n) matrix
//    double *tau;
//    ASDP_INIT(tau, double, rowNum);
//    int lwork = 3 * colNum + 2;
//    double *work;
//    ASDP_INIT(work, double, lwork);
//    int *jpvt;
//    ASDP_INIT(jpvt, int, colNum);
//    int info = 0;
//    dgeqp3(&rowNum, &colNum, UmatTran, &rowNum, jpvt, tau, work, &lwork, &info);
//    assert(info == 0);
//    double *Ltemp;
//    ASDP_INIT(Ltemp, double, rowNum * colNum);
//    int idx = 0;
//    for (int col = 0; col < colNum; ++col){
//        for (int row = 0; row < rowNum; ++row){
//            if (row > col) {
//                Ltemp[idx] = 0.0;
//            }else{
//                Ltemp[idx] = UmatTran[idx];
//            }
//
//            idx++;
//        }
//    }
//    transpose(Ltemp, rowNum, colNum, &Ltemp);
//
//    int *Ptemp;
//    ASDP_INIT(Ptemp, int, colNum * colNum);
//    for (int col = 0; col < colNum; ++col){
//        Ptemp[col * colNum + jpvt[col] - 1] = 1;
//    }
//    *P = Ptemp;
//    *L = Ltemp;
//    ASDP_FREE(UmatTran);
//    ASDP_FREE(tau);
//    ASDP_FREE(work);
//    ASDP_FREE(jpvt);
//}
//
//extern asdp_retcode rank_reduce(asdp_rk_mat_dense **RPointer, int *rankSolver, double tol){
//    asdp_retcode retcode = ASDP_RETCODE_OK;
//    // example:
//    // rank_reduce(&ASolver->R[0], &ASolver->rankDim[0], 1e-6);
//    asdp_rk_mat_dense *R = *RPointer;
//    // do rank reduce for R
//    int m = R->nRows;
//    int n = R->rank;
//    assert(n == rankSolver[0]);
//    double *L;
//    int *P;
//    LQ_decomposition(R->matElem, m, n, &L, &P);
//
//    double *nrmCol;
//    ASDP_INIT(nrmCol, double, n);
//    ASDP_MEMCHECK(nrmCol);
//    int incx = 1;
//    int nrmLen = 0;
//    for (int i = 0; i < n; ++i){
//        nrmLen = (n - i) * m;
//        nrmCol[i] = nrm2(&nrmLen, &L[m * i], &incx);
//    }
//
//
//    int allElem = m * n;
//    double nrmAll = nrmCol[0];
//
//    rscl(&n, &nrmAll, nrmCol, &incx );
//
//    int threshold = n - 1;
//    for (int i = n - 1; i >= 0; --i){
//        if (nrmCol[i] < 1e-5){
//            threshold = i;
//        }
//    }
//
//    threshold = ASDP_MAX(threshold, 1);
//
//    double *LTilde;
//    ASDP_INIT(LTilde, double, m * threshold);
//    ASDP_MEMCHECK(LTilde);
//    ASDP_MEMCPY(LTilde, L, double, m * threshold);
//
//    char transA = ACharConstantNoTrans;
//    char transB = ACharConstantNoTrans;
//    double alpha = 1.0;
//    double beta = 0.0;
//    double *RTilde;
//    ASDP_INIT(RTilde, double, m * threshold);
//    ASDP_MEMCHECK(RTilde);
//    double *doubleP;
//    ASDP_INIT(doubleP, double, m * m);
//    ASDP_MEMCHECK(doubleP);
//    for (int i = 0; i < m * m; ++i){
//        if (P[i] == 1){
//            doubleP[i] = 1.0;
//        }
//    }
//    dgemm(&transA, &transB, &m, &threshold, &m, &alpha, doubleP, &m, LTilde, &m, &beta, RTilde, &m);
//
//
//    // change R to RTilde
//    R->rank = threshold;
//    rankSolver[0] = threshold;
//    ASDP_FREE(R->matElem);
//    R->matElem = RTilde;
//
//    ASDP_FREE(nrmCol);
//    ASDP_FREE(L);
//    ASDP_FREE(P);
//    ASDP_FREE(LTilde);
//
//exit_cleanup:
//    return retcode;
//}
//
//extern asdp_retcode setASDPRkMatDense(asdp_rk_mat_dense **dstPointer, asdp_rk_mat_dense **srcPointer){
//    asdp_retcode retcode = ASDP_RETCODE_OK;
//    asdp_rk_mat_dense *dst = *dstPointer;
//    asdp_rk_mat_dense *src = *srcPointer;
//    dst->rank = src->rank;
//    dst->nRows = src->nRows;
//    ASDP_FREE(dst->matElem);
//    ASDP_INIT(dst->matElem, double, dst->rank * dst->nRows);
//    ASDP_MEMCHECK(dst->matElem);
//    ASDP_MEMCPY(dst->matElem, src->matElem, double, dst->rank * dst->nRows);
//exit_cleanup:
//    return retcode;
//}
//
//
//extern asdp_retcode ASDPRankReduce(asdp *ASolver){
//    asdp_retcode retcode = ASDP_RETCODE_OK;
//    for (int iCone = 0; iCone < ASolver->nCones; ++iCone){
//        ASDP_CALL(rank_reduce(&ASolver->V[iCone], &ASolver->rankElem[iCone], ASDP_RANK_REDUCE_TOL));
//        // copy V to R, U and S
//        ASDP_CALL(setASDPRkMatDense(&ASolver->R[iCone], &ASolver->V[iCone]));
//        ASDP_CALL(setASDPRkMatDense(&ASolver->U[iCone], &ASolver->V[iCone]));
//        ASDP_CALL(setASDPRkMatDense(&ASolver->S[iCone], &ASolver->V[iCone]));
//    }
//exit_cleanup:
//    return retcode;
//}

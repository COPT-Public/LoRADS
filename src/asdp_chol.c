//#ifdef HEADERPATH
//#include "interface/asdp.h"
//#include "interface/asdp_utils.h"
//#include "src/def_asdp_linsolver.h"
//#include "src/dense_opts.h"
//#include "src/vec_opts.h"
//#else
//#include "asdp.h"
//#include "asdp_utils.h"
//#include "def_asdp_linsolver.h"
//#include "dense_opts.h"
//#include "vec_opts.h"
//#include "asdp_cg.h"
//#include "asdp_debug.h"
//#include "def_asdp_rk_mat.h"
//#include "def_asdp_conic.h"
//#include "asdp_algo.h"
//#include "cholmod.h"
//#endif
//
//#ifdef MEMDEBUG
//#include "memwatch.h"
//#endif
//
///* The main goal is to determine the matrix is positive definite or not. */
//
//extern void ASDPPosDefDetermineSparse(void *A, double cObjNrm1, int *FLAG){
//    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *)A;
//    cholmod_sparse *cholA;
//    cholmod_common common;
//    cholmod_start(&common);
//
//    int n = sparse->nSDPCol;
//    int nnz = sparse->nTriMatElem;
//
//    cholmod_triplet T;
//    T.nrow = n;
//    T.ncol = n;
//    T.nnz = nnz;
//    T.stype = 1;// symmetric
//    T.itype = CHOLMOD_INT;
//    T.xtype = CHOLMOD_REAL;
//
//    ASDP_INIT(T.i, int, nnz);
//    ASDP_INIT(T.j, int, nnz);
//    ASDP_INIT(T.x, double, nnz);
//
//
//    int row; int col;
//    for (int i = 0; i < nnz; ++i){
//        ((int *)T.i)[i] = sparse->triMatRow[i];
//        ((int *)T.j)[i] = sparse->triMatCol[i];
//        row = sparse->triMatRow[i];
//        col = sparse->triMatCol[i];
//        if (row == col){
//            ((double *)T.x)[i] = sparse->triMatElem[i] + ASDP_TERMINATION_TOL * (1 + cObjNrm1);
//        }else{
//            ((double *)T.x)[i] = sparse->triMatElem[i];
//        }
//    }
//
//    cholA = cholmod_triplet_to_sparse(&T, nnz, &common);
//    cholmod_factor *L = cholmod_analyze(cholA, &common);
//    L->is_super = 0;
//    cholmod_factorize(cholA, L, &common);
//
//    if (L->is_super == 0){
//        // simple LDL' or simple LL'
//        double *Lx = (double *)L->x;
//        int *Li = (int *)L->i;
//        int *Lp = (int *)L->p;
//        for (int i = 0; i < n; ++i){
//            for (int j = Lp[i]; j < Lp[i+1]; ++j){
//                if (Li[j] == i){
//                    if (Lx[j] < 0){
//                        // 0 means positive definite
//                        // 1 means not positive definite
//                        *FLAG = 1;
//                        goto exit_cleanup;
//                    }
//                }
//            }
//        }
//    }else if (L->is_super == 1){
//        ;
//        // supernodal LL'
//        // result is dense
//        /*
//        double *Lx = (double *)L->x;
//        int row = 0; int idx = 0; double diagElem;
//        for (int col = 0; col < n; ++col){
//            diagElem = Lx[idx];
//            if (diagElem < 0){
//                *FLAG = 1;
//                goto exit_cleanup;
//            }
//            idx += (n - col);
//        }*/
//    }
//
//exit_cleanup:
//    cholmod_free_factor(&L, &common);
//    cholmod_free_sparse(&cholA, &common);
//    cholmod_finish(&common);
//    ASDP_FREE(T.i);
//    ASDP_FREE(T.j);
//    ASDP_FREE(T.x);
//}
//
//
//
//extern void ASDPPosDefDetermineDense(void *A, double cObjNrm1, int *FLAG){
//    TODO: not adjusted to new dense storage in dsMatElem
//    sdp_coeff_dense *dense = (sdp_coeff_dense *)A;
//    double *lowerTriMat = dense->dsMatElem;
//    int row = 0; int col = 0; int n = dense->nSDPCol;
//    double *fullMat;
//    ASDP_INIT(fullMat, double, dense->nSDPCol * dense->nSDPCol);
//    for (int idx = 0; idx < n * (n + 1)/ 2; ++idx){
//        fullMat[row * n + col] = lowerTriMat[idx];
//        fullMat[col * n + row] = lowerTriMat[idx];
//        if (row == col){
//            fullMat[row *n + col] += ASDP_TERMINATION_TOL * (cObjNrm1 + 1);
//        }
//        if (row == n - 1){
//            col++;
//            row = col;
//        }else{
//            row++;
//        }
//    }
//    char uplo = 'L';
//    int info = 0;
//    dpotrf(&uplo, &n, fullMat, &n, &info);
//
//
//    if (info != 0){
//        // 1 means not positive definite
//        *FLAG = 1;
//    }else{
//        ;
//        // 0 means positive definite
//        // no change
//        // *FLAG = 0;
//    }
//    ASDP_FREE(fullMat);
//}

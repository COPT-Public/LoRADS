#ifndef dense_opts_h
#define dense_opts_h
#include "def_lorads_sdp_data.h"



#ifdef __cplusplus
extern "C" {
#endif

extern void fds_symv( lorads_int n, double alpha, double *A, double *x, double beta, double *y );

extern void fds_symv_L(lorads_int n, double alpha, double *A, double *x, double beta, double *y);
extern void fds_syev_Min_eig (sdp_coeff *slackVar, double *minEig);
extern void sps_syev_Min_eig(sdp_coeff *slackVar, double *minEig);
extern lorads_int fds_syev( lorads_int n, double *U, double *d, double *Y,
                     double *work, lorads_int *iwork, lorads_int lwork, lorads_int liwork );
extern void fds_gemv( lorads_int m, lorads_int n, double *M, double *v, double *y );
extern void fds_gemvTran( lorads_int m, lorads_int n, double *M, double *v, double *y );
extern void fds_gemv_beta( lorads_int m, lorads_int n, double *M, double *v, double *y, double beta);
extern void fds_gemv_Trans( lorads_int m, lorads_int n, double *M, double *v, double *y );
extern void fds_symm( char side, char uplo, lorads_int m, lorads_int n, double alpha, double *a, lorads_int lda,
                      double *b, lorads_int ldb, double beta, double *c, lorads_int ldc );
extern void fds_ger( lorads_int m, lorads_int n, double alpha, double *x, lorads_int incx,
                     double *y, lorads_int incy, double *a, lorads_int lda);
extern void fds_print( lorads_int n, double *A );
extern void pds_scal( double a, lorads_int n, double *A );
extern double pds_sum_abs( lorads_int n, double *A );
extern double pds_fro_norm( lorads_int n, double *A );
extern void pds_dump( lorads_int n, double *A, double *v );
extern void pds_decompress( lorads_int nnz, lorads_int *Ci, double *Cx, double *A );
extern lorads_int pds_r1_extract( lorads_int n, double *A, double *sgn, double *a );
extern double pds_quadform( lorads_int n, double *A, double *v, double *aux );
extern void pds_spmv( char uplo, lorads_int n, double alpha, double *ap, double *x, lorads_int incx,
                      double beta, double *y, lorads_int incy );
extern void pds_syr( char uplo, lorads_int n, double alpha,
                     double *x, lorads_int incx, double *a, lorads_int lda );

extern void fds_syr2k(char uplo, char trans, lorads_int n, lorads_int k, double alpha, double *a, double *b, double beta, double *c);
extern void gemm_syr2k(lorads_int m, lorads_int n, lorads_int k, double *U, double *V, double *fullMat);

#ifdef UNDER_BLAS
extern void dpotrf_(char *uplo, lorads_int *n, double *A, lorads_int *lda, lorads_int *info);
extern void dgemv_( char *trans, lorads_int *m, lorads_int *n, double *alpha,
                   double *a, lorads_int *lda, double *x, lorads_int *incx,
                   double *beta, double *y, lorads_int *incy );
extern void dsymm_( const char *side, const char *uplo,  lorads_int *m,
                    lorads_int *n, const double *alpha, const double *a,
                    lorads_int *lda, const double *b,  lorads_int *ldb,
                   const double *beta, double *c,  lorads_int *ldc );
extern void dgemm_(char *transa, char *transb, lorads_int *m, lorads_int *n, lorads_int *k,
                  double *alpha, double *a, lorads_int *lda, double *b, lorads_int *ldb,
                  double *beta, double *c, lorads_int *ldc);
#else
extern void dpotrf(char *uplo, lorads_int *n, double *A, lorads_int *lda, lorads_int *info);
extern void dgemv( char *trans, lorads_int *m, lorads_int *n, double *alpha,
                   double *a, lorads_int *lda, double *x, lorads_int *incx,
                   double *beta, double *y, lorads_int *incy );
extern void dsymm( const char *side, const char *uplo,  lorads_int *m,
                    lorads_int *n, const double *alpha, const double *a,
                    lorads_int *lda, const double *b,  lorads_int *ldb,
                   const double *beta, double *c,  lorads_int *ldc );
extern void dgemm(char *transa, char *transb, lorads_int *m, lorads_int *n, lorads_int *k,
                  double *alpha, double *a, lorads_int *lda, double *b, lorads_int *ldb,
                  double *beta, double *c, lorads_int *ldc);
#endif


#ifdef __cplusplus
}
#endif

#endif /* dense_opts_h */

#ifndef dense_opts_h
#define dense_opts_h
#include "asdp_sdpdata.h"
#ifdef __cplusplus
extern "C" {
#endif

extern void fds_symv( int n, double alpha, double *A, double *x, double beta, double *y );

void fds_syev_Min_eig (sdp_coeff *slackVar, double *minEig);
void sps_syev_Min_eig(sdp_coeff *slackVar, double *minEig);
int fds_syev( int n, double *U, double *d, double *Y,
                     double *work, int *iwork, int lwork, int liwork );
extern void fds_gemv( int m, int n, double *M, double *v, double *y );
extern void fds_gemvTran( int m, int n, double *M, double *v, double *y );
extern void fds_gemv_beta( int m, int n, double *M, double *v, double *y, double beta);
extern void fds_gemv_Trans( int m, int n, double *M, double *v, double *y );
extern void fds_symm( char side, char uplo, int m, int n, double alpha, double *a, int lda,
                      double *b, int ldb, double beta, double *c, int ldc );
extern void fds_ger( int m, int n, double alpha, double *x, int incx,
                    double *y, int incy, double *a, int lda);
void fds_print( int n, double *A );
int fds_r1_extract( int n, double *A, double *sgn, double *a );
void fds_syr2k(char uplo, char trans, int n, int k, double alpha, double *a, double *b, double beta, double *c);


#ifdef UNDER_BLAS
extern void dpotrf_(char *uplo, int *n, double *A, int *lda, int *info);
extern void dgemv_( char *trans, int *m, int *n, double *alpha,
                   double *a, int *lda, double *x, int *incx,
                   double *beta, double *y, int *incy );
#else
extern void dpotrf(char *uplo, int *n, double *A, int *lda, int *info);
extern void dgemv( char *trans, int *m, int *n, double *alpha,
                   double *a, int *lda, double *x, int *incx,
                   double *beta, double *y, int *incy );
#endif

#ifdef __cplusplus
}
#endif

#endif /* dense_opts_h */

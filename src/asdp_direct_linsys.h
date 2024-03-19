
#ifndef ASDP_DIRECT_LINSYS_H
#define ASDP_DIRECT_LINSYS_H



#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    char uplo;
    int n;
    int lda;
    int ldb;
    int info;
    int blkDim;
    int rank;
    int nConstr;
} asdp_direct_linsys; 


extern asdp_retcode DirectLinSysSolve(void *linSys, double *a, double *b, double *res);
extern void directLinSysClear(void *linSysPtr);
extern asdp_retcode ASDPDirectLinSysCreate(asdp_direct_linsys **linsys, int blkDim, int rank, int nConstr);

/* Dense direct */
#define LAPACK_RET_OK    ( 0 )
#define LAPACK_UPLOW_LOW ('L')
#define LAPACK_NOTRANS   ('N')
#define LAPACK_TRANS     ('T')
#define LAPACK_SIDE_LEFT ('L')
#define LAPACK_DIAG_NONUNIT ('N')


#ifdef UNDER_BLAS
extern void dtrsv_( const char *uplo, const char *trans, const char *diag,
                   const int *n, const double *a, const int *lda, double *x,
                   const int *incx );
void dtrsm_( const char *side, const char *uplo, const char *transa,
            const char *diag, const int *m, const int *n, const double *alpha,
            const double *a, const int *lda, double *b, const int *ldb );
extern void dpotrf_ (char *uplo, int *n, double *a, int *lda, int *info );
void dsytrf_( const char *uplo, const int  *n, double *a, const int  *lda,
              int *ipiv, double *work, const int *lwork, int *info );
void dpotri_( const char *uplo, const int *n, double *a, const int *lda, int *info );
void dpotrs_( const char *uplo, const int *n, const int *nrhs, const double *a,
             const int *lda, double *b, const int *ldb, int *info );
void dsytrs_( const char *uplo, const int *n, const int *nrhs, const double *a,
               const int *lda, const int *ipiv, double *b, const int *ldb, int *info );
#else

extern void dtrsv( const char *uplo, const char *trans, const char *diag,
                   const int *n, const double *a, const int *lda, double *x,
                   const int *incx );
void dtrsm( const char *side, const char *uplo, const char *transa,
            const char *diag, const int *m, const int *n, const double *alpha,
            const double *a, const int *lda, double *b, const int *ldb );
extern void dpotrf (char *uplo, int *n, double *a, int *lda, int *info );
void dsytrf( const char *uplo, const int  *n, double *a, const int  *lda,
              int *ipiv, double *work, const int *lwork, int *info );
void dpotri( const char *uplo, const int *n, double *a, const int *lda, int *info );
void dpotrs( const char *uplo, const int *n, const int *nrhs, const double *a,
             const int *lda, double *b, const int *ldb, int *info );
void dsytrs( const char *uplo, const int *n, const int *nrhs, const double *a,
               const int *lda, const int *ipiv, double *b, const int *ldb, int *info );
#endif
#ifdef __cplusplus
}
#endif

#endif // ASDP_DIRECT_LINSYS_H

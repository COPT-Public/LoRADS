#ifdef HEADERPATH
#include "src/asdp_utils.h"
#include "src/dense_opts.h"
#include "src/vec_opts.h"
#include "src/asdp_sdpdata.h"
#include "src/asdp_debug.h"
#else
#include "asdp_utils.h"
#include "dense_opts.h"
#include "vec_opts.h"
#include "asdp_sdpdata.h"
#include "asdp_debug.h"
//#include "mkl.h"
#endif

#include <math.h>


#ifdef MEMDEBUG
#include "memwatch.h"
#endif

enum CBLAS_LAYOUT {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};
enum CBLAS_STORAGE {CblasPacked=151};
enum CBLAS_IDENTIFIER {CblasAMatrix=161, CblasBMatrix=162};
enum CBLAS_OFFSET {CblasRowOffset=171, CblasColOffset=172, CblasFixOffset=173};

typedef enum CBLAS_LAYOUT CBLAS_LAYOUT;
typedef enum CBLAS_TRANSPOSE CBLAS_TRANSPOSE;
typedef enum CBLAS_UPLO CBLAS_UPLO;
typedef enum CBLAS_DIAG CBLAS_DIAG;
typedef enum CBLAS_SIDE CBLAS_SIDE;
typedef enum CBLAS_STORAGE CBLAS_STORAGE;
typedef enum CBLAS_IDENTIFIER CBLAS_IDENTIFIER;
typedef enum CBLAS_OFFSET CBLAS_OFFSET;

#ifdef UNDER_BLAS
extern void dsymv_( const char *uplo, const int *n, const double *alpha,
                   const double *a, const int *lda, const double *x,
                   const int *incx, const double *beta, double *y, const int *incy );



extern void dsyevr_( const char *jobz, const char *range, const char *uplo,
                    const int  *n, double *a, const int *lda,
                    const double *vl, const double *vu, const int *il,
                    const int *iu, const double *abstol, int *m,
                    double *w, double *z, const int *ldz, int *isuppz,
                    double *work, const int *lwork, int *iwork, const int *liwork,
                    int *info );

extern void dspmv_( const char *uplo, const int *n, const double *alpha,
                   const double *ap, const double *x, const int *incx,
                   const double *beta, double *y, const int *incy );

extern void dsymm_( const char *side, const char *uplo, const int *m,
                   const int *n, const double *alpha, const double *a,
                   const int *lda, const double *b, const int *ldb,
                   const double *beta, double *c, const int *ldc );

void dsyr_( const char *uplo, const int *n, const double *alpha,
           const double *x, const int *incx, double *a, const int *lda );

void dsyr2k_(const char *uplo, const char *trans, const int *n, const int *k, const double *alpha, const double *a, const int *lda, const double *b, int *ldb, const double *beta,  double *c, const int *ldc);


extern void dger_( const int *m, const int *n, const double *alpha,
                  const double *x, const int *incx, const double *y,
                  const int *incy, double *a, const int *lda );
#else

extern void dsymv( const char *uplo, const int *n, const double *alpha,
                   const double *a, const int *lda, const double *x,
                   const int *incx, const double *beta, double *y, const int *incy );



extern void dsyevr( const char *jobz, const char *range, const char *uplo,
                    const int  *n, double *a, const int *lda,
                    const double *vl, const double *vu, const int *il,
                    const int *iu, const double *abstol, int *m,
                    double *w, double *z, const int *ldz, int *isuppz,
                    double *work, const int *lwork, int *iwork, const int *liwork,
                    int *info );

extern void dspmv( const char *uplo, const int *n, const double *alpha,
                   const double *ap, const double *x, const int *incx,
                   const double *beta, double *y, const int *incy );

extern void dsymm( const char *side, const char *uplo, const int *m,
                   const int *n, const double *alpha, const double *a,
                   const int *lda, const double *b, const int *ldb,
                   const double *beta, double *c, const int *ldc );

void dsyr( const char *uplo, const int *n, const double *alpha,
           const double *x, const int *incx, double *a, const int *lda );

void dsyr2k(const char *uplo, const char *trans, const int *n, const int *k, const double *alpha, const double *a, const int *lda, const double *b, int *ldb, const double *beta,  double *c, const int *ldc);


extern void dger( const int *m, const int *n, const double *alpha,
                  const double *x, const int *incx, const double *y,
                  const int *incy, double *a, const int *lda );
#endif

/* Full dense operations */
void fds_symv( int n, double alpha, double *A, double *x, double beta, double *y ) {
    // symmetry matrix multiplication
    // y = alpha * A * x + beta * y
    // where A is a symmetric matrix in dense format
    // n is the dimension of the matrix
    // x and y are vectors
    // alpha and beta are scalars
    // ACharConstantUploLow UPLOLOW low triangular part of A is referenced
#ifdef UNDER_BLAS
    dsymv_(&ACharConstantUploLow, &n, &alpha, A, &n, x, &AIntConstantOne, &beta, y, &AIntConstantOne);
#else
    dsymv(&ACharConstantUploLow, &n, &alpha, A, &n, x, &AIntConstantOne, &beta, y, &AIntConstantOne);
#endif

    return;
}


extern void reconstructSymmetricMatrix(sdp_coeff *slackVar, double *symmetricMatrix, int n){
    if (slackVar->dataType == SDP_COEFF_DENSE){
        sdp_coeff_dense *slackVarDense = (sdp_coeff_dense *)slackVar->dataMat;
        double *lowerTriangle = slackVarDense->dsMatElem;
        int i = 0;
        for (int col = 0; col < n; ++col) {
            i += col;
            symmetricMatrix[i] = lowerTriangle[i];
            ++i;
            for (int row = col +1; row < n; ++row) {
                symmetricMatrix[n*row + col] = symmetricMatrix[i] = lowerTriangle[i];
                ++i;
            }
        }
    }else if (slackVar->dataType == SDP_COEFF_SPARSE){
        sdp_coeff_sparse *slackVarSparse = (sdp_coeff_sparse *)slackVar->dataMat;
        int row;
        int col;
        for (int i = 0; i < slackVarSparse->nTriMatElem; ++i){
            row = slackVarSparse->triMatRow[i];
            col = slackVarSparse->triMatCol[i];
            symmetricMatrix[col * n + row] = slackVarSparse->triMatElem[i];
            symmetricMatrix[row * n + col] = slackVarSparse->triMatElem[i];
        }
    }
}

extern void fds_syev_Min_eig(sdp_coeff *slackVar, double *minEig){
    int n = slackVar->nSDPCol;
    double *completeSymMatrix;

    ASDP_INIT(completeSymMatrix, double, n * n);
    reconstructSymmetricMatrix(slackVar, completeSymMatrix, n);

    // need to check
    int retcode = ASDP_RETCODE_OK;
    char jobz = 'N'; // eigenvalues only
    char range = 'I'; // index range
    char uplo = ACharConstantUploLow;
    int il = 1;
    int iu = 1; // calculate minimal eigenvalue
    int lda = n, ldz = n;
    int m = 1;
    int info = 0;
    double *w;
    ASDP_INIT(w, double, n);
    double abstol = 1e-8;
    double *z;
    ASDP_INIT(z, double, n * n);
    int *isuppz;
    ASDP_INIT(isuppz, int, 2 * n);


    int lwork = -1;
    int liwork = -1;
    double work_query;
    int iwork_query;

#ifdef UNDER_BLAS
    dsyevr_(&jobz, &range, &uplo, &n, completeSymMatrix, &lda,
           &AblConstantZero, &AblConstantZero,
           &il, &iu, &abstol, &m, w, z,
           &ldz, isuppz, &work_query, &lwork, &iwork_query, &liwork, &info);
#else
    dsyevr(&jobz, &range, &uplo, &n, completeSymMatrix, &lda,
           &AblConstantZero, &AblConstantZero,
           &il, &iu, &abstol, &m, w, z,
           &ldz, isuppz, &work_query, &lwork, &iwork_query, &liwork, &info);
#endif

    lwork = (int)work_query;
    double *work = (double*)malloc(sizeof(double) * lwork);

    liwork = iwork_query;
    int *iwork = (int*)malloc(sizeof(int) * liwork);

#ifdef UNDER_BLAS
    dsyevr_(&jobz, &range, &uplo, &n, completeSymMatrix, &lda,
           &AblConstantZero, &AblConstantZero,
           &il, &iu, &abstol, &m, w, z,
           &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
#else
    dsyevr(&jobz, &range, &uplo, &n, completeSymMatrix, &lda,
           &AblConstantZero, &AblConstantZero,
           &il, &iu, &abstol, &m, w, z,
           &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
#endif


    if ( info != 0 ) {
        asdp_printf("info:%d", info);
        ASDP_ERROR_TRACE;
        retcode = ASDP_RETCODE_FAILED;
    }else{
        minEig[0] = w[0];
    }
    ASDP_FREE(isuppz);
    ASDP_FREE(z);
    ASDP_FREE(w);
    ASDP_FREE(work);
    ASDP_FREE(iwork);
    ASDP_FREE(completeSymMatrix);
}



extern int fds_syev( int n, double *U, double *d, double *Y,
                      double *work, int *iwork, int lwork, int liwork ) {

    int retcode = ASDP_RETCODE_OK;

    char jobz = 'V', range = 'I', uplo = ACharConstantUploUp;
    int isuppz[4] = {0};
    int il = n - 1, iu = n;
    int m = 2;
    int info = 0;

    // compute eigenvalues and eigenvectors of a symmetric matrix
#ifdef UNDER_BLAS
    dsyevr_(&jobz, &range, &uplo, &n, U, &n,
           &AblConstantZero, &AblConstantZero,
           &il, &iu, &AblConstantZero, &m, d, Y,
           &n, isuppz, work, &lwork, iwork, &liwork, &info);
#else
    dsyevr(&jobz, &range, &uplo, &n, U, &n,
           &AblConstantZero, &AblConstantZero,
           &il, &iu, &AblConstantZero, &m, d, Y,
           &n, isuppz, work, &lwork, iwork, &liwork, &info);
#endif

    if ( info != 0 ) {
        retcode = ASDP_RETCODE_FAILED;
    }

    return retcode;
}





extern void fds_gemv( int m, int n, double *M, double *v, double *y ) {
    // y = alpha * A * x + beta * y (beta = 0)

    // trans: 'N'--no transpose, 'T'--transpose, 'C'--conjugate transpose
    // m: number of rows of matrix M
    // n: number of columns of matrix M
    // alpha: scalar
    // a: matrix A
    // lda: leading dimension of a (column major means lda = m)
    // x: vector x
    // incx: the increment for the elements of x
    // beta: scalar
    // y: vector y
    // incy: the increment for the elements of y

#ifdef UNDER_BLAS
    dgemv_(&ACharConstantNoTrans, &m, &n, &AblConstantOne, M, &m,
          v, &AIntConstantOne, &AblConstantZero, y, &AIntConstantOne);
#else
    dgemv(&ACharConstantNoTrans, &m, &n, &AblConstantOne, M, &m,
          v, &AIntConstantOne, &AblConstantZero, y, &AIntConstantOne);
#endif

    return;
}

extern void fds_gemvTran( int m, int n, double *M, double *v, double *y ) {
    // y = alpha * A * x + beta * y (beta = 0)

    // trans: 'N'--no transpose, 'T'--transpose, 'C'--conjugate transpose
    // m: number of rows of matrix M
    // n: number of columns of matrix M
    // alpha: scalar
    // a: matrix A
    // lda: leading dimension of a (column major means lda = m)
    // x: vector x
    // incx: the increment for the elements of x
    // beta: scalar
    // y: vector y
    // incy: the increment for the elements of y

#ifdef UNDER_BLAS
    dgemv_(&ACharConstantTrans, &m, &n, &AblConstantOne, M, &m,
          v, &AIntConstantOne, &AblConstantZero, y, &AIntConstantOne);
#else
    dgemv(&ACharConstantTrans, &m, &n, &AblConstantOne, M, &m,
          v, &AIntConstantOne, &AblConstantZero, y, &AIntConstantOne);
#endif

    return;
}

extern void fds_gemv_beta( int m, int n, double *M, double *v, double *y , double beta) {
    // y = alpha * A * x + beta * y (beta =1.0)

    // trans: 'N'--no transpose, 'T'--transpose, 'C'--conjugate transpose
    // m: number of rows of matrix M
    // n: number of columns of matrix M
    // alpha: scalar
    // a: matrix A
    // lda: leading dimension of a (column major means lda = m)
    // x: vector x
    // incx: the increment for the elements of x
    // beta: scalar
    // y: vector y
    // incy: the increment for the elements of y

#ifdef UNDER_BLAS
    dgemv_(&ACharConstantNoTrans, &m, &n, &AblConstantOne, M, &m,
          v, &AIntConstantOne, &beta, y, &AIntConstantOne);
#else
    dgemv(&ACharConstantNoTrans, &m, &n, &AblConstantOne, M, &m,
          v, &AIntConstantOne, &beta, y, &AIntConstantOne);
#endif

    return;
}

extern void fds_gemv_Trans( int m, int n, double *M, double *v, double *y ) {
    // y = alpha * A * x + beta * y (beta = 0)

    // trans: 'N'--no transpose, 'T'--transpose, 'C'--conjugate transpose
    // m: number of rows of matrix M
    // n: number of columns of matrix M
    // alpha: scalar
    // a: matrix A
    // lda: leading dimension of a (column major means lda = m)
    // x: vector x
    // incx: the increment for the elements of x
    // beta: scalar
    // y: vector y
    // incy: the increment for the elements of y

#ifdef UNDER_BLAS
    dgemv_(&ACharConstantTrans, &m, &n, &AblConstantOne, M, &m,
          v, &AIntConstantOne, &AblConstantZero, y, &AIntConstantOne);
#else
    dgemv(&ACharConstantTrans, &m, &n, &AblConstantOne, M, &m,
          v, &AIntConstantOne, &AblConstantZero, y, &AIntConstantOne);
#endif

    return;
}

extern void fds_symm( char side, char uplo, int m, int n, double alpha, double *a, int lda,
                      double *b, int ldb, double beta, double *c, int ldc ) {
    // matrix multiply another symmetry matrix
    // side: 'L'--left, C = alpha * A * B + beta * C, 'R'--right C = alpha * B * A + beta * C
    // uplo: symmetry matrix 'U'--upper, 'L'--lower
    // m: number of rows of matrix C
    // n: number of columns of matrix C
    // alpha: scalar
    // a: symmetry matrix A, shape is (m, m)
    // lda: leading dimension of a (column major means lda = m)
    // b: symmetry matrix B, shape is (m, n) if side is 'L', shape is (n, m) if side is 'R'
#ifdef UNDER_BLAS
    dsymm_(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#else
    dsymm(&side, &uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
    return;
}

extern void fds_ger( int m, int n, double alpha, double *x, int incx,
                    double *y, int incy, double *a, int lda) {
#ifdef UNDER_BLAS
    dger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
#else
    dger(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
#endif

    return;
}

extern void fds_print( int n, double *A ) {

    for ( int i = 0; i < n; ++i ) {
        for ( int j = 0; j < n; ++j ) {
            printf("%6.3e ", A[i + n * j]);
        }
        printf("\n");
    }

    return;
}

/** @brief Check is the matrix is rank one
 *
 */
extern int fds_r1_extract( int n, double *A, double *sgn, double *a ) {

    /* Get the first nonzero */
    int i, j, k = 0;
    for ( i = 0; i < n; ++i ) {
        if ( A[k] != 0 ) {
            break;
        }
        k += n +1;
    }

    if ( i == n ) {
        return 0;
    }

    double s = ( A[k] > 0 ) ? 1.0 : -1.0;
    double v = sqrt(fabs(A[k]));
    double eps = 0.0;

    /* Extract diagonal */
    for ( k = 0; k < n; ++k ) {
        a[k] = A[n*k + i] / v;
    }

    if ( s == 1.0 ) {
        for ( i = 0; i < n; ++i ) {
            for ( j = 0; j < n - i; ++j ) {
                eps += fabs(A[j] - a[i] * a[i + j]);
            }
            A += n +1;
            if ( eps > 1e-10 ) {
                return 0;
            }
        }
    } else {
        for ( i = 0; i < n; ++i ) {
            for ( j = 0; j < n - i; ++j ) {
                eps += fabs(A[j] + a[i] * a[i + j]);
            }
            A += n +1;
            if ( eps > 1e-10 ) {
                return 0;
            }
        }
    }

    *sgn = s;
    return 1;
}

extern void fds_syr2k(char uplo, char trans, int n, int k, double alpha, double *a, double *b, double beta, double *c){
    int lda = n;
    int ldb = n;
    int ldc = n;
#ifdef UNDER_BLAS
    // C := alpha*A**T*B + alpha*B**T*A + beta*C,
    dsyr2k_(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#else
    dsyr2k(&uplo, &trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
}


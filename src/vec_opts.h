#ifndef vec_opts_h
#define vec_opts_h

#define TRANS   ('T')
#define NOTRANS ('N')
#define UPLOLOW ('L')
#define UPLOUP ('U')
#define SIDELEFT ('L')
#define SIDERIGHT ('R')

/* Constants to be used when calling BLAS. Should never be modified */
static char ACharConstantTrans = TRANS;
static char ACharConstantNoTrans = NOTRANS;
static char ACharConstantUploUp = UPLOUP;
static char ACharConstantUploLow = UPLOLOW;
static int  AIntConstantOne = 1;
static double AblConstantZero = 0.0;
static double AblConstantOne = 1.0;
static double AblConstantMinusOne = -1.0;

#ifdef __cplusplus
extern "C" {
#endif

extern double nrm1( int *n, double *x, int *incx );
extern double nrm2( int *n, double *x, int *incx );
extern void axpy( int *n, double *alpha, double *x, int *incx, double *y, int *incy );
extern void axpby( int *n, double *a, double *x, int *incx, double *b, double *y, int *incy );
extern void axpbyAddition( int *n, double *a, double *x, double *b, double *y, double *z);
extern double dot( int *n, double *x, int *incx, double *y, int *incy );
extern void scal( int *n, double *sa, double *sx, int *incx );
extern void rscl( int *n, double *sa, double *sx, int *incx );
extern void syr( char *uplo, int *n, double *alpha, double *x, int *incx, double *a, int *lda );
extern void symv( char *uplo, int *n, double *alpha, double *a, int *lda, double *x,
                  int *incx, double *beta, double *y, const int *incy );
extern int idamax( int *n, double *x, int *incx );
extern int idamin( int *n, double *x, int *incx );
extern double sumlogdet( int *n, double *x );
extern void vvscl( int *n, double *s, double *x );
extern void vvrscl( int *n, double *s, double *x );
extern double normalize( int *n, double *a );

#ifdef __cplusplus
}
#endif

#endif /* vec_opts_h */

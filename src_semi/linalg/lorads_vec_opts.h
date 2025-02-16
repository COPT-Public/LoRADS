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
static lorads_int  AIntConstantOne = 1;
static double AblConstantZero = 0.0;
static double AblConstantOne = 1.0;
static double AblConstantMinusOne = -1.0;

#ifdef __cplusplus
extern "C" {
#endif

extern double nrm1( lorads_int *n, double *x, lorads_int *incx );
extern double nrm2( lorads_int *n, double *x, lorads_int *incx );
extern void axpy( lorads_int *n, double *alpha, double *x, lorads_int *incx, double *y, lorads_int *incy );
extern void axpby( lorads_int *n, double *a, double *x, lorads_int *incx, double *b, double *y, lorads_int *incy );
extern void axpbyAddition( lorads_int *n, double *a, double *x, double *b, double *y, double *z);
extern double dot( lorads_int *n, double *x, lorads_int *incx, double *y, lorads_int *incy );
extern void scal( lorads_int *n, double *sa, double *sx, lorads_int *incx );
extern void rscl( lorads_int *n, double *sa, double *sx, lorads_int *incx );
extern void syr( char *uplo, lorads_int *n, double *alpha, double *x, lorads_int *incx, double *a, lorads_int *lda );
extern void symv( char *uplo, lorads_int *n, double *alpha, double *a, lorads_int *lda, double *x,
                  lorads_int *incx, double *beta, double *y, const lorads_int *incy );
extern double sumlogdet( lorads_int *n, double *x );
extern void vvscl( lorads_int *n, double *s, double *x );
extern void vvrscl( lorads_int *n, double *s, double *x );
extern double normalize( lorads_int *n, double *a );

#ifdef __cplusplus
}
#endif

#endif /* vec_opts_h */

#include "lorads_utils.h"

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

#include <math.h>
/* Blas functions */

#ifdef UNDER_BLAS
// ||x|| = sqrt(x[0]*x[0] + x[incx]*x[incx] + x[2*incx]*x[2*incx] + ...)
extern double dnrm2_( lorads_int *n, double *x, lorads_int *incx );

// y[i] = alpha * x[i] + y[i]，where i from 0 to n-1, incx = 1, incy = 1 generally
extern void daxpy_( lorads_int *n, double *alpha, double *x, lorads_int *incx, double *y, lorads_int *incy );

//result = x[0]*y[0] + x[incx]*y[incy] + x[2*incx]*y[2*incy] + ...
extern double ddot_( lorads_int *n, double *x, lorads_int *incx, double *y, lorads_int *incy );

//sx[i] = sa * sx[i]，where i from 0 to n-1
extern void dscal_( lorads_int *n, double *sa, double *sx, lorads_int *incx );

// sx[i] = sx[i] / sa，where i from 0 to n-1
extern void drscl_( lorads_int *n, double *sa, double *sx, lorads_int *incx );

// A = alpha * x * x^T + A
extern void dsyr_( char *uplo, lorads_int *n, double *alpha, double *x, lorads_int *incx, double *a, lorads_int *lda );

// idamax The index used to find the element of the vector with the largest absolute value. It returns the index value of the element in the vector with the largest absolute value, not the actual element value.
extern lorads_int idamax_( lorads_int *n, double *x, lorads_int *incx );

extern lorads_int idamin_( lorads_int *n, double *x, lorads_int *incx );
#else
// ||x|| = sqrt(x[0]*x[0] + x[incx]*x[incx] + x[2*incx]*x[2*incx] + ...)
extern double dnrm2( lorads_int *n, double *x, lorads_int *incx );

// y[i] = alpha * x[i] + y[i]，where i from 0 to n-1, incx = 1, incy = 1 generally
extern void daxpy( lorads_int *n, double *alpha, double *x, lorads_int *incx, double *y, lorads_int *incy );

//result = x[0]*y[0] + x[incx]*y[incy] + x[2*incx]*y[2*incy] + ...
extern double ddot( lorads_int *n, double *x, lorads_int *incx, double *y, lorads_int *incy );

//sx[i] = sa * sx[i]，where i from 0 to n-1
extern void dscal( lorads_int *n, double *sa, double *sx, lorads_int *incx );

// sx[i] = sx[i] / sa，where i from 0 to n-1
extern void drscl( lorads_int *n, double *sa, double *sx, lorads_int *incx );

// A = alpha * x * x^T + A
extern void dsyr( char *uplo, lorads_int *n, double *alpha, double *x, lorads_int *incx, double *a, lorads_int *lda );

// idamax The index used to find the element of the vector with the largest absolute value. It returns the index value of the element in the vector with the largest absolute value, not the actual element value.
extern lorads_int idamax( lorads_int *n, double *x, lorads_int *incx );

extern lorads_int idamin( lorads_int *n, double *x, lorads_int *incx );
#endif

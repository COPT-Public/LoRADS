#ifdef HEADERPATH
#include "interface/asdp_utils.h"
#include "src/vec_opts.h"
#include "src/asdp_debug.h"
#else
#include "asdp_utils.h"
#include "vec_opts.h"
#include "asdp_debug.h"
#endif

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

#include <math.h>
/* Blas functions */

#ifdef UNDER_BLAS
// ||x|| = sqrt(x[0]*x[0] + x[incx]*x[incx] + x[2*incx]*x[2*incx] + ...)
extern double dnrm2_( int *n, double *x, int *incx );

// y[i] = alpha * x[i] + y[i]，where i from 0 to n-1, incx = 1, incy = 1 generally
extern void daxpy_( int *n, double *alpha, double *x, int *incx, double *y, int *incy );

//result = x[0]*y[0] + x[incx]*y[incy] + x[2*incx]*y[2*incy] + ...
extern double ddot_( int *n, double *x, int *incx, double *y, int *incy );

//sx[i] = sa * sx[i]，where i from 0 to n-1
extern void dscal_( int *n, double *sa, double *sx, int *incx );

// sx[i] = sx[i] / sa，where i from 0 to n-1
extern void drscl_( int *n, double *sa, double *sx, int *incx );

// A = alpha * x * x^T + A
extern void dsyr_( char *uplo, int *n, double *alpha, double *x, int *incx, double *a, int *lda );

// idamax The index used to find the element of the vector with the largest absolute value. It returns the index value of the element in the vector with the largest absolute value, not the actual element value.
extern int idamax_( int *n, double *x, int *incx );

extern int idamin_( int *n, double *x, int *incx );
#else
// ||x|| = sqrt(x[0]*x[0] + x[incx]*x[incx] + x[2*incx]*x[2*incx] + ...)
extern double dnrm2( int *n, double *x, int *incx ); 

// y[i] = alpha * x[i] + y[i]，where i from 0 to n-1, incx = 1, incy = 1 generally
extern void daxpy( int *n, double *alpha, double *x, int *incx, double *y, int *incy );

//result = x[0]*y[0] + x[incx]*y[incy] + x[2*incx]*y[2*incy] + ...
extern double ddot( int *n, double *x, int *incx, double *y, int *incy );

//sx[i] = sa * sx[i]，where i from 0 to n-1
extern void dscal( int *n, double *sa, double *sx, int *incx );

// sx[i] = sx[i] / sa，where i from 0 to n-1
extern void drscl( int *n, double *sa, double *sx, int *incx );

// A = alpha * x * x^T + A
extern void dsyr( char *uplo, int *n, double *alpha, double *x, int *incx, double *a, int *lda );

// idamax The index used to find the element of the vector with the largest absolute value. It returns the index value of the element in the vector with the largest absolute value, not the actual element value.
extern int idamax( int *n, double *x, int *incx );

extern int idamin( int *n, double *x, int *incx );
#endif

extern double nrm2( int *n, double *x, int *incx ) {
// ||x|| = sqrt(x[0]*x[0] + x[incx]*x[incx] + x[2*incx]*x[2*incx] + ...)
#ifdef MYBLAS
    assert( *incx == 1 );
    
    double nrm = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        nrm += x[i] * x[i];
    }
    
    return sqrt(nrm);
#else
#ifdef UNDER_BLAS
    return dnrm2_(n, x, incx);
#else
    return dnrm2(n, x, incx);
#endif
#endif
}

extern void axpy( int *n, double *alpha, double *x, int *incx, double *y, int *incy ) {
// y[i] = alpha * x[i] + y[i]，where i from 0 to n-1, incx = 1, incy = 1 generally
#ifdef MYBLAS
    assert( *incx == 1 && *incy == 1 );
    
    for ( int i = 0; i < *n; ++i ) {
        y[i] += (*alpha) * x[i];
    }
#else
#ifdef UNDER_BLAS
    daxpy_(n, alpha, x, incx, y, incy);
#else
    daxpy(n, alpha, x, incx, y, incy);
#endif

#endif
    return;
}

extern void axpby( int *n, double *a, double *x, int *incx, double *b, double *y, int *incy ) {
// y[i] = a * x[i] + b * y[i]
    double aval = *a;
    double bval = *b;
    
    for ( int i = 0; i < *n; ++i ) {
        y[i] = aval * x[i] + bval * y[i];
    }
    
    return;
}

extern void axpbyAddition( int *n, double *a, double *x, double *b, double *y, double *z) {
    // z[i] = a * x[i] + b * y[i]
#ifdef MYBLAS
    double aval = *a;
    double bval = *b;
    int nval = *n;
    for ( int i = 0; i < nval; ++i ) {
        z[i] = aval * x[i] + bval * y[i];
    }
#else
    int incx = 1;
    ASDP_ZERO(z, double, n[0]);
    axpy(n, a, x, &incx, z, &incx);
    axpy(n, b, y, &incx, z, &incx);
#endif
    
}

extern double dot( int *n, double *x, int *incx, double *y, int *incy ) {
//result = x[0]*y[0] + x[incx]*y[incy] + x[2*incx]*y[2*incy] + ...
#ifdef MYBLAS
    
    double dres = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        dres += x[i * incx[0]] * y[i * incy[0]];
    }
    
    return dres;
#else
#ifdef UNDER_BLAS
    return ddot_(n, x, incx, y, incy);
#else
    return ddot(n, x, incx, y, incy);
#endif
#endif
}

extern void scal( int *n, double *sa, double *sx, int *incx ) {
//sx[i] = sa * sx[i]，where i from 0 to n-1
#ifdef MYBLAS
    assert( *incx == 1 );
    double a = *sa;
    
    if ( a == 1.0 ) {
        return;
    }
    
    for ( int i = 0; i < *n; ++i ) {
        sx[i] = sx[i] * a;
    }
#else
#ifdef UNDER_BLAS
    dscal_(n, sa, sx, incx);
#else
    dscal(n, sa, sx, incx);
#endif
#endif
    return;
}

/* Use standard Blas for this sensitive operation */
extern void rscl( int *n, double *sa, double *sx, int *incx ) {
// sx[i] = sx[i] / sa，where i from 0 to n-1
#if 0
    assert( *incx == 1 );
    double a = *sa;
    
    assert( a != 0.0 );
    assert( a > 0.0 );
    
    if ( a == 1.0 ) {
        return;
    }
    
    if ( fabs(a) < 1e-16 ) {
        a = (a > 0) ? 1e-16 : -1e-16;
    }
    
    for ( int i = 0; i < *n; ++i ) {
        sx[i] = sx[i] / a;
    }
#else
#ifdef UNDER_BLAS
    drscl_(n, sa, sx, incx);
#else
    drscl(n, sa, sx, incx);
#endif
#endif
    return;
}

extern void syr( char *uplo, int *n, double *alpha, double *x, int *incx, double *a, int *lda ) {
    // A = alpha * x * x^T + A
#ifdef UNDER_BLAS
    dsyr_(uplo, n, alpha, x, incx, a, lda);
#else
    dsyr(uplo, n, alpha, x, incx, a, lda);
#endif
    return;
}

extern int idamax( int *n, double *x, int *incx ) {
    
    int idmax = 0;
    double damax = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        double ax = fabs(x[i]);
        if ( ax > damax ) {
            damax = ax; idmax = i;
        }
    }
    
    return idmax;
}

extern int idamin( int *n, double *x, int *incx ) {
    
    int idmin = 0;
    double damin = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        double ax = fabs(x[i]);
        if ( ax < damin ) {
            damin = ax; idmin = i;
        }
    }
    
    return idmin;
}

extern void vvscl( int *n, double *s, double *x ) {
    // x[i] = x[i] * s[i]
    for ( int i = 0; i < *n; ++i ) {
        x[i] = x[i] * s[i];
    }
    
    return;
}

extern void vvrscl( int *n, double *s, double *x ) {
    // x[i] = x[i] / s[i];
    for ( int i = 0; i < *n; ++i ) {
        x[i] = x[i] / s[i];
    }
    
    return;
}

extern double nrm1( int *n, double *x, int *incx ) {
    // sum( abs(x[i]) )
    assert( *incx == 1 );
    
    double nrm = 0.0;
    
    for ( int i = 0; i < *n; ++i ) {
        nrm += fabs(x[i]);
    }
    
    return nrm;
}

extern double normalize( int *n, double *a ) {
    // a[i] = a[i] / norm(a)
    double norm = nrm2(n, a, &AIntConstantOne);
    
    if ( norm > 1e-16 ) {
#ifdef UNDER_BLAS
        drscl_(n, &norm, a, &AIntConstantOne);
#else
        drscl(n, &norm, a, &AIntConstantOne);
#endif
    } else {
        norm = 0.0;
        ASDP_ZERO(a, double, *n);
    }
    
    return norm;
}

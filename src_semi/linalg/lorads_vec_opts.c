#include "lorads_utils.h"
#include "lorads.h"
#include "lorads_vec_opts.h"
#include <math.h>

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

extern double dsr1_sum_abs( lorads_int n, double sign, double *factor ) {
    // return abs(sign) * ||factor||_1^2
    double nrm = 0.0;
    lorads_int incx = 1;

    nrm = nrm1(&n, factor, &incx);

    return nrm * nrm * fabs(sign);
}

extern double dsr1_fro_norm( lorads_int n, double sign, double *factor ) {
    // return abs(sign) * ||factor||_2^2
    double nrm = 0.0;
    lorads_int incx = 1;

    nrm = nrm2(&n, factor, &incx);

    return nrm * nrm * fabs(sign);
}

extern void dsr1_dump( lorads_int n, double sign, double *factor, double *v ) {
    // transform vector `factor` to matrix, and put it lorads_into `v`
    char uplolow = UPLOLOW;
    lorads_int incx = 1;
    syr(&uplolow, &n, &sign, factor, &incx, v, &n);
    LUtilMatSymmetrize(n, v);

    return;
}

extern double dsr1_quadform( lorads_int n, double sign, double *factor, double *v ) {
    // trace product of `factor` and `v` (Two rank one matrices)
    double quadform = dot(&n, factor, &AIntConstantOne, v, &AIntConstantOne);
    return sign * quadform * quadform;
}

extern double spr1_sum_abs( double sign, lorads_int nnz, double *factornz ) {

    double nrm = 0.0;
    lorads_int incx = 1;

    nrm = nrm1(&nnz, factornz, &incx);

    return nrm * nrm * fabs(sign);
}

extern void spr1_dump( lorads_int n, double sign, lorads_int nnz, lorads_int *nzidx, double *factornz ) {

    return;
}

extern double spr1_quadform( lorads_int n, double sign, lorads_int nnz, lorads_int *nzidx, double *factornzs, double *v ) {

    double quadform = 0.0;
    for ( lorads_int i = 0; i < nnz; ++i ) {
        quadform += factornzs[i] * v[nzidx[i]];
    }

    return sign * quadform * quadform;
}


extern void spr1_mat_mul(double sign, lorads_int nnz, lorads_int *nzidx, double *factornzs, double *v, double *w) {
    // w += sign * factor * factor' * v
    // factor is a sparse vector
    // v is a dense vector
    // w is a sparse vector

    double inner = 0.0;
    for ( lorads_int i = 0; i < nnz; ++i ) {
        // factor' * v
        inner += factornzs[i] * v[nzidx[i]];
    }

//    for ( lorads_int i = 0; i < nnz; ++i ) {
//        // w[nzidx[i]] += sign * inner * factornzs[i];
//        w[i] += sign * inner * factornzs[i];
//    }
    lorads_int incx = 1;
    double alpha = sign * inner;
    axpy(&nnz, &alpha, factornzs, &incx, w, &incx);

}

extern void dr1_mat_mul(lorads_int n, double sign, double *factornzs, double *v, double *w){
    // w += sign * factor * factor' * v
    // factor is a dense vector
    // v is a dense vector
    // w is a dense vector

    double inner = 0.0;
    lorads_int incx = 1;
//    for ( lorads_int i = 0; i < n; ++i ) {
//        inner += factornzs[i] * v[i];
//    }
    inner += dot(&n, factornzs, &incx, v, &incx);

//    for ( lorads_int i = 0; i < n; ++i ) {
//        w[i] += sign * inner * factornzs[i];
//    }
    double alpha = sign * inner;
    axpy(&n, &alpha, factornzs, &incx, w, &incx);
}



extern double nrm2( lorads_int *n, double *x, lorads_int *incx ) {
// ||x|| = sqrt(x[0]*x[0] + x[incx]*x[incx] + x[2*incx]*x[2*incx] + ...)
#ifdef MYBLAS
    assert( *incx == 1 );
    
    double nrm = 0.0;
    
    for ( lorads_int i = 0; i < *n; ++i ) {
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

extern void axpy( lorads_int *n, double *alpha, double *x, lorads_int *incx, double *y, lorads_int *incy ) {
// y[i] = alpha * x[i] + y[i]，where i from 0 to n-1, incx = 1, incy = 1 generally
#ifdef MYBLAS
    assert( *incx == 1 && *incy == 1 );
    
    for ( lorads_int i = 0; i < *n; ++i ) {
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

extern void axpby( lorads_int *n, double *a, double *x, lorads_int *incx, double *b, double *y, lorads_int *incy ) {
// y[i] = a * x[i] + b * y[i]
    double aval = *a;
    double bval = *b;

    for ( lorads_int i = 0; i < *n; ++i ) {
        y[i] = aval * x[i] + bval * y[i];
    }

    return;
}

extern void axpbyAddition( lorads_int *n, double *a, double *x, double *b, double *y, double *z) {
    // z[i] = a * x[i] + b * y[i]
#ifdef MYBLAS
    double aval = *a;
    double bval = *b;
    lorads_int nval = *n;
    for ( lorads_int i = 0; i < nval; ++i ) {
        z[i] = aval * x[i] + bval * y[i];
    }
#else
    lorads_int incx = 1;
    LORADS_ZERO(z, double, n[0]);
    axpy(n, a, x, &incx, z, &incx);
    axpy(n, b, y, &incx, z, &incx);
#endif

}

extern double dot( lorads_int *n, double *x, lorads_int *incx, double *y, lorads_int *incy ) {
//result = x[0]*y[0] + x[incx]*y[incy] + x[2*incx]*y[2*incy] + ...
#ifdef MYBLAS

    double dres = 0.0;
    
    for ( lorads_int i = 0; i < *n; ++i ) {
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

extern void scal( lorads_int *n, double *sa, double *sx, lorads_int *incx ) {
//sx[i] = sa * sx[i]，where i from 0 to n-1
#ifdef MYBLAS
    assert( *incx == 1 );
    double a = *sa;
    
    if ( a == 1.0 ) {
        return;
    }
    
    for ( lorads_int i = 0; i < *n; ++i ) {
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
extern void rscl( lorads_int *n, double *sa, double *sx, lorads_int *incx ) {
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
    
    for ( lorads_int i = 0; i < *n; ++i ) {
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

extern void syr( char *uplo, lorads_int *n, double *alpha, double *x, lorads_int *incx, double *a, lorads_int *lda ) {
    // A = alpha * x * x^T + A
#ifdef UNDER_BLAS
    dsyr_(uplo, n, alpha, x, incx, a, lda);
#else
    dsyr(uplo, n, alpha, x, incx, a, lda);
#endif
    return;
}

extern lorads_int idamax( lorads_int *n, double *x, lorads_int *incx ) {

    lorads_int idmax = 0;
    double damax = 0.0;

    for ( lorads_int i = 0; i < *n; ++i ) {
        double ax = fabs(x[i]);
        if ( ax > damax ) {
            damax = ax; idmax = i;
        }
    }

    return idmax;
}

extern lorads_int idamin( lorads_int *n, double *x, lorads_int *incx ) {

    lorads_int idmin = 0;
    double damin = 0.0;

    for ( lorads_int i = 0; i < *n; ++i ) {
        double ax = fabs(x[i]);
        if ( ax < damin ) {
            damin = ax; idmin = i;
        }
    }

    return idmin;
}

extern void vvscl( lorads_int *n, double *s, double *x ) {
    // x[i] = x[i] * s[i]
    for ( lorads_int i = 0; i < *n; ++i ) {
        x[i] = x[i] * s[i];
    }

    return;
}

extern void vvrscl( lorads_int *n, double *s, double *x ) {
    // x[i] = x[i] / s[i];
    for ( lorads_int i = 0; i < *n; ++i ) {
        x[i] = x[i] / s[i];
    }

    return;
}

extern double nrm1( lorads_int *n, double *x, lorads_int *incx ) {
    // sum( abs(x[i]) )
    assert( *incx == 1 );

    double nrm = 0.0;

    for ( lorads_int i = 0; i < *n; ++i ) {
        nrm += fabs(x[i]);
    }

    return nrm;
}

extern double normalize( lorads_int *n, double *a ) {
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
        LORADS_ZERO(a, double, *n);
    }

    return norm;
}


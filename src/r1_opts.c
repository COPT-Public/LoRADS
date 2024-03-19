#ifdef HEADERPATH
#include "interface/asdp_utils.h"
#include "src/r1_opts.h"
#include "src/vec_opts.h"
#include "src/asdp_debug.h"
#else
#include "asdp_utils.h"
#include "r1_opts.h"
#include "vec_opts.h"
#include "asdp_debug.h"
#endif

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

#include <math.h>
extern double dsr1_sum_abs( int n, double sign, double *factor ) {
    // return abs(sign) * ||factor||_1^2
    double nrm = 0.0;
    int incx = 1;
    
    nrm = nrm1(&n, factor, &incx);
    
    return nrm * nrm * fabs(sign);
}

extern double dsr1_fro_norm( int n, double sign, double *factor ) {
    // return abs(sign) * ||factor||_2^2
    double nrm = 0.0;
    int incx = 1;
    
    nrm = nrm2(&n, factor, &incx);
    
    return nrm * nrm * fabs(sign);
}

extern void dsr1_dump( int n, double sign, double *factor, double *v ) {
    // transform vector `factor` to matrix, and put it into `v`
    char uplolow = UPLOLOW;
    int incx = 1;
    syr(&uplolow, &n, &sign, factor, &incx, v, &n);
    AUtilMatSymmetrize(n, v);
    
    return;
}

extern double dsr1_quadform( int n, double sign, double *factor, double *v ) {
    // trace product of `factor` and `v` (Two rank one matrices)
    double quadform = dot(&n, factor, &AIntConstantOne, v, &AIntConstantOne);
    return sign * quadform * quadform;
}

extern double spr1_sum_abs( double sign, int nnz, double *factornz ) {
    
    double nrm = 0.0;
    int incx = 1;
    
    nrm = nrm1(&nnz, factornz, &incx);
    
    return nrm * nrm * fabs(sign);
}

extern void spr1_dump( int n, double sign, int nnz, int *nzidx, double *factornz ) {
    
    return;
}

extern double spr1_quadform( int n, double sign, int nnz, int *nzidx, double *factornzs, double *v ) {
    
    double quadform = 0.0;
    for ( int i = 0; i < nnz; ++i ) {
        quadform += factornzs[i] * v[nzidx[i]];
    }
    
    return sign * quadform * quadform;
}


extern void spr1_mat_mul(double sign, int nnz, int *nzidx, double *factornzs, double *v, double *w) {
    // w += sign * factor * factor' * v
    // factor is a sparse vector
    // v is a dense vector
    // w is a sparse vector

    double inner = 0.0;
    for ( int i = 0; i < nnz; ++i ) {
        // factor' * v
        inner += factornzs[i] * v[nzidx[i]];
    }

//    for ( int i = 0; i < nnz; ++i ) {
//        // w[nzidx[i]] += sign * inner * factornzs[i];
//        w[i] += sign * inner * factornzs[i];
//    }
    int incx = 1;
    double alpha = sign * inner;
    axpy(&nnz, &alpha, factornzs, &incx, w, &incx);
    
}

extern void dr1_mat_mul(int n, double sign, double *factornzs, double *v, double *w){
    // w += sign * factor * factor' * v
    // factor is a dense vector
    // v is a dense vector
    // w is a dense vector

    double inner = 0.0;
    int incx = 1;
//    for ( int i = 0; i < n; ++i ) {
//        inner += factornzs[i] * v[i];
//    }
    inner += dot(&n, factornzs, &incx, v, &incx);

//    for ( int i = 0; i < n; ++i ) {
//        w[i] += sign * inner * factornzs[i];
//    }
    double alpha = sign * inner;
    axpy(&n, &alpha, factornzs, &incx, w, &incx);
}

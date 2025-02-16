/* ========================================================================== */
/* CXSparse/Include/cs.h file */
/* ========================================================================== */

/* This is the CXSparse/Include/cs.h file.  It has the same name (cs.h) as
   the CSparse/Include/cs.h file.  The 'make install' for SuiteSparse installs
   CXSparse, and this file, instead of CSparse.  The two packages have the same
   cs.h include filename, because CXSparse is a superset of CSparse.  Any user
   program that uses CSparse can rely on CXSparse instead, with no change to the
   user code.  The #include "cs.h" line will work for both versions, in user
   code, and the function names and user-visible typedefs from CSparse all
   appear in CXSparse.  For experimenting and changing the package itself, I
   recommend using CSparse since it's simpler and easier to modify.  For
   using the package in production codes, I recommend CXSparse since it has
   more features (support for complex matrices, and both lorads_int and long
   versions).
 */

/* ========================================================================== */

#ifndef _LORADS_CS_
#define _LORADS_CS_

#ifdef __cplusplus
extern "C" {
#endif
#include "lorads.h"

typedef struct lorads_cs_sparse {
    lorads_int nzmax;
    lorads_int m;         /* number of rows */
    lorads_int n;         /* number of columns */
    lorads_int *p;        /* column pointers (size n+1) or col indices (size nzmax) */
    lorads_int *i;        /* row indices, size nzmax */
    double *x;     /* numerical values, size nzmax */
    lorads_int nz;        /* # of entries in triplet matrix, -1 for compressed-col */
} dcs ; //dcs alias

int dcs_entry (dcs *T, lorads_int i, lorads_int j, double x) ;
dcs *dcs_compress (const dcs *T) ;
double dcs_norm (const dcs *A) ;
int dcs_print (const dcs *A, int brief) ;

/* utilities */
void *dcs_calloc (lorads_int n, size_t size) ;
void *dcs_free (void *p) ;
void *dcs_realloc (void *p, lorads_int n, size_t size, int *ok) ;
dcs *dcs_spalloc (lorads_int m, lorads_int n, lorads_int nzmax, lorads_int values, lorads_int t) ;
dcs *dcs_spfree (dcs *A) ;
int dcs_sprealloc (dcs *A, lorads_int nzmax) ;
void *dcs_malloc (lorads_int n, size_t size) ;

/* utilities */
double dcs_cumsum (lorads_int *p, lorads_int *c, lorads_int n) ;
dcs *dcs_done (dcs *C, void *w, void *x, int ok) ;
lorads_int *dcs_idone (lorads_int *p, dcs *C, void *w, int ok) ;

#define IS_CSC(A) (A && (A->nz == -1))
#define IS_TRIPLET(A) (A && (A->nz >= 0))

#ifdef __cplusplus
}
#endif

#endif
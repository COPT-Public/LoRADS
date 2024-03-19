#ifndef DEF_ASDP_LANCZOS_H
#define DEF_ASDP_LANCZOS_H

#ifdef HEADERPATH
#include "interface/hdsdp.h"
#else
#include "asdp.h"
#endif

/* Define the Lanczos data structure for HDSDP */
typedef struct {
    
    /* Matrix dimension */
    int nCol;
    
    /* Maximum subspace dimension*/
    int nMaxSpaceDim;
    
    /* Matrix data */
    void *MMat;
    
    /* Auxiliary data */
    double *vVec;
    double *wVec;
    double *z1Vec;
    double *z2Vec;
    double *vaVec;
    
    double *VMat;
    double *HMat;
    double *YMat;
    double *UMat;
    
    double *dLanczosWarmStart;
    double *dArray;
    double *eigDblMat;
    int    *eigIntMat;
    
    /* Matrix vector multiplication */
    void (*Mvec) (void *, double *, double *);
    
    /* Statistics */
    int nComputed;
    
} asdp_lanczos;

#endif /* def_hdsdp_lanczos_h */

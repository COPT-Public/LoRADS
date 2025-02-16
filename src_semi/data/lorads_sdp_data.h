#ifndef LORADS_SDP_DATA_H
#define LORADS_SDP_DATA_H

#include "def_lorads_sdp_data.h"

extern void sdpDataMatCreate(sdp_coeff **psdpCoeff );
extern void sdpDataMatSetData( sdp_coeff *sdpCoeff, lorads_int nSDPCol, lorads_int dataMatNnz, lorads_int *dataMatIdx, double *dataMatElem );
extern void sdpDataMatDestroy( sdp_coeff **psdpCoeff );
extern void dataMatDenseMultiRkMat(void *A, lorads_sdp_dense *X, double *AX);
extern lorads_int hash_function(lorads_int row, lorads_int col, lorads_int size);
extern void dataMatDenseMV(void *A, double *x, double *y, lorads_int n);
extern void dataMatSparseMultiRkMat(void *A, lorads_sdp_dense *X, double *AX);
extern void dataMatSparseMV(void *A, double *x, double *y, lorads_int n);
extern void dataMatSparseZeros(void *A);
extern void dataMatSparseScale(void *A, double scaleFactor);
extern void dataMatDenseZeros(void *A);
extern void dataMatDenseScale(void *A, double scaleFactor);
#endif
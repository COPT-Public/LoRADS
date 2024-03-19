#ifndef ASDP_CHOL_H
#define ASDP_CHOL_H

#include "cholmod.h"

#ifdef __cplusplus
extern "C" {
#endif


extern void ASDPPosDefDetermineSparse(void *A, double cObjNrm1, int *FLAG);
extern void ASDPPosDefDetermineDense(void *A, double cObjNrm1, int *FLAG);


#ifdef __cplusplus
}
#endif




#endif /* ASDP_CHOL_H */

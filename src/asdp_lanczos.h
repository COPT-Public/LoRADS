#ifndef ASDP_LANCZOS_H
#define ASDP_LANCZOS_H

#ifdef HEADERPATH
#include "linalg/def_asdp_lanczos.h"
#else
#include "def_asdp_lanczos.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern asdp_retcode ALanczosCreate( asdp_lanczos **pALanczos );
extern asdp_retcode ALanczosInit( asdp_lanczos *ALanczos, int nCol, int nSpaceDim );
extern void ALanczosSetData( asdp_lanczos *ALanczos, void *MMat, void (*Mvec) (void *, double *, double *) );
extern asdp_retcode ALanczosSolve( asdp_lanczos *ALanczos, double *LanczosStart, double *dMaxStep );
extern void ALanczosClear( asdp_lanczos *ALanczos );
extern void ALanczosDestroy( asdp_lanczos **pALanczos );

#ifdef __cplusplus
}
#endif

#endif /* asdp_lanczos_h */

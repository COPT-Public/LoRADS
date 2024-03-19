#ifndef ASDP_USER_DATA_H
#define ASDP_USER_DATA_H

#ifdef HEADERPATH
#include "interface/asdp.h"
#include "interface/def_asdp_conic.h"
#else
#include "asdp.h"
#include "def_asdp_conic.h"
#endif

typedef struct asdp_user_data user_data;

#ifdef __cplusplus
extern "C" {
#endif

extern asdp_retcode AUserDataCreate( user_data **pHdata );
extern void AUserDataSetConeData( user_data *Hdata, cone_type cone, int nRow, int nCol,
                                  int *coneMatBeg, int *coneMatIdx, double *coneMatElem);
extern cone_type AUserDataChooseCone( user_data *Hdata );
extern void AUserDataClear( user_data *Hdata );
extern void AUserDataDestroy( user_data **pHdata );

#ifdef __cplusplus
}
#endif

#endif /* ASDP_USER_DATA_H */

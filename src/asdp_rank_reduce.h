
#ifndef ASDP_RANK_REDUCE_H
#define ASDP_RANK_REDUCE_H

#ifdef HEADERPATH
#include "interface/asdp.h"
#include "interface/def_asdp_conic.h"
#else
#include "asdp.h"
#include "def_asdp_conic.h"
#endif

#define ASDP_RANK_REDUCE_TOL  (1e-10)

#ifdef __cplusplus
extern "C" {
#endif

extern asdp_retcode ASDPRankReduce(asdp *ASolver);

#ifdef __cplusplus
}
#endif

#endif /* ASDP_USER_DATA_H */

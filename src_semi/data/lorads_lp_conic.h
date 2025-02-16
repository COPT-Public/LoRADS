#ifndef LORADS_LP_CONIC
#define LORADS_LP_CONIC


#include "def_lorads_lp_conic.h"
#include "lorads.h"

extern void LORADSSetLpCone(lorads_lp_cone *lp_cone, lorads_int nRows,
                            lorads_int nLpCols, lorads_int *lpMatBeg,
                            lorads_int *lpMatIdx, double *LpMatElem);

#endif
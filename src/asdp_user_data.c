/** @file ASDP\_user\_data.c
 *
 */

#ifdef HEADERPATH
#include "interface/def_asdp_user_data.h"
#include "interface/asdp_user_data.h"
#include "interface/asdp_utils.h"
#include "interface/asdp_conic.h"
#include "linalg/sparse_opts.h"
#include "src/asdp_debug.h"
#else
#include "def_asdp_user_data.h"
#include "asdp_user_data.h"
#include "asdp_utils.h"
#include "asdp_conic.h"
#include "sparse_opts.h"
#include "asdp_debug.h"
#endif

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

/** @brief Check if LP data implies bound constraint on y
 *
 */
static int AUserDataICheckLpBound( user_data *Adata ) {
    
    if ( Adata->cone != ASDP_CONETYPE_LP ) {
        return 0;
    }
    
    /* If each column holds at most one variable, then it is bound */
    for ( int i = 0; i < Adata->nConicRow; ++i ) {
        if ( Adata->coneMatBeg[i + 1] - Adata->coneMatBeg[i] >= 2 ) {
            return 0;
        }
    }
    
    return 1;
}

extern asdp_retcode AUserDataCreate( user_data **pAdata ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    if ( !pAdata ) {
        retcode = ASDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    user_data *Adata = NULL;
    ASDP_INIT(Adata, user_data, 1);
    ASDP_MEMCHECK(Adata);
    
    ASDP_ZERO(Adata, user_data, 1);
    *pAdata = Adata;
    
exit_cleanup:
    
    return retcode;
}

extern void AUserDataSetConeData( user_data *Adata, cone_type cone, int nRow, int nCol,
                                  int *coneMatBeg, int *coneMatIdx, double *coneMatElem ) {
    
    Adata->cone = cone;
    Adata->nConicRow = nRow;
    Adata->nConicCol = nCol;
    Adata->coneMatBeg = coneMatBeg;
    Adata->coneMatIdx = coneMatIdx;
    Adata->coneMatElem = coneMatElem;
    
    return;
}

extern cone_type AUserDataChooseCone( user_data *Adata ) {
        
    /* Automatic choice between different cone types*/
    if ( Adata->cone == ASDP_CONETYPE_SOCP || Adata->cone == ASDP_CONETYPE_BOUND ||
         Adata->cone == ASDP_CONETYPE_SPARSE_SDP || Adata->cone == ASDP_CONETYPE_SCALAR_BOUND ) {
        
        return Adata->cone;
        
    } else if ( Adata->cone == ASDP_CONETYPE_DENSE_SDP ) {
        
        int nzSDPCoeffs = csp_nnz_cols(Adata->nConicRow, &Adata->coneMatBeg[1]);
        return ( nzSDPCoeffs > ASDP_SPARSE_CONE_THRESHOLD * Adata->nConicRow ) ? \
                ASDP_CONETYPE_DENSE_SDP : ASDP_CONETYPE_SPARSE_SDP;
        
    } else if ( Adata->cone == ASDP_CONETYPE_LP ) {
        
        if ( AUserDataICheckLpBound(Adata) ) {
            return ASDP_CONETYPE_BOUND;
        } else {
            return ASDP_CONETYPE_LP;
        }
        
    }
    
    return ASDP_CONETYPE_UNKNOWN;
}

extern void AUserDataClear( user_data *Adata ) {
    
    if ( !Adata ) {
        return;
    }
    
    ASDP_ZERO(Adata, user_data, 1);
    
    return;
}

extern void AUserDataDestroy( user_data **pAdata ) {
    
    if ( !pAdata ) {
        return;
    }
    
    AUserDataClear(*pAdata);
    ASDP_FREE(*pAdata);
    
    return;
}

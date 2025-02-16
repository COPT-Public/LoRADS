
#include "def_lorads_user_data.h"
#include "lorads_user_data.h"
#include "lorads_utils.h"
#include "lorads_sparse_opts.h"


#ifdef MEMDEBUG
#include "memwatch.h"
#endif

/** @brief Check if LP data implies bound constraint on y
 *
 */
static lorads_int LUserDataICheckLpBound( user_data *Ldata ) {
    if ( Ldata->cone != LORADS_CONETYPE_LP ) {
        return 0;
    }
    
    /* If each column holds at most one variable, then it is bound */
    for ( lorads_int i = 0; i < Ldata->nConicRow; ++i ) {
        if ( Ldata->coneMatBeg[i + 1] - Ldata->coneMatBeg[i] >= 2 ) {
            return 0;
        }
    }
    
    return 1;
}

extern void LUserDataCreate( user_data **pLdata ) {
    

    if ( !pLdata ) {
        LORADS_ERROR_TRACE;
    }
    
    user_data *Ldata = NULL;
    LORADS_INIT(Ldata, user_data, 1);
    LORADS_MEMCHECK(Ldata);
    
    LORADS_ZERO(Ldata, user_data, 1);
    *pLdata = Ldata;
}

extern void LUserDataSetConeData( user_data *Ldata, cone_type cone, lorads_int nRow, lorads_int nCol,
                                  lorads_int *coneMatBeg, lorads_int *coneMatIdx, double *coneMatElem ) {
    
    Ldata->cone = cone;
    Ldata->nConicRow = nRow;
    Ldata->nConicCol = nCol;
    Ldata->coneMatBeg = coneMatBeg;
    Ldata->coneMatIdx = coneMatIdx;
    Ldata->coneMatElem = coneMatElem;
    
    return;
}

extern cone_type LUserDataChooseCone( user_data *Ldata ) {
        
    /* Automatic choice between different cone types*/
    if ( Ldata->cone == LORADS_CONETYPE_BOUND ||
         Ldata->cone == LORADS_CONETYPE_SPARSE_SDP || Ldata->cone == LORADS_CONETYPE_SCALAR_BOUND ) {
        
        return Ldata->cone;
        
    } else if ( Ldata->cone == LORADS_CONETYPE_DENSE_SDP ) {
        
        lorads_int nzSDPCoeffs = csp_nnz_cols(Ldata->nConicRow, &Ldata->coneMatBeg[1]);
        return ( nzSDPCoeffs > LORADS_SPARSE_CONE_THRESHOLD * Ldata->nConicRow ) ? \
                LORADS_CONETYPE_DENSE_SDP : LORADS_CONETYPE_SPARSE_SDP;
        
    } else if ( Ldata->cone == LORADS_CONETYPE_LP ) {
        
        if ( LUserDataICheckLpBound(Ldata) ) {
            return LORADS_CONETYPE_BOUND;
        } else {
            return LORADS_CONETYPE_LP;
        }
        
    }
    
    return LORADS_CONETYPE_UNKNOWN;
}

extern void LUserDataClear( user_data *Ldata ) {
    
    if ( !Ldata ) {
        return;
    }
    
    LORADS_ZERO(Ldata, user_data, 1);
    
    return;
}

extern void LUserDataDestroy( user_data **pLdata ) {
    
    if ( !pLdata ) {
        return;
    }
    
    LUserDataClear(*pLdata);
    LORADS_FREE(*pLdata);
    
    return;
}


extern void LORADSCreateSDPDatas(user_data ***SDPDatas, lorads_int nCones){
    user_data **SDPDatasTemp;
    LORADS_INIT(SDPDatasTemp, user_data *, nCones);
    *SDPDatas = SDPDatasTemp;
}

extern void LORADSClearUsrData(lorads_int **coneMatBeg, lorads_int **coneMatIdx, double **coneMatElem, lorads_int nBlks,
                               lorads_int *BlkDims, double *rowRHS, lorads_int *LpMatBeg, lorads_int *LpMatIdx, double *LpMatElem, user_data **SDPDatas){
    for (lorads_int iBlk = 0; iBlk < nBlks; ++iBlk)
    {
        LORADS_FREE(coneMatBeg[iBlk]);
        LORADS_FREE(coneMatIdx[iBlk]);
        LORADS_FREE(coneMatElem[iBlk]);
    }
    LORADS_FREE(BlkDims);
    LORADS_FREE(rowRHS);
    LORADS_FREE(LpMatBeg);
    LORADS_FREE(LpMatIdx);
    LORADS_FREE(LpMatElem);
    LORADS_FREE(coneMatBeg);
    LORADS_FREE(coneMatIdx);
    LORADS_FREE(coneMatElem);
    LORADS_FREE(SDPDatas);
}
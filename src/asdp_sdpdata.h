/** @file asdp\_sdpdata.h
 *  @brief Implenent the SDP Coeffcient data operations
 */
#ifndef asdp_sdpdata_h
#define asdp_sdpdata_h

#ifdef HEADERPATH
#include "interface/asdp.h"
#include "linalg/def_asdp_sdpdata.h"
#else
#include "asdp.h"
#include "def_asdp_sdpdata.h"
#include "asdp_utils.h"
#endif

asdp_retcode sdpDataMatCreate( sdp_coeff **psdpCoeff );
asdp_retcode sdpDataMatSetData( sdp_coeff *sdpCoeff, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem );

inline int sdpDataMatGetRank( sdp_coeff *sdpCoeff ) {

    if ( sdpCoeff->dataType == SDP_COEFF_ZERO ) {
        return 0;
    } else if ( sdpCoeff->dataType == SDP_COEFF_DSR1 || sdpCoeff->dataType == SDP_COEFF_SPR1 ) {
        return 1;
    } else if ( sdpCoeff->eigRank != -1 ) {
        return sdpCoeff->eigRank;
    }

    return sdpCoeff->nSDPCol;
}

inline void sdpDataMatScal( sdp_coeff *sdpCoeff, double scal ) {

    sdpCoeff->scal(sdpCoeff->dataMat, scal);

    return;
}

inline double sdpDataMatNorm( sdp_coeff *sdpCoeff, int type ) {

    return sdpCoeff->norm(sdpCoeff->dataMat, type);
}


inline sdp_coeff_type sdpDataMatGetType( sdp_coeff *sdpCoeff ) {

    return sdpCoeff->dataType;
}

inline void sdpDataMatClear( sdp_coeff *sdpCoeff ) {

    if ( !sdpCoeff ) {
        return;
    }

    sdpCoeff->destroy(&sdpCoeff->dataMat);

    if ( sdpCoeff->eigRank != -1 ) {
        ASDP_FREE(sdpCoeff->eigVals);
        ASDP_FREE(sdpCoeff->eigVecs);
    }

    ASDP_ZERO(sdpCoeff, sdp_coeff, 1);

    return;
}

inline void sdpDataMatDestroy( sdp_coeff **psdpCoeff ) {

    if ( !psdpCoeff ) {
        return;
    }

    sdpDataMatClear(*psdpCoeff);
    ASDP_FREE(*psdpCoeff);

    return;
}

inline void sdpDataMatView( sdp_coeff *sdpCoeff ) {

    sdpCoeff->view(sdpCoeff->dataMat);

    return;
}

inline int sdpDataMatIsEye( sdp_coeff *sdpCoeff, double *dEyeMultiple ) {

    return sdpCoeff->iseye(sdpCoeff->dataMat, dEyeMultiple);
}

inline int sdpDataMatIsUnitCol( sdp_coeff *sdpCoeff, int *iUnitCol ) {

    return sdpCoeff->isunitcol(sdpCoeff->dataMat, iUnitCol) ;
}

#ifdef UNDER_BLAS
extern void dsymm_( const char *side, const char *uplo, const int *m,
                   const int *n, const double *alpha, const double *a,
                   const int *lda, const double *b, const int *ldb,
                   const double *beta, double *c, const int *ldc );
#else
extern void dsymm( const char *side, const char *uplo, const int *m,
                   const int *n, const double *alpha, const double *a,
                   const int *lda, const double *b, const int *ldb,
                   const double *beta, double *c, const int *ldc );
#endif

#ifdef __cplusplus
}
#endif

#endif /* asdp_sdpdata_h */

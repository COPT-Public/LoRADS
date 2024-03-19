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
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern asdp_retcode sdpDataMatCreate( sdp_coeff **psdpCoeff );
extern asdp_retcode sdpDataMatSetData( sdp_coeff *sdpCoeff, int nSDPCol, int dataMatNnz, int *dataMatIdx, double *dataMatElem );
extern int sdpDataMatGetRank( sdp_coeff *sdpCoeff );
extern void sdpDataMatScal( sdp_coeff *sdpCoeff, double scal );

extern double sdpDataMatNorm( sdp_coeff *sdpCoeff, int type );
extern asdp_retcode sdpDataMatBuildUpEigs( sdp_coeff *sdpCoeff, double *dAuxFullMatrix );
extern int sdpDataMatGetNnz( sdp_coeff *sdpCoeff );
extern void sdpDataMatDump( sdp_coeff *sdpCoeff, double *dFullMatrix );
extern void sdpDataMatGetMatNz( sdp_coeff *sdpCoeff, int *iMatSpsPattern );
extern void sdpDataMatAddToBuffer( sdp_coeff *sdpCoeff, double dElem, int *iMatSpsPattern, double *dBuffer );
extern sdp_coeff_type sdpDataMatGetType( sdp_coeff *sdpCoeff );
extern void sdpDataMatClear( sdp_coeff *sdpCoeff );
extern void sdpDataMatDestroy( sdp_coeff **psdpCoeff );
extern void sdpDataMatView( sdp_coeff *sdpCoeff );
extern int sdpDataMatIsEye( sdp_coeff *sdpCoeff, double *dEyeMultiple );
extern int sdpDataMatIsUnitCol( sdp_coeff *sdpCoeff, int *iUnitCol );
extern void sdpDataCOPY(sdp_coeff *dst, sdp_coeff *src);
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

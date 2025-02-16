#ifndef DEF_LORADS_USER_DATA_H
#define DEF_LORADS_USER_DATA_H
#include "lorads_utils.h"
#include "lorads.h"
#include "def_lorads_sdp_conic.h"

/* Interface of user data */

/** @struct lorads\_user\_data
 *  @brief HDSDP user conic data for SDP and LP
 *
 * Conic data is defined differently for different cones
 *
 * SDP cone: CSC representation of an [(n + 1) \* n / 2] by [m + 1] matrix, containing SDP matrix
 * coefficients and objective coefficients in each column, only lower triangular is stored
 *
 * LP cone and bound cone: CSC representation of an [n] by [m] matrix, containg LP data
 * 
 */
struct lorads_user_data {
    
    cone_type cone;
    
    lorads_int     nConicRow; /* Number of constraints */
    lorads_int     nConicCol; /* For SDP cone: conic dimension, for LP cone: number of LP columns */
    lorads_int    *coneMatBeg;
    lorads_int    *coneMatIdx;
    double *coneMatElem;
    lorads_int nnz;
};

#endif /* def_lorads_user_data_h */

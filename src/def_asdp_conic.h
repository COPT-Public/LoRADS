#ifndef DEF_ASDP_CONIC_H
#define DEF_ASDP_CONIC_H

#ifdef HEADERPATH
#include "interface/asdp.h"
#include "linalg/def_asdp_sdpdata.h"
#else
#include "asdp.h"
#include "def_asdp_sdpdata.h"
#include "def_asdp_lpdata.h"
#endif
#include <stdint.h>

/* Define conic type */
typedef enum {
    
    ASDP_CONETYPE_UNKNOWN,    // not implement
    ASDP_CONETYPE_LP,         /* A' * y <= c */ // isolated into asdp_lp_cone
    ASDP_CONETYPE_BOUND,      /*      y <= u */
    ASDP_CONETYPE_SCALAR_BOUND,
    ASDP_CONETYPE_DENSE_SDP,
    ASDP_CONETYPE_SPARSE_SDP,
    ASDP_CONETYPE_SOCP
    
} cone_type;

/** @struct asdp\_cone
 *  @brief Define the ASDP general conic interface
 *
 * In ASDP, the general conic interface supports the following functionalities
 *
 *  Connection with interface:
 *  1. Set data
 *  2. Process data
 *  3. Destroy data
 *
 *  Algorithm:
 *  1. Conic initialize
 *  2. Conic iterate maintenance
 *  3. Conic numeric assembly of important quantities
 *  4. Conic symbolic assembly of important quantities
 *  9. Conic scal
 *
 */
typedef struct {
    
    cone_type cone;
    
    sdp_coeff *sdp_coeff_w_sum; // only two type: sparse and dense
    sdp_coeff *sdp_obj_sum; // sdp_coeff data + obj only two type: sparse and dense
    sdp_coeff *sdp_slack_var;
    sdp_coeff *UVt_w_sum;
    sdp_coeff *UVt_obj_sum;
    
    void  *usrData;
    void  *coneData;
    
    // Auxiliary variable
    double sdp_coeff_w_sum_sp_ratio;
    

    int nConstr; // when sparse is not general constraint number

    /* Conic data interface */
    asdp_retcode (*coneCreate)          ( void ** );
    asdp_retcode (*coneProcData)        ( void *, int, int, int *, int *, double * );
    asdp_retcode (*conePresolveData)    ( void * );
    void         (*coneDestroyData)     ( void ** );
    
    void         (*coneScal)            ( void *, double );
    
    // coneAUV(asdp_cone_sdp_dense (coneType) *cone, asdp_rk_mat_dense *U, asdp_rk_mat_dense *V, double *constrVal);
    void         (*coneAUV)             ( void *, asdp_rk_mat_dense *, asdp_rk_mat_dense *, double *, sdp_coeff *, int);
    
    void         (*objAUV)              ( void *, asdp_rk_mat_dense *, asdp_rk_mat_dense *, double *,  sdp_coeff *, int);
    
    void         (*coneObjNrm1)         (void *, double *);
    
    void         (*coneObjNrm2Square)   (void *, double *);
    
    void         (*coneObjNrmInf)       (void *, double *);
    
    void         (*sdpDataWSum)         (void *, double *, sdp_coeff *);
    
    void         (*addObjCoeff)         (void *, sdp_coeff *);
    
    void         (*addObjCoeffRand)     (void *, sdp_coeff *);
    
    /* Debugging */
    void         (*coneView)            ( void * ); // conView(coneData)
    
    /* Feature detection */
    void         (*getstat)             ( void *, double *, int [20], double [20] );
    void         (*nnzStat)             ( void *, int *);
    void         (*nnzStatCoeff)        ( void *, double *, int *, int *);
    void         (*dataScale)           ( void *, double, double);
} asdp_cone;


/* A dense SDP block */
typedef struct {
    
    int   nRow;
    int   nCol;
    
    sdp_coeff **sdpRow;
    sdp_coeff  *sdpObj;
    
    /* SDP block statistics */
    int sdpConeStats[5]; ///< Number of coefficients of each type
    
} asdp_cone_sdp_dense;

/* A sparse SDP block */
typedef struct {
    
    int   nRow; // all row number
    int   nCol;
    
    int nRowElem; // nnz row number
    int *rowIdx;
    sdp_coeff **sdpRow;
    sdp_coeff  *sdpObj;
    
    int sdpConeStats[5];
    
} asdp_cone_sdp_sparse;


/*
 Split an LP cone as
 ---------------------------------
          objMatElem
 ---------------------------------
 |          |          |          |
 | lp_coeff | lp_coeff | lp_coeff |
 |          |          |          |
 */
/* An LP cone */
typedef struct {
    
    int     nRow; // constraint number
    int     nCol; // dim of LP
    
    // obj coeff, full
    double  *objMatElem;
    
    // constraint coeff
    int *rowMatBeg;
    int *rowMatIdx;
    double *rowMatElem;
    
    lp_coeff **lpCol;
    double *nrm2Square;
    
} asdp_cone_lp;


typedef struct {
    void   *coneData; // asdp_cone_lp_data
    int    nCol;
    /* Conic data interface */
    asdp_retcode (*coneCreate)          ( void ** );
    asdp_retcode (*coneProcData)        ( void *, int, int, int *, int *, double * );
    asdp_retcode (*conePresolveData)    ( void * );
    void         (*coneDestroyData)     ( void ** );
    void         (*coneView)            ( void *  );
    void         (*coneAUV)             ( void *, asdp_rk_mat_lp *, asdp_rk_mat_lp *, double *, int);
    void         (*coneAUV2)            ( void *, double *, double *);
    void         (*objAUV)              ( void *, asdp_rk_mat_lp *, asdp_rk_mat_lp *, double *);
    void         (*coneObjNrm1)         ( void *, double *, int);
    void         (*coneObjNrm2Square)   ( void *, double *, int);
    void         (*coneObjNrmInf)       ( void *, double *, int);
    void         (*lpDataWSum)          ( void *, double *, double *, int);
    void         (*objCoeffSum)         ( void *, double *, int);
    void         (*scaleData)           ( void *, double  , double );
} asdpLPCone;

#endif /* def_asdp_conic_h */

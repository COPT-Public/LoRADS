#ifndef DEF_LORADS_SDP_CONIC
#define DEF_LORADS_SDP_CONIC

#include "def_lorads_elements.h"
#include "def_lorads_sdp_data.h"
#include "lorads.h"

#define INT_FEATURE_I_NULLOBJ 0
#define INT_FEATURE_I_MANYCONES 1
#define INT_FEATURE_I_NOPINTERIOR 2
#define INT_FEATURE_I_NODINTERIOR 3
#define INT_FEATURE_I_VERYDENSE 4
#define INT_FEATURE_I_IMPTRACE 5
#define INT_FEATURE_I_IMPYBOUND 6
#define INT_FEATURE_N_SUMCONEDIMS 7
#define INT_FEATURE_N_MAXCONEDIM 8
#define INT_FEATURE_N_CONES 9
#define INT_FEATURE_N_ROWS 10
#define INT_FEATURE_N_SPSDPCONES 11
#define INT_FEATURE_N_DSSDPCONES 12
#define INT_FEATURE_N_LPCONES 13
#define INT_FEATURE_N_BNDCONES 14
#define INT_FEATURE_N_ZEORMATS 15
#define INT_FEATURE_N_SPMATS 16
#define INT_FEATURE_N_DSMATS 17
#define INT_FEATURE_N_SPR1MATS 18
#define INT_FEATURE_N_DSR1MATS 19

#define DBL_FEATURE_OBJFRONORM 0
#define DBL_FEATURE_OBJONENORM 1
#define DBL_FEATURE_RHSFRONORM 2
#define DBL_FEATURE_RHSONENORM 3
#define DBL_FEATURE_RHSINFNORM 4
#define DBL_FEATURE_OBJSCALING 5
#define DBL_FEATURE_RHSSCALING 6
#define DBL_FEATURE_DATAFRONORM 7
#define DBL_FEATURE_DATAONENORM 8
#define DBL_FEATURE_IMPYBOUNDUP 9
#define DBL_FEATURE_IMPYBOUNDLOW 10
#define DBL_FEATURE_IMPTRACEX 11

/* Define conic type */
typedef enum {

    LORADS_CONETYPE_UNKNOWN,    // not implement
    LORADS_CONETYPE_LP,         /* A' * y <= c */ // isolated lorads_into lorads_lp_cone
    LORADS_CONETYPE_BOUND,      /*      y <= u */
    LORADS_CONETYPE_SCALAR_BOUND,
    LORADS_CONETYPE_DENSE_SDP,
    LORADS_CONETYPE_SPARSE_SDP,
} cone_type;

typedef struct{
    cone_type type;
    sdp_coeff *sdp_coeff_w_sum; // only two type: sparse and dense
    sdp_coeff *sdp_obj_sum; // sdp_coeff data + obj only two type: sparse and dense
    sdp_coeff *sdp_slack_var;
    sdp_coeff *UVt_w_sum;
    sdp_coeff *UVt_obj_sum;

    void  *usrData;
    void  *coneData;

    // Auxiliary variable
    double sdp_coeff_w_sum_sp_ratio;


    lorads_int nConstr; // when sparse is not general constraint number

    /* Conic data lorads_interface */
    void         (*coneCreate)          ( void ** );
    void         (*coneProcData)        ( void *, lorads_int, lorads_int, lorads_int *, lorads_int *, double * );
    void         (*conePresolveData)    ( void * );
    void         (*coneDestroyData)     ( void ** );

    // coneAUV(lorads_cone_sdp_dense (coneType) *cone, lorads_sdp_dense *U, lorads_sdp_dense *V, double *constrVal);
    void         (*coneAUV)             ( void *, lorads_sdp_dense *, lorads_sdp_dense *, double *, sdp_coeff *);
    void         (*objAUV)              ( void *, lorads_sdp_dense *, lorads_sdp_dense *, double *,  sdp_coeff *);
    void         (*coneObjNrm1)         (void *, double *);
    void         (*coneObjNrm2Square)   (void *, double *);
    void         (*coneObjNrmInf)       (void *, double *);
    void         (*sdpDataWSum)         (void *, double *, sdp_coeff *);
    void         (*addObjCoeff)         (void *, sdp_coeff *);
    void         (*addObjCoeffRand)     (void *, sdp_coeff *);
    /* Debugging */
    void         (*coneView)            ( void * ); // conView(coneData)

    /* Feature detection */
    void         (*getstat)             ( void *, double *, lorads_int [20], double [20] );
    void         (*nnzStat)             ( void *, lorads_int *);
    void         (*nnzStatCoeff)        ( void *, double *, lorads_int *, lorads_int *);
    void         (*dataScale)           ( void *, double);
    void         (*objScale)            ( void *, double);
    void         (*nnzStatAndDenseDetect)(void *, lorads_int *, bool *);
    void         (*collectNnzPos)       (void *, lorads_int *, lorads_int *);
    void         (*reConstructIndex)    (void *, Dict *);
}lorads_sdp_cone;


/* A dense SDP block */
typedef struct {

    lorads_int   nRow;
    lorads_int   nCol;

    sdp_coeff **sdpRow;
    sdp_coeff  *sdpObj;

    /* SDP block statistics */
    lorads_int sdpConeStats[5]; ///< Number of coefficients of each type

} lorads_cone_sdp_dense;

/* A sparse SDP block */
typedef struct {

    lorads_int   nRow; // all row number
    lorads_int   nCol;

    lorads_int nRowElem; // nnz row number
    lorads_int *rowIdx;
    sdp_coeff **sdpRow;
    sdp_coeff  *sdpObj;

    lorads_int sdpConeStats[5];

} lorads_cone_sdp_sparse;


typedef struct {
    lorads_int i;
    lorads_int j;
} lorads_tuple;

#endif
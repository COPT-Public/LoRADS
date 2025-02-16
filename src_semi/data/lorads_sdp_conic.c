

#include <math.h>
#include "def_lorads_sdp_conic.h"
#include "lorads_utils.h"
#include "lorads_solver.h"
#include "lorads_user_data.h"
#include "lorads_sparse_opts.h"
#include "lorads_sdp_conic.h"
#include "lorads_sdp_data.h"
#include "lorads_vec_opts.h"

extern void SDPConeDestroy(lorads_sdp_cone **pACone)
{
    LORADS_FREE(*pACone);
}

extern void sdpSparseConeCreateImpl(void **pConeIn)
{
    lorads_cone_sdp_sparse **pCone = (lorads_cone_sdp_sparse **)pConeIn;
    LORADS_NULLCHECK(pCone);
    lorads_cone_sdp_sparse *cone = NULL;
    LORADS_INIT(cone, lorads_cone_sdp_sparse, 1);
    LORADS_MEMCHECK(cone);
    LORADS_ZERO(cone, lorads_cone_sdp_sparse, 1);
    *pCone = cone;
}

extern void sdpDenseConeCreateImpl(void **pConeIn)
{
    lorads_cone_sdp_dense **pCone = (lorads_cone_sdp_dense **)pConeIn;

    LORADS_NULLCHECK(pCone);
    lorads_cone_sdp_dense *cone = NULL;

    LORADS_INIT(cone, lorads_cone_sdp_dense, 1);
    LORADS_MEMCHECK(cone);
    LORADS_ZERO(cone, lorads_cone_sdp_dense, 1);
    *pCone = cone;
}

extern void sdpDenseConeDestroyImpl(void **pConeIn)
{
    lorads_cone_sdp_dense **pCone = (lorads_cone_sdp_dense **)pConeIn;
    if (!pCone)
    {
        return;
    }
    sdpDenseConeClearImpl(*pCone);
    LORADS_FREE(*pCone);
}

extern void sdpDenseConeFeatureDetectImpl(void *coneIn, double *rowRHS,
                                          lorads_int coneIntFeatures[20], double coneDblFeatures[20])
{
    /* When there is a single SDP cone. This routine detects if the SDP has the following structures

      1. (Almost) no primal lorads_interior point
      2. (Almost) no dual lorads_interior point
      3. Implied trace bound tr(X) == Z
      4. Implied dual upperbound and lower bound l <= y <= u
      5. No objective
      6. Dense cone

     Also the number of different cones will be records as the conic features (including objective)
     */

    /* We detect the case of no primal lorads_interior by seeking constraints that look like
       trace(a * a' * X) = b \approx 0.0
     */

    /* Get statistics */
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    coneIntFeatures[INT_FEATURE_N_ZEORMATS] = cone->sdpConeStats[SDP_COEFF_ZERO];
    coneIntFeatures[INT_FEATURE_N_DSMATS] = cone->sdpConeStats[SDP_COEFF_DENSE];
    coneIntFeatures[INT_FEATURE_N_SPMATS] = cone->sdpConeStats[SDP_COEFF_SPARSE];

    lorads_int isNoPrimalInterior = 0;
    lorads_int isImpliedTraceX = 0;
    double dImpliedTraceX = 0.0;
    lorads_int iUnitCol = 0;
    lorads_int *iUnitColIdx = NULL;

    if (isNoPrimalInterior)
    {
        coneIntFeatures[INT_FEATURE_I_NOPINTERIOR] = 1;
    }

    /* Detect implied trace bound
       We detect two cases of implied bound.

       The first case is the constraint trace(I * X) == a
       The second case is the diagonal constraint diag(X) = d */

    lorads_int nSumUnit = 0;
    for (lorads_int iRow = 0; iRow < cone->nRow; ++iRow)
    {
        nSumUnit += iUnitColIdx[iRow];
    }

    if (nSumUnit == cone->nCol)
    {
        isImpliedTraceX = 1;
    }

    LORADS_FREE(iUnitColIdx);

    if (isImpliedTraceX)
    {
        coneIntFeatures[INT_FEATURE_I_IMPTRACE] = 1;
        coneDblFeatures[DBL_FEATURE_IMPTRACEX] = dImpliedTraceX;
    }

    /* Detect if there is no objective */
    if (cone->sdpObj->dataType == SDP_COEFF_ZERO)
    {
        coneIntFeatures[INT_FEATURE_I_NULLOBJ] = 1;
    }

    /* Detect dense cone */
    if (cone->sdpConeStats[SDP_COEFF_DENSE] >= 0.7 * cone->nRow)
    {
        coneIntFeatures[INT_FEATURE_I_VERYDENSE] = 1;
    }

    return;
}

extern void sdpDenseConePresolveImpl(void *coneIn)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;

    double *preConeAuxi = NULL;
    /* cone->nCol * cone->nCol for a Full matrix */
    LORADS_INIT(preConeAuxi, double, cone->nCol);
    LORADS_MEMCHECK(preConeAuxi);

    LORADS_FREE(preConeAuxi);
}

/** @brief Set data lorads_into a dense cone
 *
 */
extern void sdpDenseConeProcDataImpl(void *coneIn, lorads_int nRow, lorads_int nCol,
                                     lorads_int *coneMatBeg, lorads_int *coneMatIdx, double *coneMatElem)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;

    /* Create the cone */
    cone->nRow = nRow;
    cone->nCol = nCol;

    LORADS_INIT(cone->sdpRow, sdp_coeff *, nRow);
    LORADS_MEMCHECK(cone->sdpRow);

    LORADS_ZERO(cone->sdpConeStats, lorads_int, 5);

    /* Primal objective */
    sdpDataMatCreate(&cone->sdpObj);
    sdpDataMatSetData(cone->sdpObj, nCol, coneMatBeg[1], coneMatIdx, coneMatElem);
    cone->sdpConeStats[cone->sdpObj->dataType] += 1;

    /* Constraint data */
    for (lorads_int iRow = 0; iRow < nRow; ++iRow)
    {
        sdpDataMatCreate(&cone->sdpRow[iRow]);
        sdpDataMatSetData(cone->sdpRow[iRow], nCol, coneMatBeg[iRow + 2] - coneMatBeg[iRow + 1],
                          coneMatIdx + coneMatBeg[iRow + 1], coneMatElem + coneMatBeg[iRow + 1]);
        sdp_coeff_type coeffType = cone->sdpRow[iRow]->dataType;
        cone->sdpConeStats[coeffType] += 1;
    }
}

extern void sdpSparseConeProcDataImpl(void *coneIn, lorads_int nRow, lorads_int nCol,
                                      lorads_int *coneMatBeg, lorads_int *coneMatIdx, double *coneMatElem)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)coneIn;

    /* Create the cone */
    cone->nRow = nRow;
    cone->nCol = nCol;

    LORADS_ZERO(cone->sdpConeStats, lorads_int, 5);

    /* Primal objective */
    sdpDataMatCreate(&cone->sdpObj);
    sdpDataMatSetData(cone->sdpObj, nCol, coneMatBeg[1], coneMatIdx, coneMatElem);
    cone->sdpConeStats[cone->sdpObj->dataType] += 1;

    /* Sparse constraint data */
    lorads_int nRowElem = 0;

    /* Count #Nonzeros */
    nRowElem = csp_nnz_cols(nRow, &coneMatBeg[1]);

    LORADS_INIT(cone->sdpRow, sdp_coeff *, nRowElem);
    LORADS_INIT(cone->rowIdx, lorads_int, nRowElem);
    LORADS_MEMCHECK(cone->sdpRow);
    LORADS_MEMCHECK(cone->rowIdx);

    lorads_int nRowElemTmp = 0;
    for (lorads_int iRow = 0; iRow < nRow; ++iRow)
    {
        if (coneMatBeg[iRow + 2] - coneMatBeg[iRow + 1] > 0)
        {
            sdpDataMatCreate(&cone->sdpRow[nRowElemTmp]);
            sdpDataMatSetData(cone->sdpRow[nRowElemTmp], nCol, coneMatBeg[iRow + 2] - coneMatBeg[iRow + 1],
                              coneMatIdx + coneMatBeg[iRow + 1], coneMatElem + coneMatBeg[iRow + 1]);
            sdp_coeff_type coeffType = cone->sdpRow[nRowElemTmp]->dataType;
#ifdef LORADS_CONIC_DEBUG
            assert(coeffType != SDP_COEFF_ZERO);
#endif
            cone->sdpConeStats[coeffType] += 1;
            cone->rowIdx[nRowElemTmp] = iRow;
            nRowElemTmp += 1;
        }
        else
        {
            cone->sdpConeStats[SDP_COEFF_ZERO] += 1;
        }
    }

    cone->nRowElem = nRowElem;
#ifdef LORADS_CONIC_DEBUG
    assert(nRowElem == nRowElemTmp);
#endif
}

extern void sdpSparseConePresolveImpl(void *coneIn)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)coneIn;

    double *preConeAuxi = NULL;
    /* cone->nCol * cone->nCol for a Full matrix */
    LORADS_INIT(preConeAuxi, double, cone->nCol);
    LORADS_MEMCHECK(preConeAuxi);
    LORADS_FREE(preConeAuxi);
}
extern void sdpSparseConeFeatureDetectImpl(void *coneIn, double *rowRHS,
                                           lorads_int coneIntFeatures[20], double coneDblFeatures[20])
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    coneIntFeatures[INT_FEATURE_N_ZEORMATS] = cone->sdpConeStats[SDP_COEFF_ZERO];
    coneIntFeatures[INT_FEATURE_N_DSMATS] = cone->sdpConeStats[SDP_COEFF_DENSE];
    coneIntFeatures[INT_FEATURE_N_SPMATS] = cone->sdpConeStats[SDP_COEFF_SPARSE];
    return;
}

extern void sdpDenseConeViewImpl(void *coneIn)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
#ifdef INT32
    printf("- Dense SDP cone of %d rows. \n", cone->nRow);
#endif
#ifdef UNIX_INT64
    printf("- Dense SDP cone of %ld rows. \n", cone->nRow);
#endif
#ifdef MAC_INT64
    printf("- Dense SDP cone of %lld rows. \n", cone->nRow);
#endif
    printf("- Objective: \n");
    cone->sdpObj->view(cone->sdpObj->dataMat);

//    printf("- Constraint: \n");
//    for ( lorads_int iRow = 0; iRow < cone->nRow; ++iRow ) {
//        cone->sdpRow[iRow]->view(cone->sdpRow[iRow]->dataMat);
//    }
#ifdef UNIX_INT64
    printf("- Conic statistics: Zero %d Sp %d Ds %d SpR1 %d DsR1 %d \n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE]);
#endif

#ifdef MAC_INT64
    printf("- Conic statistics: Zero %lld Sp %lld Ds %lld \n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE]);
#endif

#ifdef LORADS_INT32
    printf("- Conic statistics: Zero %ld Sp %ld Ds %ld \n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE]);
#endif
    return;
}

extern void sdpDenseConeAUVImpl(void *coneIn, lorads_sdp_dense *U, lorads_sdp_dense *V, double *constrVal, sdp_coeff *UVt)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    for (lorads_int iRow = 0; iRow < cone->nRow; ++iRow)
    {
        cone->sdpRow[iRow]->mul_inner_rk_double(cone->sdpRow[iRow]->dataMat, U, V, &constrVal[iRow], UVt->dataMat, UVt->dataType);
    }
}

extern void sdpDenseObjAUVImpl(void *coneIn, lorads_sdp_dense *U, lorads_sdp_dense *V, double *pObj, sdp_coeff *UVt)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    double temp = 0.0;
    cone->sdpObj->mul_inner_rk_double(cone->sdpObj->dataMat, U, V, &temp, UVt->dataMat, UVt->dataType);

    pObj[0] += temp;
}

extern void sdpDenseConeObjNrm1(void *coneIn, double *objConeNrm1)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrm1(cone->sdpObj->dataMat, &temp);
    objConeNrm1[0] += temp;
}

extern void sdpDenseConeObjNrm2Square(void *coneIn, double *objConeNrm2Square)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrm2Square(cone->sdpObj->dataMat, &temp);
    objConeNrm2Square[0] += temp;
}

extern void sdpDenseConeObjNrmInf(void *coneIn, double *objConeNrmInf)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrmInf(cone->sdpObj->dataMat, &temp);
    objConeNrmInf[0] = LORADS_MAX(temp, objConeNrmInf[0]);
}

extern void sdpDenseConeAddObjCoeff(void *coneIn, sdp_coeff *w_sum)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    if (w_sum->dataType == SDP_COEFF_DENSE || w_sum->dataType == SDP_COEFF_SPARSE)
    {
        cone->sdpObj->add_sdp_coeff(cone->sdpObj->dataMat, w_sum->dataMat, 1.0, w_sum->dataType);
    }
    else
    {
        LORADS_ERROR_TRACE;
        printf("copy obj sdp coeff to other but `dense, sparse` is developing!\n");
    }
}

extern void sdpDenseConeAddObjCoeffRand(void *coneIn, sdp_coeff *w_sum)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    if (w_sum->dataType == SDP_COEFF_DENSE || w_sum->dataType == SDP_COEFF_SPARSE)
    {
        cone->sdpObj->add_sdp_coeff(cone->sdpObj->dataMat, w_sum->dataMat, rand(), w_sum->dataType);
    }
    else
    {
        LORADS_ERROR_TRACE;
        printf("copy obj sdp coeff to other but `dense, sparse` is developing!\n");
    }
}

extern void sdpDenseConeDataScale(void *coneIn, double scaleFactorSDPData)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    for (lorads_int iRow = 0; iRow < cone->nRow; ++iRow)
    {
        cone->sdpRow[iRow]->scaleData(cone->sdpRow[iRow]->dataMat, scaleFactorSDPData);
    }
}

extern void sdpDenseConeObjScale(void *coneIn, double scaleFactorObj)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    cone->sdpObj->scaleData(cone->sdpObj->dataMat, scaleFactorObj);
}

extern void sdpDenseConeNnzStat(void *coneIn, lorads_int *stat)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    for (lorads_int iRow = 0; iRow < cone->nRow; ++iRow)
    {
        cone->sdpRow[iRow]->statNnz(stat);
    }
}

extern void sdpDenseConeNnzStatCoeff(void *coneIn, double *stat, lorads_int *nnzStat, lorads_int *eleStat)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    lorads_int coeffNum = cone->sdpRow[0]->nSDPCol * (cone->sdpRow[0]->nSDPCol + 1) / 2;
    eleStat[0] += coeffNum;
    for (lorads_int iRow = 0; iRow < cone->nRow; ++iRow)
    {
        lorads_int nnzCoeff = cone->sdpRow[iRow]->getnnz(cone->sdpRow[iRow]->dataMat);
        stat[iRow] = (double)nnzCoeff / (double)coeffNum;
        nnzStat[0] += nnzCoeff;
    }
}

extern void sdpDenseNnzStatAndDenseDetect(void *coneIn, lorads_int *nnz, bool *isDense)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    lorads_int nnzCoeff = cone->sdpObj->getnnz(cone->sdpObj->dataMat);
    nnz[0] += nnzCoeff;
    if (cone->sdpObj->dataType == SDP_COEFF_DENSE)
    {
        isDense[0] = true;
    }
    for (lorads_int iRow = 0; iRow < cone->nRow; ++iRow)
    {
        nnzCoeff = cone->sdpRow[iRow]->getnnz(cone->sdpRow[iRow]->dataMat);
        nnz[0] += nnzCoeff;
        if (cone->sdpRow[iRow]->dataType == SDP_COEFF_DENSE)
        {
            isDense[0] = true;
        }
    }
}

extern void sdpDenseCollectNnzPos(void *coneIn, lorads_int *nnzPosRow, lorads_int *nnzPosCol)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    lorads_int idxStart = 0;
    cone->sdpObj->collectNnzPos(cone->sdpObj->dataMat, nnzPosRow, nnzPosCol, &idxStart);
    for (lorads_int iRow = 0; iRow < cone->nRow; ++iRow)
    {
        cone->sdpRow[iRow]->collectNnzPos(cone->sdpRow[iRow]->dataMat, nnzPosRow, nnzPosCol, &idxStart);
    }
}

extern void sdpDenseReConstructIndex(void *coneIn, Dict *dict)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    if (cone->sdpObj->dataType == SDP_COEFF_SPARSE)
    {
        cone->sdpObj->reConstructIndex(cone->sdpObj->dataMat, dict);
    }

    for (lorads_int iRow = 0; iRow < cone->nRow; ++iRow)
    {
        cone->sdpRow[iRow]->reConstructIndex(cone->sdpRow[iRow]->dataMat, dict);
    }
}

extern void sdpDenseDataWeightSumImpl(void *coneIn, double *weight, sdp_coeff *sdp_data)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    for (lorads_int iRow = 0; iRow < cone->nRow; ++iRow)
    {
        // add data to `sdp_data`;
        cone->sdpRow[iRow]->add_sdp_coeff(cone->sdpRow[iRow]->dataMat, sdp_data->dataMat, weight[iRow], sdp_data->dataType);
    }
}

extern void sdpSparseConeViewImpl(void *coneIn)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)coneIn;
#ifdef LORADS_INT32
    printf("- Sparse SDP cone of %d rows and %d nonzeros. \n", cone->nRow, cone->nRowElem);
#endif
#ifdef MAC_INT64
    printf("- Sparse SDP cone of %lld rows and %lld nonzeros. \n", cone->nRow, cone->nRowElem);
#endif

#ifdef UNIX_INT64
    printf("- Sparse SDP cone of %ld rows and %ld nonzeros. \n", cone->nRow, cone->nRowElem);
#endif
    printf("sparse cone view\n");
    printf("- Objective: \n");
    cone->sdpObj->view(cone->sdpObj->dataMat);
#ifdef LORADS_INT32
    for (lorads_int iRow = 0; iRow < cone->nRowElem; ++iRow)
    {
        printf("%d: ", cone->rowIdx[iRow]);
        sdpDataMatView(cone->sdpRow[iRow]);
    }

    printf("- Conic statistics: Zero %d Sp %d Ds %d SpR1 %d DsR1 %d \n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE],
           cone->sdpConeStats[SDP_COEFF_SPR1], cone->sdpConeStats[SDP_COEFF_DSR1]);
#endif
#ifdef MAC_INT64
    for (lorads_int iRow = 0; iRow < cone->nRowElem; ++iRow)
    {
        printf("%lld: ", cone->rowIdx[iRow]);
        cone->sdpRow[iRow]->view(cone->sdpRow[iRow]->dataMat);
    }

    printf("- Conic statistics: Zero %lld Sp %lld Ds %lld \n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE]);
#endif
#ifdef UNIX_INT64
    for (lorads_int iRow = 0; iRow < cone->nRowElem; ++iRow)
    {
        printf("%ld: ", cone->rowIdx[iRow]);
        sdpDataMatView(cone->sdpRow[iRow]);
    }

    printf("- Conic statistics: Zero %ld Sp %ld Ds %ld SpR1 %ld DsR1 %ld \n", cone->sdpConeStats[SDP_COEFF_ZERO],
           cone->sdpConeStats[SDP_COEFF_SPARSE], cone->sdpConeStats[SDP_COEFF_DENSE],
           cone->sdpConeStats[SDP_COEFF_SPR1], cone->sdpConeStats[SDP_COEFF_DSR1]);
#endif
    return;
}

extern void sdpSparseConeAUVImpl(void *coneIn, lorads_sdp_dense *U, lorads_sdp_dense *V, double *constrVal, sdp_coeff *UVt)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)coneIn;
    for (lorads_int iRow = 0; iRow < cone->nRowElem; ++iRow)
    {
        cone->sdpRow[iRow]->mul_inner_rk_double(cone->sdpRow[iRow]->dataMat, U, V, &constrVal[iRow], UVt->dataMat, UVt->dataType);
    }
}

extern void sdpSparseObjAUVImpl(void *coneIn, lorads_sdp_dense *U, lorads_sdp_dense *V, double *pObj, sdp_coeff *UVt)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)coneIn;
    double temp = 0.0;
    cone->sdpObj->mul_inner_rk_double(cone->sdpObj->dataMat, U, V, &temp, UVt->dataMat, UVt->dataType);
    pObj[0] += temp;
}

extern void sdpSparseConeObjNrm1(void *coneIn, double *objConeNrm1)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrm1(cone->sdpObj->dataMat, &temp);
    objConeNrm1[0] += temp;
}

extern void sdpSparseConeObjNrm2Square(void *coneIn, double *objConeNrm2Square)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrm2Square(cone->sdpObj->dataMat, &temp);
    objConeNrm2Square[0] += temp;
}

extern void sdpSparseConeObjNrmInf(void *coneIn, double *objConeNrmInf)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)coneIn;
    double temp = 0.0;
    cone->sdpObj->nrmInf(cone->sdpObj->dataMat, &temp);
    objConeNrmInf[0] = LORADS_MAX(temp, objConeNrmInf[0]);
}

extern void sdpSparseConeAddObjCoeff(void *pConeIn, sdp_coeff *w_sum)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)pConeIn;
    cone->sdpObj->add_sdp_coeff(cone->sdpObj->dataMat, w_sum->dataMat, 1.0, w_sum->dataType);
}

extern void sdpSparseConeAddObjCoeffRand(void *pConeIn, sdp_coeff *w_sum)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)pConeIn;
    cone->sdpObj->add_sdp_coeff(cone->sdpObj->dataMat, w_sum->dataMat, rand(), w_sum->dataType);
}

extern void sdpSparseConeNnzStat(void *pConeIn, lorads_int *stat)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)pConeIn;
    for (lorads_int iRow = 0; iRow < cone->nRowElem; ++iRow)
    {
        cone->sdpRow[iRow]->statNnz(stat);
    }
}

extern void sdpSparseConeNnzStatCoeff(void *pConeIn, double *stat, lorads_int *nnzStat, lorads_int *eleStat)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)pConeIn;
    lorads_int statSum = cone->sdpRow[0]->nSDPCol * (cone->sdpRow[0]->nSDPCol + 1) / 2;
    eleStat[0] += statSum;
    for (lorads_int iRow = 0; iRow < cone->nRowElem; ++iRow)
    {
        lorads_int nnzCoeff = cone->sdpRow[iRow]->getnnz(cone->sdpRow[iRow]->dataMat);
        stat[cone->rowIdx[iRow]] = (double)nnzCoeff / (double)statSum;
        nnzStat[0] += nnzCoeff;
    }
}

extern void sdpSparseNnzStatAndDenseDetect(void *pConeIn, lorads_int *nnz, bool *isDense)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)pConeIn;
    lorads_int nnzCoeff = cone->sdpObj->getnnz(cone->sdpObj->dataMat);
    nnz[0] += nnzCoeff;
    if (cone->sdpObj->dataType == SDP_COEFF_DENSE)
    {
        isDense[0] = true;
    }
    lorads_int nnzSum = 0;
    for (lorads_int iRow = 0; iRow < cone->nRowElem; ++iRow)
    {
        nnzSum += cone->sdpRow[iRow]->getnnz(cone->sdpRow[iRow]->dataMat);
        if (cone->sdpRow[iRow]->dataType == SDP_COEFF_DENSE)
        {
            isDense[0] = true;
        }
    }
    nnz[0] += nnzSum;
}

extern void sdpSparseCollectNnzPos(void *pConeIn, lorads_int *nnzPosRow, lorads_int *nnzPosCol)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)pConeIn;
    lorads_int idxStart = 0;
    cone->sdpObj->collectNnzPos(cone->sdpObj->dataMat, nnzPosRow, nnzPosCol, &idxStart);
    for (lorads_int iRow = 0; iRow < cone->nRowElem; ++iRow)
    {
        cone->sdpRow[iRow]->collectNnzPos(cone->sdpRow[iRow]->dataMat, nnzPosRow, nnzPosCol, &idxStart);
    }
}

extern void sdpSparseReConstructIndex(void *pConeIn, Dict *dict)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)pConeIn;
    if (cone->sdpObj->dataType == SDP_COEFF_SPARSE)
    {
        cone->sdpObj->reConstructIndex(cone->sdpObj->dataMat, dict);
    }
    for (lorads_int iRow = 0; iRow < cone->nRowElem; ++iRow)
    {
        cone->sdpRow[iRow]->reConstructIndex(cone->sdpRow[iRow]->dataMat, dict);
    }
}

extern void sdpSparseConeDataScale(void *pConeIn, double scaleFactorSDPData)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)pConeIn;
    for (lorads_int iRow = 0; iRow < cone->nRow; ++iRow)
    {
        cone->sdpRow[iRow]->scaleData(cone->sdpRow[iRow]->dataMat, scaleFactorSDPData);
    }
}

extern void sdpSparseConeObjScale(void *pConeIn, double scaleFactorObj)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)pConeIn;
    cone->sdpObj->scaleData(cone->sdpObj->dataMat, scaleFactorObj);
}

extern void sdpSparseDataWeightSumImpl(void *coneIn, double *weight, sdp_coeff *sdp_data)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)coneIn;
    for (lorads_int iRow = 0; iRow < cone->nRowElem; ++iRow)
    {
        // add data to `sdp_data`;
        cone->sdpRow[iRow]->add_sdp_coeff(cone->sdpRow[iRow]->dataMat, sdp_data->dataMat, weight[cone->rowIdx[iRow]], sdp_data->dataType);
    }
}

extern void sdpSparseConeClearImpl(void *coneIn)
{
    lorads_cone_sdp_sparse *cone = (lorads_cone_sdp_sparse *)coneIn;

    if (!cone)
    {
        return;
    }

    LORADS_FREE(cone->rowIdx);
    sdpDataMatDestroy(&cone->sdpObj);

    for (lorads_int iRow = 0; iRow < cone->nRowElem; ++iRow)
    {
        sdpDataMatDestroy(&cone->sdpRow[iRow]);
    }

    LORADS_FREE(cone->sdpRow);

    LORADS_ZERO(cone, lorads_cone_sdp_sparse, 1);

    return;
}

extern void sdpSparseConeDestroyImpl(void **coneIn)
{
    lorads_cone_sdp_sparse **pCone = (lorads_cone_sdp_sparse **)coneIn;
    if (!pCone)
    {
        return;
    }
    sdpSparseConeClearImpl(*pCone);
    LORADS_FREE(*pCone);
    return;
}

extern void sdpDenseConeClearImpl(void *coneIn)
{
    lorads_cone_sdp_dense *cone = (lorads_cone_sdp_dense *)coneIn;
    if (!cone)
    {
        return;
    }
    sdpDataMatDestroy(&cone->sdpObj);
    for (lorads_int iRow = 0; iRow < cone->nRow; ++iRow)
    {
        sdpDataMatDestroy(&cone->sdpRow[iRow]);
    }
    LORADS_FREE(cone->sdpRow);
    LORADS_ZERO(cone, lorads_cone_sdp_dense, 1);
    return;
}

extern void LConeSetData(lorads_sdp_cone *ACone, user_data *usrData)
{
    ACone->usrData = usrData;
    ACone->type = LUserDataChooseCone(usrData);
    switch (ACone->type)
    {
    case LORADS_CONETYPE_DENSE_SDP:
        ACone->coneCreate = sdpDenseConeCreateImpl;
        ACone->coneProcData = sdpDenseConeProcDataImpl;
        ACone->conePresolveData = sdpDenseConePresolveImpl;
        ACone->coneDestroyData = sdpDenseConeDestroyImpl;
        ACone->coneView = sdpDenseConeViewImpl;
        ACone->coneAUV = sdpDenseConeAUVImpl;
        ACone->objAUV = sdpDenseObjAUVImpl;
        ACone->sdpDataWSum = sdpDenseDataWeightSumImpl;
        ACone->getstat = sdpDenseConeFeatureDetectImpl;
        ACone->coneObjNrm1 = sdpDenseConeObjNrm1;
        ACone->coneObjNrm2Square = sdpDenseConeObjNrm2Square;
        ACone->coneObjNrmInf = sdpDenseConeObjNrmInf;
        ACone->addObjCoeff = sdpDenseConeAddObjCoeff;
        ACone->addObjCoeffRand = sdpDenseConeAddObjCoeffRand;
        ACone->nnzStat = sdpDenseConeNnzStat;           // high level
        ACone->nnzStatCoeff = sdpDenseConeNnzStatCoeff; // more precision
        ACone->dataScale = sdpDenseConeDataScale;
        ACone->objScale = sdpDenseConeObjScale;
        ACone->nnzStatAndDenseDetect = sdpDenseNnzStatAndDenseDetect;
        ACone->collectNnzPos = sdpDenseCollectNnzPos;
        ACone->reConstructIndex = sdpDenseReConstructIndex;
        break;
    case LORADS_CONETYPE_SPARSE_SDP:
        ACone->coneCreate = sdpSparseConeCreateImpl;
        ACone->coneProcData = sdpSparseConeProcDataImpl;
        ACone->conePresolveData = sdpSparseConePresolveImpl;
        ACone->coneDestroyData = sdpSparseConeDestroyImpl;
        ACone->coneView = sdpSparseConeViewImpl;
        ACone->coneAUV = sdpSparseConeAUVImpl;
        ACone->objAUV = sdpSparseObjAUVImpl;
        ACone->sdpDataWSum = sdpSparseDataWeightSumImpl;
        ACone->getstat = sdpSparseConeFeatureDetectImpl;
        ACone->coneObjNrm1 = sdpSparseConeObjNrm1;
        ACone->coneObjNrm2Square = sdpSparseConeObjNrm2Square;
        ACone->coneObjNrmInf = sdpSparseConeObjNrmInf;
        ACone->addObjCoeff = sdpSparseConeAddObjCoeff;
        ACone->addObjCoeffRand = sdpSparseConeAddObjCoeffRand;
        ACone->nnzStat = sdpSparseConeNnzStat;           // high level
        ACone->nnzStatCoeff = sdpSparseConeNnzStatCoeff; // more precision
        ACone->dataScale = sdpSparseConeDataScale;
        ACone->objScale = sdpSparseConeObjScale;
        ACone->nnzStatAndDenseDetect = sdpSparseNnzStatAndDenseDetect;
        ACone->collectNnzPos = sdpSparseCollectNnzPos;
        ACone->reConstructIndex = sdpSparseReConstructIndex;
        break;
    default:
        LORADS_ERROR_TRACE;
    }
}

extern void LORADSSetCone(lorads_solver *ASolver, lorads_int iCone, void *userCone)
{
    LConeSetData(ASolver->SDPCones[iCone], userCone);
}

extern void AConeProcData(lorads_sdp_cone *ACone)
{
    user_data *usrData = (user_data *)ACone->usrData;
    ACone->coneCreate(&ACone->coneData); // create a null coneData
    ACone->coneProcData(ACone->coneData, usrData->nConicRow, usrData->nConicCol,
                        usrData->coneMatBeg, usrData->coneMatIdx, usrData->coneMatElem);
}

extern void destroyForAuxiDense(void **pA)
{
    sdp_coeff_dense *dense = (sdp_coeff_dense *)pA;
    LORADS_FREE(dense->dsMatElem);
    //    for  (lorads_int i = 0; i < dense->nSDPCol; ++i){
    //        LORADS_FREE(dense->rowCol2NnzIdx[i]);
    //    }
    //    LORADS_FREE(dense->rowCol2NnzIdx);
    LORADS_FREE(dense->fullMat); // fullMat is auxiliary for calculation, not exist in data
    LORADS_FREE(pA);
}

// 比较函数，用于qsort
int cmpfunc(const void *a, const void *b)
{
    lorads_tuple *tuple_a = (lorads_tuple *)a;
    lorads_tuple *tuple_b = (lorads_tuple *)b;

    if (tuple_a->j == tuple_b->j)
    {
        return (tuple_a->i > tuple_b->i) - (tuple_a->i < tuple_b->i);
    }
    else
    {
        return (tuple_a->j > tuple_b->j) - (tuple_a->j < tuple_b->j);
    }
}

// 获取唯一元素的函数
void get_unique_tuples(lorads_int *arr1, lorads_int *arr2, lorads_int n, lorads_tuple *result, lorads_int *result_size)
{
    lorads_tuple *tuples = (lorads_tuple *)malloc(n * sizeof(lorads_tuple));

    // 将两个数组打包成元组
    for (int k = 0; k < n; k++)
    {
        tuples[k].i = arr1[k];
        tuples[k].j = arr2[k];
    }

    // 对元组数组进行排序
    qsort(tuples, n, sizeof(lorads_tuple), cmpfunc);

    // 获取唯一的元组
    *result_size = 0;
    result[*result_size] = tuples[0];
    (*result_size)++;

    for (int i = 1; i < n; i++)
    {
        if (tuples[i].i != tuples[i - 1].i || tuples[i].j != tuples[i - 1].j)
        {
            result[*result_size] = tuples[i];
            (*result_size)++;
        }
    }

    free(tuples);
}

// 初始化字典
Dict *create_dict(lorads_int size)
{
    Dict *dict = (Dict *)malloc(sizeof(Dict));
    dict->table = (DictNode **)malloc(sizeof(DictNode *) * size);
    dict->size = size;
    for (lorads_int i = 0; i < size; i++)
    {
        dict->table[i] = NULL;
    }
    return dict;
}

// 插入元素到字典
void insert_dict(Dict *dict, lorads_int row, lorads_int col, lorads_int index)
{
    lorads_int hash_index = hash_function(row, col, dict->size);
    DictNode *new_node = (DictNode *)malloc(sizeof(DictNode));
    new_node->element.row = row;
    new_node->element.col = col;
    new_node->element.index = index;
    new_node->next = dict->table[hash_index];
    dict->table[hash_index] = new_node;
}

// 释放字典的内存
void free_dict(Dict *dict)
{
    for (lorads_int i = 0; i < dict->size; i++)
    {
        DictNode *node = dict->table[i];
        while (node != NULL)
        {
            DictNode *temp = node;
            node = node->next;
            free(temp);
        }
    }
    free(dict->table);
    free(dict);
}

extern void AConePresolveData(lorads_sdp_cone *ACone, lorads_int Dim)
{
    ACone->conePresolveData(ACone->coneData);
    sdp_coeff *w_sum;
    LORADS_INIT(w_sum, sdp_coeff, 1);
    LORADS_MEMCHECK(w_sum);
    w_sum->nSDPCol = Dim;

    sdp_coeff *sdp_obj_sum;
    LORADS_INIT(sdp_obj_sum, sdp_coeff, 1);
    LORADS_MEMCHECK(sdp_obj_sum);
    sdp_obj_sum->nSDPCol = Dim;
    sdp_coeff *slackVarTemp;
    LORADS_INIT(slackVarTemp, sdp_coeff, 1);

    slackVarTemp->nSDPCol = Dim;
    if (Dim < 20)
    {
    // set dense initially
    DENSE_FLAG:
        w_sum->dataType = SDP_COEFF_DENSE;
        sdp_coeff_dense *dataMat;
        LORADS_INIT(dataMat, sdp_coeff_dense, 1);
        dataMat->nSDPCol = w_sum->nSDPCol;
        LORADS_INIT(dataMat->dsMatElem, double, Dim *(Dim + 1) / 2);
        LORADS_MEMCHECK(dataMat->dsMatElem);
        LORADS_ZERO(dataMat->dsMatElem, double, Dim *(Dim + 1) / 2);
        LORADS_INIT(dataMat->fullMat, double, Dim *Dim);
        LORADS_MEMCHECK(dataMat->fullMat);
        w_sum->dataMat = (void *)dataMat;
        w_sum->mul_rk = dataMatDenseMultiRkMat;
        w_sum->mv = dataMatDenseMV;
        w_sum->destroy = destroyForAuxiDense;
        w_sum->zeros = dataMatDenseZeros;
        w_sum->scaleData = dataMatDenseScale;
        ACone->sdp_coeff_w_sum = w_sum;

        sdp_obj_sum->dataType = SDP_COEFF_DENSE;
        sdp_coeff_dense *dataMatObj;
        LORADS_INIT(dataMatObj, sdp_coeff_dense, 1);
        dataMatObj->nSDPCol = sdp_obj_sum->nSDPCol;
        LORADS_INIT(dataMatObj->dsMatElem, double, Dim *(Dim + 1) / 2);
        LORADS_MEMCHECK(dataMatObj->dsMatElem);
        LORADS_ZERO(dataMatObj->dsMatElem, double, Dim *(Dim + 1) / 2);
        LORADS_INIT(dataMatObj->fullMat, double, Dim *Dim);
        LORADS_MEMCHECK(dataMatObj->fullMat);
        sdp_obj_sum->dataMat = (void *)dataMatObj;
        sdp_obj_sum->mul_rk = dataMatDenseMultiRkMat;
        sdp_obj_sum->mv = dataMatDenseMV;
        sdp_obj_sum->destroy = destroyForAuxiDense;
        sdp_obj_sum->zeros = dataMatDenseZeros;
        sdp_obj_sum->scaleData = dataMatDenseScale;
        ACone->sdp_obj_sum = sdp_obj_sum;

        slackVarTemp->dataType = SDP_COEFF_DENSE;
        sdp_coeff_dense *slackVarSdp_coeff;
        LORADS_INIT(slackVarSdp_coeff, sdp_coeff_dense, 1);
        slackVarSdp_coeff->nSDPCol = slackVarTemp->nSDPCol;
        LORADS_INIT(slackVarSdp_coeff->dsMatElem, double, Dim *(Dim + 1) / 2);
        LORADS_MEMCHECK(slackVarSdp_coeff->dsMatElem);
        LORADS_ZERO(slackVarSdp_coeff->dsMatElem, double, Dim *(Dim + 1) / 2);
        LORADS_INIT(slackVarSdp_coeff->fullMat, double, Dim *Dim);
        LORADS_MEMCHECK(slackVarSdp_coeff->fullMat);
        slackVarTemp->dataMat = slackVarSdp_coeff;
        slackVarTemp->mul_rk = dataMatDenseMultiRkMat;
        slackVarTemp->mv = dataMatDenseMV;
        slackVarTemp->destroy = destroyForAuxiDense;
        slackVarTemp->zeros = dataMatDenseZeros;
        slackVarTemp->scaleData = dataMatDenseScale;
        ACone->sdp_slack_var = slackVarTemp;
        lorads_int result_size = Dim * (Dim + 1) / 2;
        lorads_tuple *result;
        LORADS_INIT(result, lorads_tuple, result_size);
        lorads_int rowNum, colNum, count = 0;
        for (lorads_int colNum = 0; colNum < Dim; ++colNum)
        {
            for (lorads_int rowNum = colNum; rowNum < Dim; ++rowNum)
            {
                result[count].i = rowNum;
                result[count].j = colNum;
                count++;
            }
        }
        LORADS_MEMCHECK(result);
        Dict *dict = create_dict(result_size);
        for (lorads_int i = 0; i < result_size; ++i)
        {
            insert_dict(dict, result[i].i, result[i].j, i);
        }

        ACone->reConstructIndex(ACone->coneData, dict);

        free_dict(dict);
        LORADS_FREE(result);
        return;
    }

    // dim is too large, direct judge it sparse or dense
    lorads_int nnz[1];
    nnz[0] = 0;
    bool denseFlag = false;
    ACone->nnzStatAndDenseDetect(ACone->coneData, nnz, &denseFlag);
    if (denseFlag)
    {
        goto DENSE_FLAG;
    }
    lorads_int *nnzRow, *nnzCol;
    LORADS_INIT(nnzRow, lorads_int, nnz[0]);
    LORADS_INIT(nnzCol, lorads_int, nnz[0]);
    ACone->collectNnzPos(ACone->coneData, nnzRow, nnzCol);
    lorads_tuple *result;
    LORADS_INIT(result, lorads_tuple, nnz[0]);
    LORADS_MEMCHECK(result);
    lorads_int result_size = 0;
    get_unique_tuples(nnzRow, nnzCol, nnz[0], result, &result_size);

    //    for (lorads_int idx = 0; idx < result_size; ++idx){
    //        printf("row: %lld, col: %lld\n", result[idx].i, result[idx].j);
    //    }

    double spRatio = (double)result_size / (double)(Dim * (Dim + 1) / 2);
    if (spRatio < 0.1)
    {
        // set sparse initially
        w_sum->dataType = SDP_COEFF_SPARSE;
        sdp_coeff_sparse *dataMat;
        LORADS_INIT(dataMat, sdp_coeff_sparse, 1);
        LORADS_INIT(dataMat->triMatRow, lorads_int, result_size);
        LORADS_INIT(dataMat->triMatCol, lorads_int, result_size);
        LORADS_INIT(dataMat->triMatElem, double, result_size);
        dataMat->nTriMatElem = result_size;
        dataMat->nSDPCol = Dim;
        for (lorads_int i = 0; i < result_size; ++i)
        {
            dataMat->triMatElem[i] = 0.0;
            dataMat->triMatRow[i] = result[i].i;
            dataMat->triMatCol[i] = result[i].j;
        }
        w_sum->dataMat = (void *)dataMat;
        w_sum->mul_rk = dataMatSparseMultiRkMat;
        w_sum->mv = dataMatSparseMV;
        w_sum->destroy = destroyForAuxiSparse;
        w_sum->zeros = dataMatSparseZeros;
        w_sum->scaleData = dataMatSparseScale;
        ACone->sdp_coeff_w_sum = w_sum;

        sdp_obj_sum->dataType = SDP_COEFF_SPARSE;
        sdp_coeff_sparse *dataMatObj;
        LORADS_INIT(dataMatObj, sdp_coeff_sparse, 1);
        LORADS_INIT(dataMatObj->triMatRow, lorads_int, result_size);
        LORADS_INIT(dataMatObj->triMatCol, lorads_int, result_size);
        LORADS_INIT(dataMatObj->triMatElem, double, result_size);
        dataMatObj->nTriMatElem = result_size;
        dataMatObj->nSDPCol = Dim;
        for (lorads_int i = 0; i < result_size; ++i)
        {
            dataMatObj->triMatElem[i] = 0.0;
            dataMatObj->triMatRow[i] = result[i].i;
            dataMatObj->triMatCol[i] = result[i].j;
        }
        sdp_obj_sum->dataMat = (void *)dataMatObj;
        sdp_obj_sum->mul_rk = dataMatSparseMultiRkMat;
        sdp_obj_sum->mv = dataMatSparseMV;
        sdp_obj_sum->destroy = destroyForAuxiSparse;
        sdp_obj_sum->zeros = dataMatSparseZeros;
        sdp_obj_sum->scaleData = dataMatSparseScale;
        ACone->sdp_obj_sum = sdp_obj_sum;

        slackVarTemp->dataType = SDP_COEFF_SPARSE;
        sdp_coeff_sparse *slackVarSdp_coeff;
        LORADS_INIT(slackVarSdp_coeff, sdp_coeff_sparse, 1);
        LORADS_INIT(slackVarSdp_coeff->triMatRow, lorads_int, result_size);
        LORADS_INIT(slackVarSdp_coeff->triMatCol, lorads_int, result_size);

        LORADS_INIT(slackVarSdp_coeff->triMatElem, double, result_size);
        slackVarSdp_coeff->nTriMatElem = result_size;
        slackVarSdp_coeff->nSDPCol = Dim;
        for (lorads_int i = 0; i < result_size; ++i)
        {
            slackVarSdp_coeff->triMatElem[i] = 0.0;
            slackVarSdp_coeff->triMatRow[i] = result[i].i;
            slackVarSdp_coeff->triMatCol[i] = result[i].j;
        }
        slackVarTemp->dataMat = slackVarSdp_coeff;
        slackVarTemp->mul_rk = dataMatSparseMultiRkMat;
        slackVarTemp->mv = dataMatSparseMV;
        slackVarTemp->destroy = destroyForAuxiSparse;
        slackVarTemp->zeros = dataMatSparseZeros;
        slackVarTemp->scaleData = dataMatSparseScale;
        ACone->sdp_slack_var = slackVarTemp;

        Dict *dict = create_dict(result_size);
        for (lorads_int i = 0; i < result_size; ++i)
        {
            insert_dict(dict, result[i].i, result[i].j, i);
        }

        ACone->reConstructIndex(ACone->coneData, dict);

        free_dict(dict);
    }
    LORADS_FREE(result);
    LORADS_FREE(nnzRow);
    LORADS_FREE(nnzCol);
    if (spRatio >= 0.1)
    {
        goto DENSE_FLAG;
    }
}

extern void LORADSSumSDPData(lorads_solver *ASolver)
{
    // sum all sdp_coeff in ACone->sdp_coeff_w_sum
    // for choosing proper calculation in CG
    double *weight;
    LORADS_INIT(weight, double, ASolver->nRows);
    srand(1);
    for (lorads_int i = 0; i < ASolver->nRows; ++i)
    {
        weight[i] = (double)rand() / RAND_MAX;
        weight[i] = weight[i] * 100;
#ifdef LORADS_SDPdata_DEBUG
        weight[i] = 1.0;
#endif
    }
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        ACone->sdpDataWSum(ACone->coneData, weight, ACone->sdp_coeff_w_sum);
        ACone->sdpDataWSum(ACone->coneData, weight, ACone->sdp_obj_sum);
        ACone->addObjCoeffRand(ACone->coneData, ACone->sdp_obj_sum);
    }
    LORADS_FREE(weight);
}

extern void destroyForAuxiSparse(void **pA)
{
    sdp_coeff_sparse *sparse = (sdp_coeff_sparse *)pA;
    LORADS_FREE(sparse->triMatCol);
    LORADS_FREE(sparse->triMatRow);
    LORADS_FREE(sparse->triMatElem);
    //    for (lorads_int row = 0; row < sparse->nSDPCol; ++row){
    //        LORADS_FREE(sparse->rowCol2NnzIdx[row]);
    //    }
    //    LORADS_FREE(sparse->rowCol2NnzIdx);
    LORADS_FREE(pA);
}

extern void AConeDenseDetectSparsity(sdp_coeff **sdp_coeff_w_sum_pointer)
{
    // add modify sdp_coeff_w_sum if spRatio is smaller than 0.1
    sdp_coeff *sdp_coeff_w_sum = *sdp_coeff_w_sum_pointer;
    sdp_coeff_dense *dense = (sdp_coeff_dense *)sdp_coeff_w_sum->dataMat;
    lorads_int n = dense->nSDPCol;
    lorads_int row = 0;
    lorads_int col = 0;
    lorads_int nnz = 0;
    for (lorads_int i = 0; i < (n + 1) * n / 2; ++i)
    {
        if (fabs(dense->dsMatElem[i]) > 0.0)
        {
            nnz += 1;
        }
        row++;
        if (row == n)
        {
            col++;
            row = col;
        }
    }
    double spRatio = (double)nnz / (double)((n + 1) * n / 2);
    row = 0;
    col = 0;
    if (spRatio <= 0.1)
    {
        sdp_coeff_w_sum->dataType = SDP_COEFF_SPARSE;
        sdp_coeff_sparse *sparse;
        LORADS_INIT(sparse, sdp_coeff_sparse, 1);
        LORADS_INIT(sparse->triMatRow, lorads_int, nnz);
        LORADS_INIT(sparse->triMatCol, lorads_int, nnz);
        LORADS_INIT(sparse->triMatElem, double, nnz);
        sparse->nTriMatElem = nnz;
        sparse->nSDPCol = n;
        //        LORADS_INIT(sparse->rowCol2NnzIdx, lorads_int *, n);
        //        for (lorads_int row = 0; row < n; row++)
        //        {
        //            LORADS_INIT(sparse->rowCol2NnzIdx[row], lorads_int, row + 1);
        //            for (lorads_int i = 0; i < row + 1; ++i)
        //            {
        //                sparse->rowCol2NnzIdx[row][i] = -1;
        //            }
        //        }
        lorads_int count = 0;
        for (lorads_int i = 0; i < (n + 1) * n / 2; ++i)
        {
            if (fabs(dense->dsMatElem[i]) > 0.0)
            {
                sparse->triMatElem[count] = dense->dsMatElem[i];
                sparse->triMatRow[count] = row;
                sparse->triMatCol[count] = col;
                //                sparse->rowCol2NnzIdx[row][col] = count;
                count++;
            }
            row++;
            if (row == n)
            {
                col++;
                row = col;
            }
        }
        LORADS_FREE(dense->dsMatElem);
        LORADS_FREE(dense);
        sdp_coeff_w_sum->dataMat = (void *)sparse;
        sdp_coeff_w_sum->mul_rk = dataMatSparseMultiRkMat;
        sdp_coeff_w_sum->mv = dataMatSparseMV;
        sdp_coeff_w_sum->destroy = destroyForAuxiSparse;
        sdp_coeff_w_sum->zeros = dataMatSparseZeros;
    }
    else
    {
        sdp_coeff_w_sum->dataType = SDP_COEFF_DENSE;
        //        LORADS_INIT(dense->rowCol2NnzIdx, lorads_int *, dense->nSDPCol);
        //        for (lorads_int row = 0; row < dense->nSDPCol; ++row)
        //        {
        //            LORADS_INIT(dense->rowCol2NnzIdx[row], lorads_int, row + 1);
        //            for (lorads_int i = 0; i < row + 1; ++i)
        //            {
        //                dense->rowCol2NnzIdx[row][i] = -1;
        //            }
        //        }
        //        lorads_int count = 0;
        //        for (lorads_int col = 0; col < dense->nSDPCol; ++col)
        //        {
        //            for (lorads_int row = col; row < dense->nSDPCol; ++row)
        //            {
        //                dense->rowCol2NnzIdx[row][col] = count;
        //                count++;
        //            }
        //        }
        LORADS_INIT(dense->fullMat, double, dense->nSDPCol * dense->nSDPCol);
        LORADS_MEMCHECK(dense->fullMat);
        sdp_coeff_w_sum->destroy = destroyForAuxiDense;
    }
    *sdp_coeff_w_sum_pointer = sdp_coeff_w_sum;
}

static void modifySlackVar(sdp_coeff *w_sum, sdp_coeff **slackVarPointer)
{
    if (w_sum->dataType == SDP_COEFF_DENSE)
    {
        sdp_coeff *slackVar = *slackVarPointer;
        sdp_coeff_dense *denseSlack = (sdp_coeff_dense *)slackVar->dataMat;
        lorads_int n = w_sum->nSDPCol;
        //        LORADS_INIT(denseSlack->rowCol2NnzIdx, lorads_int *, n);
        sdp_coeff_dense *w_sumData = (sdp_coeff_dense *)w_sum->dataMat;
        //        for (lorads_int i = 0; i < n; ++i)
        //        {
        //            LORADS_INIT(denseSlack->rowCol2NnzIdx[i], lorads_int, i + 1);
        //            LORADS_MEMCPY(denseSlack->rowCol2NnzIdx[i], w_sumData->rowCol2NnzIdx[i], lorads_int, i + 1);
        //        }
        LORADS_INIT(denseSlack->fullMat, double, denseSlack->nSDPCol * denseSlack->nSDPCol);
        slackVar->destroy = destroyForAuxiDense;
        return;
    }
    else if (w_sum->dataType == SDP_COEFF_SPARSE)
    {
        sdp_coeff *slackVar = *slackVarPointer;
        slackVar->dataType = SDP_COEFF_SPARSE;
        sdp_coeff_dense *dense = (sdp_coeff_dense *)slackVar->dataMat;
        sdp_coeff_sparse *sparseSrc = (sdp_coeff_sparse *)w_sum->dataMat;
        sdp_coeff_sparse *sparse;
        lorads_int nnz = sparseSrc->nTriMatElem;
        LORADS_INIT(sparse, sdp_coeff_sparse, 1);
        sparse->nTriMatElem = nnz;
        LORADS_INIT(sparse->triMatRow, lorads_int, nnz);
        LORADS_INIT(sparse->triMatCol, lorads_int, nnz);
        LORADS_INIT(sparse->triMatElem, double, nnz);
        LORADS_MEMCPY(sparse->triMatRow, sparseSrc->triMatRow, lorads_int, nnz);
        LORADS_MEMCPY(sparse->triMatCol, sparseSrc->triMatCol, lorads_int, nnz);
        LORADS_ZERO(sparse->triMatElem, double, nnz);
        sparse->nSDPCol = sparseSrc->nSDPCol;
        //        LORADS_INIT(sparse->rowCol2NnzIdx, lorads_int *, sparse->nSDPCol);
        //        for (lorads_int row = 0; row < sparse->nSDPCol; ++row)
        //        {
        //            LORADS_INIT(sparse->rowCol2NnzIdx[row], lorads_int, row + 1);
        //            LORADS_MEMCPY(sparse->rowCol2NnzIdx[row], sparseSrc->rowCol2NnzIdx[row], lorads_int, row + 1);
        //        }
        slackVar->dataMat = (void *)sparse;
        slackVar->mul_rk = dataMatSparseMultiRkMat;
        slackVar->mv = dataMatSparseMV;
        slackVar->destroy = destroyForAuxiSparse;
        slackVar->zeros = dataMatSparseZeros;
        *slackVarPointer = slackVar;
        LORADS_FREE(dense->dsMatElem);
        LORADS_FREE(dense);
    }
}
extern void LORADSDetectSparsityOfSumSDP(lorads_solver *ASolver)
{
    for (lorads_int iCone = 0; iCone < ASolver->nCones; ++iCone)
    {
        lorads_sdp_cone *ACone = ASolver->SDPCones[iCone];
        AConeDenseDetectSparsity(&ACone->sdp_coeff_w_sum);
        AConeDenseDetectSparsity(&ACone->sdp_obj_sum);
        modifySlackVar(ACone->sdp_obj_sum, &ACone->sdp_slack_var);
    }
}

extern void AConeDestroyPresolveData(lorads_sdp_cone *ACone)
{
    ACone->sdp_coeff_w_sum->destroy(ACone->sdp_coeff_w_sum->dataMat);
    LORADS_FREE(ACone->sdp_coeff_w_sum);
    ACone->sdp_slack_var->destroy(ACone->sdp_slack_var->dataMat);
    LORADS_FREE(ACone->sdp_slack_var);
    ACone->sdp_obj_sum->destroy(ACone->sdp_obj_sum->dataMat);
    LORADS_FREE(ACone->sdp_obj_sum);
}

int dual_infeasible(void (*matvec)(void *M, double *x, double *y, lorads_int n), void *M, double *res, lorads_int n)
{
    int nev = 1;  // Number of eigenvalues to compute
    int ncv = 40; // Number of Arnoldi vectors
    if (ncv > n)
    {
        nev = 1;
        ncv = 2;
    }
    double tol = 1e-2; // Convergence tolerance
    double *resid = (double *)malloc(n * sizeof(double));
    double *v = (double *)malloc(n * ncv * sizeof(double));
    double *workd = (double *)malloc(3 * n * sizeof(double));
    double *workl = (double *)malloc((3 * ncv * ncv + 6 * ncv) * sizeof(double));
    int *iparam = (int *)malloc(11 * sizeof(int));
    int *ipntr = (int *)malloc(14 * sizeof(int));
    int *select = (int *)malloc(ncv * sizeof(int));
    double *d = (double *)malloc(nev * sizeof(double));
    int lworkl = 3 * ncv * ncv + 5 * ncv;
    int rvec = 1; // Compute eigenvectors
    char bmat = 'I';
    char which[] = "SA"; // Compute smallest eigenvalue
    int ido = 0;
    int info = 0;

    // Initialize iparam
    iparam[0] = 1;   // Use exact shift
    iparam[2] = 600; // Maximum number of iterations
    iparam[6] = 1;   // Mode 1: A*x = lambda*x

    while (ido != 99)
    {
        dsaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
        if (ido == -1 || ido == 1)
        {
            matvec(M, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1], n);
        }
    }

    if (info < 0)
    {
        printf("Error with dsaupd, info = %d\n", info);
        return 1;
    }

    dseupd_(&rvec, "A", select, d, v, &n, NULL, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
    if (info != 0)
    {
        printf("Error with dseupd, info = %d\n", info);
        return 1;
    }

    res[0] = d[0];

    free(resid);
    free(v);
    free(workd);
    free(workl);
    free(iparam);
    free(ipntr);
    free(select);
    free(d);
    return 0;
}

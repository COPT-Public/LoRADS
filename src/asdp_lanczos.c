#ifdef HEADERPATH
#include "linalg/def_ asdp_lanczos.h"
#include "linalg/ asdp_lanczos.h"
#include "linalg/vec_opts.h"
#include "linalg/dense_opts.h"
#include "interface/asdp_utils.h"
#else
#include "def_asdp_lanczos.h"
#include "asdp_lanczos.h"
#include "vec_opts.h"
#include "dense_opts.h"
#include "asdp_utils.h"
#endif

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

#include <math.h>

#define SYEV_WORK  (30)
#define SYEV_IWORK (12)

static void dArrSymmetrize( int n, double *dArray ) {
    
    double Aij, Aji;
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = i + 1; j < n; ++j ) {
            Aij = dArray[j * n + i];
            Aji = dArray[i * n + j];
            dArray[j * n + i] = dArray[i * n + j] = (Aij + Aji) * 0.5;
        }
    }
    return;
}

static void ALanczosIPrepare( int n, double *vVec ) {
    
    srand(n);
    for ( int i = 0; i < n; ++i ) {
        srand(rand());
        vVec[i] = sqrt(sqrt((rand() % 1627))) * (rand() % 2 - 0.5);
    }
    
    return;
}

static void ALanczosICleanUp(  asdp_lanczos *ALanczos ) {
    
    ASDP_ZERO(ALanczos->wVec, double, ALanczos->nCol);
    ASDP_ZERO(ALanczos->z1Vec, double, ALanczos->nCol);
    ASDP_ZERO(ALanczos->z2Vec, double, ALanczos->nCol);
    ASDP_ZERO(ALanczos->vaVec, double, ALanczos->nCol);
    ASDP_ZERO(ALanczos->VMat, double, ALanczos->nCol * (ALanczos->nMaxSpaceDim + 1));
    ASDP_ZERO(ALanczos->HMat, double, (ALanczos->nMaxSpaceDim + 1) * (ALanczos->nMaxSpaceDim + 1));
    ASDP_ZERO(ALanczos->YMat, double, ALanczos->nMaxSpaceDim * 2);
    ASDP_ZERO(ALanczos->dArray, double, ALanczos->nMaxSpaceDim);
    ASDP_ZERO(ALanczos->UMat, double, ALanczos->nMaxSpaceDim * ALanczos->nMaxSpaceDim);
    ASDP_ZERO(ALanczos->eigDblMat, double, ALanczos->nMaxSpaceDim * SYEV_WORK);
    ASDP_ZERO(ALanczos->eigIntMat, int, ALanczos->nMaxSpaceDim * SYEV_IWORK);
    
    return;
}

extern asdp_retcode ALanczosCreate(  asdp_lanczos **pALanczos ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    if ( !pALanczos ) {
        retcode = ASDP_RETCODE_FAILED;
        goto exit_cleanup;
    }
    
    asdp_lanczos *ALanczos = NULL;
    ASDP_INIT(ALanczos, asdp_lanczos, 1);
    ASDP_MEMCHECK(ALanczos);
    ASDP_ZERO(ALanczos,  asdp_lanczos, 1);
    *pALanczos = ALanczos;
    
exit_cleanup:
    return retcode;
}

extern asdp_retcode ALanczosInit(  asdp_lanczos *ALanczos, int nCol, int nSpaceDim ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    ALanczos->nCol = nCol;
    ALanczos->nMaxSpaceDim = nSpaceDim;
    
    ASDP_INIT(ALanczos->vVec, double, nCol);
    ASDP_MEMCHECK(ALanczos->vVec);
    
    ASDP_INIT(ALanczos->wVec, double, nCol);
    ASDP_MEMCHECK(ALanczos->wVec);
    
    ASDP_INIT(ALanczos->z1Vec, double, nCol);
    ASDP_MEMCHECK(ALanczos->z1Vec);
    
    ASDP_INIT(ALanczos->z2Vec, double, nCol);
    ASDP_MEMCHECK(ALanczos->z2Vec);
    
    ASDP_INIT(ALanczos->vaVec, double, nCol);
    ASDP_MEMCHECK(ALanczos->vaVec);
    
    ASDP_INIT(ALanczos->VMat, double, nCol * (ALanczos->nMaxSpaceDim + 1));
    ASDP_MEMCHECK(ALanczos->VMat);
    
    ASDP_INIT(ALanczos->HMat, double, (ALanczos->nMaxSpaceDim + 1) * (ALanczos->nMaxSpaceDim + 1));
    ASDP_MEMCHECK(ALanczos->HMat);
    
    ASDP_INIT(ALanczos->YMat, double, ALanczos->nMaxSpaceDim * 2);
    ASDP_MEMCHECK(ALanczos->YMat);
    
    ASDP_INIT(ALanczos->dArray, double, ALanczos->nMaxSpaceDim);
    ASDP_MEMCHECK(ALanczos->dArray);
    
    ASDP_INIT(ALanczos->dLanczosWarmStart, double, ALanczos->nCol);
    ASDP_MEMCHECK(ALanczos->dLanczosWarmStart);
    
    ASDP_INIT(ALanczos->UMat, double, ALanczos->nMaxSpaceDim * ALanczos->nMaxSpaceDim);
    ASDP_MEMCHECK(ALanczos->UMat);
    
    ASDP_INIT(ALanczos->eigDblMat, double, ALanczos->nMaxSpaceDim * SYEV_WORK);
    ASDP_MEMCHECK(ALanczos->eigDblMat);
    
    ASDP_INIT(ALanczos->eigIntMat, int, ALanczos->nMaxSpaceDim * SYEV_IWORK);
    ASDP_MEMCHECK(ALanczos->eigIntMat);
    
exit_cleanup:
    return retcode;
}

extern void ALanczosSetData(  asdp_lanczos *ALanczos, void *MMat, void (*Mvec) (void *, double *, double *) ) {
    
    if ( ALanczos->MMat || ALanczos->Mvec ) {
        return;
    }
    
    ALanczos->MMat = MMat;
    ALanczos->Mvec = Mvec;
    
    return;
}
#ifdef ASDP_LANCZOS_DEBUG
#undef ASDP_LANCZOS_DEBUG
#define ASDP_LANCZOS_DEBUG(format, info) printf(format, info)
#else
#define ASDP_LANCZOS_DEBUG(format, info)
#endif
#define H(i, j) ALanczos->HMat[nHRow * (j) + (i)]
#define V(i, j) ALanczos->VMat[nVRow * (j) + (i)]
extern asdp_retcode ALanczosSolve(  asdp_lanczos *ALanczos, double *LanczosStart, double *dMaxStep ) {
    
    asdp_retcode retcode = ASDP_RETCODE_OK;
    
    /* Prepare the initial point */
    if ( LanczosStart ) {
        ASDP_LANCZOS_DEBUG("Loaded user Lanczos warm-start. %s\n", "");
        ASDP_MEMCPY(ALanczos->vVec, LanczosStart, double, ALanczos->nCol);
    } else {
        if ( ALanczos->nComputed == 0 ) {
            ASDP_LANCZOS_DEBUG("Starting Lanczos from scratch. %s\n", "");
            ALanczosIPrepare(ALanczos->nCol, ALanczos->vVec);
        } else {
            ALanczosICleanUp(ALanczos);
            ASDP_LANCZOS_DEBUG("Loaded Lanczos warm start from previous iteration. %s", "");
            ASDP_MEMCPY(ALanczos->vVec, ALanczos->dLanczosWarmStart, double, ALanczos->nCol);
        }
    }
    
    /* Normalize the starting vector use it as the initial basis */
    normalize(&ALanczos->nCol, ALanczos->vVec);
    ASDP_MEMCPY(ALanczos->VMat, ALanczos->vVec, double, ALanczos->nCol);
    
    /* Configure and start Lanczos iteration */
    int k = 0;
    int LCheckFrequency = (int) ALanczos->nMaxSpaceDim / 5;
    LCheckFrequency = ASDP_MIN(LCheckFrequency, 5);
    
    int nVRow = ALanczos->nCol;
    int nHRow = ALanczos->nMaxSpaceDim + 1;
    
    int ldWork = ALanczos->nMaxSpaceDim * SYEV_WORK;
    int liWork = ALanczos->nMaxSpaceDim * SYEV_IWORK;
    
    for ( k = 0; k < ALanczos->nMaxSpaceDim; ++k ) {
        
        ALanczos->Mvec(ALanczos->MMat, ALanczos->vVec, ALanczos->wVec);
        // double normPrev = nrm2(&ALanczos->nCol, ALanczos->wVec, &AIntConstantOne);
        
        if ( k > 0 ) {
            double negHElem = - H(k, k - 1);
            axpy(&nVRow, &negHElem, &V(0, k - 1), &AIntConstantOne, ALanczos->wVec, &AIntConstantOne);
        }
        
        double vAlp = -dot(&nVRow, ALanczos->wVec, &AIntConstantOne, &V(0, k), &AIntConstantOne);
        axpy(&nVRow, &vAlp, &V(0, k), &AIntConstantOne, ALanczos->wVec, &AIntConstantOne);
        double normPres = nrm2(&nVRow, ALanczos->wVec, &AIntConstantOne);
        
        H(k, k) = - vAlp;
        ASDP_LANCZOS_DEBUG("Lanczos Alp value: %f \n", -vAlp);
        
        ASDP_MEMCPY(ALanczos->vVec, ALanczos->wVec, double, ALanczos->nCol);
        normPres = normalize(&ALanczos->nCol, ALanczos->vVec);
        ASDP_MEMCPY(&V(0, k + 1), ALanczos->vVec, double, ALanczos->nCol);
        H(k + 1, k) = H(k, k + 1) = normPres;
        
        /* Frequently check subspace */
        if ( ( k + 1 ) % LCheckFrequency == 0 || k > ALanczos->nMaxSpaceDim - 1 ) {
            
            ASDP_LANCZOS_DEBUG("Entering Lanczos internal check at iteration %d.\n", k);
            
            int kPlus1 = k + 1;
            for ( int i = 0; i < kPlus1; ++i ) {
                ASDP_MEMCPY(ALanczos->UMat + kPlus1 * i, &H(0, i), double, kPlus1);
            }
            
            dArrSymmetrize(kPlus1, ALanczos->UMat);
            ASDP_CALL(fds_syev(kPlus1, ALanczos->UMat, ALanczos->dArray, ALanczos->YMat,
                                ALanczos->eigDblMat, ALanczos->eigIntMat, ldWork, liWork));
            
            double resiVal = fabs( H(kPlus1, k) * ALanczos->YMat[kPlus1 + k] );
            ASDP_LANCZOS_DEBUG("Lanczos outer resi value: %f \n", resiVal);
            
            if ( resiVal < 1e-05 || k >= ALanczos->nMaxSpaceDim - 1 ) {
                ASDP_LANCZOS_DEBUG("Lanczos inner iteration %d \n", k);
                
                double eigMin1 = ALanczos->dArray[1];
                double eigMin2 = ALanczos->dArray[0];
                
                fds_gemv(nVRow, kPlus1, ALanczos->VMat,
                         ALanczos->YMat + kPlus1, ALanczos->z1Vec);
                ALanczos->Mvec(ALanczos->MMat, ALanczos->z1Vec, ALanczos->z2Vec);
                
                /* Record warm-start */
                ASDP_MEMCPY(ALanczos->dLanczosWarmStart, ALanczos->z2Vec, double, ALanczos->nCol);
                
                double negEig = -eigMin1;
                axpy(&nVRow, &negEig, ALanczos->z1Vec,
                     &AIntConstantOne, ALanczos->z2Vec, &AIntConstantOne);
                
                double resiVal1 = nrm2(&nVRow, ALanczos->z2Vec, &AIntConstantOne);
                fds_gemv(nVRow, kPlus1, ALanczos->VMat,
                         ALanczos->YMat, ALanczos->z2Vec);
                ALanczos->Mvec(ALanczos->MMat, ALanczos->z2Vec, ALanczos->z1Vec);
                axpy(&nVRow, &negEig, ALanczos->z2Vec,
                     &AIntConstantOne, ALanczos->z1Vec, &AIntConstantOne);
                
                /* Compute bound on the stepsize */
                double resiVal2 = nrm2(&nVRow, ALanczos->z1Vec, &AIntConstantOne);
                double resiDiff = eigMin1 - eigMin2 - resiVal2;
                double valGamma = ( resiDiff > 0 ) ? resiDiff : 1e-16;
                double resiVal1sqr = resiVal1 * resiVal1 / valGamma;
                valGamma = ASDP_MIN(resiVal1, resiVal1sqr);
                
                if ( valGamma < 1e-03 || valGamma + eigMin1 <= 0.5 ) {
                    if ( valGamma + eigMin1 <= 0.0 ) {
                        *dMaxStep = ASDP_INFINITY;
                    } else {
                        *dMaxStep = 1.0 / ( valGamma + eigMin1 );
                    }
                    break;
                    
                } else {
                    if ( normPres == 0.0 ) {
                        retcode = ASDP_RETCODE_FAILED;
                        goto exit_cleanup;
                    }
                    *dMaxStep = 1.0 / ( valGamma + eigMin1 );
                }
            }
        }
    }
    
    ALanczos->nComputed += 1;
    
exit_cleanup:
    return retcode;
}

extern void ALanczosClear(  asdp_lanczos *ALanczos ) {
    
    if ( !ALanczos ) {
        return;
    }
    
    ASDP_FREE(ALanczos->vVec);
    ASDP_FREE(ALanczos->wVec);
    ASDP_FREE(ALanczos->z1Vec);
    ASDP_FREE(ALanczos->z2Vec);
    ASDP_FREE(ALanczos->vaVec);
    
    ASDP_FREE(ALanczos->VMat);
    ASDP_FREE(ALanczos->HMat);
    ASDP_FREE(ALanczos->YMat);
    ASDP_FREE(ALanczos->UMat);
    
    ASDP_FREE(ALanczos->dArray);
    ASDP_FREE(ALanczos->dLanczosWarmStart);
    ASDP_FREE(ALanczos->eigDblMat);
    ASDP_FREE(ALanczos->eigIntMat);
    
    ASDP_ZERO(ALanczos,  asdp_lanczos, 1);
    
    return;
}

extern void ALanczosDestroy(  asdp_lanczos **pALanczos ) {
    
    if ( !pALanczos ) {
        return;
    }
    
    ALanczosClear(*pALanczos);
    ASDP_FREE(*pALanczos);
    
    return;
}

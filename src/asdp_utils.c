
#ifdef HEADERPATH
#include "interface/asdp_utils.h"
#include "interface/asdp_debug.h"
#else
#include "asdp_utils.h"
#include "asdp_debug.h"
#endif

#include <math.h>
#ifdef __WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#endif
#include <signal.h>

#ifdef MEMDEBUG
#include "memwatch.h"
#endif

#ifdef __WIN32
static BOOL monitorCtrlC(DWORD fdwCtrlType)
{
    switch (fdwCtrlType)
    {
    case CTRL_C_EVENT:
        exit(0);
        return TRUE;
        break;
    default:
        return FALSE;
        break;
    }
}
#else
static int isCtrlC = 0;
static void monitorCtrlC(int sigNum)
{
    isCtrlC = 1;
    exit(0);
    return;
}
static struct sigaction act;
#endif
/* TODO: Add compatibility for Windows platform */
#ifdef __WIN32
static double my_clock(void)
{
    LARGE_INTEGER frequency;
    LARGE_INTEGER currentTime;

    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&currentTime);

    return (double)currentTime.QuadPart / (double)frequency.QuadPart;
}
#else
static double my_clock(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (1e-06 * t.tv_usec + t.tv_sec);
}
#endif

static int dpartitioni(int *ind, double *val, int l, int h)
{

    double tmp2 = 0.0, p = val[l];
    int tmp = l, tmp3 = 0;

    while (l < h)
    {

        while (l < h && val[h] >= p)
        {
            --h;
        }
        while (l < h && val[l] <= p)
        {
            ++l;
        }

        if (l < h)
        {
            tmp2 = val[l];
            val[l] = val[h];
            val[h] = tmp2;
            tmp3 = ind[l];
            ind[l] = ind[h];
            ind[h] = tmp3;
        }
    }

    tmp2 = val[l];
    val[l] = val[tmp];
    val[tmp] = tmp2;
    tmp3 = ind[l];
    ind[l] = ind[tmp];
    ind[tmp] = tmp3;

    return l;
}

static int ipartitiond(double *ind, int *val, int l, int h)
{

    int tmp2 = 0, p = val[l], tmp = l;
    double tmp3;

    while (l < h)
    {

        while (l < h && val[h] >= p)
        {
            --h;
        }
        while (l < h && val[l] <= p)
        {
            ++l;
        }

        if (l < h)
        {
            tmp2 = val[l];
            val[l] = val[h];
            val[h] = tmp2;
            tmp3 = ind[l];
            ind[l] = ind[h];
            ind[h] = tmp3;
        }
    }

    tmp2 = val[l];
    val[l] = val[tmp];
    val[tmp] = tmp2;
    tmp3 = ind[l];
    ind[l] = ind[tmp];
    ind[tmp] = tmp3;

    return l;
}

static int ipartitioni(int *ind, int *val, int l, int h)
{

    int tmp = l, tmp2 = 0, tmp3, p = val[l];

    while (l < h)
    {

        while (l < h && val[h] <= p)
        {
            --h;
        }
        while (l < h && val[l] >= p)
        {
            ++l;
        }

        if (l < h)
        {
            tmp2 = val[l];
            val[l] = val[h];
            val[h] = tmp2;
            tmp3 = ind[l];
            ind[l] = ind[h];
            ind[h] = tmp3;
        }
    }

    tmp2 = val[l];
    val[l] = val[tmp];
    val[tmp] = tmp2;
    tmp3 = ind[l];
    ind[l] = ind[tmp];
    ind[tmp] = tmp3;

    return l;
}

static int dpartitiond(double *ind, double *val, int l, int h)
{

    int tmp = l;
    double tmp2 = 0.0, tmp3, p = val[l];

    while (l < h)
    {

        while (l < h && val[h] >= p)
        {
            --h;
        }
        while (l < h && val[l] <= p)
        {
            ++l;
        }

        if (l < h)
        {
            tmp2 = val[l];
            val[l] = val[h];
            val[h] = tmp2;
            tmp3 = ind[l];
            ind[l] = ind[h];
            ind[h] = tmp3;
        }
    }

    tmp2 = val[l];
    val[l] = val[tmp];
    val[tmp] = tmp2;
    tmp3 = ind[l];
    ind[l] = ind[tmp];
    ind[tmp] = tmp3;

    return l;
}

extern double AUtilGetTimeStamp(void)
{

    return my_clock();
}

/** @brief Symmetrize an n by n matrix whose lower triangular is filled
 *
 */
extern void AUtilMatSymmetrize(int n, double *v)
{

    for (int i = 0, j; i < n; ++i)
    {
        for (j = i + 1; j < n; ++j)
        {
            FULL_ENTRY(v, n, i, j) = FULL_ENTRY(v, n, j, i);
        }
    }

    return;
}

/* Debugging */
extern void AUtilPrintDblContent(int n, double *d)
{

    for (int i = 0; i < n; ++i)
    {
        printf("%5.3e, ", d[i]);
    }
    printf("\n");
    return;
}

extern void AUtilPrintIntContent(int n, int *d)
{

    for (int i = 0; i < n; ++i)
    {
        printf("%5d, ", d[i]);
    }
    printf("\n");
    return;
}

extern double AUtilPrintDblSum(int n, double *d)
{

    double ds = 0.0;

    for (int i = 0; i < n; ++i)
    {
        ds += d[i];
    }

    return ds;
}

extern double AUtilPrintDblAbsSum(int n, double *d)
{

    double ds = 0.0;

    for (int i = 0; i < n; ++i)
    {
        ds += fabs(d[i]);
    }

    return ds;
}

/* Sorting */
extern int AUtilCheckIfAscending(int n, int *idx)
{
    /* Check is an integer array is ascending. */

    for (int i = 0; i < n - 1; ++i)
    {
        if (idx[i] > idx[i + 1])
        {
            return 0;
        }
    }

    return 1;
}

// extern void AUtilSortIntbyDbl( int *data, double *ref, int low, int up ) {

//     if ( low < up ) {
//         int p = dpartitioni(data, ref, low, up);
//         AUtilSortIntbyDbl(data, ref, low, p - 1);
//         AUtilSortIntbyDbl(data, ref, p + 1, up);
//     }

//     return;
// }

extern void AUtilDescendSortIntByInt(int *data, int *ref, int low, int up)
{

    if (low < up)
    {
        int p = ipartitioni(data, ref, low, up);
        AUtilDescendSortIntByInt(data, ref, low, p - 1);
        AUtilDescendSortIntByInt(data, ref, p + 1, up);
    }

    return;
}

extern void AUtilAscendSortDblByInt(double *data, int *ref, int low, int up)
{

    if (low < up)
    {
        int p = ipartitiond(data, ref, low, up);
        AUtilAscendSortDblByInt(data, ref, low, p - 1);
        AUtilAscendSortDblByInt(data, ref, p + 1, up);
    }

    return;
}

extern void AUtilSortDblByDbl(double *data, double *ref, int low, int up)
{

    if (low < up)
    {
        int p = dpartitiond(data, ref, low, up);
        AUtilSortDblByDbl(data, ref, low, p - 1);
        AUtilSortDblByDbl(data, ref, p + 1, up);
    }

    return;
}
#ifdef __WIN32
extern void AUtilStartCtrlCCheck(void)
{
    SetConsoleCtrlHandler((PHANDLER_ROUTINE)monitorCtrlC, TRUE);
}
#else
extern void AUtilStartCtrlCCheck(void)
{

    act.sa_handler = monitorCtrlC;
    sigaction(SIGINT, &act, NULL);

    return;
}
extern int AUtilCheckCtrlC(void)
{

    return isCtrlC;
}

extern void AUtilResetCtrl(void)
{

    isCtrlC = 0;
}
#endif

// extern int MKL_Get_Max_Threads( void );
// extern int MKL_Set_Num_Threads( int nth );

// extern int AUtilGetGlobalMKLThreads( void ) {
//
//     return MKL_Get_Max_Threads();
// }
//
// extern void AUtilSetGlobalMKLThreads( int nTargetThreads ) {
//
//     MKL_Set_Num_Threads(nTargetThreads);
//
//     return;
// }

extern void ASDP_ONE(double *var, int size)
{
    for (int i = 0; i < size; ++i)
    {
        var[i] = (double)(i + 1) / (double)size;
        //        if (i == 0){
        //            var[i] = 1;
        //        }else{
        //            var[i] = 0;
        //        }
    }
}

// Updates the EMA value and checks for significant decrease only when updating the old EMA.
// current_ema: Pointer to the current EMA value
// old_ema: Pointer to the EMA value from a previous period
// new_value: The new observed value to update the EMA with
// alpha: The smoothing factor for the EMA, typically in the range (0,1)
// threshold: The threshold for determining significant proportion change
// update_interval: The interval at which to update the old EMA value
// counter: The counter for iterations to determine when to update the old EMA
// Returns: Boolean value indicating if there's no significant decrease. Returns true by default if not updating old EMA.
extern int AUtilUpdateCheckEma(double *current_ema, double *old_ema, double new_value, double alpha, double threshold, int update_interval, int *counter)
{
    int result = 1;                                                  // Default to true
    *current_ema = alpha * new_value + (1 - alpha) * (*current_ema); // Update the EMA value

    // Check if it's time to update the old EMA
    if (*counter >= update_interval)
    {
        // Avoid division by zero if old_ema is 0 (e.g., the first update)
        if (*old_ema != 0)
        {
            double change = (*current_ema - *old_ema) / *old_ema;     // Calculate the proportional change
            result = (change >= -threshold) && (change <= threshold); // Determine if change is within threshold bounds
            // printf("ratio: %.4f\n", result);
        }
        *old_ema = *current_ema; // Update the old EMA value
        *counter = 1;            // Reset the counter
    }
    else
    {
        (*counter)++;
    }

    return result;
}


extern void REALLOC(double **data, int nOld, int nNew){
    double *dataNewPtr;
    ASDP_INIT(dataNewPtr, double, nNew);
    ASDP_MEMCPY(dataNewPtr, *data, double, nOld);
    ASDP_FREE(*data);
    *data = dataNewPtr;
}

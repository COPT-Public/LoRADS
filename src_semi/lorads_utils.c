#include "lorads_utils.h"

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

static lorads_int dpartitioni(lorads_int *ind, double *val, lorads_int l, lorads_int h)
{

    double tmp2 = 0.0, p = val[l];
    lorads_int tmp = l, tmp3 = 0;

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

static lorads_int ipartitiond(double *ind, lorads_int *val, lorads_int l, lorads_int h)
{

    lorads_int tmp2 = 0, p = val[l], tmp = l;
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

static lorads_int ipartitioni(lorads_int *ind, lorads_int *val, lorads_int l, lorads_int h)
{

    lorads_int tmp = l, tmp2 = 0, tmp3, p = val[l];

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

static lorads_int dpartitiond(double *ind, double *val, lorads_int l, lorads_int h)
{

    lorads_int tmp = l;
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

extern double LUtilGetTimeStamp(void)
{

    return my_clock();
}

/** @brief Symmetrize an n by n matrix whose lower triangular is filled
 *
 */
extern void LUtilMatSymmetrize(lorads_int n, double *v)
{

    for (lorads_int i = 0, j; i < n; ++i)
    {
        for (j = i + 1; j < n; ++j)
        {
            FULL_ENTRY(v, n, i, j) = FULL_ENTRY(v, n, j, i);
        }
    }

    return;
}

/* Debugging */
extern void LUtilPrintDblContent(lorads_int n, double *d)
{

    for (lorads_int i = 0; i < n; ++i)
    {
        printf("%5.3e, ", d[i]);
    }
    printf("\n");
    return;
}

extern double LUtilPrintDblSum(lorads_int n, double *d)
{

    double ds = 0.0;

    for (lorads_int i = 0; i < n; ++i)
    {
        ds += d[i];
    }

    return ds;
}

extern double LUtilPrintDblAbsSum(lorads_int n, double *d)
{

    double ds = 0.0;

    for (lorads_int i = 0; i < n; ++i)
    {
        ds += fabs(d[i]);
    }

    return ds;
}

/* Sorting */
extern lorads_int LUtilCheckIfAscending(lorads_int n, lorads_int *idx)
{
    /* Check is an lorads_integer array is ascending. */

    for (lorads_int i = 0; i < n - 1; ++i)
    {
        if (idx[i] > idx[i + 1])
        {
            return 0;
        }
    }

    return 1;
}

// extern void LUtilSortIntbyDbl( lorads_int *data, double *ref, lorads_int low, lorads_int up ) {

//     if ( low < up ) {
//         lorads_int p = dpartitioni(data, ref, low, up);
//         LUtilSortIntbyDbl(data, ref, low, p - 1);
//         LUtilSortIntbyDbl(data, ref, p + 1, up);
//     }

//     return;
// }

extern void LUtilDescendSortIntByInt(lorads_int *data, lorads_int *ref, lorads_int low, lorads_int up)
{

    if (low < up)
    {
        lorads_int p = ipartitioni(data, ref, low, up);
        LUtilDescendSortIntByInt(data, ref, low, p - 1);
        LUtilDescendSortIntByInt(data, ref, p + 1, up);
    }

    return;
}

extern void LUtilAscendSortDblByInt(double *data, lorads_int *ref, lorads_int low, lorads_int up)
{

    if (low < up)
    {
        lorads_int p = ipartitiond(data, ref, low, up);
        LUtilAscendSortDblByInt(data, ref, low, p - 1);
        LUtilAscendSortDblByInt(data, ref, p + 1, up);
    }

    return;
}

extern void LUtilSortDblByDbl(double *data, double *ref, lorads_int low, lorads_int up)
{

    if (low < up)
    {
        lorads_int p = dpartitiond(data, ref, low, up);
        LUtilSortDblByDbl(data, ref, low, p - 1);
        LUtilSortDblByDbl(data, ref, p + 1, up);
    }

    return;
}
#ifdef __WIN32
extern void LUtilStartCtrlCCheck(void)
{
    SetConsoleCtrlHandler((PHANDLER_ROUTINE)monitorCtrlC, TRUE);
}
#else
extern void LUtilStartCtrlCCheck(void)
{

    act.sa_handler = monitorCtrlC;
    sigaction(SIGINT, &act, NULL);

    return;
}
extern lorads_int LUtilCheckCtrlC(void)
{

    return isCtrlC;
}

extern void LUtilResetCtrl(void)
{

    isCtrlC = 0;
}
#endif

// extern lorads_int MKL_Get_Max_Threads( void );
// extern lorads_int MKL_Set_Num_Threads( lorads_int nth );

// extern lorads_int LUtilGetGlobalMKLThreads( void ) {
//
//     return MKL_Get_Max_Threads();
// }
//
// extern void LUtilSetGlobalMKLThreads( lorads_int nTargetThreads ) {
//
//     MKL_Set_Num_Threads(nTargetThreads);
//
//     return;
// }

extern void LORADS_ONE(double *var, lorads_int size)
{
    for (lorads_int i = 0; i < size; ++i)
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
// update_interval: The lorads_interval at which to update the old EMA value
// counter: The counter for iterations to determine when to update the old EMA
// Returns: Boolean value indicating if there's no significant decrease. Returns true by default if not updating old EMA.
extern lorads_int LUtilUpdateCheckEma(double *current_ema, double *old_ema, double new_value, double alpha, double threshold, lorads_int update_interval, lorads_int *counter)
{
    lorads_int result = 1;                                                  // Default to true
    *current_ema = alpha * new_value + (1 - alpha) * (*current_ema); // Update the EMA value

    // Check if it's time to update the old EMA
    if (*counter >= update_interval)
    {
        // printf("current_ema: %.8f\n", *current_ema);
        // printf("old_ema: %.8f\n", *old_ema);
        // printf("counter: %d\n", *counter);
        // Avoid division by zero if old_ema is 0 (e.g., the first update)
        if (*old_ema != 0)
        {
            double change = (*current_ema - *old_ema) / *old_ema;     // Calculate the proportional change
            // printf("change: %.8f\n", change);
            result = (change >= -threshold) && (change <= threshold); // Determine if change is within threshold bounds
            // printf("ratio: %.4f\n", result);
        }
        *old_ema = *current_ema; // Update the old EMA value
        // printf("old_ema: %.8f\n", *old_ema);

        *counter = 1;            // Reset the counter
    }
    else
    {
        (*counter)++;
    }

    return result;
}


// extern int AUtilUpdateCheckEma(double *current_ema, double *old_ema, double new_value, double alpha, double threshold, int update_interval, int *counter)
// {
//     int result = 1;                                                  // Default to true
//     *current_ema = alpha * new_value + (1 - alpha) * (*current_ema); // Update the EMA value

//     // Check if it's time to update the old EMA
//     if (*counter >= update_interval)
//     {
//         // Avoid division by zero if old_ema is 0 (e.g., the first update)
//         if (*old_ema != 0)
//         {
//             double change = (*current_ema - *old_ema) / *old_ema;     // Calculate the proportional change
//             result = (change >= -threshold) && (change <= threshold); // Determine if change is within threshold bounds
//             // printf("ratio: %.4f\n", result);
//         }
//         *old_ema = *current_ema; // Update the old EMA value
//         *counter = 1;            // Reset the counter
//     }
//     else
//     {
//         (*counter)++;
//     }

//     return result;
// }

extern void REALLOC(double **data, lorads_int nOld, lorads_int nNew){
    double *dataNewPtr;
    LORADS_INIT(dataNewPtr, double, nNew);
    LORADS_MEMCPY(dataNewPtr, *data, double, nOld);
    LORADS_FREE(*data);
    *data = dataNewPtr;
}

#ifndef LORADS_FILE_IO_H
#define LORADS_FILE_IO_H

/** @file hdsdp\_file\_io.h
 * HDSDP input and output for SDPA format
 */
#ifndef LORADS_FILE_READER_H
#define LORADS_FILE_READER_H
#include "lorads.h"
#include "lorads_utils.h"

lorads_retcode LReadSDPA(
    char *fname, // filename
    lorads_int *pnConstrs, // constralorads_int number, row number of matrix A, the dimension of b
    lorads_int *pnBlks, // block number for semidefinite variable
    lorads_int **pblkDims, // block dimension for semidefinite variable
    double **prowRHS, // dual coefficient, b in Ax=b
    lorads_int ***pconeMatBeg, // cone matrix begin idx
    lorads_int ***pconeMatIdx, // cone matrix idx
    double ***pconeMatElem, // cone matrix element
    lorads_int *pnCols, // dimension of all variables other than LP variables
    lorads_int *pnLPCols, // dimension of LP variables
    // pnCols + pnLPCols = column number of matrix A = the dimension of x = the dimension of c
    lorads_int **pLpMatBeg,
    lorads_int **pLpMatIdx,
    double **pLpMatElem,
    lorads_int *pnElems // number of elements in A other than LP variables
    );

#endif /* LORADS_FILE_READER_H */

#endif /* LORADS_FILE_IO_H */

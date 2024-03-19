#ifndef ASDP_FILE_IO_H
#define ASDP_FILE_IO_H

/** @file hdsdp\_file\_io.h
 * HDSDP input and output for SDPA format
 */
#ifndef ASDP_FILE_READER_H
#define ASDP_FILE_READER_H
#include "asdp.h"

asdp_retcode AReadSDPA(
    char *fname, // filename
    int *pnConstrs, // constraint number, row number of matrix A, the dimension of b
    int *pnBlks, // block number for semidefinite variable
    int **pblkDims, // block dimension for semidefinite variable
    double **prowRHS, // dual coefficient, b in Ax=b
    int ***pconeMatBeg, // cone matrix begin idx
    int ***pconeMatIdx, // cone matrix idx
    double ***pconeMatElem, // cone matrix element
    int *pnCols, // dimension of all variables other than LP variables
    int *pnLPCols, // dimension of LP variables
    // pnCols + pnLPCols = column number of matrix A = the dimension of x = the dimension of c
    int **pLpMatBeg,
    int **pLpMatIdx,
    double **pLpMatElem,
    int *pnElems // number of elements in A other than LP variables
    );

#endif /* ASDP_FILE_READER_H */

#endif /* ASDP_FILE_IO_H */


#include "lorads_utils.h"
#include "lorads_file_io.h"
#include "lorads_cs.h"
#include "lorads.h"

#include <stdio.h>
#include <math.h>


#ifdef MEMDEBUG
#include "memwatch.h"
#endif

#ifndef SDPA_BUFFERSIZE
#define SDPA_BUFFERSIZE 1024
#define read_line(file, target, line) fgets(target, SDPA_BUFFERSIZE, file); \
                                      ++line;
#endif

lorads_retcode LReadSDPA( char *fname, lorads_int *pnConstrs, lorads_int *pnBlks, lorads_int **pblkDims, double **prowRHS,
                       lorads_int ***pconeMatBeg, lorads_int ***pconeMatIdx, double ***pconeMatElem, lorads_int *pnCols, lorads_int *pnLPCols,
                       lorads_int **pLpMatBeg, lorads_int **pLpMatIdx, double **pLpMatElem, lorads_int *pnElems ) {
  
  lorads_retcode retcode = LORADS_RETCODE_OK;
  FILE *file;
  
  char fileLine[SDPA_BUFFERSIZE] = "*";
  
  lorads_int whichLine = 0;
  lorads_int nConstrs = 0;
  lorads_int nCols = 0;
  lorads_int nLPCols = 0;
  lorads_int nCones = 0;
  lorads_int nElems = 0;
  lorads_int iCone = 0;
  lorads_int iCon = 0;
  lorads_int blkDim = 0;
  
  /* Allocated memory */
  lorads_int *coneDims = NULL;
  double *rowRHS = NULL;
  lorads_int **coneMatBeg = NULL;
  lorads_int **coneMatIdx = NULL;
  lorads_int *LpMatBeg = NULL;
  lorads_int *LpMatIdx = NULL;
  double *LpMatElem = NULL;
  double **coneMatElem = NULL;
  
  lorads_int warnTinyEntry = 1;
  
  /* Auxiliary memory */
  dcs *LpConeData = NULL;
  dcs **pSDPConeData = NULL;
  
  /* Open SDPA file for reading */
  file = fopen(fname, "r");
  
  if ( !file ) {
      LORADS_ERROR_TRACE;
      retcode = LORADS_RETCODE_FAILED;
      goto exit_cleanup;
  }
  
  /* Jump through comments */
  while ( !feof(file) && (fileLine[0] == '*' ||
                          fileLine[0] == '"') ) {
      read_line(file, fileLine, whichLine);
  }
  
  /* Get number of constraints */
#ifdef MAC_INT64
    if ( sscanf(fileLine, "%lld", &nConstrs) != 1 ) {
#endif
#ifdef UNIX_INT64
if ( sscanf(fileLine, "%ld", &nConstrs) != 1 ) {
#endif
#ifdef INT32
    if ( sscanf(fileLine, "%d", &nConstrs) != 1 ) {
#endif
      LORADS_ERROR_TRACE;
      retcode = LORADS_RETCODE_FAILED;
      goto exit_cleanup;
  }
  
  read_line(file, fileLine, whichLine);
  
  /* Get number of blocks */
#ifdef INT32
  if ( sscanf(fileLine, "%d", &nCones) != 1 ) {
#endif
#ifdef MAC_INT64
    if ( sscanf(fileLine, "%lld", &nCones) != 1 ) {
#endif
#ifdef UNIX_INT64
    if ( sscanf(fileLine, "%ld", &nCones) != 1 ) {
#endif
      LORADS_ERROR_TRACE;
      retcode = LORADS_RETCODE_FAILED;
      goto exit_cleanup;
  }
  
  /* Get all the dimensions */
  ++whichLine;
  LORADS_INIT(coneDims, lorads_int, nCones);
  LORADS_MEMCHECK(coneDims);
  
  for ( iCone = 0; iCone < nCones - 1; ++iCone ) {
      if ( fscanf(file, "{") == 1 || fscanf(file, "(") == 1 ||
           fscanf(file, "'") == 1 ) {
          --iCone;
#ifdef INT32
      } else if ( fscanf(file, "%d", &blkDim) == 1 ) {
#endif
#ifdef UNIX_INT64
      } else if ( fscanf(file, "%ld", &blkDim) == 1 ) {
#endif
#ifdef MAC_INT64
      } else if ( fscanf(file, "%lld", &blkDim) == 1 ) {
#endif
          if ( blkDim <= 0 ) {
              /* Reason: only one diagonal block is supported and
                         it must appear at the end
               */
              LORADS_ERROR_TRACE;
              retcode = LORADS_RETCODE_FAILED;
              goto exit_cleanup;
          }
          nCols += blkDim * blkDim;
          coneDims[iCone] = blkDim;
          
      } else {
          LORADS_ERROR_TRACE;
          retcode = LORADS_RETCODE_FAILED;
          goto exit_cleanup;
      }
  }
  
  /* Special treatment for LP */
#ifdef INT32
  if ( fscanf(file, "%d", &blkDim) == 1 ) {
#endif
#ifdef UNIX_INT64
      if ( fscanf(file, "%ld", &blkDim) == 1 ) {
#endif
#ifdef MAC_INT64
    if ( fscanf(file, "%lld", &blkDim) == 1 ) {
#endif
      if ( blkDim < 0 ) {
          --nCones;
          nLPCols = -blkDim;
      } else {
          coneDims[iCone] = blkDim;
          nCols += blkDim * blkDim;
      }
  }
  
  /* Read dual objective / primal RHS */
  read_line(file, fileLine, whichLine);
  LORADS_INIT(rowRHS, double, nConstrs);
  LORADS_MEMCHECK(rowRHS);
  
  lorads_int iCol;
  lorads_int iRow;
  double dElem;
  char charTmp;
  
  for ( iRow = 0; iRow < nConstrs; ++iRow ) {
      
      if ( fscanf(file, ",") == 1 ) {
          --iRow;
          continue;
      }
      
      while ( fscanf(file, "%lg", &dElem) != 1 ) {
          fscanf(file, "%c", &charTmp);
          if ( charTmp == '\n' ) {
              LORADS_ERROR_TRACE;
              retcode = LORADS_RETCODE_FAILED;
              goto exit_cleanup;
          }
      }
      
      rowRHS[iRow] = dElem;
  }
  
  /* Reset */
  fgets(fileLine, SDPA_BUFFERSIZE, file);
  lorads_int lineNum = whichLine;
  fseek(file, 0, SEEK_SET);
  whichLine = 0;
  
  for ( lorads_int iLine = 0; iLine < lineNum; ++iLine ) {
      charTmp = '*';
      while ( charTmp != '\n' ) {
          fscanf(file, "%c", &charTmp);
          ++whichLine;
      }
  }
  
  /* Read conic data via CSparse */
  int isMemOK = 1;
  
  if ( nLPCols > 0 ) {
      LpConeData = dcs_spalloc(nLPCols, nConstrs + 1, nConstrs, 1, 1);
      LORADS_MEMCHECK(LpConeData);
  }
  
  LORADS_INIT(pSDPConeData, dcs *, nCones);
  LORADS_MEMCHECK(pSDPConeData);
  
  for ( iCone = 0; iCone < nCones; ++iCone ) {
      pSDPConeData[iCone] = dcs_spalloc(PACK_NNZ(coneDims[iCone]), nConstrs + 1, 100, 1, 1);
      LORADS_MEMCHECK(pSDPConeData[iCone]);
  }
  
  lorads_int LpConeID = -1;
  if ( nLPCols > 0 ) {
      LpConeID = nCones;
  }
  
  while( !feof(file) ) {
      fileLine[0] = '\0';
      read_line(file, fileLine, whichLine);
#ifdef INT32
      if ( sscanf(fileLine, "%d %d %d %d %lg", &iCon, &iCone, &iRow, &iCol, &dElem) != 5 ) {
#endif
#ifdef UNIX_INT64
      if ( sscanf(fileLine, "%ld %ld %ld %ld %lg", &iCon, &iCone, &iRow, &iCol, &dElem) != 5 ) {
#endif
#ifdef MAC_INT64
      if ( sscanf(fileLine, "%lld %lld %lld %lld %lg", &iCon, &iCone, &iRow, &iCol, &dElem) != 5 ) {
#endif
          
          if ( feof(file) || strcmp(fileLine, "BEGIN.COMMENT  \n") == 0 ) {
              break;
          } else {
              LORADS_ERROR_TRACE;
              retcode = LORADS_RETCODE_FAILED;
              goto exit_cleanup;
          }
          
      } else {
      
          /* To 0-based indexing */
          iCone -= 1;
          iRow -= 1;
          iCol -= 1;
          
          if ( fabs(dElem) < 1e-12 ) {
              if ( warnTinyEntry ) {
                  printf("[Warning] Entry smaller than 1e-12 is ignored. \n");
                  warnTinyEntry = 0;
              }
              continue;
          }
          
          if ( iCone == LpConeID ) {
              
              if ( iCon == 0 ) {
                  dElem = -dElem;
              }
              
              isMemOK = dcs_entry(LpConeData, iRow, iCon, dElem);
              if ( !isMemOK ) {
                  LORADS_ERROR_TRACE;
                  retcode = LORADS_RETCODE_MEMORY;
                  return retcode;
              }
              
          } else {
              
              if ( iRow > iCol ) {
                  lorads_int iTmp = iRow;
                  iRow = iCol;
                  iCol = iTmp;
              }

              if ( iCon == 0 ) {
                  dElem = -dElem;
              }
              
              isMemOK = dcs_entry(pSDPConeData[iCone], PACK_IDX(coneDims[iCone], iCol, iRow), iCon, dElem);
              if ( !isMemOK ) {
                  LORADS_ERROR_TRACE;
                  retcode = LORADS_RETCODE_MEMORY;
                  return retcode;
              }
          }
          
          nElems += 1;
      }
  }
  
  if ( !isMemOK ) {
      LORADS_ERROR_TRACE;
      retcode = LORADS_RETCODE_MEMORY;
      return retcode;
  }
  
  /* Compress data and complete conversion */
  dcs *cscMat = NULL;
  for ( iCone = 0; iCone < nCones; ++iCone ) {
      cscMat = dcs_compress(pSDPConeData[iCone]);
      LORADS_MEMCHECK(cscMat);
      dcs_spfree(pSDPConeData[iCone]);
      pSDPConeData[iCone] = cscMat;
  }
  
  if ( nLPCols > 0 ) {
      cscMat = dcs_compress(LpConeData);
      dcs_spfree(LpConeData);
      LpConeData = cscMat;
  }
  
  /* Copy the data out */
  LORADS_INIT(coneMatBeg, lorads_int *, nCones);
  LORADS_INIT(coneMatIdx, lorads_int *, nCones);
  LORADS_INIT(coneMatElem, double *, nCones);
  
  LORADS_MEMCHECK(coneMatBeg);
  LORADS_MEMCHECK(coneMatIdx);
  LORADS_MEMCHECK(coneMatElem);
  
  for ( iCone = 0; iCone < nCones; ++iCone ) {
      
      cscMat = pSDPConeData[iCone];
      LORADS_INIT(coneMatBeg[iCone], lorads_int, cscMat->n + 1);
      LORADS_INIT(coneMatIdx[iCone], lorads_int, cscMat->p[cscMat->n]);
      LORADS_INIT(coneMatElem[iCone], double, cscMat->p[cscMat->n]);
      
      LORADS_MEMCHECK(coneMatBeg[iCone]);
      LORADS_MEMCHECK(coneMatIdx[iCone]);
      LORADS_MEMCHECK(coneMatElem[iCone]);
      
      LORADS_MEMCPY(coneMatBeg[iCone], cscMat->p, lorads_int, cscMat->n + 1);
      LORADS_MEMCPY(coneMatIdx[iCone], cscMat->i, lorads_int, cscMat->p[cscMat->n]);
      LORADS_MEMCPY(coneMatElem[iCone], cscMat->x, double, cscMat->p[cscMat->n]);
      
  }
  
  if ( nLPCols > 0 ) {
      
      LORADS_INIT(LpMatBeg, lorads_int, LpConeData->n + 1);
      LORADS_INIT(LpMatIdx, lorads_int, LpConeData->p[LpConeData->n]);
      LORADS_INIT(LpMatElem, double, LpConeData->p[LpConeData->n]);
      
      LORADS_MEMCHECK(LpMatBeg);
      LORADS_MEMCHECK(LpMatIdx);
      LORADS_MEMCHECK(LpMatElem);
      
      LORADS_MEMCPY(LpMatBeg, LpConeData->p, lorads_int, LpConeData->n + 1);
      LORADS_MEMCPY(LpMatIdx, LpConeData->i, lorads_int, LpConeData->p[cscMat->n]);
      LORADS_MEMCPY(LpMatElem, LpConeData->x, double, LpConeData->p[cscMat->n]);
  }
   
  /* Get the outputs */
  *pnConstrs = nConstrs;
  *pnBlks = nCones;
  *pblkDims = coneDims;
  *prowRHS = rowRHS;
  *pconeMatBeg = coneMatBeg;
  *pconeMatIdx = coneMatIdx;
  *pconeMatElem = coneMatElem;
  *pnCols = nCols;
  *pnLPCols = nLPCols;
  *pLpMatBeg = LpMatBeg;
  *pLpMatIdx = LpMatIdx;
  *pLpMatElem = LpMatElem;
  *pnElems = nElems;
  
exit_cleanup:
  
  /* Auxiliary memory, always freed */
  dcs_spfree(LpConeData);
  
  if ( pSDPConeData ) {
      for ( iCone = 0; iCone < nCones; ++iCone ) {
          dcs_spfree(pSDPConeData[iCone]);
      }
      LORADS_FREE(pSDPConeData);
  }
  
  if ( retcode != LORADS_RETCODE_OK ) {
      /* Free all the internal memories if failed */
      LORADS_FREE(coneDims);
      LORADS_FREE(rowRHS);
      
      LORADS_FREE(LpMatBeg);
      LORADS_FREE(LpMatIdx);
      LORADS_FREE(LpMatElem);
      
      if ( coneMatBeg ) {
          for ( iCone = 0; iCone < nCones; ++iCone ) {
              LORADS_FREE(coneMatBeg[iCone]);
          }
          LORADS_FREE(coneMatBeg);
      }
      
      if ( coneMatIdx ) {
          for ( iCone = 0; iCone < nCones; ++iCone ) {
              LORADS_FREE(coneMatIdx[iCone]);
          }
          LORADS_FREE(coneMatIdx);
      }
      
      if ( coneMatElem ) {
          for ( iCone = 0; iCone < nCones; ++iCone ) {
              LORADS_FREE(coneMatElem[iCone]);
          }
          LORADS_FREE(coneMatElem);
      }
      
  }
  
  return retcode;
}


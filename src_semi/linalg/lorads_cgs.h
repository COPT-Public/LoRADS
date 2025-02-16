#ifndef LORADS_CG_H
#define LORADS_CG_H


#include "def_lorads_elements.h"
#include "def_lorads_cgs.h"
#include "lorads_solver.h"

typedef struct
{
    lorads_sdp_cone *ACone;
    lorads_sdp_dense *noUpdateVar;
    lorads_sdp_dense *UpdateVarShell;
    double *weight;
} admmCG;

extern void LORADSCGSolverCreate(lorads_cg_linsys **pCGSolver, lorads_int blkDim, lorads_int rank, lorads_int nConstr);
extern void LORADSCGSolverReCreate(lorads_cg_linsys **pCGSolver, lorads_int blkDim, lorads_int rank, lorads_int nConstr);
extern void CGSolve(void *linSys, double *x, double *b, double cg_tol, lorads_int cg_maxIter);
extern void CGSetData (lorads_cg_linsys *cg, void *MMat, void (*Mvec) (void *, double *, double *));
extern void CGSolverClear(void *pCGSolver);
extern void LORADSCGDestroy(lorads_solver *ASolver);
#endif //LORADS_CG_H

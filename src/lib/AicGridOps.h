#ifndef AIC_GRID_OPS_H
#define AIC_GRID_OPS_H

// this actually was in AicFB.h
#include "AicCommon.h"
#include "AicDrv.h"
#include "CxFortranInt.h"

#ifndef RESUME_CLOCK
   #define RESUME_CLOCK(Index)
   #define PAUSE_CLOCK(Index)
   #define RESET_CLOCKS
#endif

namespace aic {

typedef std::size_t
   FCoMask; // largest native integer.

void EvalShellGroupOnGrid( double *pOut, uint nCompStride, uint nBfStride,
   FORTINT *pMap, FORTINT &nMap, uint iFnBase,
   double const *pGridPos, uint GridStride, uint nGridPts,
   FShellData const *pFirstBf, FShellData const *pLastBf, FVector3 const &vBfPos,
   uint DerivOrder, double ThrOrb, double LogThrOrb,
   FMemoryStack2 &Mem );

}

#endif // AIC_GRID_OPS_H


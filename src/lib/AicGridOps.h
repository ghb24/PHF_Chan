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


extern "C" {
// Find the effective range of a contracted basis function shell by explicitly
// evaluating the basis functions on a radial grid.
//
//  pOut[iCo] = (range r for which the integral Int[r,Infty] mu^2(r) d^3r <= ThrDen).
//
// i.e., ThrDen is a threshold on the neglected *electron number*.
//
// pWork must have room for at least nResolution x (2+nCo) scalars.
#define AIC_FIND_BFN_RANGE FORT_Extern(aic_find_bfn_range,AIC_FIND_BFN_RANGE)
void AIC_FIND_BFN_RANGE(double *AIC_RP pOut,
   double *AIC_RP pExp, FORTINT const &nExp, double *AIC_RP pCo, FORTINT const &nCo, FORTINT const &l,
   double const &ThrDen, FORTINT const &nResolution, double *pWork, FORTINT const &nWork);

// Evaluate the radial part of a contracted Gauss function:
//    Out = r^l * \sum_i c[i] Exp[-zeta_i r^2]
// for a number of (input) radii. Output is a nRadii x nCo matrix,
// pOut[iPt + nRadii*iCo] giving the value of the radial basis function at point
// iPt for function iCo.
#define AIC_EVAL_BFN_RADIAL FORT_Extern(aic_eval_bfn_radial,AIC_EVAL_BFN_RADIAL)
void AIC_EVAL_BFN_RADIAL(double *AIC_RP pOut,
   double *AIC_RP pExp, FORTINT const &nExp, double *AIC_RP pCo, FORTINT const &nCo, FORTINT const &l,
   double *AIC_RP pRadii, FORTINT const &nRadii);
} // extern "C"

} // namespace aic

#endif // AIC_GRID_OPS_H


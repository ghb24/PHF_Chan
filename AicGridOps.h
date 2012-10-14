#ifndef AIC_GRID_OPS_H
#define AIC_GRID_OPS_H

// this actually was in AicFB.h
#include "AicCommon.h"
#include "CxFortranInt.h"

#ifndef RESUME_CLOCK
   #define RESUME_CLOCK(Index)
   #define PAUSE_CLOCK(Index)
   #define RESET_CLOCKS
#endif

namespace aic {

typedef std::size_t
   FCoMask; // largest native integer.

struct FShellData {
   uint
      l, // angular momentum
      nCo, // number of contractions
      nPrim, // number of primitive shells (==number of exponents)
      ShOff; // primitive exponent offset in screening test density
   double const
      *pExp, // primitive exponents
      *pCo;  // nExp x nCo contraction matrix
   FVector3  // position
      vCenter;

//    FShellData(uint iGrp, FMolproBasis const &Basis);
   inline FShellData(uint l_, uint nCo_, uint nPrim_, uint ShOff_, double const *pExp_, double const *pCo_, FVector3 const &vCenter_);
   inline FShellData(uint l_, uint nCo_, uint nPrim_, uint ShOff_, double const *pExp_, double const *pCo_, double const *pvCenter_);
//    inline FShellData(FGaussBfn const &Bf, FVector3 const &vCenter_); // AicDrv.cpp
   uint nFn() const { return nCo * (2*l + 1); };
#ifdef _DEBUG
   void Print() const;
#endif // _DEGUG
};

FShellData::FShellData(uint l_, uint nCo_, uint nPrim_, uint ShOff_, double const *pExp_, double const *pCo_, FVector3 const &vCenter_)
   : l(l_), nCo(nCo_), nPrim(nPrim_), ShOff(ShOff_), pExp(pExp_), pCo(pCo_), vCenter(vCenter_)
{
}

FShellData::FShellData(uint l_, uint nCo_, uint nPrim_, uint ShOff_, double const *pExp_, double const *pCo_, double const *pvCenter_)
   : l(l_), nCo(nCo_), nPrim(nPrim_), ShOff(ShOff_), pExp(pExp_), pCo(pCo_), vCenter(pvCenter_[0], pvCenter_[1], pvCenter_[2])
{
}

void EvalShellGroupOnGrid( double *pOut, uint nCompStride, uint nBfStride,
   FORTINT *pMap, FORTINT &nMap, uint iFnBase,
   double const *pGridPos, uint GridStride, uint nGridPts,
   FShellData const *pFirstBf, FShellData const *pLastBf, FVector3 const &vBfPos,
   uint DerivOrder, double ThrOrb, double LogThrOrb,
   FMemoryStack2 &Mem );

}

#endif // AIC_GRID_OPS_H


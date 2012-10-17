/* Copyright (c) 2012  Gerald Knizia
 * 
 * This file is part of the bfint program
 * (See http://www.theochem.uni-stuttgart.de/~knizia)
 * 
 * bfint is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * bfint is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with bfint (LICENSE). If not, see http://www.gnu.org/licenses/
 */

#ifndef _AIC_DRV_H
#define _AIC_DRV_H

#include <boost/intrusive_ptr.hpp>
#include <vector>

#include "AicCommon.h"
#include "AicShells.h"
#include "AicKernels.h"

// Behold the interface for the mighty AIC integral core...
// This is a 3-index Obara-Saika core modelled after
//
//    PCCP 6, 5119 (2004)   3-index integrals
//    PCCP 8, 3072 (2006)   4-index integrals, general integral kernels
//    PCCP 10, 3390 (2008)  F12 kernels (F,F/r,[[Del, F],F])
//
// The code in AicOsrr.cpp/.h and AicSolidHarmonics.cpp/.h is generated
// and should not be modified manually. If you want to change these, please
// contact me for the python scripts generating them.
//
// The main interface for Molpro's batcher (in AicDrvFB.cpp) bypasses
// everything in AicShells.* and furthermore the EvalInt2e3c routine in
// AicDrv.cpp. It re-implementes a variant of this routine. This was
// necessary due to Fortran interfacing issues.
//
// Gerald Knizia/cgk, Jun, 2010

namespace aic {

// run FIntegralFactory functions through AicDrvFB.cpp routines instead of the
// native c++ versions implemented here. The main difference between the two is
// that the Fortran interface routines were actually optimized (due to their use
// as Molpro 3-index integral driver), while the (original) native ones were
// my first attempt at making something work and are much simpler.
// Both use the same algebraic kernel routines (OsrrA-C).
// #define RUN_THROUGH_FORTRAN_INTERFACE_CORE_ROUTINES


struct FScreeningParams {
   enum FScreenTypeAC {
      AC_Overlap = 0, // (a|c) is screened on delta(r12) (=overlap criterion). ParamAC ignored.
      AC_Slater = 1,  // (a|c) is screened on exp(-gamma * r12). ParamAC is used as gamma.
      AC_None = 0xff // no screening on (a|c) due to pointlessness (e.g., for coulomb kernel). ParamAC ignored.
   };

   FScreenTypeAC
      TypeAC;
   double
      LogThrAB, // = -log(ThrAB)
      LogThrAC, // = -log(ThrAC)
      ParamAC,  // See comments on FScreenTypeAC for meaning.
      ThrAB;

   FScreeningParams() {};
   FScreeningParams(double ThrAB_, FScreenTypeAC TypeAC_, double ThrAC_ = 1e-15, double ParamAC_ = 0) {
      Set(ThrAB_, TypeAC_, ThrAC_, ParamAC_);
   }
   void Set(double ThrAB, FScreenTypeAC TypeAC, double ThrAC = 0, double ParamAC = 0);
   void SetViaLog(double LogThrAB, FScreenTypeAC TypeAC, double LogThrAC = 0, double ParamAC = 0);
};

#ifdef RUN_THROUGH_FORTRAN_INTERFACE_CORE_ROUTINES
   struct FShEvalData {
      uint StrideA, StrideB, StrideC, StrideDeriv;
      FIntegralKernel const *pKernel;
      double mutable zcalc; // number of calculated integrals.
      double *pTestDensAB, *pTestDensAC;
      uint StrideTd;
      FScreeningParams *pScreenParams;
      double OutputFactor; // integrals will be multiplied by this.
   };
#endif


struct FIntegralFactory
{
   explicit FIntegralFactory( FIntegralKernel const *pKernel )
      : m_pKernel(pKernel)
   {
      Init(1);
   };

   explicit FIntegralFactory( FIntegralKernel const *pKernel, FScreeningParams const &ScreeningParams )
      : m_pKernel(pKernel), m_ScreeningParams(ScreeningParams)
   {
      Init(0);
   };

   // WARNING: this keeps a reference to the kernel object, it does not make a copy.
   explicit FIntegralFactory( FIntegralKernel const &rKernel )
      : m_pKernel(&rKernel)
   {
      Init(1);
   };

   // WARNING: this keeps a reference to the kernel object, it does not make a copy.
   explicit FIntegralFactory( FIntegralKernel const &rKernel, FScreeningParams const &ScreeningParams )
      : m_pKernel(&rKernel), m_ScreeningParams(ScreeningParams)
   {
      Init(0);
   };


   void SetScreeningParams(FScreeningParams const &ScreeningParams) {
      m_ScreeningParams = ScreeningParams;
   }

   // evaluates two-electron, three-center integrals (ab|ker|c),
   // multiplied by scalar OutputFactor.
   // Strides[0..2] are strides for components (a..c), respectively.
   void EvalInt2e3c( double *pOut, uint Strides[3], double OutputFactor,
      FGaussShell const &ShellA, FGaussShell const &ShellB,
      FGaussShell const &ShellC, FMemoryStack &Mem );

   void EvalInt2e2c( double *pOut, uint Strides[2], double OutputFactor,
      FGaussShell const &ShellA, FGaussShell const &ShellC, FMemoryStack &Mem );

   // evaluate two-electron, three-center gradient integrals d/d[Ax,..,Cz] (ab|ker|c),
   // multiplied by scalar OutputFactor. Stride[4] is the x/y/z component derivative stride.
   void EvalGrad2e3c( double *pOutA, double *pOutB, double *pOutC, uint Strides[4], double OutputFactor,
      FGaussShell const &ShellA, FGaussShell const &ShellB,
      FGaussShell const &ShellC, FMemoryStack &Mem );

   // evaluate two-electron, two-center gradient integrals d/d[Ax,..,Cz] (a|ker|c),
   // multiplied by scalar OutputFactor. Stride[3] is the x/y/z component derivative stride.
   void EvalGrad2e2c( double *pOutA, double *pOutC, uint Strides[3], double OutputFactor,
      FGaussShell const &ShellA, FGaussShell const &ShellC, FMemoryStack &Mem );

   FIntegralKernel const &GetKernel() const { return *m_pKernel; };
private:
   FIntegralKernel const
      *m_pKernel;
   void Init(int MakeScreeningParams);
   FScreeningParams
      m_ScreeningParams;
#ifdef RUN_THROUGH_FORTRAN_INTERFACE_CORE_ROUTINES
   FShEvalData
      m_EvalData;
#endif // RUN_THROUGH_FORTRAN_INTERFACE_CORE_ROUTINES
};


inline uint nCartX(uint l){
   // number of cartesian monomials of xyz with degree <= l.
   // maybe store this?
   return (l+1)*(l+2)*(l+3)/6;
}

inline uint nCartY(uint l){
   // number of cartesian monomials of xyz with degree == l.
   return (l+1)*(l+2)/2;
}

struct FGaussBfn;
struct FMolproBasis;

// an internal interface for low-level drivers. probably should not be here
// (all the other stuff using this is in AicFB.h, but we don't actually have that here.)
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

   FShellData(uint iGrp, FMolproBasis const &Basis);
   inline FShellData(uint l_, uint nCo_, uint nPrim_, uint ShOff_, double const *pExp_, double const *pCo_, FVector3 const &vCenter_);
   inline FShellData(uint l_, uint nCo_, uint nPrim_, uint ShOff_, double const *pExp_, double const *pCo_, double const *pvCenter_);
   inline FShellData(FGaussBfn const &Bf, FVector3 const &vCenter_); // AicDrv.cpp
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


// evaluate 2-electron 2-center integrals <a|krn|c>.
// if add is given: increment the output instead of overwriting it.
void EvalInt2e2c( double *pOut, uint *Strides,
    FShellData &ShA, FShellData &ShC, double Prefactor, bool Add,
    FIntegralKernel const *pKernel, FMemoryStack &Mem );

// evaluate 2-electron 2-center integrals <a|krn * laplace|c>
// note: to obtain the kinetic energy operator, pass an overlap kernel
//       and supply -.5 as Prefactor (ekin = -.5 laplace).
// if add is given: increment the output instead of overwriting it.
void EvalInt2e2c_LaplaceC( double *pOut, uint *Strides,
    FShellData &ShA, FShellData &ShC, double Prefactor, bool Add,
    FIntegralKernel const *pKernel, FMemoryStack &Mem );

// // r[i] += f * x[i]
// template<class FScalar>
// void Add2( FScalar *AIC_RP r, FScalar const *AIC_RP x, FScalar f, std::size_t n );
//
// // r[i] += f * x[i] * y[i]
// template<class FScalar>
// void Add2( FScalar *AIC_RP r, FScalar const *AIC_RP x, FScalar const *AIC_RP y, FScalar f, std::size_t n );
//
// ^- suncc is confused by the templates, and reports multiple definitions for
// completely unrelated symbols if these are used from a different directory
// than AicFB.h resides in; atm we only need it for doubles anyway here.

// r[i] += f * x[i]
void Add2( double *AIC_RP r, double const *AIC_RP x, double f, std::size_t n );

// r[i] += f * x[i] * y[i]
void Add2( double *AIC_RP r, double const *AIC_RP x, double const *AIC_RP y, double f, std::size_t n );



} // namespace aic

#endif // _AIC_DRV_H

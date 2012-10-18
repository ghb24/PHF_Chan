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

#include <cmath>
#include <stdexcept> // for std::runtime_error (for gauss shell interface)
#include <boost/intrusive_ptr.hpp>
#include <string.h> // for memset
#include <stdint.h> // for uint32_t vs uint64_t (used for component presence flags)

#include "AicCommon.h"
#include "AicDrv.h"
#include "AicDefs.h"
#include "AicShells.h"
#include "AicOsrr.h"
#include "AicMdrr.h"
#include "AicSolidHarmonics.h"

#include <stdio.h>

void PrintMatrix(double const *pData, uint nRows, uint RowStride, uint nCols, uint ColStride, int wd, int prec, char const *pName = 0);


namespace aic {

#include "AicIndices.h"

void OsrrBX_3cen( double *pOut, uint nStrideOut, double *pIn, double const PmQ[3], double rho,
   double InvEtaABC, double InvEtaCD, uint lab, uint la, uint lc, FMemoryStack &Mem ) AIC_NO_THROW
{
   if ( lc == 0 ) {
      // note: this case does not usually occurr, since for lc=0
      // OsrrB is not called in the first place.
      uint iBaseCompA = nCartX(la-1);
      memcpy(pOut + iBaseCompA, pIn + iBaseCompA, sizeof(double) * (nCartX(lab) - iBaseCompA));
      return;
   }

   double
      riz = rho * InvEtaCD,
      *pMemBase, *pTm0, *pTm1;
   FOsrrParamsB
      P;
   P.PmQf[0] = riz*PmQ[0];
   P.PmQf[1] = riz*PmQ[1];
   P.PmQf[2] = riz*PmQ[2];
   P.InvEtaABC2 = 0.5*InvEtaABC;
   P.AngMomAB_Max = lab;

   // note: The following special cases help mainly when lab is very small.
   //   in that case the (small!) overhead of the generic code below is
   //   still be noticeable. Since in a typical basis set there are *lots*
   //   of s and p primitives, this makes a difference.
   if ( lc == 1 ) {
      P.AngMomAB_Min = la;
      // note: this assumes that the l=1 cartesian components are ordered as x,y,z!
      assert(iCartPow[1][0] == 1 && iCartPow[2][1] == 1 && iCartPow[3][2] == 1);
      OsrrKrnB_3cen( &pOut[0*nStrideOut], &pIn[0], 0, P );
      OsrrKrnB_3cen( &pOut[1*nStrideOut], &pIn[0], 1, P );
      OsrrKrnB_3cen( &pOut[2*nStrideOut], &pIn[0], 2, P );
      return;
   }
   if ( lc == 2 ) {
      // note: this assumes the following cartesian component ordering:
      //       l=1: x,y,z
      //       l=2: x^2,y^2,z^2, xy, xz, yz
      assert(iCartPow[1][0] == 1 && iCartPow[2][1] == 1 && iCartPow[3][2] == 1);
      assert(iCartPow[4][0] == 2 && iCartPow[5][1] == 2 && iCartPow[6][2] == 2);
      assert(iCartPow[7][0] == 1 && iCartPow[8][0] == 1 && iCartPow[9][1] == 1);
      assert(iCartPow[7][1] == 1 && iCartPow[8][2] == 1 && iCartPow[9][2] == 1);
      uint
         nCompA = nCartX(P.AngMomAB_Max);
      Mem.Alloc(pMemBase, nCompA * 3);
      P.AngMomAB_Min = (la == 0)? 0 : (la-1);
      OsrrKrnB_3cen( &pMemBase[0*nCompA], &pIn[0], 0, P );
      OsrrKrnB_3cen( &pMemBase[1*nCompA], &pIn[0], 1, P );
      OsrrKrnB_3cen( &pMemBase[2*nCompA], &pIn[0], 2, P );
      P.AngMomAB_Min = la;
      OsrrKrnB_3cen( &pOut[0*nStrideOut], &pMemBase[0*nCompA], 0, P );
      OsrrKrnB_3cen( &pOut[3*nStrideOut], &pMemBase[1*nCompA], 0, P );
      OsrrKrnB_3cen( &pOut[4*nStrideOut], &pMemBase[2*nCompA], 0, P );
      OsrrKrnB_3cen( &pOut[1*nStrideOut], &pMemBase[1*nCompA], 1, P );
      OsrrKrnB_3cen( &pOut[5*nStrideOut], &pMemBase[2*nCompA], 1, P );
      OsrrKrnB_3cen( &pOut[2*nStrideOut], &pMemBase[2*nCompA], 2, P );
      Mem.Free(pMemBase);
      return;
   }
   if ( lc == 3 ) {
//  {0,0},   {0,1},   {0,2},   {1,0},   {1,1},   {1,2},   {2,0},   {2,1},   {2,2},   {0,5}
//  {3,0,0}, {1,2,0}, {1,0,2}, {2,1,0}, {0,3,0}, {0,1,2}, {2,0,1}, {0,2,1}, {0,0,3}, {1,1,1},
      uint
         nCompA = nCartX(P.AngMomAB_Max);
      Mem.Alloc(pMemBase, nCompA * (3+6));
      P.AngMomAB_Min = (la <= 1)? 0 : (la-2);
      OsrrKrnB_3cen( &pMemBase[0*nCompA], &pIn[0], 0, P );
      OsrrKrnB_3cen( &pMemBase[1*nCompA], &pIn[0], 1, P );
      OsrrKrnB_3cen( &pMemBase[2*nCompA], &pIn[0], 2, P );
      P.AngMomAB_Min = (la == 0)? 0 : (la-1);
      uint const i5 = 3;
      OsrrKrnB_3cen( &pMemBase[(3+0)*nCompA], &pMemBase[0*nCompA], 0, P );
      OsrrKrnB_3cen( &pMemBase[(3+1)*nCompA], &pMemBase[1*nCompA], 1, P );
      OsrrKrnB_3cen( &pMemBase[(3+i5)*nCompA], &pMemBase[2*nCompA], 1, P );
      OsrrKrnB_3cen( &pMemBase[(3+2)*nCompA], &pMemBase[2*nCompA], 2, P );
      P.AngMomAB_Min = la;
      OsrrKrnB_3cen( &pOut[0*nStrideOut], &pMemBase[(3+0)*nCompA], 0, P );
      OsrrKrnB_3cen( &pOut[1*nStrideOut], &pMemBase[(3+1)*nCompA], 0, P );
      OsrrKrnB_3cen( &pOut[2*nStrideOut], &pMemBase[(3+2)*nCompA], 0, P );
      OsrrKrnB_3cen( &pOut[9*nStrideOut], &pMemBase[(3+i5)*nCompA], 0, P );
      OsrrKrnB_3cen( &pOut[3*nStrideOut], &pMemBase[(3+0)*nCompA], 1, P );
      OsrrKrnB_3cen( &pOut[4*nStrideOut], &pMemBase[(3+1)*nCompA], 1, P );
      OsrrKrnB_3cen( &pOut[5*nStrideOut], &pMemBase[(3+2)*nCompA], 1, P );
      OsrrKrnB_3cen( &pOut[6*nStrideOut], &pMemBase[(3+0)*nCompA], 2, P );
      OsrrKrnB_3cen( &pOut[7*nStrideOut], &pMemBase[(3+1)*nCompA], 2, P );
      OsrrKrnB_3cen( &pOut[8*nStrideOut], &pMemBase[(3+2)*nCompA], 2, P );
      Mem.Free(pMemBase);
      return;
   }


   uint
      nCompA = nCartX(P.AngMomAB_Max),
      nCompA0 = nCompA,
      nCompC = nCartY(lc);
   Mem.Alloc(pMemBase, 2 * nCompA * nCompC);
   pTm1 = pIn;
   pTm0 = pMemBase;

   FCompFlag
      // get pointer to C recursion path mask
      *pRequiredCompsC = &OsrrB_TabRequiredCompC[OsrrB_iTabRequiredCompC[lc]]+1;
   uint
      iBaseCompC = 1;
   for ( uint ms = 1; ms <= lc; ++ ms ){ // actual m is lc-m.
      // About As: we need components la..la+lb as output, and here
      //   for each step (lc steps total) one more min la.
      P.AngMomAB_Min = ms + la - lc;
      if ( P.AngMomAB_Min < 0 )
         P.AngMomAB_Min = 0;
      if ( ms == lc ) {
         pTm0 = pOut;
         nCompA0 = nStrideOut; // to get rid of one of the strides in ShTrN.
      }

      // About Cs: For the 3-index integrals, C is finished after this. No lower Cs
      // are needed.
      FCompFlag
         RequiredCompsC = *(pRequiredCompsC++);
      char
         *pCartRD = &iCartRD[iBaseCompC][0];
      for ( uint iCompC = 0; iCompC < nCartY(ms); ++ iCompC, pCartRD += 2 ){
         if ( 0 == (RequiredCompsC & (1 << iCompC)) )
            continue;
         int iDir    = pCartRD[0];
         int iCompC1 = pCartRD[1];
         OsrrKrnB_3cen( &pTm0[iCompC * nCompA0], &pTm1[iCompC1 * nCompA], iDir, P );
      }
      if ( ms == 1 )
         pTm1 = pMemBase + nCompA * nCompC;
      std::swap(pTm0, pTm1);
      iBaseCompC += nCartY(ms);
   }

   Mem.Free(pMemBase);
}

// zero out a strided three-dimensional block of data
void WipeBlock3( double *pOut, uint nA, uint StrideA, uint nB, uint StrideB, uint nC, uint StrideC )
{
   for ( uint iC = 0; iC < nC; ++ iC )
      for ( uint iB = 0; iB < nB; ++ iB )
         for ( uint iA = 0; iA < nA; ++ iA )
            pOut[iC * StrideC + iB * StrideB + iA * StrideA] = 0;
}

void FScreeningParams::Set(double ThrAB_, FScreenTypeAC TypeAC_, double ThrAC_, double ParamAC_)
{
   TypeAC = TypeAC_;
   LogThrAB = -std::log(ThrAB_);
   ThrAB = ThrAB_;
   LogThrAC = -std::log(ThrAC_);
   ParamAC = ParamAC_;
}





} // namespace aic

#ifndef RUN_THROUGH_FORTRAN_INTERFACE_CORE_ROUTINES

#include <sstream>
void AicAssertFail( char const *pExpr, char const *pFile, int iLine )
{
   std::stringstream str;
   str << "ASSERTION FAILED: '" << pExpr << "'";
   if ( pFile != 0 ) str << " at " << pFile << ":" << iLine;
   throw std::runtime_error(str.str());
}

namespace aic {

void FIntegralFactory::Init(int iMakeScreeningParams)
{
}

struct FPrimOverlap
{
   FVector3
      vCen, // weighted center
      vDir; // direction B-A.
   double
      Eta, // ZetaA + ZetaB
      InvEta, // 1/(ZetaA + ZetaB)
      Zeta; // 1/(1/ZetaA + 1/ZetaB)
   FPrimOverlap( FVector3 const &A, double ZetaA, FVector3 const &B, double ZetaB )
   {
      Eta = ZetaA + ZetaB;
      InvEta = 1.0/Eta;
      Zeta = ZetaA * ZetaB * InvEta;
      for ( int i = 0; i < 3; ++ i ){
         vDir[i] = B[i] - A[i];
         vCen[i] = InvEta * ( ZetaA * A[i] + ZetaB * B[i] );
      }
   }
};

void FIntegralFactory::EvalInt2e3c( double *pOut_, uint Strides[3], double OutputFactor,
   FGaussShell const &ShellA_, FGaussShell const &ShellB_,
   FGaussShell const &ShellC, FMemoryStack &Mem )
{
   void
      *pBaseOfMem = Mem.Alloc(0);
   // For the OsrrC recursion we need la >= lb. Switch a and b and the
   // a/b output strides if this is not yet the case.
   bool
      NoSwapAB = (ShellA_.pFn->l >= ShellB_.pFn->l);
   uint
      StrideA = NoSwapAB? Strides[0] : Strides[1],
      StrideB = NoSwapAB? Strides[1] : Strides[0],
      StrideC = Strides[2];
   FGaussShell const
      &ShellA = *(NoSwapAB? &ShellA_ : &ShellB_),
      &ShellB = *(NoSwapAB? &ShellB_ : &ShellA_);
   FGaussBfn const
      &BfA = *ShellA.pFn, &BfB = *ShellB.pFn, &BfC = *ShellC.pFn;
   // copy some shell data into local variables.
   uint
      nShA = (2*BfA.l)+1,  nPrimA = BfA.Exponents.size(),  nCoA = BfA.Contractions.size(),
      nShB = (2*BfB.l)+1,  nPrimB = BfB.Exponents.size(),  nCoB = BfB.Contractions.size(),
      nShC = (2*BfC.l)+1,  nPrimC = BfC.Exponents.size(),  nCoC = BfC.Contractions.size(),
      nCoTotal = nCoA * nCoB * nCoC,
#ifdef _DEBUG
      nPrimTotal = nPrimA * nPrimB * nPrimC,
#endif
      TotalL = BfA.l + BfB.l + BfC.l;
   uint
      nCartY_C = nCartY(BfC.l),
      nCartX_AB = nCartX(BfA.l+BfB.l),
      iBaseCompAB = nCartX(BfA.l-1);
   FVector3 const
      &A = ShellA.vCenter, &B = ShellB.vCenter, &C = ShellC.vCenter,
      BmA = B - A;

#ifdef _DEBUG
   BfA.SanityCheck(); BfB.SanityCheck(); BfC.SanityCheck();
   assert( nCoTotal >= 1 && nPrimTotal >= 1 );
   assert( (nCoA != 1) || (BfA.Contractions.size() == 1 && BfA.Contractions[0].Coeffs.size() == nPrimA) );
   assert( (nCoB != 1) || (BfB.Contractions.size() == 1 && BfB.Contractions[0].Coeffs.size() == nPrimB) );
   assert( (nCoC != 1) || (BfC.Contractions.size() == 1 && BfC.Contractions[0].Coeffs.size() == nPrimC) );
#endif // _DEBUG

   double
       // *pGm: G(m) times Sab * Scd * (pi/(ZetaA+ZetaB+ZetaC+ZetaD))^{3/2}.
       // for all primitives directly behind each other.
      *pGm,
      // intermediate data for primitives.
      *p_A00, *p_A0C, *p_A0C_sh_,
      // intermediate contracted a0c-integrals.
      *c_A0C_sh;

   Mem.Alloc(pGm, TotalL+1);
   Mem.Alloc(p_A00, nCartX_AB);
   Mem.Alloc(p_A0C, nCartY_C * nCartX_AB );
   uint
      nA0C_sh = nCartX_AB * nShC;
   Mem.Alloc(p_A0C_sh_, nA0C_sh);

   // zero out contracted recursion intermediate. This one will
   // accumulate primitive data.
   Mem.Alloc(c_A0C_sh, nA0C_sh * nCoTotal);
   memset(c_A0C_sh, 0, sizeof(*c_A0C_sh) * nA0C_sh * nCoTotal);

   for ( uint iPrimA = 0; iPrimA < nPrimA; ++ iPrimA )
   for ( uint iPrimB = 0; iPrimB < nPrimB; ++ iPrimB )
   {
      FGaussBfn::FContraction const
         &CoA = BfA.PrimContribs[iPrimA],
         &CoB = BfB.PrimContribs[iPrimB];
      double
         ZetaA = BfA.Exponents[iPrimA],
         ZetaB = BfB.Exponents[iPrimB];
      FPrimOverlap
         OvAB(ShellA.vCenter, ZetaA, ShellB.vCenter, ZetaB);
      FVector3 const
         &P = OvAB.vCen,
         PmA = P - A,
         PmQ = P - C; // R
      double
         Sab = exp( -OvAB.Zeta * Dot(BmA,BmA) );
      for ( uint iPrimC = 0; iPrimC < nPrimC; ++ iPrimC )
      {
         double
            ZetaC = BfC.Exponents[iPrimC];
         double
            InvEtaABC = 1.0/(OvAB.Eta + ZetaC),
            rho = (OvAB.Eta * ZetaC) * InvEtaABC,
            Dummy = M_PI * InvEtaABC,
            Prefactor = sqrt(Dummy) * Dummy * Sab * OutputFactor;

         // Make I^m (00|00) = (pi/EtaABCD)^{3/2} Sab Scd G(m);
         // This is eq. (7) in PCPP 8 3072 (3073).
         m_pKernel->EvalGm(pGm, rho, rho * Dot(PmQ,PmQ), TotalL, Prefactor);

         // Make (a0|00)^m for a = 0...la+lb and m = BfC.l (inclusive, cartesian)
         // Note:
         //    - if we could absorb the rho into pGm, it may be possible to contract C
         //      already at the level of the kernel function. Check that some time.
         //    - we only need a = la...la+lb. But here this makes not enough difference
         //      to care.
         OsrrA[BfA.l+BfB.l]( p_A00, pGm + BfC.l, PmA.m, PmQ.m, rho, OvAB.InvEta );

         double
            *p_A0C_sh;
         if ( BfC.l != 0 ) {
            // move (a0|00)^m (m=lc) -> (a0|c)^(0)
            OsrrBX_3cen(p_A0C, nCartX_AB, p_A00, PmQ.m, rho, InvEtaABC, 1/ZetaC, BfA.l+BfB.l, BfA.l, BfC.l, Mem );

            // Transform c to solid harmonics. This may not be optimal in all
            // cases, since in general
            //    nCartX(2*l) >= (2*l+1)**2.
            // So here the solid transformation itself involves *more* data for c
            // than later; on the other hand the ([ab]0|->(ab| transformation becomes
            // cheaper and we don't always need all ([ab]0| components.
            p_A0C_sh = p_A0C_sh_;
            ShTrN(p_A0C_sh + iBaseCompAB, nCartX_AB, p_A0C + iBaseCompAB, nCartX_AB, BfC.l, nCartX_AB - iBaseCompAB);
         } else {
            p_A0C_sh = p_A00;
         }

         // add primitive (a0|c)-integral to all contractions in which it occurs.
         FGaussBfn::FContraction const
            &CoC = BfC.PrimContribs[iPrimC];
         double
            fCoA, fCoAB, fCoABC;
         for ( uint iCoA = CoA.nBegin; iCoA != CoA.nEnd; ++ iCoA ) {
            fCoA = CoA.Coeffs[iCoA - CoA.nBegin];
            for ( uint iCoB = CoB.nBegin; iCoB != CoB.nEnd; ++ iCoB ) {
               fCoAB = fCoA * CoB.Coeffs[iCoB - CoB.nBegin];
               for ( uint iCoC = CoC.nBegin; iCoC != CoC.nEnd; ++ iCoC ) {
                  fCoABC = fCoAB * CoC.Coeffs[iCoC - CoC.nBegin];
                  Add( c_A0C_sh + nA0C_sh * (iCoC + nCoC * (iCoB + nCoB * iCoA)),
                     p_A0C_sh, fCoABC, nA0C_sh );
               }
            }
         }
      } // C primitives
   } // AB primitives

   // wipe output data. We access each function only once, but due to the
   // important segmented case, OsrrFnC has a += characteristic. I.e., it
   // adds to the target shells.
   WipeBlock3(pOut_, nCoA*nShA, StrideA, nCoB*nShB, StrideB, nCoC*nShC, StrideC);

   // move (a0|c) (a=la..la+lb) -> (ab|c) for contracted functions
   uint
      StrideA_sh = nShA * StrideA,
      StrideB_sh = nShB * StrideB,
      StrideC_sh = nShC * StrideC;
   FOsrrFnC
      OsrrFnC = OsrrC[BfA.l][BfB.l];
   for ( uint iCoA = 0; iCoA < nCoA; ++ iCoA )
      for ( uint iCoB = 0; iCoB < nCoB; ++ iCoB )
         for ( uint iCoC = 0; iCoC < nCoC; ++ iCoC ) {
            double
               *pOut = &pOut_[iCoA * StrideA_sh + iCoB * StrideB_sh + iCoC * StrideC_sh],
               *p_A0C_sh = c_A0C_sh + nA0C_sh * (iCoC + nCoC * (iCoB + nCoB * iCoA));

            for ( uint iCompC = 0; iCompC < nShC; ++ iCompC )
               OsrrFnC( &pOut[iCompC * StrideC],
                  &p_A0C_sh[iCompC * nCartX_AB], BmA.m, StrideA, StrideB );
         }

   Mem.Free(pBaseOfMem);
};

} // namespace aic

#else

#include "AicFB.h"

namespace aic {

void FIntegralFactory::Init( int MakeScreeningParams )
{
   if ( MakeScreeningParams == 1 )
      SetScreeningParams(FScreeningParams(1e-15, FScreeningParams::AC_None));

   FShEvalData EvDat = {
      0, 0, 0, 0,
      &*m_pKernel, 0, 0, 0, 0, &m_ScreeningParams, 1.0 };
   m_EvalData = EvDat;
}

FShellData::FShellData(FGaussBfn const &Bf, FVector3 const &vCenter_)
   : l(static_cast<uint>(Bf.AngularMomentum)),
     nCo(Bf.Contractions.size()), nPrim(Bf.Exponents.size()), ShOff(static_cast<uint>(-1)),
     pExp(&Bf.Exponents[0]), pCo(&Bf.CoMatrix[0]), vCenter(vCenter_)
{
}

void FIntegralFactory::EvalInt2e3c( double *pOut, uint Strides[3], double OutputFactor,
   FGaussShell const &ShellA, FGaussShell const &ShellB, FGaussShell const &ShellC, FMemoryStack &Mem )
{
   FShellData ShA(*ShellA.pFn, ShellA.vCenter);
   FShellData ShB(*ShellB.pFn, ShellB.vCenter);
   FShellData ShC(*ShellC.pFn, ShellC.vCenter);

   m_EvalData.StrideA = Strides[0];
   m_EvalData.StrideB = Strides[1];
   m_EvalData.StrideC = Strides[2];
   m_EvalData.OutputFactor = OutputFactor;

   AicEvalShell2e3c( pOut, 0, ShA, ShB, ShC, m_EvalData, Mem );
}

void FIntegralFactory::EvalGrad2e3c( double *pOutA, double *pOutB, double *pOutC, uint Strides[4], double OutputFactor,
   FGaussShell const &ShellA, FGaussShell const &ShellB, FGaussShell const &ShellC, FMemoryStack &Mem )
{
   FShellData ShA(*ShellA.pFn, ShellA.vCenter);
   FShellData ShB(*ShellB.pFn, ShellB.vCenter);
   FShellData ShC(*ShellC.pFn, ShellC.vCenter);

   m_EvalData.StrideA = Strides[0];
   m_EvalData.StrideB = Strides[1];
   m_EvalData.StrideC = Strides[2];
   m_EvalData.StrideDeriv = Strides[3];
   m_EvalData.OutputFactor = OutputFactor;

   AicAccuContractShell3cDer1( pOutA, pOutB, pOutC,
      0, 0, 0, ShA, ShB, ShC, m_EvalData, Mem );
}

void FIntegralFactory::EvalGrad2e2c( double *pOutA, double *pOutC, uint Strides[3], double OutputFactor,
   FGaussShell const &ShellA, FGaussShell const &ShellC, FMemoryStack &Mem )
{
   FShellData ShA(*ShellA.pFn, ShellA.vCenter);
   double z0 = 0.0, z1 = 1.0;
   FShellData ShB(0, 1, 1, 0, &z0, &z1, FVector3(0,0,0));
   FShellData ShC(*ShellC.pFn, ShellC.vCenter);

   m_EvalData.StrideA = Strides[0];
   m_EvalData.StrideB = 1;
   m_EvalData.StrideC = Strides[1];
   m_EvalData.StrideDeriv = Strides[2];
   m_EvalData.OutputFactor = OutputFactor;

   AicAccuContractShell3cDer1( pOutA, 0, pOutC,
      0, 0, 0, ShA, ShB, ShC, m_EvalData, Mem );
}


} // namespace aic

#endif // RUN_THROUGH_FORTRAN_INTERFACE_CORE_ROUTINES

namespace aic {

void ShTrN_Indirect(double *AIC_RP pOut, unsigned ns, double const *AIC_RP pIn, unsigned short const *AIC_RP ii, unsigned L, unsigned Count);


// factorize [R] into [r] x S(la) with Am(R)==L, Am(r)==L-la
static void FactorizeAngularComps1(double *rl, double const *R, uint NumRSets, uint TotalL, uint la)
{
   uint
      nr1 = nCartY(TotalL - la), // number of [r]
      nrR = nCartY(TotalL);
   if (la == 0) {
      memcpy(rl, R, sizeof(*R) * nrR*NumRSets);
      return;
   }
   unsigned short const
      *ii = iCartPxx[TotalL-la][la];
   for ( uint i = 0; i < NumRSets; ++ i ) {
      aic::ShTrN_Indirect(rl, nr1, R, ii, la, nr1);
      rl += (2*la+1)*nr1;
      R  += nrR;
   }
}

// output: contracted kernels Fm(rho,T), format: (TotalL+1) x nCoA x nCoC
void Int2e2c_EvalCoKernels(double *pCoFmT, uint TotalL,
    FShellData const &ShA, FShellData const &ShC,
    double PrefactorExt, FIntegralKernel const *pKernel, FMemoryStack &Mem)
{
   double
      t = DistanceSq(ShA.vCenter, ShC.vCenter),
      *pFmT;
   Mem.Alloc(pFmT, TotalL + 1); // FmT for current primitive.

   // loop over primitives (that's all the per primitive stuff there is)
   for ( uint iExpC = 0; iExpC < ShC.nPrim; ++ iExpC )
   for ( uint iExpA = 0; iExpA < ShA.nPrim; ++ iExpA )
   {
      double
         Alpha = ShA.pExp[iExpA],
         Gamma = ShC.pExp[iExpC],
         InvEta = 1./(Alpha + Gamma),
         Rho = (Alpha * Gamma)*InvEta, // = (Alpha * Gamma)*/(Alpha + Gamma)
         Prefactor = (M_PI*InvEta)*std::sqrt(M_PI*InvEta); // = (M_PI/(Alpha+Gamma))^{3/2}

      Prefactor *= PrefactorExt;
      if(ShC.l) Prefactor *= std::pow( 1.0/(2*Gamma), (int)ShC.l); // <- use Hermites with D Ax := [1/(2 alpha)] \partial/\partial A_i.
      if(ShA.l) Prefactor *= std::pow(-1.0/(2*Alpha), (int)ShA.l); // <- -1 because \partial_A R \propto -\partial_B R!

      // calculate derivatives (D/Dt)^m exp(-rho t) with t = (A-C)^2.
      pKernel->EvalGm(pFmT, Rho, Rho*t, TotalL, Prefactor);

      // convert from Gm(rho,T) to Fm(rho,T) by absorbing powers of rho
      // (those would normally be present in the R of the MDRR)
      double
         RhoPow = 1.;
      for ( uint i = 0; i < TotalL + 1; ++ i ){
         pFmT[i] *= RhoPow;
         RhoPow *= 2*Rho;
      }

      // contract (lamely). However, normally either nCo
      // or nExp, or TotalL (or even all of them at the same time)
      // will be small, so I guess it's okay.
      for ( uint iCoC = 0; iCoC < ShC.nCo; ++ iCoC )
      for ( uint iCoA = 0; iCoA < ShA.nCo; ++ iCoA ) {
         double CoAC = ShC.pCo[iExpC + ShC.nPrim*iCoC] *
                       ShA.pCo[iExpA + ShA.nPrim*iCoA];
         Add2(&pCoFmT[(TotalL+1)*(iCoA + ShA.nCo*iCoC)],
               pFmT, CoAC, (TotalL+1));
      }
   }

   Mem.Free(pFmT);
}

// evaluate 2-electron 2-center integrals <a|krn|c>.
// if add is given: increment the output instead of overwriting it.
void EvalInt2e2c( double *pOut, uint *Strides,
    FShellData &ShA, FShellData &ShC, double Prefactor, bool Add,
    FIntegralKernel const *pKernel, FMemoryStack &Mem )
{
   FVector3
      R = ShA.vCenter - ShC.vCenter;
   uint
      lc = ShC.l, la = ShA.l,
      nShA = 2*la+1, nShC = 2*lc+1,
      iStA = Strides[0], iStC = Strides[1],
      TotalL = la + lc;
   double
      *pCoFmT, *pDataR, *pR1, *pFinal;
   Mem.ClearAlloc(pCoFmT, (TotalL+1)*ShA.nCo*ShC.nCo);
   Int2e2c_EvalCoKernels(pCoFmT, TotalL, ShA, ShC, Prefactor, pKernel, Mem);

   Mem.Alloc(pDataR, nCartY(TotalL));
   Mem.Alloc(pR1, nCartY(TotalL-lc)*(2*lc+1));
   Mem.Alloc(pFinal, (2*la+1)*(2*lc+1));

   for ( uint iCoC = 0; iCoC < ShC.nCo; ++ iCoC )
   for ( uint iCoA = 0; iCoA < ShA.nCo; ++ iCoA ) {
      double
         *pFmT = &pCoFmT[(TotalL+1)*(iCoA + ShA.nCo*iCoC)];
      aic::ShellMdrr[TotalL]( &pDataR[0], &pFmT[0], &R[0] );
      FactorizeAngularComps1(&pR1[0], &pDataR[0], 1,  TotalL, lc);
      FactorizeAngularComps1(pFinal, &pR1[0], 2*lc+1,   TotalL-lc, la);
      double
         *pOutAB = &pOut[iCoA*nShA*iStA + iCoC*nShC*iStC];
      aic::Scatter2e2c[Add](pOutAB, iStA,iStC, pFinal, la,lc);
   }

   Mem.Free(pCoFmT);
};

// evaluate 2-electron 2-center integrals <a|krn * laplace|c>
// note: to obtain the kinetic energy operator, pass an overlap kernel
//       and supply -.5 as Prefactor (ekin = -.5 laplace).
// if add is given: increment the output instead of overwriting it.
void EvalInt2e2c_LaplaceC( double *pOut, uint *Strides,
    FShellData &ShA, FShellData &ShC, double Prefactor, bool Add,
    FIntegralKernel const *pKernel, FMemoryStack &Mem )
{
   FVector3
      R = ShA.vCenter - ShC.vCenter;
   uint
      lc = ShC.l, la = ShA.l,
      nShA = 2*la+1, nShC = 2*lc+1,
      iStA = Strides[0], iStC = Strides[1],
      TotalLm2 = la + lc, TotalL = TotalLm2 + 2;
   double
      *pCoFmT, *pDataR, *pDataRm2, *pR1, *pFinal;
   Mem.ClearAlloc(pCoFmT, (TotalL+1)*ShA.nCo*ShC.nCo);
   Int2e2c_EvalCoKernels(pCoFmT, TotalL, ShA, ShC, Prefactor, pKernel, Mem);

   Mem.Alloc(pDataR, nCartY(TotalL));
   Mem.Alloc(pDataRm2, nCartY(TotalLm2));
   Mem.Alloc(pR1, nCartY(TotalLm2-lc)*(2*lc+1));
   Mem.Alloc(pFinal, (2*la+1)*(2*lc+1));

   unsigned short const
      *ic = iCartPxx[TotalLm2][2];
   for ( uint iCoC = 0; iCoC < ShC.nCo; ++ iCoC )
   for ( uint iCoA = 0; iCoA < ShA.nCo; ++ iCoA ) {
      double
         *pFmT = &pCoFmT[(TotalL+1)*(iCoA + ShA.nCo*iCoC)];
      aic::ShellMdrr[TotalL]( &pDataR[0], &pFmT[0], &R[0] );

      // form derivatives d/dxx + d/dyy + d/dzz. Due to our creative arrangement
      // of terms, that is just a matter of picking up and adding some contributions
      // from the larger shell MDRR intermediate.
      for ( uint i = 0; i < nCartY(TotalLm2); ++ i )
         pDataRm2[i] = pDataR[ic[6*i]] + pDataR[ic[6*i+1]] + pDataR[ic[6*i+2]];

      FactorizeAngularComps1(&pR1[0], &pDataRm2[0], 1,  TotalLm2, lc);
      FactorizeAngularComps1(pFinal, &pR1[0], 2*lc+1,   TotalLm2-lc, la);
      double
         *pOutAB = &pOut[iCoA*nShA*iStA + iCoC*nShC*iStC];
      aic::Scatter2e2c[Add](pOutAB, iStA,iStC, pFinal, la,lc);
   }

   Mem.Free(pCoFmT);
};

// // r[i] += f * x[i]
// template<class FScalar>
// void Add2( FScalar *AIC_RP r, FScalar const *AIC_RP x, FScalar f, std::size_t n )
// {
//    std::size_t
//       i = 0;
//    for ( ; i < (n & (~3)); i += 4 ) {
//       r[i]   += f * x[i];
//       r[i+1] += f * x[i+1];
//       r[i+2] += f * x[i+2];
//       r[i+3] += f * x[i+3];
//    }
//    for ( ; i < n; ++ i ) {
//       r[i] += f * x[i];
//    }
// };
//
// template
// void Add2<double>( double *AIC_RP r, double const *AIC_RP x, double f, std::size_t n );
//
// // r[i] += f * x[i] * y[i]
// template<class FScalar>
// void Add2( FScalar *AIC_RP r, FScalar const *AIC_RP x, FScalar const *AIC_RP y, FScalar f, std::size_t n )
// {
//    std::size_t
//       i = 0;
//    for ( ; i < (n & (~3)); i += 4 ) {
//       r[i]   += f * x[i]   * y[i];
//       r[i+1] += f * x[i+1] * y[i+1];
//       r[i+2] += f * x[i+2] * y[i+2];
//       r[i+3] += f * x[i+3] * y[i+3];
//    }
//    for ( ; i < n; ++ i ) {
//       r[i] += f * x[i] * y[i];
//    }
// };
//
// template
// void Add2<double>( double *AIC_RP r, double const *AIC_RP x, double const *AIC_RP y, double f, std::size_t n );
//
// ^- suncc doesn't get it.

// r[i] += f * x[i]
void Add2( double *AIC_RP r, double const *AIC_RP x, double f, std::size_t n )
{
   std::size_t
      i = 0;
   for ( ; i < (n & (~3)); i += 4 ) {
      r[i]   += f * x[i];
      r[i+1] += f * x[i+1];
      r[i+2] += f * x[i+2];
      r[i+3] += f * x[i+3];
   }
   for ( ; i < n; ++ i ) {
      r[i] += f * x[i];
   }
}


// r[i] += f * x[i] * y[i]
void Add2( double *AIC_RP r, double const *AIC_RP x, double const *AIC_RP y, double f, std::size_t n )
{
   std::size_t
      i = 0;
   for ( ; i < (n & (~3)); i += 4 ) {
      r[i]   += f * x[i]   * y[i];
      r[i+1] += f * x[i+1] * y[i+1];
      r[i+2] += f * x[i+2] * y[i+2];
      r[i+3] += f * x[i+3] * y[i+3];
   }
   for ( ; i < n; ++ i ) {
      r[i] += f * x[i] * y[i];
   }
}





void FIntegralFactory::EvalInt2e2c( double *pOut, uint Strides_[2], double OutputFactor,
      FGaussShell const &ShellA, FGaussShell const &ShellC, FMemoryStack &Mem )
{
   static FGaussBfnCptr
      pDummyBfn(new FGaussBfn( 0, FGaussBfn::TYPE_Unnormalized | FGaussBfn::TYPE_Spherical, 0.0 ) );
   FGaussShell
      DummyShell(-1, ShellA.vCenter, pDummyBfn);
   uint
      Strides[3] = {Strides_[0], 1, Strides_[1]};
   return EvalInt2e3c( pOut, Strides, OutputFactor, ShellA, DummyShell, ShellC, Mem );
}




} // namespace aic

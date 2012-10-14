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
#ifndef AIC_LOCAL_INCLUDES
   #define MACHINES_H_DEFINES_ONLY
   #include "util/machines.h" // for FORTINT and FORT_Extern
#else
   #include "CxFortranInt.h"
#endif
#include <stdexcept> // for std::runtime_error
#include <ostream>

#include "fmt.h"
#include "AicCommon.h"
#include "AicShells.h"

// #include "AicKernels.h" // FIXME: remove this (used for self-overlap)
// #include "AicDrv.h" // FIXME: remove this (used for self-overlap)
// #include "CtMatrix.h" // FIXME: REMOVE THIS
// #include "CtCommon.h" // FIXME: REMOVE THIS

namespace aic {

// double factorial n!! = n*(n-2)*(n-4)*...*1 (or *2 of n even)
inline uint DoubleFactT( int n )
{
   uint Result = 1;
   while( n > 1 ){
      Result *= n;
      n -= 2;
   }
   return Result;
}

static uint const
   MAX_DoubleFact = 24,
   DoubleFactTable[2 + MAX_DoubleFact]
#define D(n) DoubleFactT(n)
      = { D(-1), D(0), D(1), D(2), D(3), D(4), D(5), D(6), D(7), D(8), D(9), D(10),
         D(11), D(12), D(13), D(14), D(15), D(16), D(17), D(18), D(19), D(20), D(21),
         D(22), D(23), D(24) }; // first value: (-1)!! = 1
#undef D

inline uint DoubleFact( int n )
{
   assert( ( n >= -1 ) && ( n <= (signed)MAX_DoubleFact ) );
   return DoubleFactTable[1+n];
}

// return x^n for low n.
template<class FScalar>
inline FScalar intpow( FScalar x, int n ){
   FScalar r = 1;
   for ( ; n > 0; -- n )
      r *= x;
   return r;
}

// this is the factor which makes a unnormalized primitive spherical
// gauss orbital, i.e.   S_lm(r)*exp(-zeta r^2),'s self overlap one.
// note that the angular normalization factors are already absorbed in the S_lm
// (S_00 = 1). This function accounts for the radial part only. See purple
// book sec 6.6.4.
double GaussNormalizationSpher( double Zeta, int l )
{
   // purple book (6.6.14). (last 4*pi/(2*l+1) is for factor between
   // Ylm and Slm)
   using std::pow; using std::sqrt;
   return 2.0 * pow( 2.0 * Zeta, 0.75 ) / pow(M_PI,0.25) *
      sqrt( (1 << l) / static_cast<double>( DoubleFact(2*l+1) ) ) *
      intpow( sqrt( 2.0 * Zeta ), l ) / sqrt(4*M_PI/(2*l+1));
}

// ^- hm... Manby's program appears to use:
//       fac=(2d0*alpha/pi)**0.75d0 * sqrt(4d0*alpha)**l
// check that out. (routine basis_contraction/basis_contraction_pretty)

extern "C" {
   #define AIC_GET_GAUSS_NORMFACTOR FORT_Extern(aic_get_gauss_normfactor,AIC_GET_GAUSS_NORMFACTOR)
   double AIC_GET_GAUSS_NORMFACTOR(double const &Zeta, FORTINT const &l) {
      return GaussNormalizationSpher(Zeta, (int)l);
   }
}


FGaussBfn::FGaussBfn( uint AngularMomentum_, uint Flags_, FScalarArray const &Exponents_ )
{
   FContractionList
      Contractions_;
   Contractions_.resize(Exponents_.size());
   for ( uint iExp = 0; iExp < Exponents_.size(); ++ iExp ){
      Contractions_[iExp].nBegin = iExp;
      Contractions_[iExp].nEnd = iExp+1;
      Contractions_[iExp].Coeffs.push_back(1.0);
   }
   Initialize(AngularMomentum_, Flags_, Exponents_, Contractions_);
}

FGaussBfn::FGaussBfn( uint AngularMomentum_, uint Flags_, double Exponent_ )
{
   FContractionList
      Contractions_;
   FScalarArray
      Exponents_;
   Exponents_.push_back(Exponent_);
   Contractions_.resize(1);
   Contractions_[0].nBegin = 0;
   Contractions_[0].nEnd = 1;
   Contractions_[0].Coeffs.push_back(1.0);
   Initialize(AngularMomentum_, Flags_, Exponents_, Contractions_);
}

FGaussBfn::FGaussBfn( uint AngularMomentum_, uint Flags_, FScalarArray const &Exponents_,
   FScalarArray const &ContractionCoeffs_ )
{
   FContractionList
      Contractions_;
   Contractions_.resize(1);
   Contractions_[0].nBegin = 0;
   Contractions_[0].nEnd = Exponents_.size();
   Contractions_[0].Coeffs = ContractionCoeffs_;
   Initialize(AngularMomentum_, Flags_, Exponents_, Contractions_);
}

FGaussBfn::FGaussBfn( uint AngularMomentum_, uint Flags_, FScalarArray const &Exponents_,
   double const *pCoeffMatrix_, uint nContractions_ )
{
   FContractionList
      Contractions_;
   Contractions_.resize(nContractions_);
   // just create a set of full contractions. Initialize() will remove
   // zero contributions as needed.
   for ( uint iCo = 0; iCo < Contractions_.size(); ++ iCo ){
      FContraction
         &co = Contractions_[iCo];
      co.nBegin = 0;
      co.nEnd = Exponents_.size();
      co.Coeffs.assign( &pCoeffMatrix_[iCo * Exponents_.size()],
                        &pCoeffMatrix_[(iCo+1) * Exponents_.size()] );
   }
   Initialize(AngularMomentum_, Flags_, Exponents_, Contractions_);
}



void FGaussBfn::SanityCheck() const
{
   if ( Contractions.size() < 1 )
      throw std::runtime_error("Encountered basis-fn shell without contracted functions.");
   if ( Exponents.size() < 1 )
      throw std::runtime_error("Encountered basis-fn shell without primitive functions.");
   if ( (Flags & TYPE_Spherical) == 0 )
      throw std::runtime_error("AIC integral driver does not support Cartesian basis functions.");

   assert(nFuncs == (2*l+1) * Contractions.size());

   for ( uint i = 0; i < Exponents.size(); ++ i )
      if ( !(Exponents[i] >= 0) )
         // ^- in fact, == 0 is allowed in certain cases. We need that to emulate 2-center
         // integrals on 3-center code. Note that !(x >= y) also checks for #NaNs, #Inds, etc.
         throw std::runtime_error("Encountered basis functions with non-positive primitive exponent.");
   for ( uint i = 0; i < Contractions.size(); ++ i ) {
      FContraction const
         &co = Contractions[i];
      if ( !(( co.nBegin <= 1000 ) &&
             ( co.nEnd <= Exponents.size() ) &&
             ( co.nEnd - co.nBegin == co.Coeffs.size() ) ) )
         throw std::runtime_error("Encountered basis function with broken contraction pattern.");
   }

#ifdef _DEBUG
   // check ranges in PrimContribs
   assert( PrimContribs.size() == Exponents.size() );
   for ( uint i = 0; i < PrimContribs.size(); ++ i ) {
      FContraction const
         &pc = PrimContribs[i];
      if ( !(( pc.nBegin >= 0 ) &&
             ( pc.nEnd <= Contractions.size() ) &&
             ( pc.nEnd - pc.nBegin == pc.Coeffs.size() ) ) )
      {
         assert(!"FGaussBfn::PrimContribs broken.");
      }
   }

   // check consistency between Contractions and PrimContribs
   assert( PrimContribs.size() == Exponents.size() );
   for ( uint iExp = 0; iExp < Exponents.size(); ++ iExp ){
      FContraction const
         &pc = PrimContribs[iExp];
      assert( pc.nBegin < Contractions.size() && pc.nEnd <= Contractions.size() );
      for ( uint iCo = 0; iCo < Contractions.size(); ++ iCo ){
         FContraction const
            &co = Contractions[iCo];
         // Does iCo have contributions from iExp?
         if ( co.nBegin <= iExp && iExp < co.nEnd ){
            assert( pc.nBegin <= iCo );
            assert( pc.nEnd > iCo );
            assert( pc.Coeffs[iCo - pc.nBegin] == co.Coeffs[iExp - co.nBegin] );
         }

         // Does iExp contribute to iCo?
         if ( pc.nBegin <= iCo && iCo < pc.nEnd ){
            if ( co.nBegin <= iExp && iExp < co.nEnd ) {
               assert(pc.Coeffs[iCo - pc.nBegin] == co.Coeffs[iExp - co.nBegin]);
               assert(CoMatrix[iExp + Exponents.size() * iCo] == co.Coeffs[iExp - co.nBegin]);
            }
            else {
               assert(pc.Coeffs[iCo - pc.nBegin] == 0.0);
               assert(CoMatrix[iExp + Exponents.size() * iCo] == 0);
            }
         }
      }
   }
#endif // _DEBUG
}

void FGaussBfn::Initialize( uint AngularMomentum_, uint Flags_,
   FScalarArray const &Exponents_, FContractionList const &Contractions_ )
{
   using std::abs;
   using std::min;
   using std::max;

   AngularMomentum = AngularMomentum_;
   Flags = Flags_;
   Exponents = Exponents_;

   for ( uint i = 0; i < Exponents_.size(); ++ i )
      assert( (Exponents_[i] > 0 && 0 == (Flags & TYPE_Unnormalized)) ||
              (Exponents_[i] >= 0 && 0 != (Flags & TYPE_Unnormalized)) );
   for ( uint i = 0; i < Contractions_.size(); ++ i )
      assert( ( Contractions_[i].nBegin >= 0 ) &&
              ( Contractions_[i].nEnd <= Exponents_.size() ) );

   // copy over contractions, but omit zero contributions of primitives
   Contractions.resize(Contractions_.size());
   for ( uint iCo = 0; iCo < Contractions_.size(); ++ iCo ){
      FContraction
         &Out = Contractions[iCo];
      FContraction const
         &In = Contractions_[iCo];
      if ( In.nEnd - In.nBegin != In.Coeffs.size() )
         throw std::runtime_error("Encountered basis function with broken contraction pattern.");
      Out.nBegin = In.nBegin;
      Out.nEnd = In.nEnd;
      while ( Out.nBegin <= Out.nEnd && abs(In.Coeffs[Out.nBegin - In.nBegin]) < 1e-13 )
         Out.nBegin += 1;
      while ( Out.nEnd > Out.nBegin && abs(In.Coeffs[Out.nEnd-1 - In.nBegin]) < 1e-13 )
         Out.nEnd -= 1;
      if ( Out.nEnd - Out.nBegin == 0 || Out.nEnd - Out.nBegin > 1000 )
         throw std::runtime_error("Failed to understand basis function contraction pattern.");
      Out.Coeffs.assign( In.Coeffs.begin() + Out.nBegin - In.nBegin,
                         In.Coeffs.begin() + Out.nEnd - In.nBegin );
   }

   nFuncs = (2 * AngularMomentum + 1) * Contractions.size();

   if ( 0 == (Flags & TYPE_Unnormalized) ) {
      // absorb primitive normalization constants into contraction coefficients.
      for ( uint k = 0; k < Contractions.size(); ++ k ){
         FScalarArray
            &Coeffs = Contractions[k].Coeffs;
         double
            *pExponents = &this->Exponents[Contractions[k].nBegin];
         for ( uint i = 0; i < Coeffs.size(); ++ i )
            Coeffs[i] *= GaussNormalizationSpher( pExponents[i], AngularMomentum );
      }
   }

//    if ( 1 && Contractions.size() != 1 ) {
//       // calculate self-overlap and normalize function.
//       FOverlapKernel k;
//       FIntegralFactory
//          IntFac(k);
//       FMemoryStack2
//          Mem(2000000);
//       double
//          *pS;
//       Mem.ClearAlloc(pS, nFuncs*nFuncs);
//       uint Strides[2] = {1, nFuncs};
//       FGaussShell
//          Sh( 0, FVector3(0.,0.,0.), FGaussBfnCptr(this,false) );
//       MakeCoMatrix();
//       IntFac.EvalInt2e2c( pS, Strides, 1.0, Sh, Sh, Mem );
//       ct::PrintMatrixGen(ct::xout, pS, nFuncs,1,nFuncs,nFuncs, "self-overlap");
//       for ( uint iCo = 0; iCo < Contractions.size(); ++ iCo ) {
//          FContraction
//             &co = Contractions[iCo];
//          double f = std::sqrt(1./pS[(iCo*(2*l+1)) * (nFuncs+1)]);
//          for ( uint iExp = 0; iExp < co.Coeffs.size(); ++ iExp )
//             co.Coeffs[iExp] *= f;
//       }
//    }


   // convert .Contractions array into .PrimContribs form, which is
   // better suited for integral evaluation.
   PrimContribs.resize( Exponents.size() );

   for ( uint iExp = 0; iExp < Exponents.size(); ++ iExp ){
      FContraction
         &pc = PrimContribs[iExp];
      // find first and last contracted function to which this primitive contributes
      pc.nBegin = 0xffff;
      pc.nEnd = 0;
      for ( uint iCo = 0; iCo < Contractions.size(); ++ iCo ){
         FContraction
            &co = Contractions[iCo];
         if ( iExp >= co.nBegin && iExp < co.nEnd ) {
            pc.nBegin = min(pc.nBegin, iCo);
            pc.nEnd = max(pc.nEnd, iCo+1);
         }
      }
      if ( pc.nEnd == 0 )
         throw std::runtime_error("Encountered basis-fn shell with non-contributing primitives.");
      pc.Coeffs.resize(pc.nEnd - pc.nBegin);
      // copy over the coefficients
      for ( uint iCo = pc.nBegin; iCo < pc.nEnd; ++ iCo ){
         FContraction
            &co = Contractions[iCo];
         if ( iExp >= co.nBegin && iExp < co.nEnd )
            pc.Coeffs[iCo - pc.nBegin] = co.Coeffs[iExp - co.nBegin];
         else
            pc.Coeffs[iCo - pc.nBegin] = 0.0; // unlikely but possible.
      }
   }

   MakeCoMatrix();

   SanityCheck();
}

void FGaussBfn::MakeCoMatrix()
{
   // create CoMatrix member from Contractions member.
   uint
      nExp = Exponents.size(),
      nCo = Contractions.size();
   CoMatrix.clear();
   CoMatrix.resize(nExp * nCo, 0.0);
   for ( uint iCo = 0; iCo < nCo; ++ iCo ) {
      FContraction
         &co = Contractions[iCo];
      for ( uint iExp = co.nBegin; iExp < co.nEnd; ++ iExp )
         CoMatrix[iExp + iCo * nExp] = co.Coeffs[iExp - co.nBegin];
   }

}

FGaussBfn::~FGaussBfn()
{
}

void FGaussBfn::Print( std::ostream &out ) const
{
   using fmt::ff;
   using fmt::fi;
   using fmt::fe;
   out << "AngMom: " << AngularMomentum << "   "
/*       << "vCenter: (" << ff(vCenter[0],6,4) << " " << ff(vCenter[1],7,4) << " " << ff(vCenter[2],6,4) << ")   "
       << "iCenter: " << fi(iCenter,3) << ""*/
       << "\nExponents:\n       ";
   int wd = 10;
   for ( uint i = 0; i < Exponents.size(); ++ i )
      out << " " << fe(Exponents[i],wd,4);
   out << "\n";

   out << "Contractions:\n";
   for ( uint iCo = 0; iCo < Contractions.size(); ++ iCo ) {
      out << "   " << fi(iCo,2) << ": ";
      FContraction const &co = Contractions[iCo];
      for ( uint iExp = 0; iExp < Exponents.size(); ++ iExp ) {
         double r = 0.0;
         if ( iExp >= co.nBegin && iExp < co.nEnd )
            r = co.Coeffs[iExp - co.nBegin];
         out << " " << ff(r,wd,6);
      }
      out << "    [" << co.nBegin << "--" << co.nEnd - 1 << "]";
      out << "\n";
   }
   out << "\n";
#ifdef _DEBUG
   out << "Primitive contributions:\n";
   for ( uint iExp = 0; iExp < PrimContribs.size(); ++ iExp ) {
      out << "   " << fi(iExp,2) << ": ";
      FContraction const &pc = PrimContribs[iExp];
      for ( uint iCo = 0; iCo < Contractions.size(); ++ iCo ) {
         double r = 0.0;
         if ( iCo >= pc.nBegin && iCo < pc.nEnd )
            r = pc.Coeffs[iCo - pc.nBegin];
         out << " " << ff(r,wd,6);
      }
      out << "    [" << pc.nBegin << "--" << pc.nEnd - 1 << "]";
      out << "\n";
   }
//    out << "\n";
#endif
}


// A generally contracted gaussian orbital shell
FGaussShell::FGaussShell( FAtomIndex iCenter_, FVector3 vCenter_, FGaussBfnCptr pFn_ )
   : vCenter(vCenter_), iCenter(iCenter_), pFn(pFn_)
{
}

FGaussShell::~FGaussShell()
{
}

void FGaussShell::Print( std::ostream &out ) const
{
   using fmt::ff;
   using fmt::fi;
   using fmt::fe;
   out << "vCenter: (" << ff(vCenter[0],6,4) << " " << ff(vCenter[1],7,4) << " " << ff(vCenter[2],6,4) << ")   "
       << "iCenter: " << fi(iCenter,3) << "\n"
       << *pFn;
}

void FGaussShell::PrintAligned( std::ostream &xout, uint Indent ) const
{
   using namespace fmt;
   std::streampos
      p0 = xout.tellp(),
      p1;
//    xout << fi(iCenter, 3) << " "
   xout << fi(pFn->Contractions.size(), 3) << ":"
        << "spdfghiklm"[AngularMomentum()]
        << "   "
        << ff(vCenter[0],8,4) << " " << ff(vCenter[1],8,4) << " " << ff(vCenter[2],8,4)
        << "    ";
   p1 = xout.tellp();

   for ( uint iExp = 0; iExp < pFn->Exponents.size(); ++ iExp ){
      if ( iExp != 0 ){
         xout << "\n";
         for ( uint i = 0; i < Indent + p1 - p0; ++ i )
            xout << " ";
      }
      xout << fmt::ff(pFn->Exponents[iExp],16,7) << "  ";

      double
         fRenorm = 1.0/GaussNormalizationSpher( pFn->Exponents[iExp], pFn->AngularMomentum );
//       fRenorm = 1.;
      std::stringstream
         str;
      for ( uint iCo = 0; iCo < pFn->Contractions.size(); ++ iCo ){
         aic::FGaussBfn::FContraction const
            &Co = pFn->Contractions[iCo];
         uint
            w = 9, p = 5;
         if ( Co.nBegin <= iExp && iExp < Co.nEnd ) {
            str << " " << fmt::ff(Co.Coeffs[iExp - Co.nBegin]*fRenorm, w, p);
         } else {
            str << " " << fmt::fc("  - - - -", w);
         }
      }
      std::string
         s = str.str();
      while( !s.empty() && (s[s.size()-1] == ' ' || s[s.size()-1] == '-' ) )
         s.resize(s.size() - 1);
      xout << s;
   }
}


} // namespace aic

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

#ifndef _AIC_SHELLS_H
#define _AIC_SHELLS_H

#include <vector>
#include <iosfwd>
#include <boost/intrusive_ptr.hpp>

#include "AicCommon.h"

namespace aic {

typedef int
   FAtomIndex;
typedef std::vector<double>
   FScalarArray;

// FGaussBfn and FGaussShell: Note on intended usage
// -------------------------------------------------
// 1) FGaussBfn contains the information (exponents, contractions, AngMom)
//    of a single gaussian orbital basis function shell. It *does not* contain
//    the center coordinates. For example, one could have a FGaussBfn
//    representing 'cc-pVTZ p-shell of atom O' (without specifying where the O
//    atom is, or if there is one at all).
// 2) FGaussShell contains a center coordinate and a link to a FGaussBfn.
//    For example in the O_2 molecule, one could have two different FGaussShell
//    objects linking to the same FGaussBfn, representing the p-shells of O_2,
//    one on each O atom.
//
// The idea is that FGaussBfn objects are kept and managed by basis library
// objects, and then the FGaussShell objects are created dynamically by linking
// to them. FGaussShells are light-weight, FGaussBfns are mid-weight
// (instantiating them is not free; they may need to calcualte the
//  self-overlap, for example).
//
// Nevertheless, FGaussBfn are linked by intrusive pointers in FGaussShell;
// therefore you can also create these objects with a one-to-one relationship
// (i.e., without sharing the FGaussBfns across multiple shells). In that case
// the bfn objects would be owned by shell objects, and not by the basis set
// library object.


// A not-yet-centered generally contracted gaussian orbital basis function
// (i.e., contains information about the functions, but not the center coordinates)
// One of these may be linked to multiple FGaussShell objects.
struct FGaussBfn : public FIntrusivePtrDest
{
   union {
      int l; // 0:s, 1:p, 2:d..
      int AngularMomentum;
   };
   enum FFlags{
      TYPE_Cartesian = 0x0000,
      TYPE_Spherical = 0x0001,
      TYPE_Unnormalized = 0x0004
   };
   uint
      Flags;
   FScalarArray
      Exponents; // primitive exponents
   FScalarArray
      CoMatrix; // nExp x nCo matrix of contraction coefficients.

   struct FContraction{
      uint
         nBegin, nEnd;
      FScalarArray
         Coeffs;

      FContraction() {};
      FContraction(uint nBegin_, uint nEnd_, FScalarArray const &Coeffs_)
         : nBegin(nBegin_), nEnd(nEnd_), Coeffs(Coeffs_)
      {}
   };
   typedef std::vector<FContraction>
      FContractionList;
   FContractionList
      // [i]: Contracted function #i is composed of exponents [nBegin,nEnd) with
      // contraction coeffients .Coeffs[0..nEnd-nBegin].
      Contractions,
      // transpose of above; [i]: exponent #i contributes to contracted functions
      // [nBegin,nEnd) with prefactor Coeff.
      PrimContribs;

   // auxiliary shortcut quantities that can be calculated from the data above
   uint
      nFuncs;  // number of functions in the shell (s=1, p=3,
               // d=6 or 5, f=10 or 7,..)*nContractions

   // throws std::runtime_error if something goes wrong.
   void SanityCheck() const;

   // make a set of generally contracted functions with explicit
   // contraction specification
   FGaussBfn( uint AngularMomentum_, uint Flags_, FScalarArray const &Exponents_,
      FContractionList const &Contractions_ )
   {
      Initialize( AngularMomentum_, Flags_, Exponents_, Contractions_ );
   }

   // make a set of generally contracted functions based on contraction matrix.
   // Coefficient of function iCo for exponent iExp is
   //    pCoeffMatrix_[iExp + iCo * Exponents.size() ].
   FGaussBfn( uint AngularMomentum_, uint Flags_, FScalarArray const &Exponents_,
      double const *pCoeffMatrix_, uint nContractions_ );

   // make a single (segmented) contracted function
   FGaussBfn( uint AngularMomentum_, uint Flags_, FScalarArray const &Exponents_, FScalarArray const &ContractionCoeffs_ );

   // make a set of primitive functions
   FGaussBfn( uint AngularMomentum_, uint Flags_, FScalarArray const &Exponents_ );

   // make a single primitive function
   FGaussBfn( uint AngularMomentum_, uint Flags_, double Exponent_ );

   virtual ~FGaussBfn();

   void Print( std::ostream &out ) const;
private:
   // calculates auxiliary data and normalizes the contractions. parameters:
   //  - Exponents: primitive gaussian exponents which are used
   //    on all functions on the shell.
   //  - Coefficients: contraction coefficients in standard form, i.e. not
   //    normalized to gaussian type.
   void Initialize( uint AngularMomentum_, uint Flags_, FScalarArray const &Exponents_,
      FContractionList const &Contractions_ );

   void MakeCoMatrix();
};

typedef boost::intrusive_ptr<FGaussBfn>
   FGaussBfnPtr;
typedef boost::intrusive_ptr<FGaussBfn const>
   FGaussBfnCptr;

inline std::ostream &operator << ( std::ostream &out, FGaussBfn const &a ){
   a.Print(out);
   return out;
}

// A generally contracted gaussian orbital shell
struct FGaussShell : public FIntrusivePtrDest
{
   FVector3
      vCenter; // center position of the gaussian shell.
   FAtomIndex
      iCenter; // index of the atom this shell is centered on. -1 if on none.
   FGaussBfnCptr
      pFn; // basis function information

   FGaussShell( FAtomIndex iCenter_, FVector3 vCenter_, FGaussBfnCptr pFn_ );
   virtual ~FGaussShell();

   uint nFuncs() const { return pFn->nFuncs; };
   uint AngularMomentum() const { return pFn->AngularMomentum; }
   uint nPrimitives() const { return pFn->Exponents.size() * (2 * AngularMomentum() + 1); }
   void SanityCheck() const { assert(pFn.get() != 0); assert(pFn->m_RefCount > 0); return pFn->SanityCheck(); }

   void Print( std::ostream &out ) const;
   void PrintAligned( std::ostream &out, uint Indent = 0 ) const;
};

typedef boost::intrusive_ptr<FGaussShell>
   FGaussShellPtr;

inline std::ostream &operator << ( std::ostream &out, FGaussShell const &a ){
   a.Print(out);
   return out;
}


double GaussNormalizationSpher( double Zeta, int l );



} // namespace aic



#endif // _AIC_SHELLS_H

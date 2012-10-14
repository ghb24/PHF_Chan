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

#ifndef CT8K_BASISSET_H
#define CT8K_BASISSET_H

#include "CxPodArray.h"
#include "CtCommon.h"
#include "CtAtomSet.h"
#include "AicShells.h"
#ifdef RUN_THROUGH_FORTRAN_INTERFACE_CORE_ROUTINES
  #include "AicFB.h" // for FShellData
#endif

namespace aic {
   struct FIntegralFactory;
}

namespace ct {



// describes a single symmetry adapted basis function (SO) by giving its
// linear combination in terms of AOs.
struct FSymOrb {
   static const uint
      iNegWeight = 0x80000000,
      iAoMask = iNegWeight - 1;
   uint
      // nEqiv[iSo]:  nti = number of equivalent AOs forming SO iSo.
      nEqivAo,
      nEqivAoCase, // 0, 1, 2 or 3. nEqiv = 2^nEqivAoCase
      // iEqiv[0..nti-1, iSo]:  indices of the equivalent AOs. If >= iNegWeight,
      //   actual index is iEqiv & iAoMask but weight in stabilizer is negative.
      iEqivAo[8];
   double
      // Weight[iSo]:       1/sqrt(nti) -> weight in stabilizer (equal for all AOs).
      Weight;
};

// struct FSymGroup {
//    uint nEqivCase;
//    char SoPhase[8];
//    uint iAo[8];
//    uint iSo[8];
// };

struct FUniqueShell { // a non-redundant basis function shell
   uint iSh;
   double fWeight; // number of times it occurs symmetry-equivalently.
};


// describes the transformation between a raw basis set and a symmetry
// adapted basis set
struct FSymTrans {
   uint
      nIrreps,
      // offset of first SO transforming according to that irrep.
      // Length: nIrrep+1. Last entry contains total number.
      IrrepOffsets[9],
      // number of symmetry adapted basis functions (== number of basis functions)
      nSo,
      nUniqueShells;
   FSymOrb
      // for each SO, gives the AO linear combination forming it.
      *pSo;
   FUniqueShell
      *pUniqueShells;
};


struct FBasisSet : FIntrusivePtrDest
{
   typedef std::vector<aic::FGaussShell>
      FGaussShellArray;
   FGaussShellArray
      Shells;  // Gaussian shells centered on different atoms.
   FBasisContext
      Context; // for informative purposes only.

   void Print( std::ostream &out ) const;
   uint nPrimitiveGtos() const;
   uint nBasisFn() const;

   std::string Name;

//  FBasisSet();
   FBasisSet( FAtomSet const &AtomSet, FBasisContext Context_);
#ifdef RUN_THROUGH_FORTRAN_INTERFACE_CORE_ROUTINES
public:
   // what follows now are interfacing variants of the basis data
   // (as in this->Shells) in a format compatible with the low-level and
   // high-level driver routines (FB*) shared with Molpro.
   typedef std::vector<aic::FShellData>
      // linearized version of the basis as used in FB driver routines
      // shared with Molpro.
      FDriverShellArray;
   FDriverShellArray
      DriverShells;

   // Note: In MolproBasis, expt and cc are both set to zero.
   // The actual pointers are provided via FGrpInf(!!). That means,
   // e.g., that &cc[infg[ig][ig_ccpt]] points into the CoMatrix object
   // of the corresponding FGaussBfn object.
   FMolproBasis
      MolproBasis;
   typedef TArray<aic::FGrpInf>
      FGrpInfArray;
   FGrpInfArray
      GrpInf;
   FORTINT
      // total number of basis functions
      nFn;
   TArray<double [3]>
      Centers; // acts as rr. &Centers[0] can be passed into pCenters arrays.
#endif // RUN_THROUGH_FORTRAN_INTERFACE_CORE_ROUTINES

   // rebuilds the data for low-level interfaces (i.e., the data in this section)
   // from this->Shells
   void RebuildInterfacingData();
private:
   // makes *this a basis set corresponding to the basis descriptions in
   // AtomSet. Attempts to load the neccessary external data. Throws
   // std::runtime_error if failed.
   void LoadFromAtomSet( FAtomSet const &AtomSet, FBasisContext Context );

   void MakeAtomOffsets( uint *&pAtomShellOffsets, uint *&pAtomBfnOffsets, uint nAtoms_, FMemoryStack &Mem ) const;
   void MakeShellOffsets( uint *&pShellOffsets, FMemoryStack &Mem ) const;
   void Finalize();
public:
   void MakeSymmetryTransformation( FSymTrans &st, FAtomSet const &AtomSet, FMemoryStack &Mem );
   void PrintSoBasis( std::ostream &xout, FAtomSet const &AtomSet, FMemoryStack &Mem );
};

inline std::ostream &operator << ( std::ostream &out, FBasisSet const &BasisSet ) {
   BasisSet.Print(out);
   return out;
}

typedef boost::intrusive_ptr<FBasisSet>
   FBasisSetPtr;

struct FMatrixView;

void MakeIntMatrix( FMatrixView &Out, FBasisSet const &RowBasis,
   FBasisSet const &ColBasis, aic::FIntegralFactory &IntFactory, FMemoryStack &Mem );


enum FSymTransDir {
   ST_AoToSo = 0,
   ST_SoToAo = 1
};

/// transform, in place, the row dimension of a matrix between raw basis
/// functions (AOs) and symmetry adapted basis functions (SOs)
/// \p pData   Matrix; Input & Output
/// \p Dir     Determines if AO->SO or SO->AO transform is done
/// \p nRowSt  Stride between two consecutive rows of the matrix
/// \p nColSt  Stride between two consecutive columns of the matrix
/// \p nCols   Number of columns
/// \p st      Determines the symmetry space for the rows (by this also nRows!).
/// \p Mem     Space for temporaries
void SymTrans1( double *pData, uint nRowSt, uint nColSt, uint nCols, FSymTrans const &st, FSymTransDir Dir, FMemoryStack &Mem );
void SymTrans2( double *pData, uint nRowSt, uint nColSt, FSymTrans const &RowSy, FSymTrans const &ColSy, FSymTransDir Dir, FMemoryStack &Mem );
void SymTrans2( FMatrixView &M, FSymTrans const &RowSy, FSymTrans const &ColSy, FSymTransDir Dir, FMemoryStack &Mem );



} // namespace ct

#endif // CT8K_BASISSET_H

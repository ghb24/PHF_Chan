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

#include <stdexcept>
#include <boost/format.hpp>
using boost::format;

#include "AicDrv.h"
#include "AicShells.h"
#include "AicSolidHarmonics.h"

#include "CtCommon.h"
#include "CtBasisSet.h"
#include "CtBasisLibrary.h"
#include "CtMatrix.h"

namespace ct {


// FBasisSet::FBasisSet()
// {
// };

FBasisSet::FBasisSet( FAtomSet const &AtomSet, FBasisContext Context_ )
   : Context(Context_)
{
   LoadFromAtomSet(AtomSet, Context);
}


std::string BasisContextName( FBasisContext Context )
{
   switch( Context ) {
      case BASIS_Orbital: return "ORBITAL";
      case BASIS_JFit: return "JFIT";
      case BASIS_JkFit: return "JKFIT";
      case BASIS_Mp2Fit: return "MP2FIT";
      case BASIS_CcsdFit: return "EXTFIT";
      case BASIS_F12RI: return "OPTRI";

      default:
         assert(0);
         return "[unknown basis context]";
   };
}


void FindDefaultBasis( std::string &Out, FBasisContext Context, std::string const &In )
{
   assert( Context != BASIS_Orbital );
   // TODO: do special element/context/basis-type specific stuff here.

   Out = In + "-" + BasisContextName(Context);
}



void FBasisSet::LoadFromAtomSet( FAtomSet const &AtomSet, FBasisContext Context )
{
   std::string
      BasisName;
   bool
      AllEqual = true;

   BasisName.reserve(64);
   for ( uint iAt = 0; iAt < AtomSet.Atoms.size(); ++ iAt ){
      FAtom const
         &Atom = AtomSet.Atoms[iAt];

      // get name of basis we are supposed to be loading.
      FBasisDescs::const_iterator
         itDesc = Atom.BasisDesc.find(Context);
      if ( itDesc != Atom.BasisDesc.end() )
         // basis for context explicitly set
         BasisName = itDesc->second;
      else {
         // basis for current context not set---go find one.
         if ( Context == BASIS_Guess ) {
            // should actually check for ECP-aware guess basis. But I guess
            // this will be okay for H-Kr.
            BasisName = "cc-pVTZ";
         } else {
            itDesc = Atom.BasisDesc.find(BASIS_Orbital);
            if ( itDesc == Atom.BasisDesc.end() ) {
               std::stringstream str;
               str << "No basis set assigned to atom " << iAt << " at " << Atom.vPos << ".";
               throw std::runtime_error(str.str());
            }
            FindDefaultBasis( BasisName, Context, itDesc->second );
         }
      }

      if ( this->Name.empty() ) this->Name = BasisName;
      if ( this->Name != BasisName ) AllEqual = false;

      // we're importing molpro data... Molpro truncates basis names
      // at 32 characters. Do so here, too.
      BasisName.resize( std::min(BasisName.size(), 32ul) );

      // also, convert everything to lowercase. We do that when importing sets.
      for ( uint i = 0; i < BasisName.size(); ++ i )
         BasisName[i] = ::tolower(BasisName[i]);

      // ask basis library to add the functions.
      g_BasisSetLibrary.LoadBasisFunctions( this->Shells,
         Atom.AtomicNumber, BasisName,
         Atom.vPos, iAt );
   };

   if ( !AllEqual || this->Name.empty() )
      this->Name = BasisContextName(Context);

   // TODO:
   //   - prevent re-allocation of this->Shells during filling
   //   - make local copy of FGaussBfn objects and re-link this->Shells
   //     to these objects. This would reduce the memory spread of the data
   //     and prevent TLB misses during integration.

   Finalize();
}


void FBasisSet::Print( std::ostream &xout ) const
{
   xout << "Basis set '" << BasisContextName(Context) << "'\n" << std::endl;
/*   for ( uint iSh = 0; iSh < Shells.size(); ++ iSh )
      xout << Shells[iSh] << std::endl;*/
   xout << "  Offs   NC/AM        Center/Position        Exponents           Contractions\n";
   xout << " -----------------------------------------   -----------------   ----------------------" << std::endl;
   uint
      nOff = 0;
   for ( uint iSh = 0; iSh < Shells.size(); ++ iSh ) {
      std::streampos p0 = xout.tellp(), p1;
      xout << format("%5i   ") % nOff;
      p1 = xout.tellp();
      Shells[iSh].PrintAligned(xout, p1-p0);
      xout << std::endl;
      nOff += Shells[iSh].nFuncs();
   }

   xout << std::endl;
};


uint FBasisSet::nPrimitiveGtos() const
{
   uint r = 0;
   for ( uint iSh = 0; iSh < Shells.size(); ++ iSh )
      r += Shells[iSh].nPrimitives();
   return r;
};

uint FBasisSet::nBasisFn() const
{
   uint r = 0;
   for ( uint iSh = 0; iSh < Shells.size(); ++ iSh )
      r += Shells[iSh].nFuncs();
   return r;
};


void MakeIntMatrix( FMatrixView &Dest, FBasisSet const &RowBasis,
   FBasisSet const &ColBasis, aic::FIntegralFactory &IntFactory, FMemoryStack &Mem )
{
   bool
      MatrixSymmetric = (&RowBasis == &ColBasis);

   assert( Dest.nRows == RowBasis.nBasisFn() );
   assert( Dest.nCols == ColBasis.nBasisFn() );

   double
      *pIntResult;
   Mem.Alloc(pIntResult, Dest.nRows * Dest.nCols); // that should suffice.
   uint
      StartFnA = 0, StartFnB; // starting indices of the functions of the
                              // respective shells.
   for ( uint iShellA = 0; iShellA < RowBasis.Shells.size(); ++ iShellA ){
      FGaussShell const
         *pShellA = &RowBasis.Shells[iShellA];
      StartFnB = (!MatrixSymmetric)? 0 : StartFnA;
      for ( uint iShellB = (!MatrixSymmetric)? 0 : iShellA;
           iShellB < ColBasis.Shells.size(); ++ iShellB )
      {
         FGaussShell const
            *pShellB = &ColBasis.Shells[iShellB];
         uint
            nSizeA = pShellA->nFuncs(),
            nSizeB = pShellB->nFuncs();
         uint
            Strides[2] = {1, nSizeA};
         IntFactory.EvalInt2e2c( pIntResult, Strides, 1.0, *pShellA, *pShellB, Mem );

         assert( StartFnA + nSizeA <= RowBasis.nBasisFn() );
         assert( StartFnB + nSizeB <= ColBasis.nBasisFn() );

         // fill data we gathered into the matrices
         for ( uint j_ = 0; j_ < nSizeB; ++ j_ )
            for ( uint i_ = 0; i_ < nSizeA; ++ i_ ) {
               uint
                  i = StartFnA + i_,
                  j = StartFnB + j_;
               FScalar const
                  &r = pIntResult[ i_ + nSizeA * j_ ];
               Dest(i,j) = r;
               if ( MatrixSymmetric )
                  Dest(j,i) = r;
            }

         StartFnB += pShellB->nFuncs();
      }
      StartFnA += pShellA->nFuncs();
   }
   assert( StartFnA == RowBasis.nBasisFn() );
   Mem.Free(pIntResult);
}


void FBasisSet::MakeAtomOffsets( uint *&pAtomShellOffsets, uint *&pAtomBfnOffsets, uint nAtoms_, FMemoryStack &Mem ) const
{
   // pAtomShellOffsets: maps atom id to first basis function shell
   //    (Shells[i]) of the atom
   // pAtomBfnOffsets: maps atom id to first basis function (AO index)
   //    of the atom
   uint
      nAtoms = 0;
   for ( uint iSh = 0; iSh < Shells.size(); ++ iSh )
      if ( Shells[iSh].iCenter > 0 )
         nAtoms = std::max(nAtoms, 1 + static_cast<unsigned>(Shells[iSh].iCenter));
   assert( nAtoms <= nAtoms_ );
   // ^- might differ if there are some atoms without basis functions
   nAtoms = nAtoms_;

   Mem.Alloc(pAtomShellOffsets, nAtoms+1);
   Mem.Alloc(pAtomBfnOffsets, nAtoms+1);
   *pAtomShellOffsets = 0;
   *pAtomBfnOffsets = 0;
   // note: this code assumes that this->shells is ordered according
   // to the atom set!
   for ( uint iAt = 0; iAt < nAtoms; ++ iAt ){
      uint
         iSh = pAtomShellOffsets[iAt],
         iBf = pAtomBfnOffsets[iAt];
      while ( iSh < Shells.size() && Shells[iSh].iCenter == (signed)iAt ) {
         iBf += Shells[iSh].nFuncs();
         ++ iSh;
      }
//       _xout0("iAt: " << iAt << "   iSh: "<< iSh << "   iBf: " << iBf);
      pAtomShellOffsets[iAt+1] = iSh;
      pAtomBfnOffsets[iAt+1] = iBf;
      assert(iBf <= nBasisFn());
   };
}

void FBasisSet::MakeShellOffsets( uint *&pShellOffsets, FMemoryStack &Mem ) const
{
   Mem.Alloc(pShellOffsets, Shells.size() + 1);
   *pShellOffsets = 0;
   for ( uint iSh = 0; iSh < Shells.size(); ++ iSh )
      pShellOffsets[iSh+1] = pShellOffsets[iSh] + Shells[iSh].nFuncs();
};




void FBasisSet::RebuildInterfacingData()
{
#ifdef RUN_THROUGH_FORTRAN_INTERFACE_CORE_ROUTINES
   // copy over basis function positions. WARNING: This assumes that the
   // iCenter entries in this->Shells are all set correctly!
   uint
      nCenters = 0;
   for ( uint iShell = 0; iShell < Shells.size(); ++ iShell ) {
      FGaussShell
         &Shell = Shells[iShell];
      assert(Shell.iCenter >= 0);
      if ( static_cast<FAtomIndex>(nCenters) <= Shell.iCenter )
         nCenters = Shell.iCenter + 1;
   }

   Centers.resize(nCenters);
   for ( uint iShell = 0; iShell < Shells.size(); ++ iShell ) {
      // note that this normally copies over multiple times to the
      // same center index (because multiple shells refer to the same
      // center; namely, all with different AM on the same atom)
      // There is an assert() somewhere below that these are consistent.
      assert(Shells[iShell].iCenter >= 0 && (uint)Shells[iShell].iCenter < nCenters);
      for ( uint im = 0; im < 3; ++ im )
         Centers[Shells[iShell].iCenter][im] = Shells[iShell].vCenter.m[im];
   }

   GrpInf.resize(Shells.size());
   memset(&GrpInf[0], 0, sizeof(GrpInf[0]) * GrpInf.size());
   double
      *expt = 0, *cc = 0;
   uint
      iBfOffset = 0;
   for ( uint iShell = 0; iShell < Shells.size(); ++ iShell ) {
      FGaussShell
         &Shell = Shells[iShell];
      FGrpInf
         &infg = GrpInf[iShell];
      infg[ig_nprm] = Shell.pFn->Exponents.size();
      infg[ig_ncnt] = Shell.pFn->Contractions.size();
      infg[ig_func] = Shell.pFn->nFuncs;
      infg[ig_minl] = Shell.pFn->AngularMomentum;
      infg[ig_maxl] = infg[ig_minl];
//       infg[ig_shof] = 0xbadc0de;
//       // ^- that one is used to refer to screening densities supplied in terms
//       //    of primitives. We don't use those here.
      infg[ig_shof] = static_cast<FORTINT>(&Shell.pFn->Exponents[0] - expt);
      // ^- 0-based.
      infg[ig_ccpt] = static_cast<FORTINT>(&Shell.pFn->CoMatrix[0] - cc + 1);
      // ^- 1-based.
      infg[ig_cent] = Shell.iCenter + 1;
      // ^- 1-based.

      infg[ig_coff] = iBfOffset;
      // ^- offset of contracted function relative to start of basis. 0 based?

      infg[ig_type] = 0; // GTO. We have only GTO here (1 would be Poisson-GTO).

      iBfOffset += infg[ig_func];
      assert(aic::FVector3(Centers[Shells[iShell].iCenter]) == Shell.vCenter);
   }

   MolproBasis = FMolproBasis(&GrpInf[0], expt, cc, &Centers[0]);

   // and finally, make a set of FShellData objects for talking to the
   // low-level drivers directly.
   DriverShells.clear();
   DriverShells.reserve(GrpInf.size());
   for ( uint iGrp = 0; iGrp < GrpInf.size(); ++ iGrp ) {
      DriverShells.push_back(FShellData(iGrp, this->MolproBasis));
   }
#endif // RUN_THROUGH_FORTRAN_INTERFACE_CORE_ROUTINES
};


void FBasisSet::Finalize()
{
   RebuildInterfacingData();
};




} // namespace ct

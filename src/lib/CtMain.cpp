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
#include <iostream>
#include <string.h> // for memcpy

#include <fstream> // for importing stuff (cp2k data)
#include <stdexcept>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#define INDENTSTREAM_IMPL
#include "CxIndentStream.h"

#include "CtCommon.h"
#include "CtIo.h"
#include "CtAtomSet.h"
#include "CtBasisLibrary.h"
#include "CtBasisSet.h"
#include "CtInt1e.h"
#include "CtMatrix.h"
#include "AicDrv.h" // for making integrals

#include "CtTiming.h" // for embdedded HF timer.
#include "CtConstants.h"

#include "CxNumpyArray.h"

using boost::format;
using boost::str;

namespace ct {

static uint *MakeShellOffsets( FBasisSet::FGaussShellArray const &Shells, FMemoryStack &Mem  )
{
   uint
      *pShOff;
   Mem.Alloc(pShOff, Shells.size()+1);
   *pShOff = 0;
   for ( uint i = 0; i < Shells.size(); ++ i )
      pShOff[i+1] = pShOff[i] + Shells[i].nFuncs();
   Mem.Align(16);
   return pShOff;
}


void FormIntIJA(double *pIntFai, FMatrixView OrbA, FMatrixView OrbI, FMatrixView Jcd, FBasisSet const *pOrbBasis, FBasisSet const *pFitBasis, aic::FIntegralFactory IntFactory, FMemoryStack &Mem)
{
   void
      *pBaseOfMemory = Mem.Alloc(0);
   uint
      nAo = pOrbBasis->nBasisFn(),
      nFit = pFitBasis->nBasisFn(),
      nOrbA = OrbA.nCols,
      nOcc = OrbI.nCols;
   if ( nAo != OrbI.nRows || nAo != OrbA.nRows )
      throw std::runtime_error("MakeFittingIntegrals: Input orbitals not consistent with orbital basis set.");
   double
      *pDF_NFi;
   Mem.Alloc(pDF_NFi, nAo * nFit * nOcc);

   FBasisSet::FGaussShellArray const
      &ShellsF = pFitBasis->Shells,
      &ShellsA = pOrbBasis->Shells;
   uint
      *pShOffA = MakeShellOffsets(ShellsA, Mem),
      *pShOffF = MakeShellOffsets(ShellsF, Mem);

   for ( uint iShF = 0; iShF != ShellsF.size(); ++ iShF ){
      FGaussShell const &ShF = ShellsF[iShF];
      uint nFnF = ShF.nFuncs();

      double
         *pMNF;
      Mem.Alloc(pMNF, nAo * nAo * nFnF );

      for ( uint iShB = 0; iShB < ShellsA.size(); ++ iShB ){
         FGaussShell const &ShB = ShellsA[iShB];
         uint nFnB = ShB.nFuncs();
         for ( uint iShA = iShB; iShA < ShellsA.size(); ++ iShA ) {
            FGaussShell const &ShA = ShellsA[iShA];
            uint
               nFnA = ShA.nFuncs(),
               Strides[3] = {1, nFnA, nFnA * nFnB};
            double
               *pIntData;
            Mem.Alloc(pIntData, nFnA * nFnB * nFnF );

            IntFactory.EvalInt2e3c( pIntData, Strides, 1.0, ShA, ShB, ShF, Mem );

            for ( uint iF = 0; iF < nFnF; ++ iF )
               for ( uint iB = 0; iB < nFnB; ++ iB )
                  for ( uint iA = 0; iA < nFnA; ++ iA ) {
                     double
                        f = pIntData[iA + nFnA * (iB + nFnB * iF)];
                     // assign to (\mu\nu| and (\nu\mu|. (int has perm symmetry).
                     pMNF[(pShOffA[iShA] + iA) + nAo * (pShOffA[iShB] + iB) + nAo * nAo * iF] = f;
                     pMNF[(pShOffA[iShB] + iB) + nAo * (pShOffA[iShA] + iA) + nAo * nAo * iF] = f;
                  }

            Mem.Free(pIntData);
         }
      }

      // D[\nu A i] = (i\nu|A) = C(\mu,i) (\mu [\nu|A)]
      FMatrixView
         DF_NFI(pDF_NFi + nAo * pShOffF[iShF], nAo * nFnF, nOcc, 1, nFit * nAo);
      Mxm( DF_NFI,
           FMatrixView(pMNF, nAo * nFnF, nAo, nAo, 1),
           OrbI );

      Mem.Free(pMNF);
   }

   // Jcd version:
   //   Solve[Jcd[AB] D[\nu B i]] -> D[\nu A I]
   for ( uint i = 0; i < nOcc; ++ i ){
      FMatrixView
         D_NF(pDF_NFi + nAo * nFit * i, nAo, nFit),
         K_Fa(pIntFai + nFit * nOrbA * i, nFit, nOrbA);
//       TriangularSolve(Transpose(D_NF), Jcd);
//       Mxm(K_Fa, Transpose(D_NF), OrbA);
// ^- hm... I think I thought this was a good idea for some reason,
//          but I cannot remember why.
      Mxm(K_Fa, Transpose(D_NF), OrbA);
      TriangularSolve(K_Fa, Jcd);
   }

   Mem.Free(pDF_NFi);
   Mem.Free(pBaseOfMemory);
}




void MakeFittingIntegrals(double *pIntFai, FMatrixView OrbA, FMatrixView OrbI, FBasisSet const *pOrbBasis, FBasisSet const *pFitBasis, FMemoryStack &Mem)
{
   uint
      nAo = pOrbBasis->nBasisFn(),
      nFit = pFitBasis->nBasisFn();

   xout << boost::format("ORBITAL BASIS %s:%35t%5i FUNCTIONS") % pOrbBasis->Name % nAo << "\n"
        << boost::format("FITTING BASIS %s:%35t%5i FUNCTIONS") % pFitBasis->Name % nFit << "\n"
        << std::endl;

   aic::FCoulombKernel
      CoulombKernel;
   aic::FIntegralFactory
      IntFactory(&CoulombKernel);

   // Make fitting coefficients J^{-1/2}.
   FTimer TimerJmh;
   FStackMatrix
      Jmh(nFit, nFit, &Mem);
   MakeIntMatrix( Jmh, *pFitBasis, *pFitBasis, IntFactory, Mem );
   CalcCholeskyFactors(Jmh);
   xout << format(pTimingFmt) % "J^{-1/2} of fitting basis" % (double)TimerJmh; xout.flush();

   // make (F|ai) the lame way: assume we can keep (ai|A) in memory.
   FTimer TimerIJA;
   FormIntIJA(pIntFai, OrbA, OrbI, Jmh, pOrbBasis, pFitBasis, IntFactory, Mem);
   xout << format(pTimingFmt) % "3-index integrals" % (double)TimerIJA; xout.flush();
};



void WriteMatrixToFile_Rect(std::string const &FileName, std::string const Desc, double *pData, uint nRows, uint nCols, char const *pNumFmt)
{
    if ( pNumFmt == 0 )
        pNumFmt = "%21.14f";
    std::ofstream
        File(FileName.c_str());
    File << format("! %s, %i x %i\n") % Desc % nRows % nCols;
    if ( !File.good() )
        throw std::runtime_error("failed to open file '" + FileName + "' for writing.");
    for ( uint iRow = 0; iRow < nRows; ++ iRow ) {
        for ( uint iCol = 0; iCol < nCols; ++ iCol ) {
            File << " " << format(pNumFmt) % pData[iRow + nRows*iCol];
        };
        File << "\n";
    }
}


void WriteMatrixtoFile2(std::string const &FileName, ct::FMatrixView const &M, std::string const &MatrixName, std::string const &MatrixFormat, int Verbosity = 0)
{
   if ( MatrixFormat == "npy" ) {
      if ( M.nRowSt != 1 || M.nColSt != M.nRows )
         throw std::runtime_error("Sorry, WriteNpy can currently only export continuous matrices in default layout.");
      ct::WriteNpy(FileName, M.pData, ct::MakeShape(M.nRows, M.nCols));
   }
   else if ( MatrixFormat == "list" )
      WriteMatrixToFile(FileName, M, MatrixName);
   else if ( MatrixFormat == "rect" )
      WriteMatrixToFile_Rect(FileName, MatrixName, M.pData, M.nRows, M.nCols, " %22.15e");
   else
      throw std::runtime_error("output matrix format not recognized: '" + MatrixFormat + "'");
   if ( Verbosity >= 1 )
      ct::xout << format("* wrote %s matrix to '%s'.") % MatrixName % FileName << std::endl;
};

// Set Out := L^T In R.
void BasisChange2( FMatrixView &Out, FMatrixView const &L,
    FMatrixView const &In, FMatrixView const &R, FMemoryStack &Mem)
{
    assert( L.nRows == In.nRows && R.nRows == In.nCols );
//     Out = FMatrixView( 0, L.nCols, R.nCols );
//     Mem.Alloc(Out.pData, Out.GetStridedSize());
    assert( Out.nRows == L.nCols && Out.nCols == R.nCols );

    if ( L.nRows * L.nCols <=  R.nRows * R.nCols ) {
        // T1 := L^T In
        FStackMatrix
            T1(L.nCols, In.nCols, &Mem);
        Mxm(T1, Transpose(L), In);
        Mxm(Out, T1, R);
    } else {
        // T1 := In * R
        FStackMatrix
            T1(In.nRows, R.nCols, &Mem);
        Mxm(T1, In, R);
        Mxm(Out, Transpose(L), T1);
    }

    if ( L.pData == R.pData && In.IsSymmetric(1e-10) )
       Symmetrize(Out);
}

} // namespace ct


int main_integral_export(int argc, char *argv[])
{
   using namespace ct;
   std::cout <<
        "BASIS FUNCTIONS & INTEGRALS [v20121011]"
      "\n  -- A program for experiments in molecular electronic structure theory"
      "\n                                                (c) Gerald Knizia, 2012"
      << std::endl;

   po::options_description options_desc("Options");
   std::vector< std::string >
      FileNames_BasisLibs;
   std::string
      BasisName_Orbital,
      BasisName_Fit,
      FileName_Overlap,
      FileName_CoreH,
      FileName_2eFit,
      FileName_Atoms,
      FileName_AtomsAu,
      FileName_AtomsAng,
      Name_OrbitalTrafo,
      MatrixFormat;
   options_desc.add_options()
      ("help",
          "print this help message")
      ("basis-lib", po::value<std::vector<std::string> >(&FileNames_BasisLibs)->composing(),
          "file names of .libmol basis set libraries to load (Molpro basis lib format). Can occur multiple times. "
          "all basis sets used in the input must be defined in one of those files.")
      ("save-overlap", po::value<std::string>(&FileName_Overlap)->default_value(""),
          "if given, save the AO overlap matrix S <mu|nu> to this file.")
      ("save-coreh", po::value<std::string>(&FileName_CoreH)->default_value(""),
          "if given, save the core Hamiltonian (1e terms) to this file.")
      ("save-fint2e", po::value<std::string>(&FileName_2eFit)->default_value(""),
          "if given, save the full (F|mu nu) 2e fitting integral matrix to this file. "
          "Fitting metric is already absorbed in F; i.e., the 2e integral "
          "(mu nu|eta xi) \\approx \\sum_F (mu nu|F)(F|eta xi).")
      ("basis-orb", po::value<std::string>(&BasisName_Orbital)->default_value("MINAO"),
          "name of the orbital basis to use for r,s. (e.g., cc-pVTZ, def2-SVP)")
      ("basis-fit", po::value<std::string>(&BasisName_Fit)->default_value("def2-TZVPP-JKFIT"),
          "name of the fitting basis to use for F (e.g., def2-TZVPP-JKFIT, def2-TZVPP-MP2FI)")
      ("orb-trafo", po::value<std::string>(&Name_OrbitalTrafo)->default_value("None"),
          "can be 'None' or 'Smh'. In the latter case, all output integrals are stored in terms "
          "of symmetrically orthogonalized basis functions for mu,nu instead of raw basis "
          "functions.")
      ("matrix-format", po::value<std::string>(&MatrixFormat)->default_value("list"),
          "Can be 'npy', 'list' or 'rect'. Defines file format of output matrices. "
          "npy is a binary format that can be loaded by numpy.load(..).")
      ("atoms-au", po::value<std::string>(&FileName_AtomsAu)->default_value(""),
          "File name of input atoms (xyz format). Input coordinates are given in "
          "atomic units (1.0 == 1.0 bohr radius)")
      ("atoms-ang", po::value<std::string>(&FileName_AtomsAng)->default_value(""),
          "File name of input atoms (xyz format). Input coordinates are given in "
          "angstroms (1.0 == 10^-{10} meter)")
   ;

   po::variables_map vm;
   po::store(po::parse_command_line(argc, argv, options_desc), vm);
   po::notify(vm);

   if (vm.count("help") || argc == 1) {
      xout << "This program calculates integrals over Gaussian basis functions. See README.";
      xout << options_desc << std::endl;
      return 1;
   }

   // find file name of input atomic coordinates.
   double
      fInputAtomScale = 1.0;
   if ( FileName_AtomsAu != "" ) {
      FileName_Atoms = FileName_AtomsAu;
      fInputAtomScale = ToAng;
   } else if ( FileName_AtomsAng != "" ) {
      FileName_Atoms = FileName_AtomsAng;
      fInputAtomScale = 1.0;
   } else {
      throw std::runtime_error("either --atoms-au or --atoms-ang arguments must be given.");
   }

   // check combinations of other input options.
   if ( FileNames_BasisLibs.empty() )
      throw std::runtime_error("you must specify at least one basis set library file to load (--basis-lib missing)");
   if ( BasisName_Fit.empty() && !FileName_2eFit.empty() )
      throw std::runtime_error("if calculating 2e integrals, you must specify a fitting basis (--basis-fit)");
   bool
      TransformToSmh = false;
   if ( Name_OrbitalTrafo != "None" and Name_OrbitalTrafo != "Smh" )
      throw std::runtime_error("orbital transformation not recognized (--orb-trafo)");
   if ( Name_OrbitalTrafo == "Smh" )
      TransformToSmh = true;
   if ( MatrixFormat != "npy" && MatrixFormat != "list" && MatrixFormat != "rect" )
      throw std::runtime_error("Output matrix format must be 'npy', 'list' or 'rect' (--matrix-format)");

   xout << fmt::ind();

   MajorProgramIntro(xout, "BASIS LIBRARY");
   for ( uint iBasisLib = 0; iBasisLib < FileNames_BasisLibs.size(); ++ iBasisLib )
      g_BasisSetLibrary.ImportMolproLib(FileNames_BasisLibs[iBasisLib]);

   boost::intrusive_ptr<FAtomSet>
      pAtoms = new FAtomSet;
   FBasisDescs
      DefaultBases;
   DefaultBases[BASIS_Orbital] = BasisName_Orbital;
   DefaultBases[BASIS_JkFit] = BasisName_Fit;

   pAtoms->AddAtomsFromXyzFile(FileName_Atoms, DefaultBases);
   if ( fInputAtomScale != 1.0 ) {
      for ( uint iAtom = 0; iAtom < pAtoms->size(); ++ iAtom )
         (*pAtoms)[iAtom].vPos *= fInputAtomScale;
   }

   MajorProgramIntro(xout, "GEOMETRY & ENVIRONMENT");
   xout << "ATOMIC CONFIGURATION\n\n";
   xout << *pAtoms << std::endl;

   MajorProgramIntro(xout, "EXPORT OF INTEGRALS");
   aic::InitBoysFnTable(); // don't ask.

   // instanciate the basis sets.
   FBasisSetPtr
      pOrbBasis = new FBasisSet(*pAtoms, BASIS_Orbital),
      pFitBasis;
   if ( FileName_2eFit != "" )
      pFitBasis = new FBasisSet(*pAtoms, BASIS_JkFit);
   uint
      nAo = pOrbBasis->nBasisFn(),
      nFit = 0;
   if ( pFitBasis )
      nFit = pFitBasis->nBasisFn();

   std::size_t
      RequiredMem = 4*nAo*nAo + (500<<20);
   if ( FileName_2eFit != "" )
      RequiredMem += (2*nFit*nAo*nAo) + nFit*nFit;
   ct::FMemoryStack2
      Mem(RequiredMem);

   // we'll just do the transformation always, and push a indentity matrix
   // if raw output integrals are requested. Can simply use MakeFittingIntegrals
   // then.
   FStackMatrix
      S(nAo, nAo, &Mem),
      Orbs(nAo, nAo, &Mem);
   pAtoms->MakeOverlapMatrix(S, *pOrbBasis, *pOrbBasis, Mem);
   if ( !TransformToSmh ) {
      Orbs.SetIdentity();
   } else {
      Move(Orbs, S);
      CalcSmhMatrix(Orbs, Mem, FSmhOptions(1e-15,1e-15,("S^{-1/2} for orbital basis ("+BasisName_Orbital+")").c_str(), &ct::xout));
      S.SetIdentity();
   }
   if ( FileName_Overlap != "" )
      WriteMatrixtoFile2(FileName_Overlap, S, "Overlap", MatrixFormat, 1);

   if ( FileName_CoreH != "" ) {
      FStackMatrix
         CoreH_Ao(nAo, nAo, &Mem),
         CoreH_Mo(nAo, nAo, &Mem);
      pAtoms->MakeCoreHamiltonMatrix(CoreH_Ao, *pOrbBasis, *pOrbBasis, Mem);
      BasisChange2(CoreH_Mo, Orbs, CoreH_Ao, Orbs, Mem);
      WriteMatrixtoFile2(FileName_CoreH, CoreH_Mo, "CoreH", MatrixFormat, 1);
   }

   if ( FileName_2eFit != "" ) {
      FStackMatrix
         IntFai(nFit, nAo*nAo, &Mem);
      MakeFittingIntegrals(IntFai.pData, Orbs, Orbs, &*pOrbBasis, &*pFitBasis, Mem);
      WriteMatrixtoFile2(FileName_2eFit, IntFai, "2e-(F|rs)", MatrixFormat, 1);
   }

   xout << "\n";
   xout << format(pResultFmt) % "Nuclear repulsion energy" % pAtoms->NuclearRepulsionEnergy();

   xout.flush();
   return 0;
};


// int main(int argc, char *argv[])
// {
// //    return mainWheHeehee();
//    return main_integral_export(argc, argv);
// };

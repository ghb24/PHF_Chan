#include <memory>
#include <stdexcept>
#include <boost/format.hpp>
using boost::format;

#include "PhfBasisSet.h"

#include "lib/AicShells.h"
#include "lib/AicKernels.h"
#include "lib/AicGridOps.h"
#include "lib/CtInt1e.h" // FIXME: remove this (or better: fix it)
#include "lib/CtIo.h"

FBasisSet::FBasisSet()
{
   // super-cell translations (==periodicity of the raw basis
   // functions) go here, later.
   Data.resize(9);
   Data.clear_data();
};


// construct a basis by cloning another basis, and translating it by
// the given list of displacements
FBasisSet::FBasisSet(FBasisSet const &Source, FVector3 const *pTra, size_t nTra)
{
   // we can just copy over the actual basis data. All we need to
   // to is:
   //     (1) make new groups for the translated basis functions,
   //         by adjusting their center indices
   //     (2) add new centers for the translated functions.
   this->Data = Source.Data;
   this->Groups.reserve(nTra * Source.Groups.size());
   this->Centers.reserve(nTra * Source.Centers.size());
   for ( uint iTra = 0; iTra < nTra; ++ iTra ) {
      this->Groups.insert(this->Groups.end(), Source.Groups.begin(), Source.Groups.end());
      for ( uint iGrp = 0; iGrp < Source.Groups.size(); ++ iGrp ) {
         uint iGrpTra = iGrp + iTra*Source.Groups.size();
         this->Groups[iGrpTra].iCen += this->Centers.size();
      }

      for ( uint iCen = 0; iCen < Source.Centers.size(); ++ iCen )
         this->Centers.push_back(Source.Centers[iCen] + pTra[iTra]);
      if ( 0 ) {
         for ( uint iGrp = 0; iGrp < Source.Groups.size(); ++ iGrp ) {
            uint iGrpTra = iGrp + iTra*Source.Groups.size();
            FVector3 &v = this->Centers[this->Groups[iGrpTra].iCen];
            ct::xout << format("iTra = %i  iGrp = %i  iCen -> %i  vCen -> (%.3f %.3f %.3f)")
               % iTra % iGrpTra % Groups[iGrpTra].iCen % v[0] % v[1] % v[2] << std::endl;
         }
      }
   };
   Finalize();
};


void FBasisSet::SetPeriodicityVectors(FVector3 T[3])
{
   assert(Data.size() >= 9);
   for ( uint it = 0; it != 3; ++ it )
      for ( uint ixyz = 0; ixyz != 3; ++ ixyz )
         Data[ixyz + 3*it] = T[it][ixyz];
};



void FBasisSet::AddAicShell(aic::FGaussShell const &Shell)
{
   FGroupInfo
      gi;
   gi.nFn = Shell.pFn->nFuncs;
   gi.nExp = Shell.pFn->Exponents.size();
   gi.iExp = Data.size();
   Data.insert(Data.end(), &*Shell.pFn->Exponents.begin(), &*Shell.pFn->Exponents.end());
   gi.nCo = Shell.pFn->Contractions.size();
   gi.iCo = Data.size();
   Data.insert(Data.end(), Shell.pFn->CoMatrix.begin(), Shell.pFn->CoMatrix.end());
   gi.l = Shell.pFn->AngularMomentum;
   gi.Type = 0;  // function type (0=GTO, 1=Poisson)
   gi.iCen = Centers.size();
   assert((signed)Centers.size() == Shell.iCenter);

   // TODO: calculate effective range of the basis function.
   gi.iRange = Data.size();
   Data.push_back(1e20);
   Groups.push_back(gi);
   Centers.push_back(Shell.vCenter);
};


void FBasisSet::Finalize()
{
   // count number of basis functions.
   nFn = 0;
   for ( uint i = 0; i < Groups.size(); ++ i ) {
      assert(Groups[i].nFn == Groups[i].nCo * (2*Groups[i].l + 1));
      nFn += Groups[i].nFn;
   }
};





struct FAicIntegralContext
{
   // space for intermediates.
   ct::FMemoryStack2
      Mem;
   aic::FKernelPtr
      pKernel;
   // threshold for integral evaluation--what exactly that means
   // depends on the context. LogThr is -log(Thr/10.).
   double
      Thr, LogThr;
   uint
      iKernel;
   bool
      FreeAfterUse;
};


static FAicIntegralContext *GetContext(FORTINT iContext){
   return reinterpret_cast<FAicIntegralContext*>(iContext);
}


extern "C" {

FORTINT FD(create_integral_context)(void *pWork, FORTINT const &nWork_, double const &Thr)
{
   std::ptrdiff_t
      nWork = nWork_;
   // well... anyone not convinced this is a fantastic idea?
   if ( nWork <= 0 ) {
      nWork = 20000000/sizeof(double); // ~20mb
      pWork = static_cast<double*>(::malloc(sizeof(double)*nWork));
   }

   FAicIntegralContext
      *pContext;
   pContext = new((void*)pWork) FAicIntegralContext;
   pContext->FreeAfterUse = (nWork != nWork_);
   pContext->Mem.AssignMemory(reinterpret_cast<char*>(pWork) + sizeof(FAicIntegralContext),
      sizeof(double)*nWork);
   pContext->Mem.Align(16);
   pContext->Thr = std::min(Thr, (double)0.);
   if ( pContext->Thr > 1e-42 )
      pContext->LogThr = -log(Thr/10.);
   else
      pContext->LogThr = 99.;
   pContext->iKernel = INTKERNEL_None;
   if ( sizeof(FORTINT) < sizeof(FAicIntegralContext*) )
      throw std::runtime_error("sorry, our hacky C++/Fortran interface assumes that Fortran integers can hold C pointers. Please use 8-byte integers on 64bit systems!");
   return reinterpret_cast<FORTINT>(pContext);
};

void FD(destroy_integral_context)(FORTINT &iContext)
{
   FAicIntegralContext
      *ic = GetContext(iContext);
   bool Free = ic->FreeAfterUse;
   ic->~FAicIntegralContext();
   if (Free)
      ::free(ic);
};

void FD(assign_integral_kernel)(FORTINT &iContext, FORTINT const &iKernel, FORTINT *pParamsI, double *pParamsF)
{
   FAicIntegralContext
      *ic = GetContext(iContext);
   ic->iKernel = iKernel;
   switch ( iKernel ) {
      case INTKERNEL_Overlap:
         ic->pKernel = 0;
         break;
      case INTKERNEL_Kinetic:
         ic->pKernel = 0;
         break;
      case INTKERNEL_Coulomb:
         ic->pKernel = new aic::FCoulombKernel();
         break;
      // note: the screened kernels are not yet implemented. It is, however,
      //       simple to do so by c/p'ing FCoulombKernel and modifying it.
      //       The formulas for Gm(T,rho) are noted at the end
      //       of http://dx.doi.org/10.1039/b605188j
//       INTKERNEL_Coulomb_LongRange_Erf = 4,  // pParamsF[0]: Screening length
//       INTKERNEL_Coulomb_ShortRange_Erfc = 5 // pParamsF[0]: Screening length
      default:
         throw std::runtime_error("integral kernel not recognized.");
   }
}

void FD(eval_basis_int1e)(double *pOut, FORTINT *Strides, double const &Factor,
   FBasisSet const &BasisA, FBasisSet const &BasisB, FORTINT &iContext)
{
   size_t
      iFnB = 0;
   for ( size_t iGrpB = 0; iGrpB < BasisB.size(); ++ iGrpB )
   {
      size_t
         iFnA = 0;
      for ( size_t iGrpA = 0; iGrpA < BasisA.size(); ++ iGrpA )
      {
         FD(eval_group_int1e)(&pOut[iFnA * Strides[0] + iFnB * Strides[1]],
            Strides, Factor, iGrpA, BasisA, iGrpB, BasisB, iContext);
         iFnA += BasisA[iGrpA].nFn;
      };
      iFnB += BasisB[iGrpB].nFn;
   }
};

void FD(eval_group_int1e)(double *pOut, FORTINT *Strides, double const &Factor,
   FORTINT const &iGrpA, FBasisSet const &BasisA,
   FORTINT const &iGrpB, FBasisSet const &BasisB, FORTINT &iContext)
{
   FAicIntegralContext
      *ic = GetContext(iContext);
   using namespace aic;
   using namespace ct;
   // FIXME: - do a proper implementation, using the xyz-factorization in the
   //          reciprocal lattice directions.
   //        - this actually needs to be fast! (not only ``not totally
   //          embarrassing as it is now'')
   //        - super-cell range and screening range needs to be
   //          determined from GroupInfo.Range (*requires GroupInfo.Range to actually be set)
   assert(ic->iKernel == INTKERNEL_Kinetic || ic->iKernel == INTKERNEL_Overlap);
   FGroupInfo const &giA = BasisA[iGrpA];
   FGaussShell
      ShA(giA.iCen, BasisA.Centers[giA.iCen],
         new FGaussBfn(giA.l, FGaussBfn::TYPE_Unnormalized | FGaussBfn::TYPE_Spherical,
            FScalarArray(&BasisA.Data[giA.iExp],&BasisA.Data[giA.iExp+giA.nExp]),
            &BasisA.Data[giA.iCo], giA.nCo ));
   FGroupInfo const &giB = BasisB[iGrpB];
   FGaussShell
      ShB(giB.iCen, BasisB.Centers[giB.iCen],
         new FGaussBfn(giB.l, FGaussBfn::TYPE_Unnormalized | FGaussBfn::TYPE_Spherical,
            FScalarArray(&BasisB.Data[giB.iExp],&BasisB.Data[giB.iExp+giB.nExp]),
            &BasisB.Data[giB.iCo], giB.nCo ));

   FDoublettIntegralFactory
      *pFactory = 0;
   FDoublettIntegralFactoryOverlap      OverlapFactory;
   FDoublettIntegralFactoryKineticTerm  KineticFactory;

   switch ( ic->iKernel ) {
      case INTKERNEL_Kinetic:
         pFactory = &KineticFactory;
         break;
      case INTKERNEL_Overlap:
         pFactory = &OverlapFactory;
         break;
      default:
         assert(0);
   }

   FShellDoublettIntegral
      r;
   for ( int iFnA = 0; iFnA < giA.nFn; ++ iFnA )
      for ( int iFnB = 0; iFnB < giB.nFn; ++ iFnB )
         pOut[Strides[0]*iFnA + Strides[1]*iFnB] = 0;

   FVector3 const
      *pSuperCellT = reinterpret_cast<FVector3 const*>(&BasisB.Data[0]);

//    int const N = 10;
   int const N = 5;
   for ( int iT0 = -N; iT0 <= N; ++ iT0 ) {
      for ( int iT1 = -N; iT1 <= N; ++ iT1 ) {
         for ( int iT2 = -N; iT2 <= N; ++ iT2 ) {
            ShB.vCenter = BasisB.Centers[giB.iCen] +
               (double)iT0 * pSuperCellT[0] +
               (double)iT1 * pSuperCellT[1] +
               (double)iT2 * pSuperCellT[2];
            pFactory->EvalDoublett(r, &ShA, &ShB, ic->Mem);
            assert(r.nSizeA == giA.nFn && r.nSizeB == giB.nFn);
            for ( int iFnB = 0; iFnB < giB.nFn; ++ iFnB )
               for ( int iFnA = 0; iFnA < giA.nFn; ++ iFnA )
                  pOut[Strides[0]*iFnA + Strides[1]*iFnB] += Factor * r.Int(iFnA,iFnB);
         }
      }
   }
   // ^- just seeing that already hurts. I can't believe I actually wrote it...
};

static aic::FShellData MakeAicShellData(FBasisSet const &Basis, FORTINT iGrp)
{
   FGroupInfo const &gi = Basis[iGrp];
   return aic::FShellData(gi.l, gi.nCo, gi.nExp, -1, &Basis.Data[gi.iExp], &Basis.Data[gi.iCo], Basis.Centers[gi.iCen]);
};



// c/p'd from AIC_EVAL_BFN_ON_GRID
void FD(eval_basis_fn_on_grid)(double *pOut, FORTINT const &nCompSt,
   FORTINT *pCentersOut, FORTINT *pMap, FORTINT &nMap,
   FBasisSet const &Basis, double  (*pGridPt)[3], FORTINT const &nGridPt,
   FORTINT const &DerivOrder, FORTINT &iContext)
{
   // FIXME:
   //   - this evaluates *raw* gaussians only, not periodically symmetrized ones.
   //   - not sure how to best handle screening in that case... probably
   //     need to modify EvalShellGroupOnGrid.
   FAicIntegralContext
      *ic = GetContext(iContext);
   void *pBase = ic->Mem.Alloc(0);

   nMap = 0;
   size_t
      iFn = 0; // orbital offset of current group, in terms of basis functions
   for ( size_t iGrp = 0; iGrp < Basis.size(); ) {
      // find last group sitting on the same center as iGrp.
      FORTINT
         iCenter = Basis[iGrp].iCen;
      size_t
         iGrpEnd = iGrp + 1;
      while ( iGrpEnd < Basis.size() && Basis[iGrpEnd].iCen == iCenter )
         ++ iGrpEnd;

      aic::FShellData
         *pBfs;
      ic->Mem.Alloc(pBfs, iGrpEnd - iGrp);
      for ( uint iBf = 0; iBf < iGrpEnd - iGrp; ++ iBf )
         pBfs[iBf] = MakeAicShellData(Basis, iGrp+iBf);

      FORTINT
         nMap0 = nMap;

      RESUME_CLOCK(1)
      EvalShellGroupOnGrid( pOut + nGridPt * nCompSt * nMap, nGridPt, nGridPt * nCompSt,
         pMap, nMap, 1+iFn,
         &pGridPt[0][0], 3, nGridPt,
         pBfs, pBfs + iGrpEnd - iGrp, pBfs[0].vCenter,
         DerivOrder, ic->Thr, ic->LogThr, ic->Mem );
      PAUSE_CLOCK(1)

      // note: exact post-screening code deleted here to make things simpler.

      for ( FORTINT iMap = nMap0; iMap < nMap; ++ iMap )
         pCentersOut[iMap] = iCenter;

      iGrp = iGrpEnd;
      iFn += Basis[iGrp].nFn;
      ic->Mem.Free(pBfs);
   }

   ic->Mem.Free(pBase);
   PAUSE_CLOCK(0)
}


//    double const
//       // just a random large number. Use this for setting up basis sets of multipoles.
//       fPointMultipoleExponent = 1.84391e20;


void FD(eval_basis_int2e_contract_point_charges)(double *pOut, FORTINT *Strides, double const &Factor,
   FBasisSet const &BasisA, FBasisSet const &BasisB,
   double const *pCoeffC, FVector3 const *pCentersC, size_t nCentersC, FORTINT &iContext)
{
   assert(0);
}








} // extern "C"



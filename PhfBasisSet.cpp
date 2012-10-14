#include <memory>
#include <stdexcept>
#include "CxMemoryStack.h"
#include "AicShells.h"
#include "AicKernels.h"
#include "PhfBasisSet.h"

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
      for ( uint iGrp = 0; iGrp < Source.Groups.size(); ++ iGrp )
         this->Groups[iGrp].iCen += this->Centers.size();

      for ( uint iCen = 0; iCen < Source.Centers.size(); ++ iCen )
         this->Centers.push_back(Source.Centers[iCen] + pTra[iTra]);
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
   gi.Type = 0;  //! function type (0=GTO, 1=Poisson)
   gi.iCen = Shell.iCenter;  //! center index

   // TODO: calculate effective range of the basis function.
   gi.iRange = Data.size();
   Data.push_back(1e20);
   Groups.push_back(gi);
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
   // depends on the context.
   double
      Thr;
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
   pContext->Thr = Thr;
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




};




} // extern "C"
































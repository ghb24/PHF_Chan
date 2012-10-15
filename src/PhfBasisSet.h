#ifndef PHF_BASIS_SET_H
#define PHF_BASIS_SET_H

#include "PhfTypes.h"

namespace aic {
   struct FGaussShell;
};


// A note on Gaussian basis sets (general):
//  - A single Gaussian basis basis function (e.g., a 2px function) cannot be
//    evaluated efficiently. Normally, there are many functions which share the
//    same evaluation intermediates (e.g., p functions on a center); these have
//    to be evaluated together
//
//  - For this reason a basis set is subdivided into /shell groups/. A shell
//    group bundles all functions using shared intermediates. Normally a shell
//    group consists of:
//
//      o Functions on a single center, with the same angular
//        momentum
//      o One set of primitive exponents used by all of them
//      o A contraction matrix defining how the primitives are
//        linearly combined to form the actual basis functions.
//
//    In rare cases (the one we're about to implement being one) it is sensible
//    to share exponents between the s- and p functions (and thus create
//    dual-angular momentum shell groups). For simplicitly reasons this
//    is not done here.
//
//  - A shell group has (2*l + 1)*nCo functions, where 'l' is the angular
//    momentum and nCo is the number of contractions.
//
// In the solid case:
//  - For solids, technically the basis functions *themselves* are periodic,
//    with the period of the super-cell size (see supercells_and_kspace.pdf).
//
//  - I(cgk)'ll try to make the integral/basis-function-on-grid routines etc
//    take that into account automatically. For this reason the basis set
//    objects carry super-cell periodicity vectors.


// note:
//  - the objects here and on the Fortran side need to be consistent!
//
//  - Within this program, all indices are 0-based!

struct FGroupInfo {
   FORTINT nFn;   //! number of functions (contractions x components)
   FORTINT nExp;  //! number of primitive exponents
   FORTINT iExp;  //! start index in exponent array, in externally supplied *pExp / *pBasisData
   FORTINT nCo;   //! number of contractions
   FORTINT iCo;   //! start index of nExp x nCo contraction matrix, in externally supplied *pCo / *pBasisData
   FORTINT l;     //! angular momentum
   FORTINT Type;  //! function type (0=GTO, 1=Poisson)
   FORTINT iCen;  //! center index
   FORTINT iRange; //! pBasisData[i]: effective range at which we can consider the basis function to be zero.
};

// /// Data on a group of shells we can evaluate in one round.
// /// Actual integral routines work with this.
// struct FShellGroup
// {
//    FORTINT
//       nExp, nCo, Am;
//    double const
//       /// pointer to nExp x nCo contraction matrix.
//       /// pCo[iExp + nExp*iCo] is the contribution of the
//       /// unnormalized primitive Gaussian iExp to contraction iCo.
//       *pCo,
//       /// Pointer to exponent array. pExp[i] is the exponent of primitive #i
//       *pExp,
//       /// pointer to center of Gaussian. pCen[0],[1],[2] are the x,y,z coordinates.
//       *pCen,
//       /// pointer to 3 x 3 matrix defining the (super-cell) lattice translations
//       /// which define the periodicity of the basis function.
//       *pTra;
// };

// Holds information about a gaussian basis set.
struct FBasisSet
{
   TArray<FGroupInfo>
      /// shell groups: defines the actual basis functions
      Groups;
   TArray<double>
      /// data on sequences of (flattened) contraction matrices
      /// and exponent arrays. Data in *pGroupInfo refers to these.
      Data;
   TArray<FVector3>
      /// centers of the raw Gaussians
      Centers;
   FORTINT
      /// total number of basis functions.
      nFn;

   FGroupInfo const& operator [] ( size_t iGroup ) const { return Groups[iGroup]; }
   size_t size() const { return Groups.size(); }

//    // combine the information for a shell group and return
//    // it as a single object, referring only to the shell group.
//    FShellGroup GetGroup( size_t iGroup ) {
//       FGroupInfo &ig = Groups[iGroup];
//       FShellGroup r = { ig.nExp, ig.nCo, ig.l, &Data[ig.iCo], &Data[ig.iExp], &Centers[ig.iCen][0], &Data[0] };
//       return r;
//    }

   FBasisSet();

   /// construct a basis by cloning another basis, and translating it by
   /// the given list of displacements
   explicit FBasisSet( FBasisSet const &Source, FVector3 const *pTra, size_t nTra);

   // convert an AIC shell object into this basis format and add it to the
   // controlled objects.
   void AddAicShell(aic::FGaussShell const &Shell);
   // store the lattice translations defining the periodicity of the
   // basis functions (== SuperCell.T).
   void SetPeriodicityVectors(FVector3 T[3]);

   void Finalize(); // set derived properties from Groups array.
};


extern "C" {
   enum FPhfIntgralKernel {
      INTKERNEL_None = 0,
      INTKERNEL_Overlap = 1,
      INTKERNEL_Kinetic = 2,
      INTKERNEL_Coulomb = 3,
      INTKERNEL_Coulomb_LongRange_Erf = 4,  // pParamsF[0]: Screening length
      INTKERNEL_Coulomb_ShortRange_Erfc = 5 // pParamsF[0]: Screening length
   };

   /// create an internal integral evaluation context object, and return a handle to it.
   /// The handle must be passed to subsequent integral routines and be freed (via DestroyIntegralContext)
   /// after use.
   ///   \p pWork     Work space -- temporary storage where the integral driver can store intermediates
   ///   \p nWork     Size of work space, in doubles. If <= 0, pWork is not used but memory is allocated form the global heap.
   ///   \p Thr       Threshold for integral/basis function evaluation
   FORTINT FD(create_integral_context)(void *pWork, FORTINT const &nWork, double const &Thr);
//    #define CreateIntegralContext FD(create_integral_context)

   /// initialize an integral kernel for evaluating integrals of a given type.
   ///   \p iKernel   INTKERNEL_* constant defining the kernel
   ///   \p pParamsI  Integer parameters to the kernel (use depends on kernel. Can be 0 if unused)
   ///   \p pParamsF  Float parameters to the kernel (use depends on kernel, can be 0 if unused)
   void FD(assign_integral_kernel)(FORTINT &iContext, FORTINT const &iKernel, FORTINT *pParamsI, double *pParamsF);
//    #define AssignIntegralKernel FD(assign_integral_kernel)

   /// deallocate a kernel object.
   void FD(destroy_integral_context)(FORTINT &iContext);
//    #define DestroyIntegralContext FD(destroy_integral_context)

   /// evaluate 1-electron integrals < a|krn|b>. where a is the shell-group BasisA[iGrpA]
   /// and b is the shell-group BasisB[iGrpB].
   ///   \p pOut      Output will be written to pOut[iFnA*Strides[0] + iFnB*Strides[1]]
   ///   \p Strides   row/column strides of output data. [0] is for A, [1] is for B
   ///   \p Factor    prefactor: integral data will be multiplied by this.
   ///   \p iContext  Integral context defining the parameters (kernel, thresholds, etc).
   ///                \see CreateIntegralContext
   void FD(eval_group_int1e)(double *pOut, FORTINT *Strides, double const &Factor,
      FORTINT const &iGrpA, FBasisSet const &BasisA,
      FORTINT const &iGrpB, FBasisSet const &BasisB, FORTINT &iContext);

   /// evaluate 2-electron integrals (ab|krn|cd). where x=a/b/c/d are the shell-groups
   /// BasisABCD[iGrpX], but their centers are translated by the vectors vTraX
   /// (note: in Fortran, read this as 'vTraX :: double(3)'.
   /// All groups are taken from the same basis set, BasisABCD.
   /// This function *INCREMENTS* the output buffer! It does not clear it.
   ///   \p pOut      Output will be written to pOut[iFnA*Strides[0] + iFnB*Strides[1] + ...]
   ///   \p Strides   row/column strides of output data. [0] is for A, [1] is for B
   ///   \p Factor    prefactor: integral data will be multiplied by this.
   ///   \p iContext  Integral context defining the parameters (kernel, thresholds, etc).
   ///                \see CreateIntegralContext
   void FD(eval_group_int2e_tra_incr)(double *pOut, FORTINT *Strides, double const &Factor,
      FORTINT const &iGrpA, FVector3 const &vTraA,
      FORTINT const &iGrpB, FVector3 const &vTraB,
      FORTINT const &iGrpC, FVector3 const &vTraC,
      FORTINT const &iGrpD, FVector3 const &vTraD,
      FBasisSet const &BasisABCD, FORTINT &iContext);

   /// evaluate 1-electron integrals < a|krn|b>. where a and b run over the entire basis A/B.
   ///   \p pOut      Output will be written to pOut[iFnA*Strides[0] + iFnB*Strides[1]]
   ///   \p Strides   row/column strides of output data. [0] is for A, [1] is for B
   ///   \p Factor    prefactor: integral data will be multiplied by this.
   ///   \p iContext  Integral context defining the parameters (kernel, thresholds, etc).
   ///                \see CreateIntegralContext
   void FD(eval_basis_int1e)(double *pOut, FORTINT *Strides, double const &Factor,
      FBasisSet const &BasisA, FBasisSet const &BasisB, FORTINT &iContext);

   /// evaluate contracted 2-electron integrals \sum_c (ab|krn|c) pCoeffC[c],
   /// where a/b are the shell-group BasisA[iGrpA]/BasisB[iGrpB], and
   /// c runs over the entire basis C. The contraction vector pCoeffC[iFnC]
   /// gives the coefficients of the basis functions in C.
   ///
   /// Note: Can be used in conjunction with with MakePointCharge basis on C
   /// for short-range n/e integrals. In that case pCoeffC are the point charges.
   ///   \p pOut      Output will be written to pOut[iFnA*Strides[0] + iFnB*Strides[1]]
   ///   \p Strides   row/column strides of output data. [0] is for A, [1] is for B
   ///   \p Factor    prefactor: integral data will be multiplied by this.
   ///   \p pCoeffC   contraction coefficients of basis functions in C.
   ///   \p iContext  Integral context defining the parameters (kernel, thresholds, etc).
   ///                \see CreateIntegralContext
   void FD(eval_group_int2e_contract_v)(double *pOut, FORTINT *Strides, double const &Factor,
      FORTINT const &iGrpA, FBasisSet const &BasisA,
      FORTINT const &iGrpB, FBasisSet const &BasisB,
      double const *pCoeffC, FBasisSet const &BasisC, FORTINT &iContext);

   /// Same as eval_group_int2e_contract_v, but for the entire basis.
   void FD(eval_basis_int2e_contract_v)(double *pOut, FORTINT *Strides, double const &Factor,
      FBasisSet const &BasisA, FBasisSet const &BasisB,
      double const *pCoeffC, FBasisSet const &BasisC, FORTINT &iContext);

   /// evaluate contracted 2-electron integrals \sum_c (ab|krn|c) pCoeffC[c],
   /// where c runs over the point charges specified by pCentersC
   ///
   /// Note: BasisB will be periodically symmetrized with the super-cell translation
   /// length, but the centers on C will be taken as they are. Therefore, in a short-range
   /// case, all images to consider must be explicitly unpacked and supplied.
   void FD(eval_basis_int2e_contract_point_charges)(double *pOut, FORTINT *Strides, double const &Factor,
      FBasisSet const &BasisA, FBasisSet const &BasisB,
      double const *pCoeffC, FVector3 const *pCentersC, size_t nCentersC, FORTINT &iContext);

   /// evaluate basis functions (of the entire basis) on a grid.
   ///   \p pOut        Output will be written to pOut[iGridPt + nGridPt * (iComp + nCompSt * iMap)]. Must hold room for nGridPt x nCompSt x nBasisfn entries.
   ///   \p nCompSt     Derivative component strides: Normally 1 for densities, 4 for densities+gradients. Can be larger.
   ///   \p pCentersOut [iMap]: center index of basis function Map[iMap]. Note: 0-based.
   ///   \p pMap        Output: Indices of basis functions retained.
   ///   \p nMap        Output: number of basis functions retained
   ///   \p Basis       Basis functions to evaluate
   ///   \p nGridPt     Number of points supplied in pGridPt
   ///   \p pGridPt     Array of 3d grid points on which the functions will be evaluated.
   ///   \p DerivOrder  0: make just densities. 1: make density and d/dx..d/dz, 2: make up to 2nd derivatives.
   ///   \p iContext    Integral context defining the parameters. Threshold is taken from here.
   ///                  \see CreateIntegralContext
   /// FIXME: doesn't add periodic super-cell images yet!
   void FD(eval_basis_fn_on_grid)(double *pOut, FORTINT const &nCompSt,
      FORTINT *pCentersOut, FORTINT *pMap, FORTINT &nMap,
      FBasisSet const &Basis, double  (*pGridPt)[3], FORTINT const &nGridPt,
      FORTINT const &DerivOrder, FORTINT &iContext);
}








#endif // PHF_BASIS_SET_H

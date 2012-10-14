#include <boost/format.hpp>

#include "CtBasisLibrary.h"

#include "PhfSolidDef.h"
#include "PhfBasisSet.h"
#include "PhfScf.h"

#include "CtIo.h"
using namespace ct;
using boost::format;



int main(int argc, char *argv[])
{
   // load basis set libraries
   g_BasisSetLibrary.ImportMolproLib("libmol/cp2k-gth.libmol");


   // setup some dummy objects. These things should later be handled
   // by some input/output infrastructure.
//    FScfOptions
//       ScfOptions = {1e-8, 1e-8};
   FSolidModel
      Solid;

   // begin with a slightly skewed lattice. Should catch
   // various problems with orthogonality.
   double f = 1.4;
   Solid.Lattice = FLattice(
      f*FVector3(1.0,0.08,0.00),
      f*FVector3(0.0,1.00,0.12),
      f*FVector3(0.0,0.05,1.00) );

   FBasisDescs DefaultBases;
   DefaultBases[BASIS_Orbital] = "SZV-GTH";
   char const *
      pUnitCellAtomsXyz =
         "2                                    "
       "\n(this is xyz format--but a.u input!) "
       "\nH   0.0000   0.100   -0.710          "
       "\nH   0.1200   0.000   +0.710          "
       "\n";
   Solid.UnitCell.AddAtomsFromXyzData(pUnitCellAtomsXyz, 1.0, DefaultBases);
   Solid.UnitCell.Volume = Solid.Lattice.UnitCellVolume;
   Solid.SuperCell.Init(FVector3i(4,4,4), Solid.Lattice, Solid.UnitCell);
   Solid.UnitCell.OrbBasis.SetPeriodicityVectors(Solid.SuperCell.T);

   if ( 1 ) {
      FOpMatrix
         Overlap(Solid),
         Kinetic(Solid);
      FORTINT
         ic = FD(create_integral_context)(0,0, 1e-10),
         Strides[2] = {1, Overlap.nRows};
      FD(assign_integral_kernel)(ic, INTKERNEL_Overlap, 0, 0);
      FD(eval_basis_int1e)(&Overlap[0], Strides, 1.0, Solid.UnitCell.OrbBasis, Solid.SuperCell.OrbBasis, ic);
      FD(assign_integral_kernel)(ic, INTKERNEL_Kinetic, 0, 0);
      FD(eval_basis_int1e)(&Kinetic[0], Strides, 1.0, Solid.UnitCell.OrbBasis, Solid.SuperCell.OrbBasis, ic);
      FD(destroy_integral_context)(ic);

      Overlap.Print(xout, "OVERLAP UnitCell x SuperCell");
      Kinetic.Print(xout, "KINETIC UnitCell x SuperCell");
   }

   xout << format("wheee!!") << std::endl;
};





// define scf:
//    overlap & kinetic energy build in nAo x nAo x nSuperCell.nSize
//
//    (transformation to k-space basis for smh and orbital update)
//
//    back-transformation.
//
//    fock build in nAo x nAo x nSuperCell.nSize
//
//



// can I just use \mu(r) = \bar \sum_T exp(ikT) \mu(r+T)
// as symmetry adapted basis functions?
// (instead of exp(ikr) \sum_T \mu(r+T))

// for Gaussian poission solver: I can also do the SR/LR split and
// run the short range part through BAS_ACT/BAS_VCA. Combined with
// the core distribution.




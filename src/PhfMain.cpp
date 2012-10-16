#include <boost/format.hpp>

#include "lib/CtBasisLibrary.h"

#include "PhfSolidDef.h"
#include "PhfBasisSet.h"
#include "PhfScf.h"

#include "lib/CtIo.h"
using namespace ct;
using boost::format;

#include "Ewald_pt1.h"


extern "C" {
    void environment_report_();
    void print_unit_cell_(FUnitCell&);
    void exchangesum_(double&,FOpMatrix&,FLattice&,FUnitCell&,FSuperCell&,FOpMatrix&); 
}

int main(int argc, char *argv[])
{

   environment_report_(); 

   // load basis set libraries
   g_BasisSetLibrary.ImportMolproLib("../libmol/cp2k-gth.libmol");


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
   Solid.SuperCell.Init(FVector3i(2,2,2), Solid.Lattice, Solid.UnitCell);
   Solid.UnitCell.OrbBasis.SetPeriodicityVectors(Solid.SuperCell.T);

   // allocate memory for density matrix, exchange matrix and coulomb matrix
   FOpMatrix Density(Solid), Exchange(Solid), Coulomb(Solid);

   if ( 0 ) {
      FOpMatrix
         Overlap(Solid),
         Kinetic(Solid);
      FORTINT
         ic = create_integral_context_(0,0, 1e-10),
         Strides[2] = {1, Overlap.nRows};
      assign_integral_kernel_(ic, INTKERNEL_Overlap, 0, 0);
      eval_basis_int1e_(&Overlap[0], Strides, 1.0, Solid.UnitCell.OrbBasis, Solid.SuperCell.OrbBasis, ic);
      assign_integral_kernel_(ic, INTKERNEL_Kinetic, 0, 0);
      eval_basis_int1e_(&Kinetic[0], Strides, 1.0, Solid.UnitCell.OrbBasis, Solid.SuperCell.OrbBasis, ic);
      destroy_integral_context_(ic);

      Overlap.Print(xout, "OVERLAP UnitCell x SuperCell");
      Kinetic.Print(xout, "KINETIC UnitCell x SuperCell");
   }

   if ( 0 ) {
      // evaluating short-range part of point-charge lattice:
      FOpMatrix
         Nuclear_ShortRange(Solid);
      FORTINT
         ic = create_integral_context_(0,0, 1e-10),
         Strides[2] = {1, Nuclear_ShortRange.nRows};
      double
         // the omega of  g(r) = erfc(omega r)/r
         Omega = 5.;
      TArray<double>
         PointCharges;
      TArray<FVector3>
         PointCenters;
      // fill up PointCharges and PointCenters with all the nuclear images
      // within screening range of erfc(omega*r)/r.

      //  [ ... ]

      // evaluate the integrals.
      assign_integral_kernel_(ic, INTKERNEL_Coulomb_ShortRange_Erfc, 0, &Omega);
      eval_basis_int2e_contract_point_charges_(&Nuclear_ShortRange[0],
         Strides, 1.0, Solid.UnitCell.OrbBasis, Solid.SuperCell.OrbBasis,
         &PointCharges[0], &PointCenters[0], PointCenters.size(), ic);
      destroy_integral_context_(ic);
   }

   xout << format("wheee!!") << std::endl;


   if (0) {

      double eta = 1.2;
      Ewald_pt1 MrEwald = Ewald_pt1(eta, Solid);
      int lmax = 1;
      int dimlin = (lmax+1)*(lmax+2)/2;
      double * result = new double[dimlin*dimlin];
      MrEwald.eN(result, 1.1, lmax, 0.0, 0.0, 0.0, 1.2, lmax, 2.0, 1.9, 0.7); //Random numbers to do compile and run test.
      delete [] result;
      MrEwald.NN();

   }


   // test call of fortran
   double ExchangeEnergy;
   exchangesum_(ExchangeEnergy,Exchange,Solid.Lattice,Solid.UnitCell,Solid.SuperCell,Density);

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




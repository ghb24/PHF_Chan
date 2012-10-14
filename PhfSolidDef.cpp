#include <sstream>
#include <stdexcept>
#include "CtBasisLibrary.h"
#include "AicShells.h"   // <- that's what the basis library stores.

#include "PhfSolidDef.h"

namespace ct {
   uint ElementNumberFromName( std::string const &Name );
   std::string ElementNameFromNumber( uint AtomicNumber );
   std::string tolower(std::string const&);
}
uint ElementNumberFromName( std::string const &Name ) { return ct::ElementNumberFromName(Name); }
std::string ElementNameFromNumber(uint AtomicNumber) { return ct::ElementNameFromNumber(AtomicNumber); }


void FUnitCell::AddAtom(FVector3 const &vPos, int Element, FBasisDescs const &DefaultBases)
{
   // get orbital basis functions from basis library.
   // Other contexts ignored for now.
   std::vector<aic::FGaussShell>
      Shells;
   ct::g_BasisSetLibrary.LoadBasisFunctions(Shells, Element,
      ct::tolower(DefaultBases.find(BASIS_Orbital)->second), vPos, Coords.size());

   // convert from AIC shells to our flattened basis format.
   for ( uint i = 0; i < Shells.size(); ++ i )
      OrbBasis.AddAicShell(Shells[i]);
   OrbBasis.Finalize();

   Coords.push_back(vPos);
   Elements.push_back(Element);
   EcpCharges.push_back(0);
}


void FUnitCell::AddAtomsFromXyzData(char const *pFileText, double InputToAuFactor, FBasisDescs const &DefaultBases)
{
   // c/p'd from ct::FAtomSet.
   std::stringstream
      str(pFileText);
   uint
      nAtoms;
   str >> nAtoms;
   Coords.reserve(nAtoms);
   Elements.reserve(nAtoms);
   EcpCharges.reserve(nAtoms);

   str.ignore( 0xffff, '\n');
   // now in second line (comment line)
   std::string
      CurLine;
   // read comment line.
   std::getline( str, CurLine );
   for ( uint nAtom = 0; nAtom < nAtoms; ++ nAtom ){
      std::string
         Element;
      double
         x, y, z;
      str >> Element >> x >> y >> z;
      if ( str.bad() )
         throw std::runtime_error("AddAtomsFromXyzFile: FAILED to load .xyz File '" + std::string(pFileText) + "'. Actual format not understood.");

      Element = ct::tolower(Element);
      double
         f = InputToAuFactor;
      FVector3
         vPos( f*x, f*y, f*z );
      AddAtom(vPos, ElementNumberFromName(Element), DefaultBases);
      str.ignore( 0xffff, '\n');
   }
}

FLattice::FLattice(FVector3 const &T0, FVector3 const &T1, FVector3 const &T2)
{
   T[0] = T0;
   T[1] = T1;
   T[2] = T2;
   Init();
};

FLattice::FLattice(FVector3 const (&T_)[3])
   : T(T_)
{
   Init();
}

void FLattice::Init()
{
   // construct the bi-orthogonal basis of the k-space. 2pi factors are absorbed.
   Cross(K[0], T[1], T[2]);
   UnitCellVolume = Dot(K[0], T[0]);
   double
      f = 2*M_PI/UnitCellVolume;
   K[0] *= f;
   K[1]  = f*Cross(T[2], T[0]);
   K[2]  = f*Cross(T[0], T[1]);
};

void FSuperCell::Init(FVector3i const &Size_, FLattice const &Lattice, FUnitCell const &UnitCell)
{
   Size = Size_;
   // setup super-cell translation vectors.
   for ( uint i = 0; i < 3; ++ i )
      T[i] = (double)Size[i] * Lattice.T[i];
   // ...and unit-cell translation vectors.
   Ts.clear();
   Ts.reserve(Size[0] * Size[1] * Size[2]);
   for ( uint iz = 0; iz < Size[2]; ++ iz )
      for ( uint iy = 0; iy < Size[1]; ++ iy )
         for ( uint ix = 0; ix < Size[0]; ++ ix )
            Ts.push_back(double(ix)*Lattice.T[0] +
                         double(iy)*Lattice.T[1] +
                         double(iz)*Lattice.T[2]);

   // make a basis set object for the super-cell. All we have to do is translate
   // the unit-cell basis by all our Ts vectors.
   OrbBasis = FBasisSet(UnitCell.OrbBasis, &Ts[0], Ts.size());
   OrbBasis.SetPeriodicityVectors(T);
};


FSuperCell::FSuperCell(FVector3i const &Size_, FLattice const &Lattice, FUnitCell const &UnitCell)
{
   Init(Size_, Lattice, UnitCell);
};




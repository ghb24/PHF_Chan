#ifndef PHF_SOLID_DEF_H
#define PHF_SOLID_DEF_H

#include <map>
#include <string>

#include "PhfTypes.h"
#include "PhfBasisSet.h"


// conversion element name <-> atomic number
std::string ElementNameFromNumber( uint AtomicNumber );
uint ElementNumberFromName( std::string const &Name );

enum FBasisContext {
   BASIS_Orbital,
   BASIS_Ecp,
   BASIS_JFit,
   BASIS_JkFit,
   BASIS_Mp2Fit
};
typedef std::map<FBasisContext, std::string>
   FBasisDescs;




// note:
//  - The type definitions here are designed with the primary
//    purpose of being binary compatible with something we can pass
//    to Fortran (that is: integer arrays and double arrays instead
//    of compound types).
//  - For that reason I am also reluctant to put in constructors.
//    Having a constructor will technically make the classes non-POD
//    types and by that allow the compiler to re-align their memory
//    layouts in arbitrary ways (although they might not actually do that)


/// defines atoms in the unit-cell and their properties
struct FUnitCell
{
   TArray<FVector3>
      /// Coord[i][ixyz]: xyz coordinates of the atoms
      Coords;
   TArray<FORTINT>
      /// Identifies atom types. E.g., 2 == He
      Elements;
   TArray<double>
      /// [i]: number of electrons of atom #i which are absorbed
      /// in its effective core potential/pseudo potential
      EcpCharges;

   FBasisSet
      /// Gaussian basis used to expand the molecular orbitals
      /// ...with centers on the unit-cell.
      /// Note: Actual basis functions are periodically repeated
      /// with the super-cell frequency!
      OrbBasis;

   double
      /// 3d volume (=det(Lattice.T))
      /// (should this be here? technically it is a lattice property, not
      ///  a unit-cell property)
      Volume;

   void AddAtom(FVector3 const &Pos, int Element, FBasisDescs const &DefaultBases);
   /// add an atom from xyz-file; content of file is provided in pFileText.
   /// Basis sets and atomic properties are assigned via DefaultBases.
   void AddAtomsFromXyzData(char const *pFileText, double InputToAuFactor, FBasisDescs const &DefaultBases);
};

/// defines the physical lattice: allowed periodic translations of the unit-cell.
// (note: Lattice symmetry stuff, classical lattice summations(*), and lattice
//  reduction algorithms should probably go here.
//  (*): That is a routine which takes a number of point charges in the unit
//       cell (neutral in total), and a (large) number of real-space grid
//       points, and evaluates the classical Madelung potential on the grid.)
struct FLattice
{
   FVector3
      /// real-space lattice vectors (translations between two unit-cells)
      T[3],
      /// k-space lattice vectors (2pi*bi-orth basis of UnitCellT).
      K[3];
   // ^- all of those should be binary compatible with a 3x3 double matrix
   //    at &xxx[0][0], but need to look that up to make sure.
   double
      UnitCellVolume;

   FLattice() {};
   FLattice(FVector3 const &T0, FVector3 const &T1, FVector3 const &T2);
   explicit FLattice(FVector3 const (&T_)[3]);
private:
   void Init(); // initialize derived members from T.
};

/// defines our calculation model: number of unit-cells explicitly treated
/// before being periodically repeated.
/// Notes:
///   - Unit-cell placement: First unit-cell is at T = (0,0,0), and other
///     unit-cells are displaced in *positive* Lattice.T direction (i.e.,
///     the first unit-cell is at the super-cell boundary).
///     This might simplify some implementation aspects, and theoretically
///     it should not matter since the calculation is periodic)
///   - As the k-point grid is effectively a different representation of
///     the super-cell information, it should probably also go in here
///     once we get to that point.
struct FSuperCell
{
   FVector3
      /// real-space displacements between two super-cells in the
      /// directions Lattice.T[i].
      /// This is the same as the translation between a raw Gaussian
      /// basis function and its periodic images.
      /// (note: SuperCell.T[i] == Size[i] * Lattice.T[i])
      T[3];
   TVector3<FORTINT>
      /// [i]: number of times the unit-cell is repeated in Lattice.T[i]
      /// direction to form a super-cell. In total we obtain a
      ///    Size[0] x Size[1] x Size[2]
      /// grid of super-cells.
      ///
      /// Note: Total number of super-cells is SuperCell.Ts.size().
      Size;
   // ^- FIXME: name is bad. Come up with something better.
   TArray<FVector3>
      /// set of all translation vectors between the first unit-cell
      /// and the other unit-cells in the super-cell (i.e., T x Size unpacked)
      /// Note: The 0-vector (first to first) is included in the set.
      Ts;

   FBasisSet
      /// Gaussian basis used to expand the molecular orbitals
      /// ...with centers on the *entire* super-cell.
      OrbBasis;

   FSuperCell() {};
   FSuperCell(FVector3i const &Size, FLattice const &Lattice, FUnitCell const &UnitCell);
   void Init(FVector3i const &Size, FLattice const &Lattice, FUnitCell const &UnitCell);
};

/// defines the actual calculation model we will be dealing with
// (maybe skip that and add a FSolidHfContext directly?
//  the difference is that the latter would contain temporary data
//  used only during the optimization, while the solid model actually
//  specifies what is supposed to happen (in particular, it could be
//  run with multiple different paramters (e.g., HF and KS)...
struct FSolidModel
{
   FUnitCell
      UnitCell;
   FLattice
      Lattice;
   FSuperCell
      SuperCell;
};







#endif // PHF_SOLID_DEF_H

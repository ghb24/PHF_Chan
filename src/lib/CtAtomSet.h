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

#ifndef CT8K_ATOMSET_H
#define CT8K_ATOMSET_H

#include <map>
#include <string>

#include "CtCommon.h"

namespace ct {

// conversion element name <-> atomic number
std::string ElementNameFromNumber( uint AtomicNumber );
uint ElementNumberFromName( std::string const &Name );

enum FBasisContext {
   BASIS_Orbital,
   BASIS_JFit,
   BASIS_JkFit,
   BASIS_Mp2Fit,
   BASIS_CcsdFit,
   BASIS_F12RI,
   BASIS_Guess
};

std::string BasisContextName( FBasisContext Context );


typedef aic::FVector3
   FVector3;

typedef std::string
   FBasisDesc;
typedef std::string
   FEcpDesc;
typedef std::map<FBasisContext, FBasisDesc>
   FBasisDescs;

struct FDoublettIntegralFactory;
struct FMatrixView;

struct FAtom
{
   FVector3
      vPos; // center of atom, in a_bohr (always in a_bohr)
   uint
      AtomicNumber;
   uint
      // number of electrons implicit due to effective
      // core potential ("pseudopotential")
      nEcpElectrons;
   int
      Charge; // might differ from AtomicNumber if ECPs are used.
   FBasisDescs
      BasisDesc;
   FEcpDesc
      EcpDesc;

   std::string const &GetOrbBasis() const;
   std::string GetElementName() const;

   FAtom() {};
   FAtom( FVector3 const &vPos_, std::string const &ElementName_, std::string
            const &BasisDesc_ = "def2-SVP", std::string const &Ecp_ = "" );
   FAtom( FVector3 const &vPos_, std::string const &ElementName_, FBasisDescs
            const &BasisDesc_, std::string const &Ecp_ = ""  );

private:
   void Init( FVector3 const &vPos_, std::string const &ElementName_, std::string const &EcpName_ );
};


// describes global effects acting on the current calculation. For example,
// external fields.
struct FEnvironment
{

};

struct FBasisSet;

struct FAtomSymGroup {
	// number of equivalnt atoms
	uint nEqiv;
	// index into Atoms[]
	uint iEqiv[8];
	// symmetry op moving iEqiv[0] to iEqiv[i]
	uint SymOp[8];
	// bit pattern of XYZ symmetry op not affecting the atom positions,
	// since the relevant coordinates are zero.
	uint Stabilizer;
};

// An atom set specifies which atoms are where and which properties these
// atoms do have in calculations (i.e. which basis/ecp to use and so on).
// Additionally,
// So an instance of this class basically describes the problem which
// we currently handle.
class FAtomSet : public FIntrusivePtrDest
{
public:
	typedef std::vector<FAtom>
		FAtomList;
	FAtomList
		Atoms;
	FEnvironment
		Environment;


	// ecp-reduced nuclear repulsion energy.
	FScalar NuclearRepulsionEnergy() const;
	void 	AddAtom( FAtom const &Atom );
	// ecp-reduced total nuclear charge.
	int 	NuclearCharge() const;

	FVector3 NuclearDipoleMoment( FVector3 const &ExpansionPoint = FVector3(0,0,0) ) const;

	// returns ecp-reduced number of electrons in rare gas configurations.
	uint	nCoreElectrons() const;

	struct FXyzLoadOptions{
		typedef std::map<std::string,std::string>
			FElementBasesMap;
		FElementBasesMap
			// defines basis strings for elements which are loaded. If an element
			// is not present here, the standard basis string (e.g. "dzp") is
			// assigned for it. Case of element names (keys) don't matter.
			ElementBases,
			ElementEcps;
		FScalar
			// if the file coordinates are not specified in angstroms (e.g.
			// if they are in a.u. (=a_bohr) or something), this entry allows
			// you to specify a factor which converts them to A. Default
			// value is 1.0.
			InputCoordsToAngstromFactor;
		FXyzLoadOptions();
	};
	// Adds the molecules out of the Rasmol XYZ file to the current atom set.
	// Atoms will be added with basis string StandardBasis unless some other
	// basis is specified in *pOtherOptions.
	void AddAtomsFromXyzFile( std::string const &FileName,
		FBasisDescs const &DefaultBases,
		FXyzLoadOptions const *pOtherOptions = 0 );



	// get up to 8 symmetry operators of the subgroups of D2h.
	// Each op: bit 0 -> x, bit 1->y, bit 2->z; E.g., Op 011b means invariance
	// to simultaneous mirroring around x and y.
	// Gen[i] indexes into Ops. It designates which of the operators form the minimal
	// generator set on which the symmetry adapted basis functions are supposed to
	// be based.
	void FindMirrorSymmetryOps(uint Ops[8], uint &nOps, uint Gen[3], uint &nGen) const;
	void FindEquivalentAtoms(FAtomSymGroup *&pGroups, uint &nGroups, FMemoryStack &Mem) const;


	// makes a matrix representing a given one electron operator, and adds it
	// to Dest (which must either be 0x0 or have compatible dimensions).
	// The operator here can be specified generally using a fitting integral
	// driver object for it. For some operators of interest, direct evaluation
	// functions are provided (see next functions).
	// FDoublettIntegralFactory and derived classes are described in integrals.h
	void Add1eIntegralMatrix( FMatrixView &Dest,
		FBasisSet const &RowBasis, FBasisSet const &ColBasis,
		FDoublettIntegralFactory &IntFactory, double Factor, FMemoryStack &Mem ) const;
	// the following functions are more or less dummy functions that just call
	// Make1eIntegralMatrix.
	void MakeCoreHamiltonMatrix( FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem ) const;
	void MakeOverlapMatrix( FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem ) const;
	void AddKineticMatrix( FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem ) const;
	void AddNuclearAttractionMatrix( FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem ) const;
	void AddDipoleMatrix( FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis,
		FVector3 const &Direction, FVector3 const &ExpansionPoint, FMemoryStack &Mem ) const;

public:
	uint size() const { return Atoms.size(); }
	FAtom const &operator [] (uint iAt) const { return Atoms[iAt]; }
	FAtom &operator [] (uint iAt) { return Atoms[iAt]; }
	void clear() { Atoms.clear(); }
	bool empty() const { return Atoms.empty(); }
};

std::ostream &operator << ( std::ostream &out, FAtomSet const& );
std::ostream &operator << ( std::ostream &out, FVector3 const &Pos );

} // namespace ct


#endif // CT8K_ATOMSET_H

// kate: space-indent off; tab-indent on; indent-width 4; mixedindent off; indent-mode normal;

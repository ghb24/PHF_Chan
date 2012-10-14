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
#include <cmath>

#include "AicShells.h"

#include "CtCommon.h"
#include "CtAtomSet.h"
#include "CtIo.h"
#include "CtMatrix.h"
#include "CtBasisSet.h"
#include "CtInt1e.h"
#include "CtConstants.h"

namespace ct {

double
	SymmetryTolerance = 1e-10;


FScalar
	// atomic positions will be multiplied with this value when FAtom(pos)
	// is called. This may be useful to convert A to a.u. or reverse.
	g_AtomPositionInputFactor = 1.0;


static uint const
	s_nElementNames = 87;
static char const * // element names in periodic system order.
	s_ElementNames[s_nElementNames] = { "H", "He", "Li", "Be", "B", "C", "N",
			"O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
			"K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
			"Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
			"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
				"Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
			"XX" // special kind of dummy atom without charge.
	};


// conversion element name <-> atomic number
std::string ElementNameFromNumber( uint AtomicNumber )
{
	if ( AtomicNumber == 0 )
		return "XX";
	assert( ( AtomicNumber >= 1 ) && ( AtomicNumber <= s_nElementNames ) );
	if ( false == ( AtomicNumber >= 1 ) && ( AtomicNumber <= s_nElementNames ) )
		throw std::range_error((format("No name for atom with number '%i' present.") % AtomicNumber).str());
	return std::string( s_ElementNames[AtomicNumber-1] );
};


uint ElementNumberFromName( std::string const &Name )
{
	if ( strcasecmp( Name.c_str(), "XX" ) == 0 )
		return 0;
	for ( uint i = 0; i < s_nElementNames; ++i )
		if ( 0 == strcasecmp( Name.c_str(), s_ElementNames[i] ) )
			return 1+i;
	throw std::runtime_error( "Failed to look up atomic number for element symbol '"+Name+"'" );
};


std::ostream &operator << ( std::ostream &out, FVector3 const &Pos ){
	out << format("(%7.5f %7.5f %7.5f)") %  Pos[0] % Pos[1] % Pos[2];
	return out;
}


FScalar FAtomSet::NuclearRepulsionEnergy() const
{
	FScalar
		fEnergy = 0;
	FAtomList::const_iterator
		itAtom1, itAtom2;
	for ( itAtom1 = Atoms.begin(); itAtom1 != Atoms.end(); ++itAtom1 )
		for ( itAtom2 = itAtom1, ++itAtom2; itAtom2 != Atoms.end(); ++itAtom2 )
			fEnergy += itAtom1->Charge * itAtom2->Charge *
						1.0/Distance( itAtom1->vPos, itAtom2->vPos );
	return fEnergy;
};

FVector3 FAtomSet::NuclearDipoleMoment( FVector3 const &ExpansionPoint ) const
{
	FVector3
		Result(0,0,0);
	FAtomList::const_iterator
		itAtom;
	for ( itAtom = Atoms.begin(); itAtom != Atoms.end(); ++itAtom ){
		Result += static_cast<FScalar>( itAtom->Charge ) *
			( itAtom->vPos - ExpansionPoint );
	}
	return Result;
};



// ecp-reduced total nuclear charge.
int FAtomSet::NuclearCharge() const
{
	int
		Result = 0;
	FAtomList::const_iterator
		itAtom;
	for ( itAtom = Atoms.begin(); itAtom != Atoms.end(); ++itAtom )
		Result += itAtom->Charge;
	return Result;
};


void FAtomSet::AddAtom( FAtom const &Atom )
{
	Atoms.push_back( Atom );
};

// converts a map with string keys into another map of the same format but
// which's keys have been lowercased.
template<class FValue>
std::map<std::string,FValue> MakeStringKeysLowercaseInMap(
	std::map<std::string,FValue> In )
{
	std::map<std::string,FValue>
		Out;
	typename std::map<std::string,FValue>::iterator
		it;
	_for_each( it, In )
		Out[tolower(it->first)] = it->second;
	return Out;
};


// Adds the molecules out of the Rasmol XYZ file to the current atom set.
// Atoms will be added with basis string StandardBasis unless some other
// basis is specified in *pOtherOptions.
void FAtomSet::AddAtomsFromXyzFile( std::string const &FileName,
	FBasisDescs const &DefaultBases,
	FXyzLoadOptions const *pOtherOptions )
{
	// the .xyz file format is quite popular because it is extremely simple.
	// format is as follows:
	// - first line: number of atoms =: N
	// - second line: some title or information
	// - N x [Element Name Xcoord Ycoord Zcoord]. Elements may be separated
	//   by whitespace. Following the lines some colons may stand or don't
	//   stand.
	TArray<char>
		pFileContent;
	if ( false == LoadFileIntoMemory( pFileContent, FileName ) )
		throw std::runtime_error("AddAtomsFromXyzFile: FAILED to open/read input file: '" + FileName + "'" );
	std::stringstream
		str( &pFileContent[0], std::stringstream::in );

	FXyzLoadOptions
		LoadOptions;
	if ( 0 != pOtherOptions ){
		LoadOptions = *pOtherOptions;
		LoadOptions.ElementBases = MakeStringKeysLowercaseInMap( pOtherOptions->ElementBases );
		LoadOptions.ElementEcps = MakeStringKeysLowercaseInMap( pOtherOptions->ElementEcps );
	};

	uint
		nAtoms;
	str >> nAtoms;
	str.ignore( 0xffff, '\n');
	// now in second line (comment line)
	std::string
		CurLine;
	// read comment line.
	std::getline( str, CurLine );
	for ( uint nAtom = 0; nAtom < nAtoms; ++ nAtom ){
		std::string
			Element;
		FScalar
			x, y, z;
		str >> Element >> x >> y >> z;
		if ( str.bad() )
			throw std::runtime_error("AddAtomsFromXyzFile: FAILED to load .xyz File '" + FileName + "'. Actual format not understood.");

		Element = tolower( Element );
/*		FXyzLoadOptions::FElementBasesMap::iterator
			itBasis, itEcp;
		itBasis = LoadOptions.ElementBases.find( Element );
		if ( itBasis != LoadOptions.ElementBases.end() )
			Basis = itBasis->second;
		itEcp = LoadOptions.ElementBases.find( Element );
		if ( itEcp != LoadOptions.ElementBases.end() )
			Ecp = itEcp->second;*/

		FScalar
			// our internal data structures are actually input in a_bohrs.
			f = LoadOptions.InputCoordsToAngstromFactor/ToAng;
		FVector3
			vPos( f*x, f*y, f*z );

		AddAtom( FAtom(vPos, Element, DefaultBases) );
		str.ignore( 0xffff, '\n');
	}


};


// returns ecp-reduced number of orbitals in rare gas configurations.
uint FAtomSet::nCoreElectrons() const
{
	uint
		Result = 0;
	uint
		nRareGasAtoms = 7,
		RareGasAtomicNumbers[] = { 0, 2, 10, 18, 36, 54, 86 };
	FAtomList::const_iterator
		itAtom;
	_for_each( itAtom, Atoms ){
		uint
			nAtomCoreOrbitals = 0;
		for ( uint i = 0; i < nRareGasAtoms; ++ i ){
			if ( itAtom->AtomicNumber > RareGasAtomicNumbers[i] ){
				nAtomCoreOrbitals = RareGasAtomicNumbers[i];
				nAtomCoreOrbitals -= itAtom->nEcpElectrons;
				assert_rt( nAtomCoreOrbitals >= 0 );
			} else
				break;
		}
		Result += nAtomCoreOrbitals;
	}
	return Result;
};



FAtomSet::FXyzLoadOptions::FXyzLoadOptions()
	 : InputCoordsToAngstromFactor( 1.0 )
{};


std::string const &FAtom::GetOrbBasis() const
{
	FBasisDescs::const_iterator
		itDesc = BasisDesc.find(BASIS_Orbital);
	assert(itDesc != BasisDesc.end());
	return itDesc->second;
}

FAtom::FAtom( FVector3 const &vPos_, std::string const &ElementName_, std::string
	const &BasisDesc_, std::string const &Ecp_ )
{
	Init(vPos_, ElementName_, Ecp_);
	BasisDesc[BASIS_Orbital] = BasisDesc_;
};

FAtom::FAtom( FVector3 const &vPos_, std::string const &ElementName_, FBasisDescs
			const &BasisDesc_, std::string const &Ecp_ )
	: BasisDesc(BasisDesc_)
{
	Init(vPos_, ElementName_, Ecp_);
};

void FAtom::Init( FVector3 const &vPos_, std::string const &ElementName_, std::string const &EcpName_ )
{
	vPos = vPos_;
	AtomicNumber = ElementNumberFromName(ElementName_);
	nEcpElectrons = 0;
	Charge = AtomicNumber;
	if ( EcpName_.empty() )
		EcpDesc = "(none)";
	else
		EcpDesc = EcpName_;
};



std::ostream& operator << ( std::ostream &out, FAtom const &a )
{
	out << boost::format( "%+2s:   (%10.6f,%10.6f,%10.6f)    %s" )
		% ElementNameFromNumber( a.AtomicNumber )
		% Distance( a.vPos.x(), false )
		% Distance( a.vPos.y(), false )
		% Distance( a.vPos.z(), false )
		% a.GetOrbBasis();
	return out;
}

// std::ostream& operator << ( std::ostream &out, FAtom const &a )
// {
// 	out << boost::format( "%+2s: (%+.5f,%+.5f,%+.5f) basis: '%s' ecp: '%s'" )
// 		% ElementNameFromNumber( a.AtomicNumber )
// 		% Distance( a.vPos.x(), false )
// 		% Distance( a.vPos.y(), false )
// 		% Distance( a.vPos.z(), false )
// 		% a.GetOrbBasis() % a.EcpDesc;
// 	return out;
// }

// std::ostream& operator << ( std::ostream &out, FAtomSet const &a )
// {
// 	out << "ATOMS: " << a.Atoms.size() << " atoms. (distance unit: "
// 		<< DistanceUnit() << ")\n";
// 	for( uint i = 0; i < a.Atoms.size(); ++ i ){
//
// 		out << boost::format("%+6d  ") % i << a.Atoms[i] << std::endl;
// 	}
// 	return out;
// }

std::ostream& operator << ( std::ostream &xout, FAtomSet const &a )
{
	xout << "   ATOM  ELM       POS/X      POS/Y      POS/Z        BASIS\n";
	for( uint i = 0; i < a.Atoms.size(); ++ i ){
		xout << boost::format("%+6d   ") % i << a.Atoms[i] << std::endl;
	}
	xout << "\n   (" << a.Atoms.size() << " atoms;"
	     << " " << a.NuclearCharge() - a.nCoreElectrons() << " valence electrons;"
	     << " distance unit: " << DistanceUnit() << ")\n";


	uint Ops[8], nOps,  Gen[3], nGen;
	a.FindMirrorSymmetryOps(Ops, nOps, Gen, nGen);

	xout << "\n";
	if ( 1 ) {
		// Generators 	Point group
		// (null card) 	$C_1$ (i.e. no point group symmetry)
		// X (or Y or Z) 	$C_s$
		// XY 	$C_2$
		// XYZ 	$C_i$
		// X,Y 	$C_{2v}$
		// XY,Z 	$C_{2h}$
		// XZ,YZ 	$D_2$
		// X,Y,Z 	$D_{2h}$
		//
		// 1 X Y Z XY XZ YZ XYZ
		char const
			*pPointGroupName = "?!";
		if ( nOps == 1 ) pPointGroupName = "C1";
		else if ( nOps == 8 ) pPointGroupName = "D2h";
		else if ( nOps == 2 && Ops[0] == 7 ) pPointGroupName = "Ci";
		else if ( nOps == 2 && (Ops[0] == 1 || Ops[0] == 2 || Ops[0] == 4 ) ) pPointGroupName = "Cs";
		else if ( nOps == 2 && (Ops[0] == 3 || Ops[0] == 6 || Ops[0] == 5 ) ) pPointGroupName = "C2";
		else {
			uint
				n1 = 0,
				n2 = 0;
			for ( uint i = 0; i < nOps; ++ i ) {
				uint
					nBits = ((Ops[i] & 1)>>0) + ((Ops[i] & 2)>>1) + ((Ops[i] & 4)>>2);
				if (nBits == 1) n1 += 1;
				if (nBits == 2) n2 += 1;
			}
			if ( n1 == 2 && n2 == 1 ) pPointGroupName = "C2v"; // id X Y XY
			if ( n1 == 1 && n2 == 1 ) pPointGroupName = "C2h"; // id XY Z XYZ
			if ( n1 == 0 && n2 == 3 ) pPointGroupName = "D2"; // id XZ YZ XY
		}

		if ( Verbosity >= 0 ) {
			_xout0( format("%31s%s") % "POINT GROUP:" % pPointGroupName );
			xout << format("%31s") % "SYMMETRY OPS:";
			for ( uint iOp = 0; iOp < nOps; ++ iOp ) {
				if ( iOp != 0 )
					xout << "  ";
				for ( uint i = 0; i < 3; ++ i )
					if ( (Ops[iOp] & (1 << i)) != 0 )
						xout << "XYZ"[i];
				for ( uint i = 0; i < nGen; ++ i )
					if ( iOp == Gen[i] )
						xout << "*";
				if ( Ops[iOp] == 0 )
					xout << "id";
			}
			xout << std::endl;
		}
	}

	return xout;
}










void FAtomSet::Add1eIntegralMatrix( FMatrixView &Dest,
		FBasisSet const &RowBasis, FBasisSet const &ColBasis,
		FDoublettIntegralFactory &IntFactory, double Factor, FMemoryStack &Mem ) const
{
	bool
		MatrixSymmetric = (&RowBasis == &ColBasis);

	assert( Dest.nRows == RowBasis.nBasisFn() );
	assert( Dest.nCols == ColBasis.nBasisFn() );

	FShellDoublettIntegral
		IntResult;
	uint
		StartFnA = 0, StartFnB; // starting indices of the functions of the
								// respective shells.
	for ( uint nShellA = 0; nShellA < RowBasis.Shells.size(); ++ nShellA ){
		FGaussShell const
			*pShellA = &RowBasis.Shells[nShellA];
		StartFnB = (!MatrixSymmetric)? 0 : StartFnA;
		for ( uint nShellB = (!MatrixSymmetric)? 0 : nShellA;
			nShellB < ColBasis.Shells.size(); ++ nShellB )
		{
			FGaussShell const
				*pShellB = &ColBasis.Shells[nShellB];
			IntFactory.EvalDoublett( IntResult, pShellA, pShellB, Mem );

			assert( IntResult.nSizeA == pShellA->nFuncs() );
			assert( IntResult.nSizeB == pShellB->nFuncs() );
			assert( StartFnA + IntResult.nSizeA <= RowBasis.nBasisFn() );
			assert( StartFnB + IntResult.nSizeB <= ColBasis.nBasisFn() );

			// fill data we gathered into the matrices
			for ( uint i_ = 0; i_ < IntResult.nSizeA; ++ i_ )
				for ( uint j_ = 0; j_ < IntResult.nSizeB; ++ j_ )
				{
					uint
						i = StartFnA + i_,
						j = StartFnB + j_;
					FScalar const
						&r = IntResult.Int( i_, j_ );
					Dest(i,j) += Factor * r;
					if ( MatrixSymmetric && nShellA != nShellB )
						Dest(j,i) += Factor * r;
				}

			StartFnB += pShellB->nFuncs();
		}
		StartFnA += pShellA->nFuncs();
	}
	assert( StartFnA == RowBasis.nBasisFn() );
}

void FAtomSet::MakeCoreHamiltonMatrix( FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem ) const
{
	Out.Clear();
	AddKineticMatrix(Out, RowBasis, ColBasis, Mem);
	AddNuclearAttractionMatrix(Out, RowBasis, ColBasis, Mem);

	// TODO: apply other CoreH patches from Environment object etc.

// 	AddDipoleMatrix( Out, RowBasis, ColBasis, FVector3(0.0,0.0,1e-5), FVector3(0,0,0), Mem );
};

void FAtomSet::MakeOverlapMatrix( FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem ) const
{
	FDoublettIntegralFactoryOverlap
		IntFactory; // can't bind temporary to non-const reference.
	Out.Clear();
	Add1eIntegralMatrix( Out, RowBasis, ColBasis, IntFactory, 1.0, Mem );
};

void FAtomSet::AddKineticMatrix( FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem ) const
{
	FDoublettIntegralFactoryKineticTerm
		IntFactory;
	Out.Clear();
	Add1eIntegralMatrix( Out, RowBasis, ColBasis, IntFactory, 1.0, Mem );
}

void FAtomSet::AddNuclearAttractionMatrix( FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis, FMemoryStack &Mem ) const
{
	typedef FDoublettIntegralFactoryFieldTerms
		FNucFactory;
	FNucFactory::FPointChargeList
		PointCharges;
	FAtomSet::FAtomList::const_iterator
		itAtom;
	// make one point charge for every atom. charge: ecp-reduced core charge.
	_for_each( itAtom, Atoms )
		PointCharges.push_back( FNucFactory::FPointCharge(
				itAtom->vPos, static_cast<FScalar>( itAtom->Charge ) ) );
	FNucFactory
		NuclearAttractionIntegrator( TVector3<uint>( 0, 0, 0 ), PointCharges );
	Add1eIntegralMatrix( Out, RowBasis, ColBasis, NuclearAttractionIntegrator, 1.0, Mem );
}


void FAtomSet::AddDipoleMatrix( FMatrixView &Out, FBasisSet const &RowBasis, FBasisSet const &ColBasis,
	FVector3 const &Direction, FVector3 const &ExpansionPoint, FMemoryStack &Mem ) const
{
/*	FVector3
		DirNormalized = (-1.0/std::sqrt(Dot(Direction,Direction))) * Direction;*/
		// ^- -1.0: Electrons are very negatively charged. this is actually
		// the charge factor with the elementary charge (which is 1 in a.u.)
	//xout << boost::format( "Dipole Direction: %f %f %f" )
	//	% DirNormalized.x() % DirNormalized.y() % DirNormalized.z() << std::endl;

	//  make dipole matrices into the three cartesian directions.
	for ( uint CartComp = 0; CartComp < 3; ++ CartComp ){
		if ( fabs( Direction[CartComp] ) < 1e-15 )
			continue;

		FDoublettIntegralFactoryMultipoleMoment::FMoment
			CartMoment( 0, 0, 0 );
		// this specifies the derivative power in cart direction Cart.
		CartMoment[CartComp] = 1;

		FDoublettIntegralFactoryMultipoleMoment
			IntFactory( CartMoment, ExpansionPoint );
		// add dot product component of current direction times dipole matrix.
		Add1eIntegralMatrix( Out, RowBasis, ColBasis, IntFactory,
			Direction[CartComp], Mem );
		//xout << "DipoleMomentOp" << CartComp << ":"
		//	 << ResultDir[CartComp].View( 0, 7, 0, std::min( 40u, ResultDir[CartComp].SizeY() ) );
	}
};


bool IsEquivalent(FAtom const &A, FVector3 const &vA, FAtom const &B)
{
	if ( A.AtomicNumber != B.AtomicNumber )
		return false;
	if ( A.Charge != B.Charge )
		return false;
	if ( DistanceSq(vA, B.vPos) > sqr(SymmetryTolerance) )
		return false;
	if ( A.BasisDesc != B.BasisDesc )
		return false;
	if ( A.EcpDesc != B.EcpDesc )
		return false;
	return true;
}

FVector3 ApplySymOp( FVector3 const &In, uint Op ){
	FVector3
		MirrorPos = In;
	for ( uint iXyz = 0; iXyz < 3; ++ iXyz )
		if ( Op & (1<<iXyz) )
			MirrorPos[iXyz] = -MirrorPos[iXyz];
	return MirrorPos;
}


void FAtomSet::FindEquivalentAtoms(FAtomSymGroup *&pGroups, uint &nGroups, FMemoryStack &Mem) const
{
	Mem.Alloc(pGroups, Atoms.size()); // actual data will be less.
	bool
		*pDone;
	Mem.ClearAlloc(pDone, Atoms.size());

	uint Ops[8], nOps, Gen[3], nGen;
	FindMirrorSymmetryOps(Ops, nOps, Gen, nGen);


	FAtomSymGroup
		*pCurGroup = pGroups;
	for ( uint iAt0 = 0; iAt0 < Atoms.size(); ++ iAt0 ){
		if ( pDone[iAt0] )
			continue; // already seen as part of another atom symmetry group.

		pCurGroup->nEqiv = 0;

		// look up other equivalent atoms and add them to the group.
		for ( uint iOp = 0; iOp < nOps; ++ iOp ){
			FVector3
				MirrorPos = ApplySymOp(Atoms[iAt0].vPos, Ops[iOp]);
			for ( uint iAt1 = iAt0; iAt1 < Atoms.size(); ++ iAt1 )
				if ( !pDone[iAt1] &&
						IsEquivalent(Atoms[iAt0], MirrorPos, Atoms[iAt1]) ) {
					pCurGroup->iEqiv[pCurGroup->nEqiv] = iAt1;
					pCurGroup->SymOp[pCurGroup->nEqiv] = Ops[iOp];
					pCurGroup->nEqiv += 1;
					assert(pCurGroup->nEqiv <= 8);
					pDone[iAt1] = true;
				}
		}

		// find out directions in which current atoms lie
		// on the mirror plane.
		pCurGroup->Stabilizer = 0;
		for ( uint iXyz = 0; iXyz < 3; ++ iXyz ){
			if ( sqr(Atoms[iAt0].vPos[iXyz]) <= SymmetryTolerance )
				pCurGroup->Stabilizer |= (1 << iXyz);
		}

		++ pCurGroup;
	};
	nGroups = pCurGroup - pGroups;
	assert(nGroups <= Atoms.size());

	Mem.Free(pDone);
	// still on stack: pGroups[1..Atoms.size()].
};



void FAtomSet::FindMirrorSymmetryOps(uint Ops[8], uint &nOps, uint Gen[3], uint &nGen) const
{
	nOps = 0;

	if ( 1 ) {
		// This disables symmetry everywhere.
		nOps = 1;
		nGen = 0;
		Ops[0] = 0;
		nOps = 1;
		return;
	}

	uint
		AllOps[8] = { 0, 1, 2, 4, 3, 5, 6, 7 }; // simple ones first: id X Y Z XY XZ YZ XYZ

	for ( uint iOp_ = 0; iOp_ < 8; ++ iOp_ ){
		uint
			iOp = AllOps[iOp_];
		bool
			Keep = true;

		for ( uint iAt0 = 0; iAt0 < Atoms.size(); ++ iAt0 ){
			bool
				FoundEquiv = false;
			FVector3
				MirrorPos = ApplySymOp(Atoms[iAt0].vPos, iOp);

			for ( uint iAt1 = 0; iAt1 < Atoms.size(); ++ iAt1 )
				if ( IsEquivalent(Atoms[iAt0], MirrorPos, Atoms[iAt1]) ) {
					FoundEquiv = true;
					break;
				}

			// no atom equivalent under symmetry xform found for iAt0.
			if ( !FoundEquiv ) {
				Keep = false;
				break;
			}
		}

		if ( Keep ) {
			Ops[nOps] = iOp;
			nOps += 1;
		}
	};

	nGen = 0;
	for ( uint i = 0; i < 3; ++ i )
		Gen[i] = 0;

	for ( uint iBaseOp = 0; iBaseOp < nOps; ++ iBaseOp ){
		uint
			BaseOp = Ops[iBaseOp];
		bool
			Redundant = false || BaseOp == 0;
		// Can the current operator can be formed as some combination
		// of the previous generators?...
		for ( uint iPat = 0; iPat < (1u<<nGen); ++ iPat ){
			uint
				TrialOp = 0;
			for ( uint iGen = 0; iGen < nGen; ++ iGen )
				if ( iPat & (1u << iGen) ) {
					// generator iGen activated in trial pattern.
					TrialOp ^= Gen[iGen];
				}
			if ( TrialOp == BaseOp ) {
				Redundant = true;
				break;
			}
		}
		if ( Redundant == true )
			continue;

		// ...nope. Include it as new generator.
		Gen[nGen] = BaseOp;
		nGen += 1;
	}
};

std::string FAtom::GetElementName() const {
	return ElementNameFromNumber(AtomicNumber);
}


} // namespace ct


// kate: space-indent off; tab-indent on; indent-width 4; mixedindent off; indent-mode normal;

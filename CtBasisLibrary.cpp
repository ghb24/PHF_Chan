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

// This file is released under the GNU General Public License ("GPL", version 2)
// as part of the CT8K program. Copying, modification, creating derivative works
// and redistribution of the file is allowed, but _only_ subject to the terms
// of that GPL. You should have received a version of this license along
// with this source code. The program comes "as is", without any kind of
// warranty.
//
// Authors/Copyright holders:
//  - Gerald Knizia, 2006 (tag: cgk, contact: cgk.d@gmx.net)


// #include <iostream> // FIXME: REMOVE THIS

#include <fstream>
#include <stdexcept>
#include <memory>
#include <set>
#include <sstream>
#include <boost/format.hpp>
using boost::str;
using boost::format;

#include "CtCommon.h"
#include "CtBasisLibrary.h"
#include "CtIo.h"
#include "CtAtomSet.h" // for ElementNumberFromName

namespace ct {

FBasisSetLibrary
	g_BasisSetLibrary;


FBasisSetLibrary::FBasisKey FBasisSetLibrary::MakeKey(std::string const &Name, uint iElement, bool AssertExists) const
{
	FBasisNameSet::const_iterator
		itName = m_BasisNames.find(Name);
	FBasisKey
		r;
	if ( itName == m_BasisNames.end() ) {
		if ( AssertExists ) {
			xerr << "Recognized names:\n";
			_for_each( itName, m_BasisNames )
				xerr << "   " << *itName << "\n";
			xerr << "\n";
			throw std::runtime_error("FBasisSetLibrary: Basis/ECP entry '" + Name + "' not found.");
		}
		r.pName = 0;
		r.iElement = 0;
		return r;
	}
	r.pName = &*itName;
	r.iElement = iElement;
	return r;
};


void FBasisSetLibrary::ImportMolproLib( std::string const &FileName )
{
	// this is a kind of ugly file formatm but MOLPRO seems to have nearly
	// every kind of basis which was ever conceived useful for anything
	// distributed with it, which makes it attractive. We try to parse
	// it here and maybe later export it as XML or something.
	// Also the EMSL gaussian basis set library can output libmol files.
	// see http://www.emsl.pnl.gov/forms/basisform.html.
	std::ifstream
		File( FileName.c_str() );
	if ( false == File.good() )
		throw std::runtime_error( "FBasisSetLibrary: Failed to open file '"+FileName+"' for basis library import." );

	// how this files look:
	//	Note: Lines with * are comments.
		// He s STO-3G STO3G : 3 1 1.3
		// STO-3G
		// 6.3624214 1.158923 0.31364979 0.15432897 0.53532814 0.44463454
		// Li s STO-3G STO3G : 6 2 1.3 4.6
		// STO-3G
		// 16.119575 2.9362007 0.7946505 0.6362897 0.1478601 0.0480887 0.15432897
		// 0.53532814 0.44463454 -0.09996723 0.39951283 0.70011547
		// Li p STO-3G STO3G : 3 1 1.3
		// STO-3G
		// 0.6362897 0.1478601 0.0480887 0.15591627 0.60768372 0.39195739
		// B s cc-pVDZ VDZ : 9 3 1.9 1.9 9.9
		// cc-pVDZ
		// 4570 685.9 156.5 44.47 14.48 5.131 1.898 0.3329 0.1043 0.000696
		// 0.005353 0.027134 0.10138 0.272055 0.448403 0.290123 0.014322 -0.003486
		// -0.000139 -0.001097 -0.005444 -0.021916 -0.059751 -0.138732 -0.131482
		// 0.539526 0.580774 1
	// interpretation: First Line:
	//		ElementName OrbitalType N*[AlternativeNameI] :
	//				#PrimOrbitals #OrbitalsTheyAreContractedTo N*[ContractFrom.ContractTo]
	//		Any string, usually reference to be cited for the basis set
	//		M*[PrimitiveOrbitalExponents] L*[PrimitiveOrbitalCoefficients]
	//	for each contraction contraction length coefficients are supplied.
	// and data sets like this are also seen:
		// H  P 6-311G2P :   2  0
		// G92 6-311G polarization
		//    .15000000D+01   .37500000D+00
		// H  P 6-311G3P :   3  0
		// G92 6-311G polarization
		//    .30000000D+01   .75000000D+00   .18750000D+00
		// H  S 6-31G** 6-31G* 6-31G :   4  1      1.3
		// G92 6-31G**
 		// .18731137D+02   .28253944D+01   .64012169D+00   .16127776D+00   .33494604D-01
		// .23472695D+00   .81375733D+00
	// all/some primitive orbitals are probably left uncontracted. Therefore no
	// contraction coefficients are supplied for these (only exponents).

	// Also I originally tought that the names of the basis set after the first
	// one were alternative Names of the set (like cc-pVDZ = VDZ). This is
	// however not how it works, at least not in all cases: It seems that names
	// which are supplied actually mean that the currently described elements
	// have to be inserted into the basis sets of ALL names provided, which
	// may or may not be identical basis sets (for example, for the s and p
	// orbitals of the 931G basis set, also the starred names are listed, which
	// means that these basis sets share these orbitals, altough in general they
	// are different).

	// so.. let's begin the mess.
	// read the entire file into an stringstream object.
	std::auto_ptr<char>
		pFileContent;
	std::size_t
		FileLength;
	File.seekg( 0, std::ios::end );
	FileLength = File.tellg();
	pFileContent.reset( static_cast<char*>( calloc( 1, 2 + FileLength ) ) );
	File.seekg( 0, std::ios::beg );
	File.read( pFileContent.get(), FileLength );
	// ^- not quite sure how/if this is terminated by read, therefore
	// calloc instead of malloc.  pFileContent[FileLength] = 0; <- might
	// not work as excpected if e.g. 0x0A0D are converted to \n or something

	// okay, now this is very bad. This FORTRAN stuff denotes exponents
	// in scientific notation not as "23423e-3" but as "23423D-03" and
	// similiar nonsense. We do a lame attempt to convert that. This might break
	// in some situations. // FIXME: correct this somehow.
	for ( uint i = 0; i < FileLength - 2; ++i ){
		if ( (pFileContent.get()[i]=='D') &&
			 ((pFileContent.get()[i+1]=='+')||(pFileContent.get()[i+1]=='-')) &&
			 (pFileContent.get()[i+2]=='0') )
			pFileContent.get()[i] = 'E';
	}

	FBasisNameSet
		AllBasisNames;

	std::stringstream
		str( pFileContent.get(), std::stringstream::in );
		// NOTE: "in" means from stream, into local data, not
		// into stream (file naming convention)
	pFileContent.reset(0);

	try{
		while( str.good() )
		{
			// clear exception mask, now stream will not throw() when something
			// unexpected happens.
			str.exceptions(std::ios::goodbit);
			std::string
				s;
			std::getline( str, s );
			if ( ( s.size() == 0 ) || ( s[0] == '*' ) || ( s[0] == '!' ) )
				// empty or comment line, throw it away and go on with the next
				continue;
			if ( !str.good() ) // eof, bad etc.
				break;
			str.exceptions( std::ios::failbit );
			// ^- when something fails, throw an exception. this will happen
			// if the actual file format does not match the one I had in mind
			// when coding this.

			// expected format: ElementName[w]OrbitalType[w]AlternativeNames[w] :
			// using namespace std;
			// cout << "** read line: '" << s << "'" << endl;
			std::stringstream
				line( s, std::stringstream::in );
			line.exceptions( std::ios::badbit | std::ios::failbit );
			std::string
				Element, Type;
			line >> Element >> Type;
			if ( Type.size() != 1 )
				throw std::runtime_error( "Parsing error, cannot interpret orbital type '" + Type + "'" );


			int
				AngMom;
			aic::FGaussBfn::FContractionList
				Contractions;
			aic::FScalarArray
				Exponents;

			char
				cAngMom = ::tolower(Type[0]);
			for ( AngMom = 0; AngMom < 9; ++ AngMom )
				if ( cAngMom == "spdfghikl"[AngMom] )
					break;
			if ( AngMom == 9 )
				throw std::runtime_error((format("Failed to understand angular momentum '%s'.") % cAngMom).str());


			// cout << "Element " << Element << " Type " << Type << endl;

			FStringList
				BasisNames; // all names of basis sets in which the
							// current entry is to be inserted.
			for ( line >> s; s != ":"; line >> s ) {
				// cout << "Alternative Name:" << s << std::endl;
				BasisNames.push_back( tolower(stripwhitespace(s)) );
				AllBasisNames.insert(stripwhitespace(s));
			}

			// expected format: #prim orbitals #contractions (#contr.)*[a.b]
			// denoting indices of begin and end of a contraction with the
			// following exponents/contraction coefficients.
			int
				nPrimOrbitals,
				nContractions,
				nCoefficients(0),
				nHighestPrimOrbitalInContraction(0); // 1-based index.
			line >> nPrimOrbitals >> nContractions;
			// cout << "#Prim " << nPrimOrbitals << " #Contr " << nContractions << endl;
			Contractions.reserve(nContractions);
			for ( int i = 0; i < nContractions; ++ i ){
				Contractions.push_back(aic::FGaussBfn::FContraction());
				char Dot;
				int nContrBegin, nContrEnd; // [begin,end]
				line >> nContrBegin >> Dot >> nContrEnd;
				if( Dot != '.' )
					throw std::runtime_error("GTO-Contraction read format error.");
				// cout << "  Contr: #" << nContrBegin << "-#" << nContrEnd << endl;
				Contractions.back().nBegin = nContrBegin - 1; // idx 0-based.
				Contractions.back().nEnd = nContrEnd; // end is excl.
				if( ( nContrBegin > nContrEnd ) || ( nContrEnd > nPrimOrbitals ) )
					throw std::runtime_error("GTO-Contraction logical error.");
				nCoefficients += nContrEnd - nContrBegin + 1;
				nHighestPrimOrbitalInContraction = std::max( nHighestPrimOrbitalInContraction, nContrEnd );
			}

			std::string
				EntryComment;
			do{ // read name, maybe skip comments.
				getline( str, EntryComment );
			} while ( ( EntryComment.size() != 0 ) && ( EntryComment[0] == '*' ) );
			// cout << "Entry Comment: " << EntryComment << endl;


			Exponents.resize(nPrimOrbitals);
			aic::FScalarArray
				Coefficients(nCoefficients, 0);
			// now read exponents and contraction coefficients;
			// (this will break if comments are present inbetween)
			for ( int i = 0; i < nPrimOrbitals; ++ i ){
				double Exponent;
				str >> Exponent;
				// cout << "Exponent: " << Exponent << endl;
				Exponents[i] = Exponent;
			}

			// read in coefficients.
			if ( nContractions != 0 )
				for ( int i = 0; i < nCoefficients; ++ i ){
					double Coefficient;
					str >> Coefficient;
					// cout << "Coefficient: " << Coefficient << endl;
					Coefficients[i] = Coefficient;
				}
			// copy over the contraction coefficients to the contractions.
			int
				iCoefficient = 0;
			for ( int i = 0; i < nContractions; ++ i ){
				int iCoefficientEnd = iCoefficient + (Contractions[i].nEnd - Contractions[i].nBegin);
				Contractions[i].Coeffs.assign( Coefficients.begin() + iCoefficient,
					Coefficients.begin() + iCoefficientEnd );
				iCoefficient = iCoefficientEnd;
			}


			// in some files some primitive orbitals are left uncontracted.
			// but these are not stored as 1-GTO contractions but the
			// coefficients for these are just not present in the file.
			int
				nAdditionalContr = nPrimOrbitals - nHighestPrimOrbitalInContraction;
			if ( 0 != nAdditionalContr )
			{	// generate 1-GTO-each contractions manually.
				int
					nStart = nHighestPrimOrbitalInContraction;
					// ^- 0 based index, the rhs one is 1-based.
				for ( int i = 0; i < nAdditionalContr; ++ i ){
					Contractions.push_back( aic::FGaussBfn::FContraction() );
					Contractions.back().nBegin = i + nStart;
					Contractions.back().nEnd = i + nStart + 1;
					Contractions.back().Coeffs.push_back(1.0);
				}
			}

			// import all names of the basis function
			FStringList::const_iterator
				itName;
			_for_each(itName, BasisNames)
				m_BasisNames.insert(*itName);

			// make the actual basis function and link it to all the names.
			FGaussBfnPtr
				pBfn( new FGaussBfn(AngMom, FGaussBfn::TYPE_Spherical, Exponents, Contractions) );
			int
				iElement = ElementNumberFromName(Element);
			_for_each(itName, BasisNames)
				m_BasisFns.insert( FBasisFnMap::value_type(MakeKey(*itName, iElement), pBfn) );

			// chew the EOL marker if present, leave loop otherwise.
			str.exceptions(std::ios::goodbit);
			str.ignore(0xbad, '\n');
		};
	}
	catch( std::ios_base::failure &e){
		// this is not exactly something i would usually
		// call "error handling" but i really hate this string-
		// fiddling stuff and we can't really do anything better
		// about it anyway.
		xerr << "PARSER EXCEPTION:" << e.what() << std::endl;
		throw std::runtime_error( "Parsing of LIBMOL file FAILED because the actual syntax did not match the expected one. Last entry successfully processed: ");
	}
	catch( std::exception &e ){
		xerr << "Exception during LibmolFile parsing: " << e.what() << std::endl;
		throw;
	};


	if ( 1 ) {
		std::size_t
			iDirSep = FileName.rfind('/');
		if ( iDirSep == std::string::npos )
			iDirSep = 0;
		else
			iDirSep += 1;
		xout << format("LOADED %25s") % FileName.substr(iDirSep);
		if ( 1 ) {
			xout << "[";
			FBasisNameSet::const_iterator
				itSet;
			uint
				nLen = 0;
			_for_each(itSet, AllBasisNames) {
				if ( nLen >= 40 ) {
					xout << ",...";
					break;
				}

				if ( itSet != AllBasisNames.begin() )
					xout << ", ";
				xout << *itSet;
				nLen += itSet->size();
			}
			xout << "]";
		}
		xout << std::endl;
	}
};


// converts something like "sto3g>>spd" into "sto3g" and ['s','p','d'].
// First two parameters are Out parameters.
void ParseBasisDesc( std::string &BasisName, std::set<char> &OrbitalTypes,
	FBasisDesc const &BasisDesc )
{
	OrbitalTypes.clear();
	std::size_t
		nSeparatorIdx = BasisDesc.find(">>");
	// separator present in string?
	if ( nSeparatorIdx == std::string::npos ){
		// no->primitive declaration, "import all orbitals from [...]".
		BasisName = tolower(stripwhitespace(BasisDesc));
	} else { // yes. find out which orbitals to import.
		BasisName = tolower(stripwhitespace(BasisDesc.substr(0,nSeparatorIdx)));
		for ( std::size_t i = nSeparatorIdx+2; i < BasisDesc.size(); ++i ){
			char c = ::tolower( BasisDesc[i] );
			if ( ( c >= 'a' ) && ( c <= 'z' ) )
				// ^- we should consider using more b- and q-orbitals.
				OrbitalTypes.insert(c);
		}
	}
};

bool starts_with(std::string const &s, char const *p)
{
	std::size_t
		N = s.size();
	for ( std::size_t i = 0; i < N && *p && s[i] == *p; ++ i, ++p ){
	}
	return *p == 0;
};

int AngMomFromChar(char c) {
	switch (c) {
		case 's': return 0;
		case 'p': return 1;
		case 'd': return 2;
		case 'f': return 3;
		case 'g': return 4;
		case 'h': return 5;
		case 'i': return 6;
		case 'k': return 7;
		default: {
			std::stringstream str;
			str << "Angular momentum not recognized: '" << c << "'. Should be one of spdfghik.";
			throw std::runtime_error(str.str());
		}
	};
}

void FBasisSetLibrary::TryMakeBasis(std::string const &BasisDesc, int iElement)
{
	if ( starts_with(BasisDesc, ".et[") ) {
		// syntax: .ET[spdf,<center>, <ratio>, <powmin>, <powmax>]
		char
			AngMoms[20];
		double
			fCenter, fRatio;
		int
			iPowMin, iPowMax, iScan;
		iScan = std::sscanf(BasisDesc.c_str(), ".et[%10[spdfghik],%lf,%lf,%i,%i]",
			&AngMoms[0], &fCenter, &fRatio, &iPowMax, &iPowMin);
		if ( !((iScan == 4 && iPowMin > 0) || iScan == 5 )) {
			std::stringstream str;
			str << "Failed to understand even tempered basis declaration '" << BasisDesc << "'."
				<< " Expected: .et[<spd...>, <center>, <ratio>, <powmax>(, <powmin>)].";
			throw std::runtime_error(str.str());
		}
		if ( iScan == 4 )
			iPowMin = iPowMax;
		iPowMin *= -1;

		// make the basis function and store it.
		aic::FScalarArray
			Exps;
		Exps.reserve(iPowMax - iPowMin + 1);
		for ( int i = iPowMax; i >= iPowMin; -- i )
			Exps.push_back( fCenter * std::pow(fRatio, (double)i) );

		m_BasisNames.insert(BasisDesc);

		for ( char *pAm = &AngMoms[0]; *pAm != 0; ++ pAm ) {
			FGaussBfnPtr
				pBfn( new FGaussBfn(AngMomFromChar(*pAm), FGaussBfn::TYPE_Spherical, Exps) );
			m_BasisFns.insert( FBasisFnMap::value_type(MakeKey(BasisDesc, iElement), pBfn) );
		}
	} else if ( starts_with(BasisDesc, "aug-" ) ||
				starts_with(BasisDesc, "daug-" ) ||
				starts_with(BasisDesc, "taug-" ) ) {
// 		xout << "  in .aug generation code." << std::endl;
		uint
			iAugLevel = 1;
		if ( starts_with(BasisDesc, "daug-" ) )
			iAugLevel = 2;
		if ( starts_with(BasisDesc, "taug-" ) )
			iAugLevel = 3;
		// look up the parent basis.
		typedef FBasisFnMap::const_iterator
			FBfnIt;
		std::string
			ParentDesc(BasisDesc.substr((iAugLevel==1)? 4 : 5)); // includes the '-'.
		std::pair<FBfnIt,FBfnIt>
			itBfs = m_BasisFns.equal_range(MakeKey(ParentDesc, iElement));
		if ( itBfs.first == m_BasisFns.end() )
			// not there.
			return;

		m_BasisNames.insert(BasisDesc);

		uint
			AmMask = 0; // bitmask of angular momenta already handled.
		for ( FBfnIt it = itBfs.first; it != itBfs.second; ++ it ){
			// note: this code assumes that all functions within a angular
			// momentum are stored inside one single contracted shell!
			FGaussBfnCptr
				pFn = it->second;
			if ( (AmMask & (1 << pFn->AngularMomentum)) != 0 ) {
				std::stringstream str;
				str << "sorry, diffuse augmentation code got confused. Encountered multiple shells with angular momentum l=" << pFn->AngularMomentum << ".";
				throw std::runtime_error(str.str());
			}
			AmMask |= 1 << pFn->AngularMomentum;

			aic::FScalarArray
				Exps = pFn->Exponents;
			FGaussBfn::FContractionList
				Cos = pFn->Contractions;
			if ( Exps.empty() )
				throw std::runtime_error("Encountered a basis function without exponents during augmentation.");
			double
				fCen = Exps.back(),
				fRatio = 1./2.5;
			if ( Exps.size() > 1 ) {
				fRatio = Exps[Exps.size()-1]/Exps[Exps.size()-2];
			}
			for ( uint iAug = 1; iAug <= iAugLevel; ++ iAug ) {
				Exps.push_back(fCen * std::pow(fRatio, (double)iAug));
				aic::FScalarArray
					Coeffs(1, GaussNormalizationSpher(Exps.back(), pFn->AngularMomentum));
				Cos.push_back( FGaussBfn::FContraction(Exps.size()-1, Exps.size(), Coeffs) );
			}

			FGaussBfnPtr
				pBfn( new FGaussBfn(pFn->AngularMomentum, FGaussBfn::TYPE_Spherical |
					FGaussBfn::TYPE_Unnormalized, Exps, Cos) );
			m_BasisFns.insert( FBasisFnMap::value_type(MakeKey(BasisDesc, iElement), pBfn) );
		}
	}
};


void FBasisSetLibrary::LoadBasisFunctions( std::vector<aic::FGaussShell> &Shells,
	int iElement, std::string const &BasisDesc_,
	ct::FVector3 const &vAtomPos, FAtomIndex iAtomIdx ) const
{
	// todo: tokenize BasisDesc at ';'s here (and strip whitespace).
	// do the same this as he have below for the individual parts.

	std::string BasisDesc = BasisDesc_;
	if ( BasisDesc == "cc-pvtz-ig" ) BasisDesc = "cc-pvtz";
	if ( BasisDesc == "def2-tzvp-ig" ) BasisDesc = "def2-tzvp";
	if ( BasisDesc == "def2-sv(p)-jfit-ig" ) BasisDesc = "def2-sv(p)-jfit";

	typedef FBasisFnMap::const_iterator
		FBfnIt;
	FBasisKey const
		&Key = MakeKey(BasisDesc, iElement, false);
	std::pair<FBfnIt,FBfnIt>
		itBfs = m_BasisFns.equal_range(Key);
	if ( Key.pName == 0 || itBfs.first == itBfs.second ) {
		// basis not yet officially registered, but maybe it is one we
		// can make by modifying another one (e.g., aug-cc-pVTZ from cc-pVTZ)
		const_cast<FBasisSetLibrary*>(this)->TryMakeBasis(BasisDesc, iElement);
		itBfs = m_BasisFns.equal_range(MakeKey(BasisDesc, iElement));
	}

	if ( itBfs.first == m_BasisFns.end() ) {
		std::stringstream str;
		str << "Requested basis '" << BasisDesc << "' for element '"
		    << ElementNameFromNumber(iElement) << "' was not found.";
		throw std::runtime_error(str.str());
	}

	// fixme: remove this.
#if 1

	if ( BasisDesc_ == "cc-pvtz-ig" ) {
		for ( FBfnIt it = itBfs.first; it != itBfs.second; ++ it ){
			// extract only the functions of the basis which are contracted.
			// No polarization functions, no higher AM.
			FGaussBfnCptr
				pFn = it->second;

			FGaussBfn::FContractionList
				CoFns;
			for ( uint i = 0; i < pFn->Contractions.size(); ++ i )
				if ( pFn->Contractions[i].Coeffs.size() != 1 )
					CoFns.push_back( pFn->Contractions[i] );
			if ( !CoFns.empty() ) {
				FGaussBfnCptr
					pFn2( new FGaussBfn(pFn->AngularMomentum, pFn->Flags | FGaussBfn::TYPE_Unnormalized, pFn->Exponents, CoFns ) );
				Shells.push_back(FGaussShell(iAtomIdx, vAtomPos, pFn2));
			}
			// add sv(p) polarization for higher AM...
			if ( CoFns.empty() ) {
				FBasisKey const
					&Key = MakeKey("def2-sv(p)", iElement, false);
				std::pair<FBfnIt,FBfnIt>
					itBfs = m_BasisFns.equal_range(Key);
				for ( FBfnIt it = itBfs.first; it != itBfs.second; ++ it )
					if ( it->second->AngularMomentum == pFn->AngularMomentum )
						Shells.push_back(FGaussShell(iAtomIdx, vAtomPos, it->second));
			}
		}
	} else {
		int
			MaxL = 99,
			ReduceL = 0;
		if ( BasisDesc_ == "def2-sv(p)-jfit-ig" ) ReduceL = 2;
		if ( BasisDesc_ == "def2-tzvp-ig" ) ReduceL = 1;

		if ( ReduceL != 0 ) {
			MaxL = 0;
			for ( FBfnIt it = itBfs.first; it != itBfs.second; ++ it )
				MaxL = std::max(MaxL, it->second->AngularMomentum);
			MaxL -= ReduceL;
		}

		// note that we do NOT clear Shells but rather append only.
		for ( FBfnIt it = itBfs.first; it != itBfs.second; ++ it ){
			FGaussBfnCptr
				pFn = it->second;
			if ( pFn->AngularMomentum > MaxL )
				continue;
			Shells.push_back(FGaussShell(iAtomIdx, vAtomPos, pFn));
		}
	}

#else
	// note that we do NOT clear Shells but rather append only.
	for ( FBfnIt it = itBfs.first; it != itBfs.second; ++ it ){
		FGaussBfnCptr
			pFn = it->second;
		Shells.push_back(FGaussShell(iAtomIdx, vAtomPos, pFn));
	}
#endif
};







// void FBasisSetLibrary::LoadBasisFunctions( std::vector<aic::FGaussShell> &Shells,
// 	int iElement, std::string const &BasisDesc,
// 	ct::FVector3 const &vAtomPos, FAtomIndex iAtomIdx ) const
// {
// 	typedef FBasisFnMap::const_iterator
// 		FBfnIt;
//
// 	// note that we do NOT clear Shells but rather append only.
// 	std::pair<FBfnIt,FBfnIt>
// 		itBfs = m_BasisFns.equal_range(MakeKey(BasisDesc, iElement));
// 	if ( itBfs.first == m_BasisFns.end() ) {
// 		std::stringstream str;
// 		str << "Requested basis '" << BasisDesc << "' for element '"
// 		    << ElementNameFromNumber(iElement) << "' was not found.";
// 		throw std::runtime_error(str.str());
// 	}
//
// 	for ( FBfnIt it = itBfs.first; it != itBfs.second; ++ it ){
// 		FGaussBfnCptr
// 			pFn = it->second;
// 		Shells.push_back(FGaussShell(iAtomIdx, vAtomPos, pFn));
// 	}
// };













// // syntax: "sto3g" for "all orbitals from there" or "sto3g|spd" for limited basis sets.
// // further: "sto3g>>spd|cc-vdz>>f" for multiple orbitals from multiple basis sets?
// // we leave that open for the future.
// bool FBasisSetLibrary::LoadBasisFunctions( FBasisSet::FGaussShellArray &Shells,
// 		bool MakeSolidHarmonicGtos, uint nElement, FBasisDesc const &BasisDesc,
// 		FVector3 const &vAtomPos, FAtomIndex nAtomIdx )
// {
// 	// note that we do NOT clear Shells but rather append only.
//
// 	std::set<char>
// 		RequestedOrbitalTypes, ActualOrbitalTypes;
// 	std::string
// 		BasisLookupName;
// 	ParseBasisDesc( BasisLookupName, RequestedOrbitalTypes, BasisDesc );
//
// 	// get a ref to the basis set entry for the element we are to load.
// 	FBasisSetLookupMap::iterator
// 		itSet = m_LookupMap.find( BasisLookupName );
// 	if ( itSet == m_LookupMap.end() ){
// 		xerr << "BasisSetLibrary: A basis matching '" << BasisLookupName << "' is not present in the library." << std::endl;
// 		return false;
// 	}
//
// 	FBasisSetEntry
// 		*pBasisSetEntry = itSet->second;
//
// 	FBasisSetEntry::FElementOrbitalListMap::iterator
// 		itElement = pBasisSetEntry->ElementOrbitals.find(nElement);
// 	if ( itElement == pBasisSetEntry->ElementOrbitals.end() ){
// 		xerr << "BasisSetLibrary: Element " << ElementNameFromNumber(nElement)
// 			 << " is not contained in the basis set '" << BasisDesc << "' as i've loaded it." << std::endl;
// 		return false;
// 	}
//
// 	// if no orbital types were supplied, load all orbitals contained in the
// 	// basis set for the given element into the shell array. Otherwisde, load
// 	// only the requested shell types.
// 	FBasisSetOrbitalList
// 		&BasisSetOrbitalList = itElement->second;
// 	FBasisSetOrbitalList::iterator
// 		itBasisOrbital;
// 	_for_each( itBasisOrbital, BasisSetOrbitalList )
// 	{
// 		//assert( itBasisOrbital->Coefficients.size() == itBasisOrbital->Exponents.size() );
// 		assert( itBasisOrbital->Contractions[0].size() == itBasisOrbital->Contractions[1].size() );
//
// 		// if orbital types to load from this set were supplied, check if
// 		// the current basis function belongs to the supplied types.
// 		if ( ( !RequestedOrbitalTypes.empty() ) &&
// 			 ( RequestedOrbitalTypes.find( itBasisOrbital->Type ) == RequestedOrbitalTypes.end() ) )
// 			// nop, not present
// 			continue;
//
// 		// each concrete element in the basis set library (e.g. 's' for some
// 		// basis set) may result several contracted gaussian shells, which in
// 		// turn may correspond to multiple basis functions (i.e. library p ->
// 		// fn p_x, fn p_y).
//
// 		FAngularMomentum
// 			AngularMomentum;
// 		try{
// 			AngularMomentum = AngularNumberNumberFromName( itBasisOrbital->Type );
// 		}
// 		catch( std::exception &e ){
// 			xerr << "BasisSetLibrary: FAILED to select orbital type '" << itBasisOrbital->Type
// 					<< "' because relevant program structure is not yet present. Ignoring corresponding basis function." << std::endl;
// 			continue;
// 		}
// 		ActualOrbitalTypes.insert( itBasisOrbital->Type );
//
// 		uint
// 			LastContractionCoefficientOffset = 0;
//
// 		// generate a shell structure for the current basis set entry. after
// 		// adjusting the code to general contractions, this is actually a
// 		// one to one correspondence.
// 		FGaussShell
// 			NewShell;
// 		FGaussShell::FContractionList
// 			Contractions;
//
// 		for ( int nContraction = 0; nContraction <
// 				itBasisOrbital->Contractions[0].size(); ++nContraction )
// 		{
// 			FGaussShell::FContraction
// 				Contraction;
// 			FScalarArray
// 				&Coefficients = Contraction.Coeffs;
//
// 			Contraction.nBegin = itBasisOrbital->Contractions[0][nContraction],
// 			Contraction.nEnd = itBasisOrbital->Contractions[1][nContraction];
// 				// ^- indices into the exponent array of the *itBasisOrb
//
//
// 			// sort out exponents and coefficients for the current contraction
// 			Coefficients.resize( Contraction.nEnd - Contraction.nBegin );
// 			for ( int nPrimOrbital = 0; nPrimOrbital < Contraction.nEnd - Contraction.nBegin; ++ nPrimOrbital )
// 				Coefficients[nPrimOrbital] = itBasisOrbital->Coefficients[
// 					LastContractionCoefficientOffset + nPrimOrbital];
// 			LastContractionCoefficientOffset += Contraction.nEnd - Contraction.nBegin;
//
// 			Contractions.push_back( Contraction );
// 		};
// 		// normalize the shell, generate additional data and append it to
// 		// the list.
// 		NewShell.Initialize( nAtomIdx, vAtomPos, AngularMomentum,
// 			MakeSolidHarmonicGtos, Contractions, itBasisOrbital->Exponents );
//
// 		Shells.push_back( NewShell );
// 	};
//
// 	if ( ( !RequestedOrbitalTypes.empty() ) &&
// 		 ( RequestedOrbitalTypes != ActualOrbitalTypes ) )
// 	{
// 		xerr << "BasisSetLibrary: WARNING: Of basis desc. '" << BasisDesc
// 			 <<"' for element '" << ElementNameFromNumber( nElement )
// 			 << "' only the following orbital types were present in the library: ";
// 		std::set<char>::iterator
// 			itType;
// 		_for_each( itType, ActualOrbitalTypes )
// 			xerr << *itType << " ";
// 		xerr << "\b. The other requested types were IGNORED." << std::endl;
// 	};
//
// 	return true;
// };


/*
void ImportTest()
{
	g_BasisSetLibrary.ImportMolproLib( "/home/knizia/Programs/n82/basis_lib/sto3g.libmol" );
	g_BasisSetLibrary.ImportMolproLib( "/home/knizia/Programs/n82/basis_lib/sto6g.libmol" );

	FBasisSet::FGaussShellArray
		TestBasis;
	g_BasisSetLibrary.LoadBasisFunctions( TestBasis, ElementNumberFromName("O"),
			"sto3g", FVector3(1.0,0.0,0.0), 0 );
	FBasisSet::FBasisFnArray::iterator
		itBasisFn;
	_for_each( itBasisFn, TestBasis )
	{
		xout << "Contracted Basis Fn: #Prims: " << itBasisFn->Gaussians.size() << "\n";
		FGaussianArray::iterator
			itGaussFn;
		_for_each( itGaussFn, itBasisFn->Gaussians )
			xout << "   " << *itGaussFn << "\n";
	};
};*/

// returns g_BasisSetLibrary. Used for python interface because global vars
// can't really be exported as it seems.
FBasisSetLibrary& GetBasisSetLibrary()
{
	return g_BasisSetLibrary;
}

} // namespace ct

// kate: space-indent off; tab-indent on; indent-width 4; mixedindent off; indent-mode normal;

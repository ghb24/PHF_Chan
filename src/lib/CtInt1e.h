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

#ifndef CT_INT1E_H
#define CT_INT1E_H

#include <vector>
#include "CtCommon.h"

namespace aic {
	struct FGaussShell;
}

namespace ct {

typedef double
	FScalar;
typedef std::vector<FScalar>
	FScalarArray;
typedef uint
	FAngularMomentum;

using aic::FGaussShell;


struct FAngularComp : public TVector3<signed char> // struct size: one uint32
{
	typedef TVector3<signed char>
		base_type;
	FAngularComp()
	{};

	FAngularComp( value_type nx, value_type ny, value_type nz )
		: TVector3<value_type>( nx, ny, nz )
	{};

	// the actual l-angular momentum number this momentum object belongs to
	FAngularMomentum AngMom() const {
		return this->x() + this->y() + this->z();
	};

#ifdef _DEBUG
	value_type operator [] ( uint idx ){
		return base_type::operator []( idx );
	}
	value_type operator [] ( uint idx ) const {
		return base_type::operator []( idx );
	}
#endif // _DEBUG
};

typedef std::vector<FAngularComp>
	FAngularCompList;
FAngularCompList AngularComponentsCart( FAngularMomentum const &l );







// ####################################################################
// ********************************************************************
//
//		Now starting: interface for 1e^- integrals (i.e. overlap,
//		multipoles, kinetic terms, nuclear attraction and similiar stuff)
//
// --------------------------------------------------------------------


struct FShellDoublettIntegral
{
	FScalarArray
		Data;
	uint
		// number of functions in the respective shells.
		nSizeA, nSizeB;
	// parameters: shell indices, i.e. indices of functions in their shells.
	// return value: <a|op|b>
	FScalar& Int( uint a, uint b ){
		return Data[ a + nSizeA * b ];
	}
};

struct FDoublettIntegralFactory
{
	// use: instantiate one of the derived classes to specify which integrals
	// are to be evaluated. Then do something like
	//		FDoublettIntegralFactoryOverlap IntFactory;
	//		FMemoryStack2 tmp1(2000000);
	//		FShellDoublettIntegral Res;
	//		IntFactory.EvalDoublett( Res, a, b, tmp1 );
	// Notes:
	//	Neither this class nor the derived classes have any significant
	//  member variables, therefore class instantiation (to change the integral
	//  type) is very cheap. To preserve this behavior and simplify writing
	//  thread safe code, the temporary arrays are allocated outside of the
	//  class.
	void EvalDoublett( FShellDoublettIntegral &Result,
		FGaussShell const *a, FGaussShell const *b,
		FMemoryStack &Mem );

	virtual ~FDoublettIntegralFactory();
protected:
	struct FPrimitiveIntegralInfo{
		FPrimitiveIntegralInfo( FAngularCompList const& AngCompsA_,
			FAngularCompList const &AngCompsB_, FMemoryStack &TempStorage_ )
			: AngCompsA( AngCompsA_ ), AngCompsB( AngCompsB_ ), Mem( TempStorage_ )
		{};

		FAngularMomentum
			AngMomA,
			AngMomB;
		FAngularCompList const
			&AngCompsA,
			&AngCompsB;
		FVector3
			vCenterA,
			vCenterB,
			vAtoB;
		FMemoryStack
			&Mem;

		FScalar
			ZetaA,
			ZetaB,
			ZetaSum,
			ZetaGeomMean;
		FVector3
			vWeightedCenter, vAtoWC, vBtoWC, vMtoWC;
			// ^- center of the overlap distribution of the two primitives
			// (from gaussian product theorem)

	};

	// sets contribution of the given primitive cartesian shell to integrals at
	// Result[AngCompA + AngCompB*NumAngCompA]. Implemented for concrete
	// integral types in derived classes.
	virtual void EvalPrimitiveIntegral( FScalar *Result, FPrimitiveIntegralInfo const &info ) = 0;

	// calculates zeta-dependent data in *info and calls EvalPrimitiveIntegral
	void CalcDataAndEvalPrimitiveIntegral( FScalar *Result, FScalar ZetaA, FScalar ZetaB,
		FPrimitiveIntegralInfo &info );
};


// this one can generate cartesian multipole moments, i.e. the expectation
// values
//		<a| (x-mx)^m0 (y-my)^m1 (z-mz)^m2 |b>.
// where (mx,my,mz) is the point around which the multipole expansion takes
// place and m0 .. m2 are the cartesian powers (total moment: m0+m1+m2).
// note: actual evaluation function is FDoublettIntegralFactory::EvalDoublett()
struct FDoublettIntegralFactoryMultipoleMoment : public FDoublettIntegralFactory
{
	typedef TVector3<uint>
		FMoment;
	FDoublettIntegralFactoryMultipoleMoment( FMoment CartesianMoment, FVector3 MomentCenter )
		: vMoment( CartesianMoment ), vMomentCenter( MomentCenter )
	{};

protected:
	void EvalPrimitiveIntegral( FScalar *Result, FPrimitiveIntegralInfo const &info ); // override
private:
	FMoment
		vMoment;
	FVector3
		vMomentCenter;
};

// makes overlap integrals: <a|b>
//
// note: actual evaluation function is FDoublettIntegralFactory::EvalDoublett()
struct FDoublettIntegralFactoryOverlap : public FDoublettIntegralFactoryMultipoleMoment
{
	FDoublettIntegralFactoryOverlap()
		: FDoublettIntegralFactoryMultipoleMoment( FMoment(0,0,0), FVector3(0.0,0.0,0.0) )
	{};
};



// makes non-relativistic kinetic energy integrals:
//		<a|T|b> = <a| -1/2 (\del_x^2 + \del_y^2 + \del_z^2) |b>
//
// note: actual evaluation function is FDoublettIntegralFactory::EvalDoublett()
struct FDoublettIntegralFactoryKineticTerm : public FDoublettIntegralFactory
{
	// note: actual evaluation function is FDoublettIntegralFactory::EvalDoublett()
protected:
	void EvalPrimitiveIntegral( FScalar *Result, FPrimitiveIntegralInfo const &info ); // override
};

// used to evaluate point charge potential (->nuclear attraction) integrals,
// integrals over electric fields of point charges, field gradients and similar:
//		\sum_i	Charge_i *
//			<a|(\del_Ri[x])^e (\del_Ri[y])^f (\del_Ri[z])^g) 1/(r - R_i)|b>
// with (R_i, Charge_i) a set of point charges and (e,f,g) = FieldType.
struct FDoublettIntegralFactoryFieldTerms : public FDoublettIntegralFactory
{
	struct FPointCharge{
		FVector3
			Pos;
		FScalar
			Charge;
		FPointCharge();
		FPointCharge( FVector3 Pos_, FScalar Charge_ )
			: Pos( Pos_ ), Charge( Charge_ )
		{}
	};
	typedef std::vector<FPointCharge>
		FPointChargeList;

	typedef TVector3<uint>
		FFieldType; // (0,0,0): potential  (1,0,0): x-comp of electric field,
					// (2,0,0): \del_x^2 potential -> \del_x E

	// initializes the integrator with a given set of point charges
	FDoublettIntegralFactoryFieldTerms( FFieldType FieldType, FPointChargeList PointCharges )
		: m_FieldType( FieldType ), m_PointCharges( PointCharges )
	{};

	// initializes the integrator with a single point charge given by
	// Pos and Charge.
	FDoublettIntegralFactoryFieldTerms( FFieldType FieldType, FVector3 Pos, FScalar Charge )
		: m_FieldType( FieldType )
	{
		m_PointCharges.push_back( FPointCharge( Pos, Charge ) );
	};

protected:
	void EvalPrimitiveIntegral( FScalar *Result, FPrimitiveIntegralInfo const &info ); // override

private:
	FFieldType
		m_FieldType;
	FPointChargeList
		m_PointCharges;
};

// used to evaluate electron densities
//			<a| \delta(r - R) |b>
// with R provided.
struct FDoublettIntegralFactoryDensity : public FDoublettIntegralFactory
{
	FDoublettIntegralFactoryDensity( FVector3 Pos )
		: m_Pos( Pos )
	{}
protected:
	void EvalPrimitiveIntegral( FScalar *Result, FPrimitiveIntegralInfo const &info ); // override
private:
	FVector3
		m_Pos;
};


} // namespace ct

#endif // CT_INT1E_H

// kate: space-indent off; tab-indent on; indent-width 4; mixedindent off; indent-mode normal;

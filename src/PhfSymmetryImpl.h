/*
 * Oct. 2012, Written by Naoki Nakatani
 * ----------------------------------------------------------------------------------------------------
 * Implementation of symmetry objects which is independent from main skeleton
 */
#ifndef PHF_SYMMETRY_IMPL_H
#define PHF_SYMMETRY_IMPL_H

#include "PhfTypes.h"
#include "PhfSolidDef.h"
#include "PhfSymmetry.h"

enum FSymOperator
{
    SOp_E   = 0x0001,
    SOp_X   = 0x0002,
    SOp_Y   = 0x0004,
    SOp_Z   = 0x0008,
    SOp_XY  = 0x0010,
    SOp_XZ  = 0x0020,
    SOp_YZ  = 0x0040,
    SOp_XYZ = 0x0080
};

void GetTSymIndex(const FSuperCell& SuperCell, TArray<FORTINT>& TSymIndex);

FVector3 SymOp_E  (const FVector3& q);
FVector3 SymOp_X  (const FVector3& q);
FVector3 SymOp_Y  (const FVector3& q);
FVector3 SymOp_Z  (const FVector3& q);
FVector3 SymOp_XY (const FVector3& q);
FVector3 SymOp_XZ (const FVector3& q);
FVector3 SymOp_YZ (const FVector3& q);
FVector3 SymOp_XYZ(const FVector3& q);

const double ThreSymmetric = 1e-8;
bool IsSymmetric(const FVector3& q1, const FVector3& q2);
// compute table of symmetry-operation
uint GetSymTable(const FVector3& q);

// determine lattice point group
uint GetSymLattice(const FLattice& lattice);

#endif

/*
 * Oct. 2012, Written by Naoki Nakatani
 * ----------------------------------------------------------------------------------------------------
 * Implementation of symmetry objects which is independent from main skeleton
 */
#ifndef PHF_SYMMETRY_IMPL_H
#define PHF_SYMMETRY_IMPL_H

#include "PhfSolidDef.h"
#include "PhfSymmetry.h"

enum FSymOperator
{
    SOp_E   = 0x01,
    SOp_X   = 0x02,
    SOp_Y   = 0x04,
    SOp_Z   = 0x08,
    SOp_XY  = 0x10,
    SOp_XZ  = 0x20,
    SOp_YZ  = 0x40,
    SOp_XYZ = 0x80
};

void GetTransSymMatrix(const FSuperCell& SuperCell, TArray<FORTINT> TSymMatrix);

#endif

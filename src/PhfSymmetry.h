/*
 * Oct. 2012, Written by Naoki Nakatani
 * ----------------------------------------------------------------------------------------------------
 * Skeleton interfaces to symmetry operation for unit-cell, super-cell, and orbitals
 * They can be called from both Fortran and C/C++ but will be changed to fit global data structure
 */
#ifndef PHF_SYMMETRY_H
#define PHF_SYMMETRY_H

#include "PhfTypes.h"

/*
 * Simple symmetry class which might be usable from FORTRAN
 * other information is utilized from FSolidModel
 *
 * Otherwise, can UnitCell be defined as Symmetry-Unique ?
 * And, can symmetry-operation be implemented in construction of SuperCell, too ?
 * In such case, loop over symmetry-unique atoms is just considered loop over atoms in UnitCell
 */
struct FSymmetry
{
    // point group : controlled by unsigned integer number
    FORTINT
        PointGroup;

    // # of irreducible representation : 1 ~ 8
    FORTINT
        nRep;

    // # of symmetry-adapted AO for each irr-rep
    FORTINT
        nSO[8];

    // index list of symmetry-unique atom center in unit-cell
    TArray<FORTINT>
        SymAtoms;

    // index matrix of translational lattice vector
    // kTs = TSymIndex[iTs, jTs] satisfies [0, kTs] == [iTs, jTs]
    // iTs, jTs, and kTs \in super-cell
    TArray<FORTINT>
        TSymIndex;
};

#endif

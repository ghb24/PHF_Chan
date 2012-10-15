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
 * enumeration of point groups (only Abelian group)
 * numbers are temporary
 */
enum FPointGroup
{
    Pg_C1  = 0,
    Pg_CS  = 1,
    Pg_C2  = 2,
    Pg_D2  = 3,
    Pg_C2V = 4,
    Pg_C2H = 5,
    Pg_D2H = 6
};

/*
 * Simple symmetry class which might be usable from FORTRAN
 * other information is utilized from FSolidModel
 */
struct FSymmetry
{
    // point group : controlled by integer number
    FPointGroup
        PGroup;

    // # of irreducible representation : 1 ~ 8
    FORTINT
        nRep;

    // # of symmetry-adapted AO for each irr-rep
    FORTINT
        nSO[8];

    // index list of symmetry-unique atom center in unit-cell
    TArray<FORTINT>
        SymAtoms;
};

#endif

/*
 * Oct. 2012, Written by Naoki Nakatani
 * ----------------------------------------------------------------------------------------------------
 * Skeleton interfaces to symmetry operation for unit-cell, super-cell, and orbitals
 * They can be called from both Fortran and C/C++ but will be changed to fit global data structure
 */
#ifndef PHF_SYM_ITERFACE_H
#define PHF_SYM_ITERFACE_H

#include "PhfTypes.h"
#include "PhfSolidDef.h"
#include "PhfSymmetry.h"

/*
 * Compute Symmetry-Unique Atoms
 * + Solid  : SolidModel object in "PhfSolidDef.h"
 * + sym    : Symmetry object in "PhfSymmetry.h"
 */
void GetSymmetry(const FSolidModel& Solid, FSymmetry& sym);

/*
 * Matrix transformation between AO and SO
 * SymTrans1 : [ Packed-AO(UNIT) x AO(UNIT) ] -> [ SO(UNIT) x SO(UNIT) ]
 * SymTrans2 : [ Packed-AO(UNIT) x AO(SUPER) ] -> [ SO(SUPER) x SO(SUPER) ]
 * SymTrans3 : [ AO(SUPER) x AO(SUPER) ] -> [ SO(SUPER) x SO(SUPER) ]
 * SymTrans4 : [ SO(SUPER) x SO(SUPER) ] -> [ AO(SUPER) x AO(SUPER) ]
 *
 * + DataAO : data matrix spanned by AO (Real(nAO*nAO), in/out)
 * + DataSO : data matrix spanned by SO (Real(nSO*nSO), in/out)
 * + Solid  : SolidModel object in "PhfSolidDef.h"
 * + sym    : Symmetry object in "PhfSymmetry.h"
 */
void SymTrans1(TArray<double>& DataAO, TArray<double>& DataSO, const FSolidModel& Solid, const FSymmetry& sym);
void SymTrans2(TArray<double>& DataAO, TArray<double>& DataSO, const FSolidModel& Solid, const FSymmetry& sym);
void SymTrans3(TArray<double>& DataAO, TArray<double>& DataSO, const FSolidModel& Solid, const FSymmetry& sym);
void SymTrans4(TArray<double>& DataAO, TArray<double>& DataSO, const FSolidModel& Solid, const FSymmetry& sym);

/*
 * Taking whole information about super-cell (FSolidModel) as an input looks less efficient
 * but, probably better to avoid any confliction an/or complicated dependency (?)
 */
#endif

/*
 * Oct. 2012, Written by Naoki Nakatani
 * ----------------------------------------------------------------------------------------------------
 * Skeleton interfaces to symmetry operation for unit-cell, super-cell, and orbitals
 * They can be called from both Fortran and C/C++ but will be changed to fit global data structure
 */

#include "PhfSymInterface.h"
#include "PhfSymmetryImpl.h"

/*
 * Compute Symmetry-Unique Atoms
 * + Solid  : SolidModel object in "PhfSolidDef.h"
 * + sym    : Symmetry object in "PhfSymmetry.h"
 */
void GetSymmetry(const FSolidModel& Solid, FSymmetry& sym)
{
    uint nAtoms = Solid.UnitCell.Coords.size();
    uint nBasis = Solid.UnitCell.OrbBasis.nFn;
    // default
    sym.PointGroup = 0x0001; // suppose C1, only (SOP_E & 0x01) gives true
    sym.nRep   = 1;
    sym.nSO[0] = nBasis;
    sym.SymAtoms.resize(nAtoms);
    for(uint iAtom = 0; iAtom < nAtoms; ++iAtom) {
        sym.SymAtoms[iAtom] = iAtom;
    }
    GetTSymIndex(Solid.SuperCell, sym.TSymIndex);
};

/*
 * Transformation between AO and SO
 *
 * + DataAO : data matrix spanned by AO (Real(nAO*nAO), in/out)
 * + DataSO : data matrix spanned by SO (Real(nSO*nSO), in/out)
 * + Solid  : SolidModel object in "PhfSolidDef.h"
 * + sym    : Symmetry object in "PhfSymmetry.h"
 */

// SymTrans1 : Packed-AO(UNIT) -> SO(UNIT)
void SymTrans1(TArray<double>& DataAO, TArray<double>& DataSO, const FSolidModel& Solid, const FSymmetry& sym)
{
    DataSO = DataAO; // suppose C1
}

// SymTrans2 : Packed-AO(UNIT) -> SO(SUPER)
void SymTrans2(TArray<double>& DataAO, TArray<double>& DataSO, const FSolidModel& Solid, const FSymmetry& sym)
{
    uint SCellSize = Solid.SuperCell.Ts.size();
    uint nUCellBfn = Solid.UnitCell.OrbBasis.nFn;
    uint nSCellBfn = Solid.SuperCell.OrbBasis.nFn;

    DataSO.clear();
    DataSO.reserve(nSCellBfn * nSCellBfn);

    DataSO.insert(DataSO.end(), DataAO.begin(), DataAO.end()); // suppose C1
    for(uint iTs = 1; iTs < SCellSize; ++iTs) {
      for(uint jTs = 0; jTs < SCellSize; ++jTs) {
          uint iOffBfn = nUCellBfn * sym.TSymIndex[iTs * SCellSize + jTs];
          DataSO.insert(DataSO.end(), DataAO.begin() + iOffBfn, DataAO.begin() + iOffBfn + nUCellBfn);
      }
    }
}

// SymTrans3 : AO(SUPER) -> SO(SUPER)
void SymTrans3(TArray<double>& DataAO, TArray<double>& DataSO, const FSolidModel& Solid, const FSymmetry& sym)
{
    DataSO = DataAO; // suppose C1
}

// SymTrans4 : SO(SUPER) -> AO(SUPER)
void SymTrans4(TArray<double>& DataAO, TArray<double>& DataSO, const FSolidModel& Solid, const FSymmetry& sym)
{
    DataAO = DataSO; // suppose C1
}

// Compute symmetry-unique atom list from input data and return # of symmetry-unique atoms
uint GetSymUniqueAtoms(const TArray<uint>& Elements, const TArray<FVector3>& Coords, TArray<uint>& SymUniqueList)
{
    uint nAtoms = Elements.size();
    uint nSymAtoms = 0;
    SymUniqueList.clear();
    SymUniqueList.resize(nAtoms);
    // initialize: every atom is unique
    for(uint iAtom = 0; iAtom < nAtoms; ++iAtom) SymUniqueList[iAtom] = 1;

    for(uint iAtom = 0; iAtom < nAtoms; ++iAtom) {
      if(SymUniqueList[iAtom] == 0) continue;
      nSymAtoms++;
//    for(uint iSOp = 0; iSOp < nSOps; ++iSOp) {
        // transforming by point group symmetry-operator
//      FVector3 iEqCoord = SOp[iSOp](Coords[iAtom]);
        FVector3 iEqCoord = Coords[iAtom]; // suppose C1, SymOp_E

        for(uint jAtom = iAtom + 1; jAtom < nAtoms; ++jAtom) {
          if(SymUniqueList[jAtom] == 0) continue;

          if(Elements[iAtom] == Elements[jAtom] && IsSymmetric(Coords[jAtom], iEqCoord)) {
            SymUniqueList[iAtom]++;
            SymUniqueList[jAtom]--;
          }
        }
//    }
    }
    return nSymAtoms;
}



/*
 * Oct. 2012, Written by Naoki Nakatani
 * ----------------------------------------------------------------------------------------------------
 * Implementation of symmetry objects which is independent from main skeleton
 */

#include <cmath>
#include "PhfSymmetryImpl.h"

/*
 * Getting an index matrix for translational symmetry
 * + TSymIndex[ iTs, jTs ] returns index of identical translational vector from first unit-cell
 *   (i.e. [ iTs, jTs ] = [ 0, kTs ]
 */
void GetTSymIndex(const FSuperCell& SuperCell, TArray<FORTINT>& TSymIndex)
{
    const TVector3<FORTINT> &SCellSize = SuperCell.Size;
    uint nSCells = SCellSize[0] * SCellSize[1] * SCellSize[2];

    TSymIndex.clear();
    TSymIndex.reserve(nSCells * nSCells);

    for(uint iz = 0; iz < SCellSize[2]; ++iz) {
      for(uint iy = 0; iy < SCellSize[1]; ++iy) {
        for(uint ix = 0; ix < SCellSize[0]; ++ix) {

          for(uint jz = 0; jz < SCellSize[2]; ++jz) {
            FORTINT zOffSet = jz - ix;
            if(zOffSet < 0) zOffSet += SCellSize[2];

            for(uint jy = 0; jy < SCellSize[1]; ++jy) {
              FORTINT yOffSet = jy - iy;
              if(yOffSet < 0) yOffSet += SCellSize[1];

              for(uint jx = 0; jx < SCellSize[0]; ++jx) {
                FORTINT xOffSet = jx - ix;
                if(xOffSet < 0) xOffSet += SCellSize[0];

                TSymIndex.push_back((xOffSet * SCellSize[1] + yOffSet) * SCellSize[2] + zOffSet);
              }
            }
          }
        }
      }
    }
}

FVector3 SymOp_E  (const FVector3& q)
{
    return FVector3( q);
}

FVector3 SymOp_X  (const FVector3& q)
{
    return FVector3(-q[0], q[1], q[2]);
}

FVector3 SymOp_Y  (const FVector3& q)
{
    return FVector3( q[0],-q[1], q[2]);
}

FVector3 SymOp_Z  (const FVector3& q)
{
    return FVector3( q[0], q[1],-q[2]);
}

FVector3 SymOp_XY (const FVector3& q)
{
    return FVector3(-q[0],-q[1], q[2]);
}

FVector3 SymOp_XZ (const FVector3& q)
{
    return FVector3(-q[0], q[1],-q[2]);
}

FVector3 SymOp_YZ (const FVector3& q)
{
    return FVector3( q[0],-q[1],-q[2]);
}

FVector3 SymOp_XYZ(const FVector3& q)
{
    return FVector3(-q[0],-q[1],-q[2]);
}

bool IsSymmetric(const FVector3& q1, const FVector3& q2)
{
    FVector3 qDiff = q2 - q1;
    if(qDiff.LengthSq() < ThreSymmetric)
      return true;
    else
      return false;
}

uint GetSymTable(const FVector3& q1, const FVector3& q2)
{
    uint table = 0x0001;
    return table;
}

uint GetSymLattice(const FLattice& lattice)
{
    uint table = 0x0001;
    return table;
}

/*

void GetSymUniqueAtoms(const TArray<FVector3>& Coords, const TArray<std::string>& Elements, Symmetry& sym)
{
    sym.SymAtoms.clear();
    sym.SymAtoms.reserve(nAtoms);

    TArray<FORTINT> SymUniqueList;
    SymUniqueList.resize(nAtoms);
    // initialize: every atom is unique
    for(uint iAtom = 0; iAtom < nAtoms; ++iAtom) SymUniqueList[iAtom] = 1;

    for(uint iAtom = 0; iAtom < nAtoms; ++iAtom) {
      if(SymUniqueList[iAtom] == 0) continue;
      sym.SymAtoms.push_back(iAtom);
      for(uint iSOp = 0; iSOp < nSOps; ++iSOp) {
        // transforming by point group symmetry-operator
        FVector3 iEqCoord = SOp[iSOp](Coords[iAtom]);

        for(uint jAtom = iAtom + 1; jAtom < nAtoms; ++jAtom) {
          if(SymUniqueList[jAtom] == 0) continue;

          if(std::strcmp(Elements[iAtom], Elements[jAtom]) == 0 && IsSymmetric(Coords[jAtom], iEqCoord))
            SymUniqueList[iAtom]++;
            SymUniqueList[jAtom]--;
          }
        }
      }
    }
}

*/

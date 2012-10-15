/*
 * Oct. 2012, Written by Naoki Nakatani
 * ----------------------------------------------------------------------------------------------------
 * Implementation of symmetry objects which is independent from main skeleton
 */

#include "PhfSymmetryImpl.h"

/*
 * Getting an index matrix for translational symmetry
 * + TSymMatrix[ iTs, jTs ] returns index of identical translational vector from first unit-cell
 *   (i.e. [ iTs, jTs ] = [ 0, kTs ]
 */
void GetTransSymMatrix(const FSuperCell& SuperCell, TArray<FORTINT> TSymMatrix)
{
    const TVector3<FORTINT> &SCellSize = SuperCell.Size;
    uint nSCells = SCellSize[0] * SCellSize[1] * SCellSize[2];

    TSymMatrix.clear();
    TSymMatrix.reserve(nSCells * nSCells);

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

                TSymMatrix.push_back((xOffSet * SCellSize[1] + yOffSet) * SCellSize[2] + zOffSet);
              }
            }
          }
        }
      }
    }
}

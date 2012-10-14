#include "AicGridOps.h"
#include "AicSolidHarmonics.h"


namespace aic {
inline unsigned iSlmX(int l, int m) {
   assert( -l <= m && m <= l );
   return (unsigned)( l*l + l + m );
}

inline unsigned nSlmX(int l) {
   int l1 = l + 1;
   return (unsigned)(l1*l1);
}


const uint g_nGridDerivComp[3] = {1,4,10};

struct FEvalGridBfParams
{
   FVector3 vGridBf; // vGridPoint - vBfPos
   double fDistSq;   // norm(vGridBf)^2
   double *pPowR;    // norm(vGridBf)^l for l = 0.. pBf->l.
   uint MaxPowL;

   uint DerivOrder;
   uint nComp; // 1,4,10 for DerivOrder 0,1,2.
   uint nScalComp; // scalar components. 1 + DerivOrder.
   double LogThrOrb;
   double ThrOrb;
   double *pExpV; // temporary

   enum {
      GEOM_ScreenStage = 1
   };

   FEvalGridBfParams(uint nExpMax, uint nAmMax, uint nCoMax, uint nCompMax, double ThrOrb, double LogThrOrb, FMemoryStack &Mem);

   void SetGeometry(double const *pGridPos, FVector3 const &vBfPos, uint Flags);
   void SetDerivOrder(uint DerivOrder_) {
      DerivOrder = DerivOrder_;
      nComp = g_nGridDerivComp[DerivOrder_];
      nScalComp = DerivOrder + 1;
      assert(DerivOrder_ <= 2);
   };

   void EvalScalar( double *AIC_RP pValues, FShellData const *pBf );
};


FEvalGridBfParams::FEvalGridBfParams(uint nExpMax_, uint nAmMax_, uint nCoMax_, uint nCompMax_, double ThrOrb_, double LogThrOrb_, FMemoryStack &Mem)
{
   MaxPowL = nAmMax_;
   Mem.Alloc(pExpV, nExpMax_ + (nAmMax_+1));
   pPowR = pExpV + nExpMax_;
   ThrOrb = ThrOrb_;
   LogThrOrb = LogThrOrb_;
}


void FEvalGridBfParams::SetGeometry(double const *pGridPos, FVector3 const &vBfPos, uint Flags)
{
   vGridBf[0] = pGridPos[0] - vBfPos[0];
   vGridBf[1] = pGridPos[1] - vBfPos[1];
   vGridBf[2] = pGridPos[2] - vBfPos[2];
   fDistSq = vGridBf.LengthSq();

   // evaluate r^l. This is the leading term of the largest solid harmonic
   // centered on vBfPos. Need that for screening.
   double
      fDist = std::sqrt(fDistSq);
   pPowR[0] = fDist;
   for ( uint l = 1; l <= MaxPowL; ++ l )
      pPowR[l] = pPowR[l-1] * fDist;

   if ( 0 != (GEOM_ScreenStage & Flags) ) {
      // these are geometry parameters for the closest grid point, as used for screening.
      // in this case we may /NOT/ remove basis functions due to small r^l values, since
      // r^l will be larger for other grid points!
      if ( pPowR[0] < 1.0 )
         for ( uint l = 0; l <= MaxPowL; ++ l )
            pPowR[l] = 1.0;
   };
}



void FEvalGridBfParams::EvalScalar( double *AIC_RP pValues, FShellData const *pBf )
{
//    xout << boost::format("   !EvalScalar: fDistSq=%8.3f  LogThr=%8.3f") % fDistSq % LogThrOrb << std::endl;
   bool
      NothingThere = true;
   uint const
      nExp = pBf->nPrim,
      nScalComp_ = this->nScalComp;
   // evaluate the scalar part of the function. First the exponents.
   // here we assume that pPowR can only make the function larger
   // (i.e., the inner parts of high-am functions are not screened out)
   RESUME_CLOCK(23)
   for ( uint iExp = 0; iExp < nExp; ++ iExp ){
      double fExp = pBf->pExp[iExp];
      if ( fDistSq * fExp > LogThrOrb )
         pExpV[iExp] = 0;
      else {
         pExpV[iExp] = std::exp(-fDistSq * fExp);
         NothingThere = false;
      }
   };
   PAUSE_CLOCK(23)

//    xout << boost::format("   !EvalScalar: NothingThere? %i") % (int)NothingThere << std::endl;
   RESUME_CLOCK(24)
   if ( NothingThere ) {
      // all primitives screened out
      for ( uint i = 0; i < nScalComp_ * pBf->nCo; ++ i )
         pValues[i] = 0;
   } else {
      double const
         fPowR = pPowR[pBf->l];
      double const
         *AIC_RP pCo = pBf->pCo,
         *AIC_RP pExp = pBf->pExp;
      for ( uint iCo = 0; iCo < pBf->nCo; ++ iCo, pCo += nExp) {
         // now the total contracted function.
         if ( DerivOrder == 0) {
            assert(nScalComp_ == 1);
            double
               v0 = 0;
            for ( uint iExp = 0; iExp < nExp; ++ iExp )
               v0 += pCo[iExp] * pExpV[iExp];
            if ( std::abs(v0 * fPowR) < ThrOrb )
               v0 = 0;
            pValues[iCo] = v0;
         } else if ( DerivOrder == 1 ) {
            assert(nScalComp_ == 2);
            double v0, v1;
            v0 = 0; // \sum C[iZ] exp(-Z r^2)
            v1 = 0; // \sum C[iZ] (-2Z) exp(-Z r^2)
            // only the scalar parts here.
            for ( uint iExp = 0; iExp != nExp; ++ iExp ) {
               double
                  e = -2 * pExp[iExp],
                  g = pCo[iExp] * pExpV[iExp];
               v0 += g;
               v1 += e * g;
            }

//             if ( std::abs(v0 * fPowR) < ThrOrb && std::abs(v1 * fPowR) < ThrOrb ){
            if ( std::abs(v0 * fPowR) < ThrOrb ){
               v0 = 0;
               v1 = 0;
            }
            pValues[2*iCo    ] = v0; // should be "iCo*nScalComp_"
            pValues[2*iCo + 1] = v1;
         } else {
            assert(DerivOrder == 2);
            assert(nScalComp_ == 3);
            double v0, v1, v2;
            v0 = 0; // \sum C[iZ] exp(-Z r^2)
            v1 = 0; // \sum C[iZ] (-2Z) exp(-Z r^2)
            v2 = 0; // \sum C[iZ] 4 Z^2 exp(-Z r^2)
            // only the scalar parts here. see below.
            for ( uint iExp = 0; iExp != nExp; ++ iExp ) {
               double
                  e = -2 * pExp[iExp],
                  g = pCo[iExp] * pExpV[iExp];
               v0 += g;
               v1 += e * g;
               v2 += e * e * g;
            }

            if ( std::abs(v0 * fPowR) < ThrOrb ){
               v0 = 0;
               v1 = 0;
               v2 = 0;
            }
            pValues[3*iCo    ] = v0;
            pValues[3*iCo + 1] = v1;
            pValues[3*iCo + 2] = v2;
         };
      }
   }
   PAUSE_CLOCK(24)
}


void FindClosestGridPoint( uint &iGridPtMin, double const *pGridPos, uint GridStride, uint nGridPts, FVector3 const &vBfPos )
{
   double
      fDistSqMin = 1e99;
   iGridPtMin = 0xffffffff;
   for ( uint iGridPt = 0; iGridPt < nGridPts; ++ iGridPt ) {
      FVector3
         vGridPos(pGridPos[iGridPt * GridStride],
                  pGridPos[iGridPt * GridStride + 1],
                  pGridPos[iGridPt * GridStride + 2]);
      double
         fDistSq = (vGridPos - vBfPos).LengthSq();
      if ( fDistSq < fDistSqMin ) {
         fDistSqMin = fDistSq;
         iGridPtMin = iGridPt;
      }
   }
}

void FindBfRangeProperties(uint &nCoSum, uint &nFnSum, uint &nCoMax, uint &nExpMax, uint &nAmMax, FShellData const *pFirstBf, FShellData const *pLastBf )
{
   nCoSum = 0;
   nAmMax = 0;
   nExpMax = 0;
   nCoMax = 0;
   nFnSum = 0; // note: not actually required apart from addressing hack. should be removed!
   for ( FShellData const *pBf = pFirstBf; pBf != pLastBf; ++ pBf ) {
      nCoSum += pBf->nCo;
      nFnSum += pBf->nCo * (2*pBf->l + 1);
      if ( pBf->nPrim > nExpMax )
         nExpMax = pBf->nPrim;
      if ( pBf->l > nAmMax )
         nAmMax = pBf->l;
      if ( pBf->nCo > nCoMax )
         nCoMax = pBf->nCo;
   }
}

struct FGridPtBfn
{
   FShellData const
      *pBf;
   uint
      iCo,     // contraction index in pBf
      iFnOff;  // offset of *function* (i.e., including spherical components) regarding iFnBase. //regarding pFirstBf
   FCoMask // an uint.
      iCoMask;
};

// increment the nGridPts x #bf matrix of basis functions centered at vBfPos
// on the given grid points.
//    pOut[iPt + nCompStride * iDiff + nBfStride * iBf] = (value).
//        iDiff: derivative component.
//        iBf: basis function index starting at pFirstBf
void EvalShellGroupOnGrid( double *pOut, uint nCompStride, uint nBfStride,
   FORTINT *pMap, FORTINT &nMap, uint iFnBase,
   double const *pGridPos, uint GridStride, uint nGridPts,
   aic::FShellData const *pFirstBf, aic::FShellData const *pLastBf, aic::FVector3 const &vBfPos, // <- namespace names for doxygen only
   uint DerivOrder, double ThrOrb, double LogThrOrb,
   FMemoryStack2 &Mem )
{
   void
      *pBaseOfStorage = Mem.Alloc(0);

   // determine number of derivative components.
   assert(DerivOrder <= 2);
   uint nComp = g_nGridDerivComp[DerivOrder];

   RESUME_CLOCK(10)

   // find grid point closest to the current basis function center
   uint
      iGridPtMin;
   FindClosestGridPoint( iGridPtMin, pGridPos, GridStride, nGridPts, vBfPos);


   // initialize data structure for evaluating scalar parts of basis functions
   uint
      nCoSum, nFnSum, nCoMax, nAmMax, nExpMax;
   FindBfRangeProperties(nCoSum, nFnSum, nCoMax, nExpMax, nAmMax, pFirstBf, pLastBf);


   FEvalGridBfParams
      P(nExpMax, nAmMax, nCoMax, nComp, ThrOrb, LogThrOrb, Mem);
   double
      *pValues;
   Mem.Alloc(pValues, nCoSum * (1+DerivOrder));

   // make list of basis functions we need to keep
   FGridPtBfn
      *pFirstKeptBf, *pLastKeptBf;
   Mem.Alloc(pFirstKeptBf, nCoSum);
   pLastKeptBf = pFirstKeptBf;
   P.SetDerivOrder(0); // only actual function values for screening.
   P.SetGeometry(&pGridPos[iGridPtMin * GridStride], vBfPos, FEvalGridBfParams::GEOM_ScreenStage);

   uint
      iFnOff = iFnBase;
   int
      MaxAmKept = -1;
   for ( FShellData const *pBf = pFirstBf; pBf != pLastBf; iFnOff += pBf->nCo * (2*pBf->l + 1), ++ pBf ) {
      P.EvalScalar(pValues, pBf);
      FCoMask
         iCoMask = 0;
      uint
         nSh = 2*pBf->l + 1;
      for ( uint iCo = 0; iCo < pBf->nCo; ++ iCo )
         if ( pValues[iCo] != 0 ) {
            iCoMask |= (1ul << iCo);

            for ( uint iSh = 0; iSh < nSh; ++ iSh ){
               pMap[nMap] = iFnOff + iSh + iCo * nSh;
               ++ nMap;
            }
         }

      if ( iCoMask != 0 ) { // at least one non-vanishing contraction. this function must be kept.
         pLastKeptBf->pBf = pBf;
         pLastKeptBf->iFnOff = iFnOff;
         pLastKeptBf->iCoMask = iCoMask;
         if ( (int)pBf->l > MaxAmKept )
            MaxAmKept = (int)pBf->l;
         ++ pLastKeptBf;
      }
   };
   PAUSE_CLOCK(10)

   uint
      nSlmA = nSlmX(MaxAmKept),
      nSlmA1 = nComp * nSlmA;
   double
      *pSlmA;
   Mem.Alloc(pSlmA, nSlmA1);
   P.SetDerivOrder(DerivOrder);
   assert( P.nScalComp == (1+DerivOrder));

   for ( uint iGridPt = 0; iGridPt < nGridPts; ++ iGridPt ) {

      RESUME_CLOCK(20);
      P.SetGeometry(&pGridPos[iGridPt * GridStride], vBfPos, 0);

      int
         nMaxAmOnPt = -1;
      // evaluate the scalar parts of the functions.
      double
         *pVal1 = pValues;
      for ( FGridPtBfn const *pKeptBf = pFirstKeptBf; pKeptBf != pLastKeptBf; ++ pKeptBf ){
         FShellData const *pBf = pKeptBf->pBf;
         P.EvalScalar(pVal1, pBf);
         if ( (int)pBf->l > nMaxAmOnPt )
            nMaxAmOnPt = (int)pBf->l;
         pVal1 += P.nScalComp * pBf->nCo;
      };
      PAUSE_CLOCK(20)

      if ( nMaxAmOnPt == -1 ) {
         // everything screened out.
         continue;
      }

      RESUME_CLOCK(21)
      // evaluate the angular part of the function
      if ( DerivOrder == 0 )
         EvalSlmX( pSlmA, &P.vGridBf[0], nMaxAmOnPt );
      else if ( DerivOrder == 1 )
         EvalSlmX_Grd1( pSlmA, &P.vGridBf[0], nMaxAmOnPt );
      else if ( DerivOrder == 2 )
         EvalSlmX_Grd2( pSlmA, &P.vGridBf[0], nMaxAmOnPt );
      else
         assert(0);
      PAUSE_CLOCK(21)


      // multiply scalar and angular parts together to get the final grid values.
      RESUME_CLOCK(22)
      pVal1 = pValues;
      double
         *AIC_RP pOut0 = &pOut[iGridPt];
      double const
         x = P.vGridBf.m[0],
         y = P.vGridBf.m[1],
         z = P.vGridBf.m[2];
      for ( FGridPtBfn const *pKeptBf = pFirstKeptBf; pKeptBf != pLastKeptBf; ++ pKeptBf ){
         FShellData const
            *pBf = pKeptBf->pBf;
         int
            la = (int)pBf->l;
         uint
            nSh = 2*la + 1;
         double const
            // base of current basis function components in pSlmA
            *AIC_RP pSlm = &pSlmA[nComp * nSlmX(la-1)];

         // TODO: here set everything to exactly zero if it lies below
         // the target threshold.
         if ( DerivOrder == 0 ) {
            for ( uint iCo = 0; iCo < pBf->nCo; ++ iCo ) {
               if ( 0 == (pKeptBf->iCoMask & (1ul<<iCo)))
                  continue;
               for ( uint iSh = 0; iSh < nSh; ++ iSh )
                  pOut0[nBfStride * iSh] = pSlm[iSh] * pVal1[iCo];
               pOut0 += nBfStride * nSh;
            }
         } else if ( DerivOrder == 1 ) {
            //       phi(r) = Slm(r - A)      * \sum C[iZ] exp(-Z (r-A)^2)
            //  d/dx phi(r) = (d/dx Slm(r-A)) * \sum C[iZ] exp(-Z (r-A)^2) +
            //                Slm(r - A) * \sum C[iZ] (-2*Z*(rx-Ax))* exp(-Z (r-A)^2)
            for ( uint iCo = 0; iCo < pBf->nCo; ++ iCo ) {
               if ( 0 == (pKeptBf->iCoMask & (1ul<<iCo)))
                  continue;
               for ( uint iSh = 0; iSh < nSh; ++ iSh ) {
                  double
                     *AIC_RP pOut1 = &pOut0[iSh*nBfStride];
                  double const
                     *AIC_RP pSlm1 = &pSlm[4*iSh];
                  assert(nComp==4);
                  double
                     fSlm = pSlm1[0],
                     fValue = pVal1[2*iCo],
                     fSlmD1 = pVal1[2*iCo+1] * fSlm;
                  pOut1[0*nCompStride] = fSlm * fValue;
                  pOut1[1*nCompStride] = pSlm1[1] * fValue + x * fSlmD1;
                  pOut1[2*nCompStride] = pSlm1[2] * fValue + y * fSlmD1;
                  pOut1[3*nCompStride] = pSlm1[3] * fValue + z * fSlmD1;
               }
               pOut0 += nBfStride * nSh;
            }
         } else if ( DerivOrder == 2 ) {
            //      phi(r)    = Slm(r - A)      * \sum C[iZ] exp(-Z (r-A)^2)
            // d/dx phi(r)    = (d/dx Slm(r-A)) * \sum C[iZ] exp(-Z (r-A)^2) +
            //                  Slm(r - A) * {\sum C[iZ] (-2*Z*(rx-Ax))* exp(-Z (r-A)^2)}
            // d^2/dxx phi(r) = (d^2/dxx Slm(r-A)) * \sum C[iZ] exp(-Z (r-A)^2) +
            //                  2 (d/dx Slm(r-A)) * \sum (-2*Z*(rx-Ax)) C[iZ] exp(-Z (r-A)^2) +
            //                   Slm(r - A) * {\sum C[iZ] ((-2*Z*(rx-Ax))^2 - 2 Z) * exp(-Z (r-A)^2)}
            // d^2/dxy phi(r) = (d^2/dxy Slm(r-A)) * \sum C[iZ] exp(-Z (r-A)^2) +
            //                  (d/dx Slm(r-A)) * \sum (-2*Z*(ry-Ay)) C[iZ] exp(-Z (r-A)^2) +
            //                  (d/dy Slm(r-A)) * \sum (-2*Z*(rx-Ax)) C[iZ] exp(-Z (r-A)^2) +
            //                   Slm(r - A) * {\sum C[iZ] ((-2*Z*(rx-Ax)) (-2*Z*(ry-Ay))) * exp(-Z (r-A)^2)}
            // again, only the scalar parts here.
            for ( uint iCo = 0; iCo < pBf->nCo; ++ iCo ) {
               if ( 0 == (pKeptBf->iCoMask & (1ul<<iCo)))
                  continue;
               for ( uint iSh = 0; iSh < nSh; ++ iSh ) {
                  double
                     *AIC_RP pOut1 = &pOut0[iSh*nBfStride];
                  double const
                     *AIC_RP pSlm1 = &pSlm[10*iSh];
                  assert(nComp==10);
                  double
                     fSlm = pSlm1[0],
                     fSlmDX = pSlm1[1],
                     fSlmDY = pSlm1[2],
                     fSlmDZ = pSlm1[3],
                     fValue = pVal1[3*iCo],
                     fDer1 = pVal1[3*iCo+1],
                     fDer2 = pVal1[3*iCo+2];
                  pOut1[0*nCompStride] = fSlm * fValue;

                  pOut1[1*nCompStride] = fSlmDX * fValue + x * fSlm * fDer1; // d/dx
                  pOut1[2*nCompStride] = fSlmDY * fValue + y * fSlm * fDer1; // d/dy
                  pOut1[3*nCompStride] = fSlmDZ * fValue + z * fSlm * fDer1; // d/dz

                  // hmpf. apparently there are two different sets of orders for second derivatives.
                  // The code was generated for xx/yy/zz/xy/xz/yz, but this function is supposed
                  // to return xx/xy/xz/yy/yz/zz. this is the reason for this re-ordering business here.
                  //     4  5  6  7  8  9
                  //    xx/xy/xz/yy/yz/zz
                  pOut1[4*nCompStride] = pSlm1[4] * fValue + 2 * fSlmDX * x * fDer1 + fSlm * (x * x * fDer2 + fDer1); // d^2/d[xx]
                  pOut1[7*nCompStride] = pSlm1[5] * fValue + 2 * fSlmDY * y * fDer1 + fSlm * (y * y * fDer2 + fDer1); // d^2/d[yy]
                  pOut1[9*nCompStride] = pSlm1[6] * fValue + 2 * fSlmDZ * z * fDer1 + fSlm * (z * z * fDer2 + fDer1); // d^2/d[zz]

                  pOut1[5*nCompStride] = pSlm1[7] * fValue + (x * fSlmDY + y * fSlmDX) * fDer1 + fSlm * y * x * fDer2; // d^2/d[xy]
                  pOut1[6*nCompStride] = pSlm1[8] * fValue + (x * fSlmDZ + z * fSlmDX) * fDer1 + fSlm * z * x * fDer2; // d^2/d[xz]
                  pOut1[8*nCompStride] = pSlm1[9] * fValue + (y * fSlmDZ + z * fSlmDY) * fDer1 + fSlm * y * z * fDer2; // d^2/d[yz]
               }
               pOut0 += nBfStride * nSh;
            }
         }

         pVal1 += P.nScalComp * pBf->nCo;
      }
      PAUSE_CLOCK(22)
   }

   Mem.Free(pBaseOfStorage);
}


} // namespace aic

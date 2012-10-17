#include "AicMdrr.h"

namespace aic{

// Scatter 2-shell result: write a (2la + 1) x (2lc + 1) block of data at pIn
// into memory at pOut[ia*sa + ic*sc], overwriting the previous result.
void Scatter2e2c_O(FReal0 *AIC_RP pOut, unsigned sa, unsigned sc, FReal2 const *AIC_RP pIn, unsigned la, unsigned lc)
{
   // overkill? yes.
   if ( la < 2 ) {
      switch(la + lc*2)
      {
         case 0 + 6*2:
            pOut[0*sa+12*sc] = pIn[0+1*12];
            pOut[0*sa+11*sc] = pIn[0+1*11];
         case 0 + 5*2:
            pOut[0*sa+10*sc] = pIn[0+1*10];
            pOut[0*sa+9*sc] = pIn[0+1*9];
         case 0 + 4*2:
            pOut[0*sa+8*sc] = pIn[0+1*8];
            pOut[0*sa+7*sc] = pIn[0+1*7];
         case 0 + 3*2:
            pOut[0*sa+6*sc] = pIn[0+1*6];
            pOut[0*sa+5*sc] = pIn[0+1*5];
         case 0 + 2*2:
            pOut[0*sa+4*sc] = pIn[0+1*4];
            pOut[0*sa+3*sc] = pIn[0+1*3];
         case 0 + 1*2:
            pOut[0*sa+2*sc] = pIn[0+1*2];
            pOut[0*sa+1*sc] = pIn[0+1*1];
         case 0 + 0*2:
            pOut[0*sa+0*sc] = pIn[0+1*0];
            return;
         case 1 + 6*2:
            pOut[0*sa+12*sc] = pIn[0+3*12];
            pOut[1*sa+12*sc] = pIn[1+3*12];
            pOut[2*sa+12*sc] = pIn[2+3*12];
            pOut[0*sa+11*sc] = pIn[0+3*11];
            pOut[1*sa+11*sc] = pIn[1+3*11];
            pOut[2*sa+11*sc] = pIn[2+3*11];
         case 1 + 5*2:
            pOut[0*sa+10*sc] = pIn[0+3*10];
            pOut[1*sa+10*sc] = pIn[1+3*10];
            pOut[2*sa+10*sc] = pIn[2+3*10];
            pOut[0*sa+9*sc] = pIn[0+3*9];
            pOut[1*sa+9*sc] = pIn[1+3*9];
            pOut[2*sa+9*sc] = pIn[2+3*9];
         case 1 + 4*2:
            pOut[0*sa+8*sc] = pIn[0+3*8];
            pOut[1*sa+8*sc] = pIn[1+3*8];
            pOut[2*sa+8*sc] = pIn[2+3*8];
            pOut[0*sa+7*sc] = pIn[0+3*7];
            pOut[1*sa+7*sc] = pIn[1+3*7];
            pOut[2*sa+7*sc] = pIn[2+3*7];
         case 1 + 3*2:
            pOut[0*sa+6*sc] = pIn[0+3*6];
            pOut[1*sa+6*sc] = pIn[1+3*6];
            pOut[2*sa+6*sc] = pIn[2+3*6];
            pOut[0*sa+5*sc] = pIn[0+3*5];
            pOut[1*sa+5*sc] = pIn[1+3*5];
            pOut[2*sa+5*sc] = pIn[2+3*5];
         case 1 + 2*2:
            pOut[0*sa+4*sc] = pIn[0+3*4];
            pOut[1*sa+4*sc] = pIn[1+3*4];
            pOut[2*sa+4*sc] = pIn[2+3*4];
            pOut[0*sa+3*sc] = pIn[0+3*3];
            pOut[1*sa+3*sc] = pIn[1+3*3];
            pOut[2*sa+3*sc] = pIn[2+3*3];
         case 1 + 1*2:
            pOut[0*sa+2*sc] = pIn[0+3*2];
            pOut[1*sa+2*sc] = pIn[1+3*2];
            pOut[2*sa+2*sc] = pIn[2+3*2];
            pOut[0*sa+1*sc] = pIn[0+3*1];
            pOut[1*sa+1*sc] = pIn[1+3*1];
            pOut[2*sa+1*sc] = pIn[2+3*1];
         case 1 + 0*2:
            pOut[0*sa+0*sc] = pIn[0+3*0];
            pOut[1*sa+0*sc] = pIn[1+3*0];
            pOut[2*sa+0*sc] = pIn[2+3*0];
            return;
      }
   }
   unsigned ic;
   switch(la)
   {
      case 2:
         for ( ic = 0; ic < 2*lc+1; ++ ic ) {
            pOut[0*sa] = pIn[0];
            pOut[1*sa] = pIn[1];
            pOut[2*sa] = pIn[2];
            pOut[3*sa] = pIn[3];
            pOut[4*sa] = pIn[4];
            pOut += sc;
            pIn += 5;
         }
         return;
      case 3:
         for ( ic = 0; ic < 2*lc+1; ++ ic ) {
            pOut[0*sa] = pIn[0];
            pOut[1*sa] = pIn[1];
            pOut[2*sa] = pIn[2];
            pOut[3*sa] = pIn[3];
            pOut[4*sa] = pIn[4];
            pOut[5*sa] = pIn[5];
            pOut[6*sa] = pIn[6];
            pOut += sc;
            pIn += 7;
         }
         return;
      case 4:
         for ( ic = 0; ic < 2*lc+1; ++ ic ) {
            pOut[0*sa] = pIn[0];
            pOut[1*sa] = pIn[1];
            pOut[2*sa] = pIn[2];
            pOut[3*sa] = pIn[3];
            pOut[4*sa] = pIn[4];
            pOut[5*sa] = pIn[5];
            pOut[6*sa] = pIn[6];
            pOut[7*sa] = pIn[7];
            pOut[8*sa] = pIn[8];
            pOut += sc;
            pIn += 9;
         }
         return;
      case 5:
         for ( ic = 0; ic < 2*lc+1; ++ ic ) {
            pOut[0*sa] = pIn[0];
            pOut[1*sa] = pIn[1];
            pOut[2*sa] = pIn[2];
            pOut[3*sa] = pIn[3];
            pOut[4*sa] = pIn[4];
            pOut[5*sa] = pIn[5];
            pOut[6*sa] = pIn[6];
            pOut[7*sa] = pIn[7];
            pOut[8*sa] = pIn[8];
            pOut[9*sa] = pIn[9];
            pOut[10*sa] = pIn[10];
            pOut += sc;
            pIn += 11;
         }
         return;
      case 6:
         for ( ic = 0; ic < 2*lc+1; ++ ic ) {
            pOut[0*sa] = pIn[0];
            pOut[1*sa] = pIn[1];
            pOut[2*sa] = pIn[2];
            pOut[3*sa] = pIn[3];
            pOut[4*sa] = pIn[4];
            pOut[5*sa] = pIn[5];
            pOut[6*sa] = pIn[6];
            pOut[7*sa] = pIn[7];
            pOut[8*sa] = pIn[8];
            pOut[9*sa] = pIn[9];
            pOut[10*sa] = pIn[10];
            pOut[11*sa] = pIn[11];
            pOut[12*sa] = pIn[12];
            pOut += sc;
            pIn += 13;
         }
         return;
   }
}

// Scatter 2-shell result: write a (2la + 1) x (2lc + 1) block of data at pIn
// into memory at pOut[ia*sa + ic*sc], adding to the previous result.
void Scatter2e2c_A(FReal0 *AIC_RP pOut, unsigned sa, unsigned sc, FReal2 const *AIC_RP pIn, unsigned la, unsigned lc)
{
   // did you know that the Zuse Z3 lacked conditional branches?
   // I fully adhere to this philosophy. Begone, evil loops!
   if ( la < 2 ) {
      switch(la + lc*2)
      {
         case 0 + 6*2:
            pOut[0*sa+12*sc] += pIn[0+1*12];
            pOut[0*sa+11*sc] += pIn[0+1*11];
         case 0 + 5*2:
            pOut[0*sa+10*sc] += pIn[0+1*10];
            pOut[0*sa+9*sc] += pIn[0+1*9];
         case 0 + 4*2:
            pOut[0*sa+8*sc] += pIn[0+1*8];
            pOut[0*sa+7*sc] += pIn[0+1*7];
         case 0 + 3*2:
            pOut[0*sa+6*sc] += pIn[0+1*6];
            pOut[0*sa+5*sc] += pIn[0+1*5];
         case 0 + 2*2:
            pOut[0*sa+4*sc] += pIn[0+1*4];
            pOut[0*sa+3*sc] += pIn[0+1*3];
         case 0 + 1*2:
            pOut[0*sa+2*sc] += pIn[0+1*2];
            pOut[0*sa+1*sc] += pIn[0+1*1];
         case 0 + 0*2:
            pOut[0*sa+0*sc] += pIn[0+1*0];
            return;
         case 1 + 6*2:
            pOut[0*sa+12*sc] += pIn[0+3*12];
            pOut[1*sa+12*sc] += pIn[1+3*12];
            pOut[2*sa+12*sc] += pIn[2+3*12];
            pOut[0*sa+11*sc] += pIn[0+3*11];
            pOut[1*sa+11*sc] += pIn[1+3*11];
            pOut[2*sa+11*sc] += pIn[2+3*11];
         case 1 + 5*2:
            pOut[0*sa+10*sc] += pIn[0+3*10];
            pOut[1*sa+10*sc] += pIn[1+3*10];
            pOut[2*sa+10*sc] += pIn[2+3*10];
            pOut[0*sa+9*sc] += pIn[0+3*9];
            pOut[1*sa+9*sc] += pIn[1+3*9];
            pOut[2*sa+9*sc] += pIn[2+3*9];
         case 1 + 4*2:
            pOut[0*sa+8*sc] += pIn[0+3*8];
            pOut[1*sa+8*sc] += pIn[1+3*8];
            pOut[2*sa+8*sc] += pIn[2+3*8];
            pOut[0*sa+7*sc] += pIn[0+3*7];
            pOut[1*sa+7*sc] += pIn[1+3*7];
            pOut[2*sa+7*sc] += pIn[2+3*7];
         case 1 + 3*2:
            pOut[0*sa+6*sc] += pIn[0+3*6];
            pOut[1*sa+6*sc] += pIn[1+3*6];
            pOut[2*sa+6*sc] += pIn[2+3*6];
            pOut[0*sa+5*sc] += pIn[0+3*5];
            pOut[1*sa+5*sc] += pIn[1+3*5];
            pOut[2*sa+5*sc] += pIn[2+3*5];
         case 1 + 2*2:
            pOut[0*sa+4*sc] += pIn[0+3*4];
            pOut[1*sa+4*sc] += pIn[1+3*4];
            pOut[2*sa+4*sc] += pIn[2+3*4];
            pOut[0*sa+3*sc] += pIn[0+3*3];
            pOut[1*sa+3*sc] += pIn[1+3*3];
            pOut[2*sa+3*sc] += pIn[2+3*3];
         case 1 + 1*2:
            pOut[0*sa+2*sc] += pIn[0+3*2];
            pOut[1*sa+2*sc] += pIn[1+3*2];
            pOut[2*sa+2*sc] += pIn[2+3*2];
            pOut[0*sa+1*sc] += pIn[0+3*1];
            pOut[1*sa+1*sc] += pIn[1+3*1];
            pOut[2*sa+1*sc] += pIn[2+3*1];
         case 1 + 0*2:
            pOut[0*sa+0*sc] += pIn[0+3*0];
            pOut[1*sa+0*sc] += pIn[1+3*0];
            pOut[2*sa+0*sc] += pIn[2+3*0];
            return;
      }
   }
   unsigned ic;
   switch(la)
   {
      case 2:
         for ( ic = 0; ic < 2*lc+1; ++ ic ) {
            pOut[0*sa] += pIn[0];
            pOut[1*sa] += pIn[1];
            pOut[2*sa] += pIn[2];
            pOut[3*sa] += pIn[3];
            pOut[4*sa] += pIn[4];
            pOut += sc;
            pIn += 5;
         }
         return;
      case 3:
         for ( ic = 0; ic < 2*lc+1; ++ ic ) {
            pOut[0*sa] += pIn[0];
            pOut[1*sa] += pIn[1];
            pOut[2*sa] += pIn[2];
            pOut[3*sa] += pIn[3];
            pOut[4*sa] += pIn[4];
            pOut[5*sa] += pIn[5];
            pOut[6*sa] += pIn[6];
            pOut += sc;
            pIn += 7;
         }
         return;
      case 4:
         for ( ic = 0; ic < 2*lc+1; ++ ic ) {
            pOut[0*sa] += pIn[0];
            pOut[1*sa] += pIn[1];
            pOut[2*sa] += pIn[2];
            pOut[3*sa] += pIn[3];
            pOut[4*sa] += pIn[4];
            pOut[5*sa] += pIn[5];
            pOut[6*sa] += pIn[6];
            pOut[7*sa] += pIn[7];
            pOut[8*sa] += pIn[8];
            pOut += sc;
            pIn += 9;
         }
         return;
      case 5:
         for ( ic = 0; ic < 2*lc+1; ++ ic ) {
            pOut[0*sa] += pIn[0];
            pOut[1*sa] += pIn[1];
            pOut[2*sa] += pIn[2];
            pOut[3*sa] += pIn[3];
            pOut[4*sa] += pIn[4];
            pOut[5*sa] += pIn[5];
            pOut[6*sa] += pIn[6];
            pOut[7*sa] += pIn[7];
            pOut[8*sa] += pIn[8];
            pOut[9*sa] += pIn[9];
            pOut[10*sa] += pIn[10];
            pOut += sc;
            pIn += 11;
         }
         return;
      case 6:
         for ( ic = 0; ic < 2*lc+1; ++ ic ) {
            pOut[0*sa] += pIn[0];
            pOut[1*sa] += pIn[1];
            pOut[2*sa] += pIn[2];
            pOut[3*sa] += pIn[3];
            pOut[4*sa] += pIn[4];
            pOut[5*sa] += pIn[5];
            pOut[6*sa] += pIn[6];
            pOut[7*sa] += pIn[7];
            pOut[8*sa] += pIn[8];
            pOut[9*sa] += pIn[9];
            pOut[10*sa] += pIn[10];
            pOut[11*sa] += pIn[11];
            pOut[12*sa] += pIn[12];
            pOut += sc;
            pIn += 13;
         }
         return;
   }
}

// In: Gm[0 .. 0] (inclusive), Out: [r]^0, ordered as AngularComps(0)
// (moment 0, 1 entries)
void Mdrr0(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   pOut[0] = pIn[0];
   // 0 flops, 2 mops   (0.00 flops/integral, 2.00 mops/integral)
}

// In: Gm[0 .. 1] (inclusive), Out: [r]^0, ordered as AngularComps(1)
// (moment 1, 3 entries)
void Mdrr1(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   // m = 0
   pOut[0] = R[0] * pIn[1];
   pOut[1] = R[1] * pIn[1];
   pOut[2] = R[2] * pIn[1];
   // 3 flops, 9 mops   (1.00 flops/integral, 3.00 mops/integral)
}

// In: Gm[0 .. 2] (inclusive), Out: [r]^0, ordered as AngularComps(2)
// (moment 2, 6 entries)
void Mdrr2(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 1
   FReal2 r_100_1 = Rx * pIn[2];
   FReal2 r_010_1 = Ry * pIn[2];
   FReal2 r_001_1 = Rz * pIn[2];

   // m = 0
   pOut[0] = Rx * r_100_1 - pIn[1];
   pOut[1] = Ry * r_010_1 - pIn[1];
   pOut[2] = Rz * r_001_1 - pIn[1];
   pOut[3] = Rx * r_010_1;
   pOut[4] = Rz * r_100_1;
   pOut[5] = Ry * r_001_1;
   // 9 flops, 24 mops   (1.50 flops/integral, 4.00 mops/integral)
}

// In: Gm[0 .. 3] (inclusive), Out: [r]^0, ordered as AngularComps(3)
// (moment 3, 10 entries)
void Mdrr3(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 2
   FReal2 r_100_2 = Rx * pIn[3];
   FReal2 r_010_2 = Ry * pIn[3];
   FReal2 r_001_2 = Rz * pIn[3];

   // m = 1
   FReal2 r_200_1 = Rx * r_100_2 - pIn[2];
   FReal2 r_020_1 = Ry * r_010_2 - pIn[2];
   FReal2 r_002_1 = Rz * r_001_2 - pIn[2];
   FReal2 r_100_1 = Rx * pIn[2];
   FReal2 r_010_1 = Ry * pIn[2];
   FReal2 r_001_1 = Rz * pIn[2];
   FReal2 r_101_1 = Rz * r_100_2;

   // m = 0
   pOut[0] = Rx * r_200_1 - 2 * r_100_1;
   pOut[1] = Rx * r_020_1;
   pOut[2] = Rx * r_002_1;
   pOut[3] = Ry * r_200_1;
   pOut[4] = Ry * r_020_1 - 2 * r_010_1;
   pOut[5] = Ry * r_002_1;
   pOut[6] = Rz * r_200_1;
   pOut[7] = Rz * r_020_1;
   pOut[8] = Rz * r_002_1 - 2 * r_001_1;
   pOut[9] = Ry * r_101_1;
   // 23 flops, 49 mops   (2.30 flops/integral, 4.90 mops/integral)
}

// In: Gm[0 .. 4] (inclusive), Out: [r]^0, ordered as AngularComps(4)
// (moment 4, 15 entries)
void Mdrr4(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 3
   FReal2 r_100_3 = Rx * pIn[4];
   FReal2 r_010_3 = Ry * pIn[4];
   FReal2 r_001_3 = Rz * pIn[4];

   // m = 2
   FReal2 r_200_2 = Rx * r_100_3 - pIn[3];
   FReal2 r_020_2 = Ry * r_010_3 - pIn[3];
   FReal2 r_002_2 = Rz * r_001_3 - pIn[3];
   FReal2 r_100_2 = Rx * pIn[3];
   FReal2 r_010_2 = Ry * pIn[3];
   FReal2 r_001_2 = Rz * pIn[3];

   // m = 1
   FReal2 r_200_1 = Rx * r_100_2 - pIn[2];
   FReal2 r_020_1 = Ry * r_010_2 - pIn[2];
   FReal2 r_002_1 = Rz * r_001_2 - pIn[2];
   FReal2 r_300_1 = Rx * r_200_2 - 2 * r_100_2;
   FReal2 r_120_1 = Rx * r_020_2;
   FReal2 r_030_1 = Ry * r_020_2 - 2 * r_010_2;
   FReal2 r_012_1 = Ry * r_002_2;
   FReal2 r_201_1 = Rz * r_200_2;
   FReal2 r_003_1 = Rz * r_002_2 - 2 * r_001_2;

   // m = 0
   pOut[0] = Rx * r_300_1 - 3 * r_200_1;
   pOut[1] = Rx * r_120_1 - r_020_1;
   pOut[2] = Rz * r_201_1 - r_200_1;
   pOut[3] = Ry * r_030_1 - 3 * r_020_1;
   pOut[4] = Ry * r_012_1 - r_002_1;
   pOut[5] = Rz * r_003_1 - 3 * r_002_1;
   pOut[6] = Ry * r_300_1;
   pOut[7] = Rx * r_030_1;
   pOut[8] = Rx * r_012_1;
   pOut[9] = Rz * r_300_1;
   pOut[10] = Rz * r_120_1;
   pOut[11] = Rx * r_003_1;
   pOut[12] = Ry * r_201_1;
   pOut[13] = Rz * r_030_1;
   pOut[14] = Ry * r_003_1;
   // 39 flops, 84 mops   (2.60 flops/integral, 5.60 mops/integral)
}

// In: Gm[0 .. 5] (inclusive), Out: [r]^0, ordered as AngularComps(5)
// (moment 5, 21 entries)
void Mdrr5(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 4
   FReal2 r_100_4 = Rx * pIn[5];
   FReal2 r_010_4 = Ry * pIn[5];
   FReal2 r_001_4 = Rz * pIn[5];

   // m = 3
   FReal2 r_200_3 = Rx * r_100_4 - pIn[4];
   FReal2 r_020_3 = Ry * r_010_4 - pIn[4];
   FReal2 r_002_3 = Rz * r_001_4 - pIn[4];
   FReal2 r_100_3 = Rx * pIn[4];
   FReal2 r_010_3 = Ry * pIn[4];
   FReal2 r_001_3 = Rz * pIn[4];

   // m = 2
   FReal2 r_200_2 = Rx * r_100_3 - pIn[3];
   FReal2 r_020_2 = Ry * r_010_3 - pIn[3];
   FReal2 r_002_2 = Rz * r_001_3 - pIn[3];
   FReal2 r_300_2 = Rx * r_200_3 - 2 * r_100_3;
   FReal2 r_120_2 = Rx * r_020_3;
   FReal2 r_100_2 = Rx * pIn[3];
   FReal2 r_030_2 = Ry * r_020_3 - 2 * r_010_3;
   FReal2 r_012_2 = Ry * r_002_3;
   FReal2 r_010_2 = Ry * pIn[3];
   FReal2 r_201_2 = Rz * r_200_3;
   FReal2 r_003_2 = Rz * r_002_3 - 2 * r_001_3;
   FReal2 r_001_2 = Rz * pIn[3];

   // m = 1
   FReal2 r_400_1 = Rx * r_300_2 - 3 * r_200_2;
   FReal2 r_220_1 = Rx * r_120_2 - r_020_2;
   FReal2 r_202_1 = Rz * r_201_2 - r_200_2;
   FReal2 r_040_1 = Ry * r_030_2 - 3 * r_020_2;
   FReal2 r_022_1 = Ry * r_012_2 - r_002_2;
   FReal2 r_004_1 = Rz * r_003_2 - 3 * r_002_2;
   FReal2 r_300_1 = Rx * r_200_2 - 2 * r_100_2;
   FReal2 r_030_1 = Ry * r_020_2 - 2 * r_010_2;
   FReal2 r_003_1 = Rz * r_002_2 - 2 * r_001_2;
   FReal2 r_310_1 = Ry * r_300_2;
   FReal2 r_130_1 = Rx * r_030_2;
   FReal2 r_301_1 = Rz * r_300_2;
   FReal2 r_103_1 = Rx * r_003_2;
   FReal2 r_031_1 = Rz * r_030_2;
   FReal2 r_013_1 = Ry * r_003_2;

   // m = 0
   pOut[0] = Rx * r_400_1 - 4 * r_300_1;
   pOut[1] = Ry * r_310_1 - r_300_1;
   pOut[2] = Rz * r_301_1 - r_300_1;
   pOut[3] = Rx * r_040_1;
   pOut[4] = Rx * r_022_1;
   pOut[5] = Rx * r_004_1;
   pOut[6] = Ry * r_400_1;
   pOut[7] = Rx * r_130_1 - r_030_1;
   pOut[8] = Ry * r_202_1;
   pOut[9] = Ry * r_040_1 - 4 * r_030_1;
   pOut[10] = Rz * r_031_1 - r_030_1;
   pOut[11] = Ry * r_004_1;
   pOut[12] = Rz * r_400_1;
   pOut[13] = Rz * r_220_1;
   pOut[14] = Rx * r_103_1 - r_003_1;
   pOut[15] = Rz * r_040_1;
   pOut[16] = Ry * r_013_1 - r_003_1;
   pOut[17] = Rz * r_004_1 - 4 * r_003_1;
   pOut[18] = Ry * r_301_1;
   pOut[19] = Rz * r_130_1;
   pOut[20] = Rx * r_013_1;
   // 69 flops, 144 mops   (3.29 flops/integral, 6.86 mops/integral)
}

// In: Gm[0 .. 6] (inclusive), Out: [r]^0, ordered as AngularComps(6)
// (moment 6, 28 entries)
void Mdrr6(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 5
   FReal2 r_100_5 = Rx * pIn[6];
   FReal2 r_010_5 = Ry * pIn[6];
   FReal2 r_001_5 = Rz * pIn[6];

   // m = 4
   FReal2 r_200_4 = Rx * r_100_5 - pIn[5];
   FReal2 r_020_4 = Ry * r_010_5 - pIn[5];
   FReal2 r_002_4 = Rz * r_001_5 - pIn[5];
   FReal2 r_100_4 = Rx * pIn[5];
   FReal2 r_010_4 = Ry * pIn[5];
   FReal2 r_001_4 = Rz * pIn[5];

   // m = 3
   FReal2 r_200_3 = Rx * r_100_4 - pIn[4];
   FReal2 r_020_3 = Ry * r_010_4 - pIn[4];
   FReal2 r_002_3 = Rz * r_001_4 - pIn[4];
   FReal2 r_300_3 = Rx * r_200_4 - 2 * r_100_4;
   FReal2 r_100_3 = Rx * pIn[4];
   FReal2 r_030_3 = Ry * r_020_4 - 2 * r_010_4;
   FReal2 r_010_3 = Ry * pIn[4];
   FReal2 r_201_3 = Rz * r_200_4;
   FReal2 r_003_3 = Rz * r_002_4 - 2 * r_001_4;
   FReal2 r_001_3 = Rz * pIn[4];

   // m = 2
   FReal2 r_400_2 = Rx * r_300_3 - 3 * r_200_3;
   FReal2 r_202_2 = Rz * r_201_3 - r_200_3;
   FReal2 r_040_2 = Ry * r_030_3 - 3 * r_020_3;
   FReal2 r_004_2 = Rz * r_003_3 - 3 * r_002_3;
   FReal2 r_200_2 = Rx * r_100_3 - pIn[3];
   FReal2 r_020_2 = Ry * r_010_3 - pIn[3];
   FReal2 r_002_2 = Rz * r_001_3 - pIn[3];
   FReal2 r_300_2 = Rx * r_200_3 - 2 * r_100_3;
   FReal2 r_030_2 = Ry * r_020_3 - 2 * r_010_3;
   FReal2 r_201_2 = Rz * r_200_3;
   FReal2 r_003_2 = Rz * r_002_3 - 2 * r_001_3;
   FReal2 r_310_2 = Ry * r_300_3;
   FReal2 r_130_2 = Rx * r_030_3;
   FReal2 r_301_2 = Rz * r_300_3;
   FReal2 r_031_2 = Rz * r_030_3;
   FReal2 r_013_2 = Ry * r_003_3;

   // m = 1
   FReal2 r_400_1 = Rx * r_300_2 - 3 * r_200_2;
   FReal2 r_202_1 = Rz * r_201_2 - r_200_2;
   FReal2 r_040_1 = Ry * r_030_2 - 3 * r_020_2;
   FReal2 r_004_1 = Rz * r_003_2 - 3 * r_002_2;
   FReal2 r_500_1 = Rx * r_400_2 - 4 * r_300_2;
   FReal2 r_320_1 = Ry * r_310_2 - r_300_2;
   FReal2 r_302_1 = Rz * r_301_2 - r_300_2;
   FReal2 r_140_1 = Rx * r_040_2;
   FReal2 r_410_1 = Ry * r_400_2;
   FReal2 r_230_1 = Rx * r_130_2 - r_030_2;
   FReal2 r_212_1 = Ry * r_202_2;
   FReal2 r_050_1 = Ry * r_040_2 - 4 * r_030_2;
   FReal2 r_032_1 = Rz * r_031_2 - r_030_2;
   FReal2 r_014_1 = Ry * r_004_2;
   FReal2 r_401_1 = Rz * r_400_2;
   FReal2 r_203_1 = Rz * r_202_2 - 2 * r_201_2;
   FReal2 r_041_1 = Rz * r_040_2;
   FReal2 r_023_1 = Ry * r_013_2 - r_003_2;
   FReal2 r_005_1 = Rz * r_004_2 - 4 * r_003_2;
   FReal2 r_130_1 = Rx * r_030_2;
   FReal2 r_301_1 = Rz * r_300_2;
   FReal2 r_013_1 = Ry * r_003_2;

   // m = 0
   pOut[0] = Rx * r_500_1 - 5 * r_400_1;
   pOut[1] = Ry * r_410_1 - r_400_1;
   pOut[2] = Rx * r_302_1 - 3 * r_202_1;
   pOut[3] = Rx * r_140_1 - r_040_1;
   pOut[4] = Ry * r_212_1 - r_202_1;
   pOut[5] = Rz * r_203_1 - 3 * r_202_1;
   pOut[6] = Ry * r_050_1 - 5 * r_040_1;
   pOut[7] = Rz * r_041_1 - r_040_1;
   pOut[8] = Ry * r_014_1 - r_004_1;
   pOut[9] = Rz * r_005_1 - 5 * r_004_1;
   pOut[10] = Ry * r_500_1;
   pOut[11] = Rx * r_230_1 - 2 * r_130_1;
   pOut[12] = Ry * r_302_1;
   pOut[13] = Rx * r_050_1;
   pOut[14] = Rx * r_032_1;
   pOut[15] = Rx * r_014_1;
   pOut[16] = Rx * r_401_1 - 4 * r_301_1;
   pOut[17] = Rz * r_320_1;
   pOut[18] = Rz * r_302_1 - 2 * r_301_1;
   pOut[19] = Rz * r_140_1;
   pOut[20] = Rx * r_023_1;
   pOut[21] = Rx * r_005_1;
   pOut[22] = Ry * r_401_1;
   pOut[23] = Rz * r_230_1;
   pOut[24] = Ry * r_203_1;
   pOut[25] = Rz * r_050_1;
   pOut[26] = Ry * r_023_1 - 2 * r_013_1;
   pOut[27] = Rz * r_014_1 - 4 * r_013_1;
   // 111 flops, 220 mops   (3.96 flops/integral, 7.86 mops/integral)
}

// In: Gm[0 .. 7] (inclusive), Out: [r]^0, ordered as AngularComps(7)
// (moment 7, 36 entries)
void Mdrr7(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 6
   FReal2 r_100_6 = Rx * pIn[7];
   FReal2 r_010_6 = Ry * pIn[7];
   FReal2 r_001_6 = Rz * pIn[7];

   // m = 5
   FReal2 r_200_5 = Rx * r_100_6 - pIn[6];
   FReal2 r_020_5 = Ry * r_010_6 - pIn[6];
   FReal2 r_002_5 = Rz * r_001_6 - pIn[6];
   FReal2 r_100_5 = Rx * pIn[6];
   FReal2 r_010_5 = Ry * pIn[6];
   FReal2 r_001_5 = Rz * pIn[6];

   // m = 4
   FReal2 r_200_4 = Rx * r_100_5 - pIn[5];
   FReal2 r_020_4 = Ry * r_010_5 - pIn[5];
   FReal2 r_002_4 = Rz * r_001_5 - pIn[5];
   FReal2 r_300_4 = Rx * r_200_5 - 2 * r_100_5;
   FReal2 r_100_4 = Rx * pIn[5];
   FReal2 r_030_4 = Ry * r_020_5 - 2 * r_010_5;
   FReal2 r_010_4 = Ry * pIn[5];
   FReal2 r_003_4 = Rz * r_002_5 - 2 * r_001_5;
   FReal2 r_001_4 = Rz * pIn[5];

   // m = 3
   FReal2 r_400_3 = Rx * r_300_4 - 3 * r_200_4;
   FReal2 r_040_3 = Ry * r_030_4 - 3 * r_020_4;
   FReal2 r_004_3 = Rz * r_003_4 - 3 * r_002_4;
   FReal2 r_200_3 = Rx * r_100_4 - pIn[4];
   FReal2 r_020_3 = Ry * r_010_4 - pIn[4];
   FReal2 r_002_3 = Rz * r_001_4 - pIn[4];
   FReal2 r_300_3 = Rx * r_200_4 - 2 * r_100_4;
   FReal2 r_100_3 = Rx * pIn[4];
   FReal2 r_030_3 = Ry * r_020_4 - 2 * r_010_4;
   FReal2 r_010_3 = Ry * pIn[4];
   FReal2 r_003_3 = Rz * r_002_4 - 2 * r_001_4;
   FReal2 r_001_3 = Rz * pIn[4];
   FReal2 r_130_3 = Rx * r_030_4;
   FReal2 r_301_3 = Rz * r_300_4;
   FReal2 r_013_3 = Ry * r_003_4;

   // m = 2
   FReal2 r_400_2 = Rx * r_300_3 - 3 * r_200_3;
   FReal2 r_040_2 = Ry * r_030_3 - 3 * r_020_3;
   FReal2 r_004_2 = Rz * r_003_3 - 3 * r_002_3;
   FReal2 r_500_2 = Rx * r_400_3 - 4 * r_300_3;
   FReal2 r_302_2 = Rz * r_301_3 - r_300_3;
   FReal2 r_140_2 = Rx * r_040_3;
   FReal2 r_104_2 = Rx * r_004_3;
   FReal2 r_300_2 = Rx * r_200_3 - 2 * r_100_3;
   FReal2 r_410_2 = Ry * r_400_3;
   FReal2 r_230_2 = Rx * r_130_3 - r_030_3;
   FReal2 r_050_2 = Ry * r_040_3 - 4 * r_030_3;
   FReal2 r_014_2 = Ry * r_004_3;
   FReal2 r_030_2 = Ry * r_020_3 - 2 * r_010_3;
   FReal2 r_401_2 = Rz * r_400_3;
   FReal2 r_041_2 = Rz * r_040_3;
   FReal2 r_023_2 = Ry * r_013_3 - r_003_3;
   FReal2 r_005_2 = Rz * r_004_3 - 4 * r_003_3;
   FReal2 r_003_2 = Rz * r_002_3 - 2 * r_001_3;
   FReal2 r_130_2 = Rx * r_030_3;
   FReal2 r_301_2 = Rz * r_300_3;
   FReal2 r_013_2 = Ry * r_003_3;

   // m = 1
   FReal2 r_600_1 = Rx * r_500_2 - 5 * r_400_2;
   FReal2 r_420_1 = Ry * r_410_2 - r_400_2;
   FReal2 r_402_1 = Rz * r_401_2 - r_400_2;
   FReal2 r_240_1 = Rx * r_140_2 - r_040_2;
   FReal2 r_204_1 = Rx * r_104_2 - r_004_2;
   FReal2 r_060_1 = Ry * r_050_2 - 5 * r_040_2;
   FReal2 r_042_1 = Rz * r_041_2 - r_040_2;
   FReal2 r_024_1 = Ry * r_014_2 - r_004_2;
   FReal2 r_006_1 = Rz * r_005_2 - 5 * r_004_2;
   FReal2 r_500_1 = Rx * r_400_2 - 4 * r_300_2;
   FReal2 r_302_1 = Rz * r_301_2 - r_300_2;
   FReal2 r_140_1 = Rx * r_040_2;
   FReal2 r_230_1 = Rx * r_130_2 - r_030_2;
   FReal2 r_050_1 = Ry * r_040_2 - 4 * r_030_2;
   FReal2 r_014_1 = Ry * r_004_2;
   FReal2 r_401_1 = Rz * r_400_2;
   FReal2 r_023_1 = Ry * r_013_2 - r_003_2;
   FReal2 r_005_1 = Rz * r_004_2 - 4 * r_003_2;
   FReal2 r_510_1 = Ry * r_500_2;
   FReal2 r_330_1 = Rx * r_230_2 - 2 * r_130_2;
   FReal2 r_312_1 = Ry * r_302_2;
   FReal2 r_150_1 = Rx * r_050_2;
   FReal2 r_501_1 = Rz * r_500_2;
   FReal2 r_303_1 = Rz * r_302_2 - 2 * r_301_2;
   FReal2 r_123_1 = Rx * r_023_2;
   FReal2 r_105_1 = Rx * r_005_2;
   FReal2 r_231_1 = Rz * r_230_2;
   FReal2 r_051_1 = Rz * r_050_2;
   FReal2 r_033_1 = Ry * r_023_2 - 2 * r_013_2;
   FReal2 r_015_1 = Ry * r_005_2;

   // m = 0
   pOut[0] = Rx * r_600_1 - 6 * r_500_1;
   pOut[1] = Ry * r_510_1 - r_500_1;
   pOut[2] = Rx * r_402_1 - 4 * r_302_1;
   pOut[3] = Rx * r_240_1 - 2 * r_140_1;
   pOut[4] = Ry * r_312_1 - r_302_1;
   pOut[5] = Rz * r_303_1 - 3 * r_302_1;
   pOut[6] = Rx * r_060_1;
   pOut[7] = Rx * r_042_1;
   pOut[8] = Rx * r_024_1;
   pOut[9] = Rx * r_006_1;
   pOut[10] = Ry * r_600_1;
   pOut[11] = Rx * r_330_1 - 3 * r_230_1;
   pOut[12] = Ry * r_402_1;
   pOut[13] = Rx * r_150_1 - r_050_1;
   pOut[14] = Rz * r_231_1 - r_230_1;
   pOut[15] = Ry * r_204_1;
   pOut[16] = Ry * r_060_1 - 6 * r_050_1;
   pOut[17] = Rz * r_051_1 - r_050_1;
   pOut[18] = Ry * r_024_1 - 2 * r_014_1;
   pOut[19] = Rz * r_015_1 - 5 * r_014_1;
   pOut[20] = Rx * r_501_1 - 5 * r_401_1;
   pOut[21] = Rz * r_420_1;
   pOut[22] = Rz * r_402_1 - 2 * r_401_1;
   pOut[23] = Rz * r_240_1;
   pOut[24] = Rx * r_123_1 - r_023_1;
   pOut[25] = Rx * r_105_1 - r_005_1;
   pOut[26] = Rz * r_060_1;
   pOut[27] = Ry * r_033_1 - 3 * r_023_1;
   pOut[28] = Rz * r_024_1 - 4 * r_023_1;
   pOut[29] = Rz * r_006_1 - 6 * r_005_1;
   pOut[30] = Ry * r_501_1;
   pOut[31] = Rz * r_330_1;
   pOut[32] = Ry * r_303_1;
   pOut[33] = Rz * r_150_1;
   pOut[34] = Rx * r_033_1;
   pOut[35] = Rx * r_015_1;
   // 160 flops, 311 mops   (4.44 flops/integral, 8.64 mops/integral)
}

// In: Gm[0 .. 8] (inclusive), Out: [r]^0, ordered as AngularComps(8)
// (moment 8, 45 entries)
void Mdrr8(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 7
   FReal2 r_100_7 = Rx * pIn[8];
   FReal2 r_010_7 = Ry * pIn[8];
   FReal2 r_001_7 = Rz * pIn[8];

   // m = 6
   FReal2 r_200_6 = Rx * r_100_7 - pIn[7];
   FReal2 r_020_6 = Ry * r_010_7 - pIn[7];
   FReal2 r_002_6 = Rz * r_001_7 - pIn[7];
   FReal2 r_100_6 = Rx * pIn[7];
   FReal2 r_010_6 = Ry * pIn[7];
   FReal2 r_001_6 = Rz * pIn[7];

   // m = 5
   FReal2 r_200_5 = Rx * r_100_6 - pIn[6];
   FReal2 r_020_5 = Ry * r_010_6 - pIn[6];
   FReal2 r_002_5 = Rz * r_001_6 - pIn[6];
   FReal2 r_300_5 = Rx * r_200_6 - 2 * r_100_6;
   FReal2 r_100_5 = Rx * pIn[6];
   FReal2 r_030_5 = Ry * r_020_6 - 2 * r_010_6;
   FReal2 r_010_5 = Ry * pIn[6];
   FReal2 r_003_5 = Rz * r_002_6 - 2 * r_001_6;
   FReal2 r_001_5 = Rz * pIn[6];

   // m = 4
   FReal2 r_400_4 = Rx * r_300_5 - 3 * r_200_5;
   FReal2 r_040_4 = Ry * r_030_5 - 3 * r_020_5;
   FReal2 r_004_4 = Rz * r_003_5 - 3 * r_002_5;
   FReal2 r_200_4 = Rx * r_100_5 - pIn[5];
   FReal2 r_020_4 = Ry * r_010_5 - pIn[5];
   FReal2 r_002_4 = Rz * r_001_5 - pIn[5];
   FReal2 r_300_4 = Rx * r_200_5 - 2 * r_100_5;
   FReal2 r_100_4 = Rx * pIn[5];
   FReal2 r_030_4 = Ry * r_020_5 - 2 * r_010_5;
   FReal2 r_010_4 = Ry * pIn[5];
   FReal2 r_003_4 = Rz * r_002_5 - 2 * r_001_5;
   FReal2 r_001_4 = Rz * pIn[5];
   FReal2 r_130_4 = Rx * r_030_5;
   FReal2 r_301_4 = Rz * r_300_5;
   FReal2 r_013_4 = Ry * r_003_5;

   // m = 3
   FReal2 r_400_3 = Rx * r_300_4 - 3 * r_200_4;
   FReal2 r_040_3 = Ry * r_030_4 - 3 * r_020_4;
   FReal2 r_004_3 = Rz * r_003_4 - 3 * r_002_4;
   FReal2 r_200_3 = Rx * r_100_4 - pIn[4];
   FReal2 r_020_3 = Ry * r_010_4 - pIn[4];
   FReal2 r_002_3 = Rz * r_001_4 - pIn[4];
   FReal2 r_500_3 = Rx * r_400_4 - 4 * r_300_4;
   FReal2 r_302_3 = Rz * r_301_4 - r_300_4;
   FReal2 r_140_3 = Rx * r_040_4;
   FReal2 r_300_3 = Rx * r_200_4 - 2 * r_100_4;
   FReal2 r_230_3 = Rx * r_130_4 - r_030_4;
   FReal2 r_050_3 = Ry * r_040_4 - 4 * r_030_4;
   FReal2 r_014_3 = Ry * r_004_4;
   FReal2 r_030_3 = Ry * r_020_4 - 2 * r_010_4;
   FReal2 r_401_3 = Rz * r_400_4;
   FReal2 r_023_3 = Ry * r_013_4 - r_003_4;
   FReal2 r_005_3 = Rz * r_004_4 - 4 * r_003_4;
   FReal2 r_003_3 = Rz * r_002_4 - 2 * r_001_4;
   FReal2 r_130_3 = Rx * r_030_4;
   FReal2 r_301_3 = Rz * r_300_4;
   FReal2 r_013_3 = Ry * r_003_4;

   // m = 2
   FReal2 r_600_2 = Rx * r_500_3 - 5 * r_400_3;
   FReal2 r_402_2 = Rz * r_401_3 - r_400_3;
   FReal2 r_240_2 = Rx * r_140_3 - r_040_3;
   FReal2 r_060_2 = Ry * r_050_3 - 5 * r_040_3;
   FReal2 r_024_2 = Ry * r_014_3 - r_004_3;
   FReal2 r_006_2 = Rz * r_005_3 - 5 * r_004_3;
   FReal2 r_400_2 = Rx * r_300_3 - 3 * r_200_3;
   FReal2 r_040_2 = Ry * r_030_3 - 3 * r_020_3;
   FReal2 r_004_2 = Rz * r_003_3 - 3 * r_002_3;
   FReal2 r_500_2 = Rx * r_400_3 - 4 * r_300_3;
   FReal2 r_302_2 = Rz * r_301_3 - r_300_3;
   FReal2 r_140_2 = Rx * r_040_3;
   FReal2 r_230_2 = Rx * r_130_3 - r_030_3;
   FReal2 r_050_2 = Ry * r_040_3 - 4 * r_030_3;
   FReal2 r_014_2 = Ry * r_004_3;
   FReal2 r_401_2 = Rz * r_400_3;
   FReal2 r_023_2 = Ry * r_013_3 - r_003_3;
   FReal2 r_005_2 = Rz * r_004_3 - 4 * r_003_3;
   FReal2 r_510_2 = Ry * r_500_3;
   FReal2 r_330_2 = Rx * r_230_3 - 2 * r_130_3;
   FReal2 r_130_2 = Rx * r_030_3;
   FReal2 r_303_2 = Rz * r_302_3 - 2 * r_301_3;
   FReal2 r_105_2 = Rx * r_005_3;
   FReal2 r_301_2 = Rz * r_300_3;
   FReal2 r_051_2 = Rz * r_050_3;
   FReal2 r_033_2 = Ry * r_023_3 - 2 * r_013_3;
   FReal2 r_013_2 = Ry * r_003_3;

   // m = 1
   FReal2 r_600_1 = Rx * r_500_2 - 5 * r_400_2;
   FReal2 r_402_1 = Rz * r_401_2 - r_400_2;
   FReal2 r_240_1 = Rx * r_140_2 - r_040_2;
   FReal2 r_060_1 = Ry * r_050_2 - 5 * r_040_2;
   FReal2 r_024_1 = Ry * r_014_2 - r_004_2;
   FReal2 r_006_1 = Rz * r_005_2 - 5 * r_004_2;
   FReal2 r_700_1 = Rx * r_600_2 - 6 * r_500_2;
   FReal2 r_520_1 = Ry * r_510_2 - r_500_2;
   FReal2 r_502_1 = Rx * r_402_2 - 4 * r_302_2;
   FReal2 r_340_1 = Rx * r_240_2 - 2 * r_140_2;
   FReal2 r_304_1 = Rz * r_303_2 - 3 * r_302_2;
   FReal2 r_160_1 = Rx * r_060_2;
   FReal2 r_124_1 = Rx * r_024_2;
   FReal2 r_106_1 = Rx * r_006_2;
   FReal2 r_610_1 = Ry * r_600_2;
   FReal2 r_430_1 = Rx * r_330_2 - 3 * r_230_2;
   FReal2 r_412_1 = Ry * r_402_2;
   FReal2 r_250_1 = Ry * r_240_2 - 4 * r_230_2;
   FReal2 r_070_1 = Ry * r_060_2 - 6 * r_050_2;
   FReal2 r_052_1 = Rz * r_051_2 - r_050_2;
   FReal2 r_034_1 = Ry * r_024_2 - 2 * r_014_2;
   FReal2 r_016_1 = Ry * r_006_2;
   FReal2 r_601_1 = Rz * r_600_2;
   FReal2 r_403_1 = Rz * r_402_2 - 2 * r_401_2;
   FReal2 r_241_1 = Rz * r_240_2;
   FReal2 r_205_1 = Rx * r_105_2 - r_005_2;
   FReal2 r_061_1 = Rz * r_060_2;
   FReal2 r_043_1 = Ry * r_033_2 - 3 * r_023_2;
   FReal2 r_025_1 = Rz * r_024_2 - 4 * r_023_2;
   FReal2 r_007_1 = Rz * r_006_2 - 6 * r_005_2;
   FReal2 r_330_1 = Rx * r_230_2 - 2 * r_130_2;
   FReal2 r_303_1 = Rz * r_302_2 - 2 * r_301_2;
   FReal2 r_033_1 = Ry * r_023_2 - 2 * r_013_2;
   FReal2 r_331_1 = Rz * r_330_2;
   FReal2 r_313_1 = Ry * r_303_2;
   FReal2 r_133_1 = Rx * r_033_2;

   // m = 0
   pOut[0] = Rx * r_700_1 - 7 * r_600_1;
   pOut[1] = Ry * r_610_1 - r_600_1;
   pOut[2] = Rx * r_502_1 - 5 * r_402_1;
   pOut[3] = Rx * r_340_1 - 3 * r_240_1;
   pOut[4] = Ry * r_412_1 - r_402_1;
   pOut[5] = Rz * r_403_1 - 3 * r_402_1;
   pOut[6] = Rx * r_160_1 - r_060_1;
   pOut[7] = Rz * r_241_1 - r_240_1;
   pOut[8] = Rx * r_124_1 - r_024_1;
   pOut[9] = Rx * r_106_1 - r_006_1;
   pOut[10] = Ry * r_070_1 - 7 * r_060_1;
   pOut[11] = Rz * r_061_1 - r_060_1;
   pOut[12] = Ry * r_034_1 - 3 * r_024_1;
   pOut[13] = Rz * r_025_1 - 5 * r_024_1;
   pOut[14] = Rz * r_007_1 - 7 * r_006_1;
   pOut[15] = Ry * r_700_1;
   pOut[16] = Rx * r_430_1 - 4 * r_330_1;
   pOut[17] = Ry * r_502_1;
   pOut[18] = Ry * r_340_1 - 4 * r_330_1;
   pOut[19] = Rz * r_331_1 - r_330_1;
   pOut[20] = Ry * r_304_1;
   pOut[21] = Rx * r_070_1;
   pOut[22] = Rx * r_052_1;
   pOut[23] = Rx * r_034_1;
   pOut[24] = Rx * r_016_1;
   pOut[25] = Rz * r_700_1;
   pOut[26] = Rz * r_520_1;
   pOut[27] = Rx * r_403_1 - 4 * r_303_1;
   pOut[28] = Rz * r_340_1;
   pOut[29] = Ry * r_313_1 - r_303_1;
   pOut[30] = Rz * r_304_1 - 4 * r_303_1;
   pOut[31] = Rz * r_160_1;
   pOut[32] = Rx * r_043_1;
   pOut[33] = Rx * r_025_1;
   pOut[34] = Rx * r_007_1;
   pOut[35] = Ry * r_601_1;
   pOut[36] = Rz * r_430_1;
   pOut[37] = Ry * r_403_1;
   pOut[38] = Rz * r_250_1;
   pOut[39] = Rx * r_133_1 - r_033_1;
   pOut[40] = Ry * r_205_1;
   pOut[41] = Rz * r_070_1;
   pOut[42] = Ry * r_043_1 - 4 * r_033_1;
   pOut[43] = Rz * r_034_1 - 4 * r_033_1;
   pOut[44] = Ry * r_007_1;
   // 224 flops, 426 mops   (4.98 flops/integral, 9.47 mops/integral)
}

// In: Gm[0 .. 9] (inclusive), Out: [r]^0, ordered as AngularComps(9)
// (moment 9, 55 entries)
void Mdrr9(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 8
   FReal2 r_100_8 = Rx * pIn[9];
   FReal2 r_010_8 = Ry * pIn[9];
   FReal2 r_001_8 = Rz * pIn[9];

   // m = 7
   FReal2 r_200_7 = Rx * r_100_8 - pIn[8];
   FReal2 r_020_7 = Ry * r_010_8 - pIn[8];
   FReal2 r_002_7 = Rz * r_001_8 - pIn[8];
   FReal2 r_100_7 = Rx * pIn[8];
   FReal2 r_010_7 = Ry * pIn[8];
   FReal2 r_001_7 = Rz * pIn[8];

   // m = 6
   FReal2 r_200_6 = Rx * r_100_7 - pIn[7];
   FReal2 r_020_6 = Ry * r_010_7 - pIn[7];
   FReal2 r_002_6 = Rz * r_001_7 - pIn[7];
   FReal2 r_300_6 = Rx * r_200_7 - 2 * r_100_7;
   FReal2 r_100_6 = Rx * pIn[7];
   FReal2 r_030_6 = Ry * r_020_7 - 2 * r_010_7;
   FReal2 r_010_6 = Ry * pIn[7];
   FReal2 r_003_6 = Rz * r_002_7 - 2 * r_001_7;
   FReal2 r_001_6 = Rz * pIn[7];

   // m = 5
   FReal2 r_400_5 = Rx * r_300_6 - 3 * r_200_6;
   FReal2 r_040_5 = Ry * r_030_6 - 3 * r_020_6;
   FReal2 r_004_5 = Rz * r_003_6 - 3 * r_002_6;
   FReal2 r_200_5 = Rx * r_100_6 - pIn[6];
   FReal2 r_020_5 = Ry * r_010_6 - pIn[6];
   FReal2 r_002_5 = Rz * r_001_6 - pIn[6];
   FReal2 r_300_5 = Rx * r_200_6 - 2 * r_100_6;
   FReal2 r_100_5 = Rx * pIn[6];
   FReal2 r_030_5 = Ry * r_020_6 - 2 * r_010_6;
   FReal2 r_010_5 = Ry * pIn[6];
   FReal2 r_003_5 = Rz * r_002_6 - 2 * r_001_6;
   FReal2 r_001_5 = Rz * pIn[6];
   FReal2 r_301_5 = Rz * r_300_6;

   // m = 4
   FReal2 r_400_4 = Rx * r_300_5 - 3 * r_200_5;
   FReal2 r_040_4 = Ry * r_030_5 - 3 * r_020_5;
   FReal2 r_004_4 = Rz * r_003_5 - 3 * r_002_5;
   FReal2 r_200_4 = Rx * r_100_5 - pIn[5];
   FReal2 r_020_4 = Ry * r_010_5 - pIn[5];
   FReal2 r_002_4 = Rz * r_001_5 - pIn[5];
   FReal2 r_500_4 = Rx * r_400_5 - 4 * r_300_5;
   FReal2 r_302_4 = Rz * r_301_5 - r_300_5;
   FReal2 r_140_4 = Rx * r_040_5;
   FReal2 r_300_4 = Rx * r_200_5 - 2 * r_100_5;
   FReal2 r_100_4 = Rx * pIn[5];
   FReal2 r_410_4 = Ry * r_400_5;
   FReal2 r_050_4 = Ry * r_040_5 - 4 * r_030_5;
   FReal2 r_014_4 = Ry * r_004_5;
   FReal2 r_030_4 = Ry * r_020_5 - 2 * r_010_5;
   FReal2 r_010_4 = Ry * pIn[5];
   FReal2 r_401_4 = Rz * r_400_5;
   FReal2 r_041_4 = Rz * r_040_5;
   FReal2 r_005_4 = Rz * r_004_5 - 4 * r_003_5;
   FReal2 r_003_4 = Rz * r_002_5 - 2 * r_001_5;
   FReal2 r_001_4 = Rz * pIn[5];
   FReal2 r_301_4 = Rz * r_300_5;

   // m = 3
   FReal2 r_600_3 = Rx * r_500_4 - 5 * r_400_4;
   FReal2 r_420_3 = Ry * r_410_4 - r_400_4;
   FReal2 r_402_3 = Rz * r_401_4 - r_400_4;
   FReal2 r_240_3 = Rx * r_140_4 - r_040_4;
   FReal2 r_060_3 = Ry * r_050_4 - 5 * r_040_4;
   FReal2 r_042_3 = Rz * r_041_4 - r_040_4;
   FReal2 r_024_3 = Ry * r_014_4 - r_004_4;
   FReal2 r_006_3 = Rz * r_005_4 - 5 * r_004_4;
   FReal2 r_400_3 = Rx * r_300_4 - 3 * r_200_4;
   FReal2 r_040_3 = Ry * r_030_4 - 3 * r_020_4;
   FReal2 r_004_3 = Rz * r_003_4 - 3 * r_002_4;
   FReal2 r_500_3 = Rx * r_400_4 - 4 * r_300_4;
   FReal2 r_302_3 = Rz * r_301_4 - r_300_4;
   FReal2 r_140_3 = Rx * r_040_4;
   FReal2 r_104_3 = Rx * r_004_4;
   FReal2 r_300_3 = Rx * r_200_4 - 2 * r_100_4;
   FReal2 r_410_3 = Ry * r_400_4;
   FReal2 r_050_3 = Ry * r_040_4 - 4 * r_030_4;
   FReal2 r_014_3 = Ry * r_004_4;
   FReal2 r_030_3 = Ry * r_020_4 - 2 * r_010_4;
   FReal2 r_401_3 = Rz * r_400_4;
   FReal2 r_041_3 = Rz * r_040_4;
   FReal2 r_005_3 = Rz * r_004_4 - 4 * r_003_4;
   FReal2 r_003_3 = Rz * r_002_4 - 2 * r_001_4;
   FReal2 r_510_3 = Ry * r_500_4;
   FReal2 r_150_3 = Rx * r_050_4;
   FReal2 r_303_3 = Rz * r_302_4 - 2 * r_301_4;
   FReal2 r_105_3 = Rx * r_005_4;
   FReal2 r_301_3 = Rz * r_300_4;
   FReal2 r_051_3 = Rz * r_050_4;
   FReal2 r_015_3 = Ry * r_005_4;

   // m = 2
   FReal2 r_600_2 = Rx * r_500_3 - 5 * r_400_3;
   FReal2 r_420_2 = Ry * r_410_3 - r_400_3;
   FReal2 r_402_2 = Rz * r_401_3 - r_400_3;
   FReal2 r_240_2 = Rx * r_140_3 - r_040_3;
   FReal2 r_204_2 = Rx * r_104_3 - r_004_3;
   FReal2 r_060_2 = Ry * r_050_3 - 5 * r_040_3;
   FReal2 r_042_2 = Rz * r_041_3 - r_040_3;
   FReal2 r_024_2 = Ry * r_014_3 - r_004_3;
   FReal2 r_006_2 = Rz * r_005_3 - 5 * r_004_3;
   FReal2 r_700_2 = Rx * r_600_3 - 6 * r_500_3;
   FReal2 r_520_2 = Ry * r_510_3 - r_500_3;
   FReal2 r_502_2 = Rx * r_402_3 - 4 * r_302_3;
   FReal2 r_340_2 = Rx * r_240_3 - 2 * r_140_3;
   FReal2 r_304_2 = Rz * r_303_3 - 3 * r_302_3;
   FReal2 r_500_2 = Rx * r_400_3 - 4 * r_300_3;
   FReal2 r_140_2 = Rx * r_040_3;
   FReal2 r_104_2 = Rx * r_004_3;
   FReal2 r_430_2 = Ry * r_420_3 - 2 * r_410_3;
   FReal2 r_250_2 = Rx * r_150_3 - r_050_3;
   FReal2 r_070_2 = Ry * r_060_3 - 6 * r_050_3;
   FReal2 r_052_2 = Rz * r_051_3 - r_050_3;
   FReal2 r_034_2 = Ry * r_024_3 - 2 * r_014_3;
   FReal2 r_410_2 = Ry * r_400_3;
   FReal2 r_050_2 = Ry * r_040_3 - 4 * r_030_3;
   FReal2 r_014_2 = Ry * r_004_3;
   FReal2 r_403_2 = Rz * r_402_3 - 2 * r_401_3;
   FReal2 r_205_2 = Rx * r_105_3 - r_005_3;
   FReal2 r_043_2 = Rz * r_042_3 - 2 * r_041_3;
   FReal2 r_025_2 = Ry * r_015_3 - r_005_3;
   FReal2 r_007_2 = Rz * r_006_3 - 6 * r_005_3;
   FReal2 r_401_2 = Rz * r_400_3;
   FReal2 r_041_2 = Rz * r_040_3;
   FReal2 r_005_2 = Rz * r_004_3 - 4 * r_003_3;
   FReal2 r_510_2 = Ry * r_500_3;
   FReal2 r_150_2 = Rx * r_050_3;
   FReal2 r_501_2 = Rz * r_500_3;
   FReal2 r_303_2 = Rz * r_302_3 - 2 * r_301_3;
   FReal2 r_051_2 = Rz * r_050_3;
   FReal2 r_015_2 = Ry * r_005_3;
   FReal2 r_313_2 = Ry * r_303_3;

   // m = 1
   FReal2 r_800_1 = Rx * r_700_2 - 7 * r_600_2;
   FReal2 r_620_1 = Rx * r_520_2 - 5 * r_420_2;
   FReal2 r_602_1 = Rx * r_502_2 - 5 * r_402_2;
   FReal2 r_440_1 = Rx * r_340_2 - 3 * r_240_2;
   FReal2 r_404_1 = Rz * r_403_2 - 3 * r_402_2;
   FReal2 r_260_1 = Ry * r_250_2 - 5 * r_240_2;
   FReal2 r_206_1 = Rz * r_205_2 - 5 * r_204_2;
   FReal2 r_080_1 = Ry * r_070_2 - 7 * r_060_2;
   FReal2 r_062_1 = Ry * r_052_2 - 5 * r_042_2;
   FReal2 r_044_1 = Ry * r_034_2 - 3 * r_024_2;
   FReal2 r_026_1 = Rz * r_025_2 - 5 * r_024_2;
   FReal2 r_008_1 = Rz * r_007_2 - 7 * r_006_2;
   FReal2 r_700_1 = Rx * r_600_2 - 6 * r_500_2;
   FReal2 r_502_1 = Rz * r_501_2 - r_500_2;
   FReal2 r_340_1 = Rx * r_240_2 - 2 * r_140_2;
   FReal2 r_304_1 = Rx * r_204_2 - 2 * r_104_2;
   FReal2 r_430_1 = Ry * r_420_2 - 2 * r_410_2;
   FReal2 r_250_1 = Rx * r_150_2 - r_050_2;
   FReal2 r_070_1 = Ry * r_060_2 - 6 * r_050_2;
   FReal2 r_034_1 = Ry * r_024_2 - 2 * r_014_2;
   FReal2 r_403_1 = Rz * r_402_2 - 2 * r_401_2;
   FReal2 r_043_1 = Rz * r_042_2 - 2 * r_041_2;
   FReal2 r_025_1 = Ry * r_015_2 - r_005_2;
   FReal2 r_007_1 = Rz * r_006_2 - 6 * r_005_2;
   FReal2 r_710_1 = Ry * r_700_2;
   FReal2 r_530_1 = Ry * r_520_2 - 2 * r_510_2;
   FReal2 r_512_1 = Ry * r_502_2;
   FReal2 r_350_1 = Rx * r_250_2 - 2 * r_150_2;
   FReal2 r_314_1 = Ry * r_304_2;
   FReal2 r_170_1 = Rx * r_070_2;
   FReal2 r_134_1 = Rx * r_034_2;
   FReal2 r_701_1 = Rz * r_700_2;
   FReal2 r_503_1 = Rz * r_502_2 - 2 * r_501_2;
   FReal2 r_341_1 = Rz * r_340_2;
   FReal2 r_323_1 = Ry * r_313_2 - r_303_2;
   FReal2 r_305_1 = Rz * r_304_2 - 4 * r_303_2;
   FReal2 r_143_1 = Rx * r_043_2;
   FReal2 r_125_1 = Rx * r_025_2;
   FReal2 r_107_1 = Rx * r_007_2;
   FReal2 r_431_1 = Rz * r_430_2;
   FReal2 r_413_1 = Ry * r_403_2;
   FReal2 r_251_1 = Rz * r_250_2;
   FReal2 r_071_1 = Rz * r_070_2;
   FReal2 r_053_1 = Rz * r_052_2 - 2 * r_051_2;
   FReal2 r_035_1 = Ry * r_025_2 - 2 * r_015_2;
   FReal2 r_017_1 = Ry * r_007_2;
   FReal2 r_313_1 = Ry * r_303_2;

   // m = 0
   pOut[0] = Rx * r_800_1 - 8 * r_700_1;
   pOut[1] = Ry * r_710_1 - r_700_1;
   pOut[2] = Rx * r_602_1 - 6 * r_502_1;
   pOut[3] = Rx * r_440_1 - 4 * r_340_1;
   pOut[4] = Ry * r_512_1 - r_502_1;
   pOut[5] = Rx * r_404_1 - 4 * r_304_1;
   pOut[6] = Ry * r_350_1 - 5 * r_340_1;
   pOut[7] = Rz * r_341_1 - r_340_1;
   pOut[8] = Ry * r_314_1 - r_304_1;
   pOut[9] = Rz * r_305_1 - 5 * r_304_1;
   pOut[10] = Rx * r_080_1;
   pOut[11] = Rx * r_062_1;
   pOut[12] = Rx * r_044_1;
   pOut[13] = Rx * r_026_1;
   pOut[14] = Rx * r_008_1;
   pOut[15] = Ry * r_800_1;
   pOut[16] = Rx * r_530_1 - 5 * r_430_1;
   pOut[17] = Ry * r_602_1;
   pOut[18] = Rx * r_350_1 - 3 * r_250_1;
   pOut[19] = Rz * r_431_1 - r_430_1;
   pOut[20] = Ry * r_404_1;
   pOut[21] = Rx * r_170_1 - r_070_1;
   pOut[22] = Rz * r_251_1 - r_250_1;
   pOut[23] = Rx * r_134_1 - r_034_1;
   pOut[24] = Ry * r_206_1;
   pOut[25] = Ry * r_080_1 - 8 * r_070_1;
   pOut[26] = Rz * r_071_1 - r_070_1;
   pOut[27] = Ry * r_044_1 - 4 * r_034_1;
   pOut[28] = Rz * r_035_1 - 5 * r_034_1;
   pOut[29] = Ry * r_008_1;
   pOut[30] = Rz * r_800_1;
   pOut[31] = Rz * r_620_1;
   pOut[32] = Rx * r_503_1 - 5 * r_403_1;
   pOut[33] = Rz * r_440_1;
   pOut[34] = Ry * r_413_1 - r_403_1;
   pOut[35] = Rz * r_404_1 - 4 * r_403_1;
   pOut[36] = Rz * r_260_1;
   pOut[37] = Rx * r_143_1 - r_043_1;
   pOut[38] = Rx * r_125_1 - r_025_1;
   pOut[39] = Rx * r_107_1 - r_007_1;
   pOut[40] = Rz * r_080_1;
   pOut[41] = Ry * r_053_1 - 5 * r_043_1;
   pOut[42] = Ry * r_035_1 - 3 * r_025_1;
   pOut[43] = Rz * r_026_1 - 6 * r_025_1;
   pOut[44] = Rz * r_008_1 - 8 * r_007_1;
   pOut[45] = Ry * r_701_1;
   pOut[46] = Rz * r_530_1;
   pOut[47] = Rx * r_413_1 - 4 * r_313_1;
   pOut[48] = Rz * r_350_1;
   pOut[49] = Ry * r_323_1 - 2 * r_313_1;
   pOut[50] = Rz * r_314_1 - 4 * r_313_1;
   pOut[51] = Rz * r_170_1;
   pOut[52] = Rx * r_053_1;
   pOut[53] = Rx * r_035_1;
   pOut[54] = Rx * r_017_1;
   // 321 flops, 597 mops   (5.84 flops/integral, 10.85 mops/integral)
}

// In: Gm[0 .. 10] (inclusive), Out: [r]^0, ordered as AngularComps(10)
// (moment 10, 66 entries)
void Mdrr10(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 9
   FReal2 r_100_9 = Rx * pIn[10];
   FReal2 r_010_9 = Ry * pIn[10];
   FReal2 r_001_9 = Rz * pIn[10];

   // m = 8
   FReal2 r_200_8 = Rx * r_100_9 - pIn[9];
   FReal2 r_020_8 = Ry * r_010_9 - pIn[9];
   FReal2 r_002_8 = Rz * r_001_9 - pIn[9];
   FReal2 r_100_8 = Rx * pIn[9];
   FReal2 r_010_8 = Ry * pIn[9];
   FReal2 r_001_8 = Rz * pIn[9];

   // m = 7
   FReal2 r_200_7 = Rx * r_100_8 - pIn[8];
   FReal2 r_020_7 = Ry * r_010_8 - pIn[8];
   FReal2 r_002_7 = Rz * r_001_8 - pIn[8];
   FReal2 r_300_7 = Rx * r_200_8 - 2 * r_100_8;
   FReal2 r_100_7 = Rx * pIn[8];
   FReal2 r_030_7 = Ry * r_020_8 - 2 * r_010_8;
   FReal2 r_010_7 = Ry * pIn[8];
   FReal2 r_003_7 = Rz * r_002_8 - 2 * r_001_8;
   FReal2 r_001_7 = Rz * pIn[8];

   // m = 6
   FReal2 r_400_6 = Rx * r_300_7 - 3 * r_200_7;
   FReal2 r_040_6 = Ry * r_030_7 - 3 * r_020_7;
   FReal2 r_004_6 = Rz * r_003_7 - 3 * r_002_7;
   FReal2 r_200_6 = Rx * r_100_7 - pIn[7];
   FReal2 r_020_6 = Ry * r_010_7 - pIn[7];
   FReal2 r_002_6 = Rz * r_001_7 - pIn[7];
   FReal2 r_300_6 = Rx * r_200_7 - 2 * r_100_7;
   FReal2 r_100_6 = Rx * pIn[7];
   FReal2 r_030_6 = Ry * r_020_7 - 2 * r_010_7;
   FReal2 r_010_6 = Ry * pIn[7];
   FReal2 r_003_6 = Rz * r_002_7 - 2 * r_001_7;
   FReal2 r_001_6 = Rz * pIn[7];

   // m = 5
   FReal2 r_400_5 = Rx * r_300_6 - 3 * r_200_6;
   FReal2 r_040_5 = Ry * r_030_6 - 3 * r_020_6;
   FReal2 r_004_5 = Rz * r_003_6 - 3 * r_002_6;
   FReal2 r_200_5 = Rx * r_100_6 - pIn[6];
   FReal2 r_020_5 = Ry * r_010_6 - pIn[6];
   FReal2 r_002_5 = Rz * r_001_6 - pIn[6];
   FReal2 r_500_5 = Rx * r_400_6 - 4 * r_300_6;
   FReal2 r_140_5 = Rx * r_040_6;
   FReal2 r_300_5 = Rx * r_200_6 - 2 * r_100_6;
   FReal2 r_100_5 = Rx * pIn[6];
   FReal2 r_050_5 = Ry * r_040_6 - 4 * r_030_6;
   FReal2 r_014_5 = Ry * r_004_6;
   FReal2 r_030_5 = Ry * r_020_6 - 2 * r_010_6;
   FReal2 r_010_5 = Ry * pIn[6];
   FReal2 r_401_5 = Rz * r_400_6;
   FReal2 r_005_5 = Rz * r_004_6 - 4 * r_003_6;
   FReal2 r_003_5 = Rz * r_002_6 - 2 * r_001_6;
   FReal2 r_001_5 = Rz * pIn[6];

   // m = 4
   FReal2 r_600_4 = Rx * r_500_5 - 5 * r_400_5;
   FReal2 r_402_4 = Rz * r_401_5 - r_400_5;
   FReal2 r_240_4 = Rx * r_140_5 - r_040_5;
   FReal2 r_060_4 = Ry * r_050_5 - 5 * r_040_5;
   FReal2 r_024_4 = Ry * r_014_5 - r_004_5;
   FReal2 r_006_4 = Rz * r_005_5 - 5 * r_004_5;
   FReal2 r_400_4 = Rx * r_300_5 - 3 * r_200_5;
   FReal2 r_040_4 = Ry * r_030_5 - 3 * r_020_5;
   FReal2 r_004_4 = Rz * r_003_5 - 3 * r_002_5;
   FReal2 r_200_4 = Rx * r_100_5 - pIn[5];
   FReal2 r_020_4 = Ry * r_010_5 - pIn[5];
   FReal2 r_002_4 = Rz * r_001_5 - pIn[5];
   FReal2 r_500_4 = Rx * r_400_5 - 4 * r_300_5;
   FReal2 r_140_4 = Rx * r_040_5;
   FReal2 r_300_4 = Rx * r_200_5 - 2 * r_100_5;
   FReal2 r_050_4 = Ry * r_040_5 - 4 * r_030_5;
   FReal2 r_014_4 = Ry * r_004_5;
   FReal2 r_030_4 = Ry * r_020_5 - 2 * r_010_5;
   FReal2 r_401_4 = Rz * r_400_5;
   FReal2 r_005_4 = Rz * r_004_5 - 4 * r_003_5;
   FReal2 r_003_4 = Rz * r_002_5 - 2 * r_001_5;
   FReal2 r_510_4 = Ry * r_500_5;
   FReal2 r_150_4 = Rx * r_050_5;
   FReal2 r_501_4 = Rz * r_500_5;
   FReal2 r_105_4 = Rx * r_005_5;
   FReal2 r_051_4 = Rz * r_050_5;
   FReal2 r_015_4 = Ry * r_005_5;

   // m = 3
   FReal2 r_600_3 = Rx * r_500_4 - 5 * r_400_4;
   FReal2 r_402_3 = Rz * r_401_4 - r_400_4;
   FReal2 r_240_3 = Rx * r_140_4 - r_040_4;
   FReal2 r_060_3 = Ry * r_050_4 - 5 * r_040_4;
   FReal2 r_024_3 = Ry * r_014_4 - r_004_4;
   FReal2 r_006_3 = Rz * r_005_4 - 5 * r_004_4;
   FReal2 r_400_3 = Rx * r_300_4 - 3 * r_200_4;
   FReal2 r_040_3 = Ry * r_030_4 - 3 * r_020_4;
   FReal2 r_004_3 = Rz * r_003_4 - 3 * r_002_4;
   FReal2 r_700_3 = Rx * r_600_4 - 6 * r_500_4;
   FReal2 r_520_3 = Ry * r_510_4 - r_500_4;
   FReal2 r_502_3 = Rz * r_501_4 - r_500_4;
   FReal2 r_340_3 = Rx * r_240_4 - 2 * r_140_4;
   FReal2 r_106_3 = Rx * r_006_4;
   FReal2 r_500_3 = Rx * r_400_4 - 4 * r_300_4;
   FReal2 r_140_3 = Rx * r_040_4;
   FReal2 r_610_3 = Ry * r_600_4;
   FReal2 r_250_3 = Rx * r_150_4 - r_050_4;
   FReal2 r_070_3 = Ry * r_060_4 - 6 * r_050_4;
   FReal2 r_052_3 = Rz * r_051_4 - r_050_4;
   FReal2 r_034_3 = Ry * r_024_4 - 2 * r_014_4;
   FReal2 r_050_3 = Ry * r_040_4 - 4 * r_030_4;
   FReal2 r_014_3 = Ry * r_004_4;
   FReal2 r_403_3 = Rz * r_402_4 - 2 * r_401_4;
   FReal2 r_205_3 = Rx * r_105_4 - r_005_4;
   FReal2 r_061_3 = Rz * r_060_4;
   FReal2 r_025_3 = Ry * r_015_4 - r_005_4;
   FReal2 r_007_3 = Rz * r_006_4 - 6 * r_005_4;
   FReal2 r_401_3 = Rz * r_400_4;
   FReal2 r_005_3 = Rz * r_004_4 - 4 * r_003_4;
   FReal2 r_510_3 = Ry * r_500_4;
   FReal2 r_150_3 = Rx * r_050_4;
   FReal2 r_501_3 = Rz * r_500_4;
   FReal2 r_105_3 = Rx * r_005_4;
   FReal2 r_051_3 = Rz * r_050_4;
   FReal2 r_015_3 = Ry * r_005_4;

   // m = 2
   FReal2 r_800_2 = Rx * r_700_3 - 7 * r_600_3;
   FReal2 r_620_2 = Ry * r_610_3 - r_600_3;
   FReal2 r_602_2 = Rx * r_502_3 - 5 * r_402_3;
   FReal2 r_440_2 = Rx * r_340_3 - 3 * r_240_3;
   FReal2 r_404_2 = Rz * r_403_3 - 3 * r_402_3;
   FReal2 r_260_2 = Ry * r_250_3 - 5 * r_240_3;
   FReal2 r_206_2 = Rx * r_106_3 - r_006_3;
   FReal2 r_080_2 = Ry * r_070_3 - 7 * r_060_3;
   FReal2 r_062_2 = Rz * r_061_3 - r_060_3;
   FReal2 r_044_2 = Ry * r_034_3 - 3 * r_024_3;
   FReal2 r_026_2 = Rz * r_025_3 - 5 * r_024_3;
   FReal2 r_008_2 = Rz * r_007_3 - 7 * r_006_3;
   FReal2 r_600_2 = Rx * r_500_3 - 5 * r_400_3;
   FReal2 r_402_2 = Rz * r_401_3 - r_400_3;
   FReal2 r_240_2 = Rx * r_140_3 - r_040_3;
   FReal2 r_060_2 = Ry * r_050_3 - 5 * r_040_3;
   FReal2 r_024_2 = Ry * r_014_3 - r_004_3;
   FReal2 r_006_2 = Rz * r_005_3 - 5 * r_004_3;
   FReal2 r_700_2 = Rx * r_600_3 - 6 * r_500_3;
   FReal2 r_520_2 = Ry * r_510_3 - r_500_3;
   FReal2 r_502_2 = Rz * r_501_3 - r_500_3;
   FReal2 r_340_2 = Rx * r_240_3 - 2 * r_140_3;
   FReal2 r_106_2 = Rx * r_006_3;
   FReal2 r_610_2 = Ry * r_600_3;
   FReal2 r_250_2 = Rx * r_150_3 - r_050_3;
   FReal2 r_070_2 = Ry * r_060_3 - 6 * r_050_3;
   FReal2 r_052_2 = Rz * r_051_3 - r_050_3;
   FReal2 r_034_2 = Ry * r_024_3 - 2 * r_014_3;
   FReal2 r_403_2 = Rz * r_402_3 - 2 * r_401_3;
   FReal2 r_205_2 = Rx * r_105_3 - r_005_3;
   FReal2 r_061_2 = Rz * r_060_3;
   FReal2 r_025_2 = Ry * r_015_3 - r_005_3;
   FReal2 r_007_2 = Rz * r_006_3 - 6 * r_005_3;
   FReal2 r_530_2 = Ry * r_520_3 - 2 * r_510_3;
   FReal2 r_350_2 = Rx * r_250_3 - 2 * r_150_3;
   FReal2 r_134_2 = Rx * r_034_3;
   FReal2 r_510_2 = Ry * r_500_3;
   FReal2 r_150_2 = Rx * r_050_3;
   FReal2 r_503_2 = Rz * r_502_3 - 2 * r_501_3;
   FReal2 r_341_2 = Rz * r_340_3;
   FReal2 r_305_2 = Rx * r_205_3 - 2 * r_105_3;
   FReal2 r_501_2 = Rz * r_500_3;
   FReal2 r_105_2 = Rx * r_005_3;
   FReal2 r_413_2 = Ry * r_403_3;
   FReal2 r_053_2 = Rz * r_052_3 - 2 * r_051_3;
   FReal2 r_035_2 = Ry * r_025_3 - 2 * r_015_3;
   FReal2 r_051_2 = Rz * r_050_3;
   FReal2 r_015_2 = Ry * r_005_3;

   // m = 1
   FReal2 r_800_1 = Rx * r_700_2 - 7 * r_600_2;
   FReal2 r_602_1 = Rx * r_502_2 - 5 * r_402_2;
   FReal2 r_440_1 = Rx * r_340_2 - 3 * r_240_2;
   FReal2 r_404_1 = Rz * r_403_2 - 3 * r_402_2;
   FReal2 r_260_1 = Ry * r_250_2 - 5 * r_240_2;
   FReal2 r_080_1 = Ry * r_070_2 - 7 * r_060_2;
   FReal2 r_044_1 = Ry * r_034_2 - 3 * r_024_2;
   FReal2 r_026_1 = Rz * r_025_2 - 5 * r_024_2;
   FReal2 r_008_1 = Rz * r_007_2 - 7 * r_006_2;
   FReal2 r_900_1 = Rx * r_800_2 - 8 * r_700_2;
   FReal2 r_720_1 = Rx * r_620_2 - 6 * r_520_2;
   FReal2 r_702_1 = Rx * r_602_2 - 6 * r_502_2;
   FReal2 r_540_1 = Ry * r_530_2 - 3 * r_520_2;
   FReal2 r_504_1 = Rz * r_503_2 - 3 * r_502_2;
   FReal2 r_360_1 = Ry * r_350_2 - 5 * r_340_2;
   FReal2 r_342_1 = Rz * r_341_2 - r_340_2;
   FReal2 r_306_1 = Rx * r_206_2 - 2 * r_106_2;
   FReal2 r_180_1 = Rx * r_080_2;
   FReal2 r_144_1 = Rx * r_044_2;
   FReal2 r_126_1 = Rx * r_026_2;
   FReal2 r_108_1 = Rx * r_008_2;
   FReal2 r_810_1 = Ry * r_800_2;
   FReal2 r_630_1 = Ry * r_620_2 - 2 * r_610_2;
   FReal2 r_612_1 = Ry * r_602_2;
   FReal2 r_450_1 = Rx * r_350_2 - 3 * r_250_2;
   FReal2 r_414_1 = Ry * r_404_2;
   FReal2 r_270_1 = Ry * r_260_2 - 6 * r_250_2;
   FReal2 r_234_1 = Rx * r_134_2 - r_034_2;
   FReal2 r_090_1 = Ry * r_080_2 - 8 * r_070_2;
   FReal2 r_072_1 = Ry * r_062_2 - 6 * r_052_2;
   FReal2 r_054_1 = Rz * r_053_2 - 3 * r_052_2;
   FReal2 r_036_1 = Rz * r_035_2 - 5 * r_034_2;
   FReal2 r_018_1 = Ry * r_008_2;
   FReal2 r_801_1 = Rz * r_800_2;
   FReal2 r_603_1 = Rx * r_503_2 - 5 * r_403_2;
   FReal2 r_441_1 = Rz * r_440_2;
   FReal2 r_423_1 = Ry * r_413_2 - r_403_2;
   FReal2 r_405_1 = Rx * r_305_2 - 3 * r_205_2;
   FReal2 r_261_1 = Rz * r_260_2;
   FReal2 r_207_1 = Rz * r_206_2 - 6 * r_205_2;
   FReal2 r_081_1 = Rz * r_080_2;
   FReal2 r_063_1 = Rz * r_062_2 - 2 * r_061_2;
   FReal2 r_045_1 = Ry * r_035_2 - 3 * r_025_2;
   FReal2 r_027_1 = Rz * r_026_2 - 6 * r_025_2;
   FReal2 r_009_1 = Rz * r_008_2 - 8 * r_007_2;
   FReal2 r_530_1 = Ry * r_520_2 - 2 * r_510_2;
   FReal2 r_350_1 = Rx * r_250_2 - 2 * r_150_2;
   FReal2 r_134_1 = Rx * r_034_2;
   FReal2 r_503_1 = Rz * r_502_2 - 2 * r_501_2;
   FReal2 r_341_1 = Rz * r_340_2;
   FReal2 r_305_1 = Rx * r_205_2 - 2 * r_105_2;
   FReal2 r_413_1 = Ry * r_403_2;
   FReal2 r_053_1 = Rz * r_052_2 - 2 * r_051_2;
   FReal2 r_035_1 = Ry * r_025_2 - 2 * r_015_2;
   FReal2 r_531_1 = Rz * r_530_2;
   FReal2 r_513_1 = Ry * r_503_2;
   FReal2 r_351_1 = Rz * r_350_2;
   FReal2 r_315_1 = Ry * r_305_2;
   FReal2 r_153_1 = Rx * r_053_2;
   FReal2 r_135_1 = Rx * r_035_2;

   // m = 0
   pOut[0] = Rx * r_900_1 - 9 * r_800_1;
   pOut[1] = Ry * r_810_1 - r_800_1;
   pOut[2] = Rx * r_702_1 - 7 * r_602_1;
   pOut[3] = Rx * r_540_1 - 5 * r_440_1;
   pOut[4] = Ry * r_612_1 - r_602_1;
   pOut[5] = Rx * r_504_1 - 5 * r_404_1;
   pOut[6] = Ry * r_450_1 - 5 * r_440_1;
   pOut[7] = Rz * r_441_1 - r_440_1;
   pOut[8] = Ry * r_414_1 - r_404_1;
   pOut[9] = Rz * r_405_1 - 5 * r_404_1;
   pOut[10] = Rx * r_180_1 - r_080_1;
   pOut[11] = Rz * r_261_1 - r_260_1;
   pOut[12] = Rx * r_144_1 - r_044_1;
   pOut[13] = Rx * r_126_1 - r_026_1;
   pOut[14] = Rx * r_108_1 - r_008_1;
   pOut[15] = Ry * r_090_1 - 9 * r_080_1;
   pOut[16] = Rz * r_081_1 - r_080_1;
   pOut[17] = Ry * r_054_1 - 5 * r_044_1;
   pOut[18] = Rz * r_045_1 - 5 * r_044_1;
   pOut[19] = Rz * r_027_1 - 7 * r_026_1;
   pOut[20] = Rz * r_009_1 - 9 * r_008_1;
   pOut[21] = Ry * r_900_1;
   pOut[22] = Rx * r_630_1 - 6 * r_530_1;
   pOut[23] = Ry * r_702_1;
   pOut[24] = Rx * r_450_1 - 4 * r_350_1;
   pOut[25] = Rz * r_531_1 - r_530_1;
   pOut[26] = Ry * r_504_1;
   pOut[27] = Ry * r_360_1 - 6 * r_350_1;
   pOut[28] = Rz * r_351_1 - r_350_1;
   pOut[29] = Rx * r_234_1 - 2 * r_134_1;
   pOut[30] = Ry * r_306_1;
   pOut[31] = Rx * r_090_1;
   pOut[32] = Rx * r_072_1;
   pOut[33] = Ry * r_144_1 - 4 * r_134_1;
   pOut[34] = Rz * r_135_1 - 5 * r_134_1;
   pOut[35] = Rx * r_018_1;
   pOut[36] = Rz * r_900_1;
   pOut[37] = Rz * r_720_1;
   pOut[38] = Rx * r_603_1 - 6 * r_503_1;
   pOut[39] = Rx * r_441_1 - 4 * r_341_1;
   pOut[40] = Ry * r_513_1 - r_503_1;
   pOut[41] = Rz * r_504_1 - 4 * r_503_1;
   pOut[42] = Ry * r_351_1 - 5 * r_341_1;
   pOut[43] = Rz * r_342_1 - 2 * r_341_1;
   pOut[44] = Ry * r_315_1 - r_305_1;
   pOut[45] = Rz * r_306_1 - 6 * r_305_1;
   pOut[46] = Rz * r_180_1;
   pOut[47] = Rx * r_063_1;
   pOut[48] = Rx * r_045_1;
   pOut[49] = Rx * r_027_1;
   pOut[50] = Rx * r_009_1;
   pOut[51] = Ry * r_801_1;
   pOut[52] = Rz * r_630_1;
   pOut[53] = Rx * r_513_1 - 5 * r_413_1;
   pOut[54] = Rz * r_450_1;
   pOut[55] = Ry * r_423_1 - 2 * r_413_1;
   pOut[56] = Rz * r_414_1 - 4 * r_413_1;
   pOut[57] = Rz * r_270_1;
   pOut[58] = Rx * r_153_1 - r_053_1;
   pOut[59] = Rx * r_135_1 - r_035_1;
   pOut[60] = Ry * r_207_1;
   pOut[61] = Rz * r_090_1;
   pOut[62] = Ry * r_063_1 - 6 * r_053_1;
   pOut[63] = Ry * r_045_1 - 4 * r_035_1;
   pOut[64] = Rz * r_036_1 - 6 * r_035_1;
   pOut[65] = Ry * r_009_1;
   // 419 flops, 765 mops   (6.35 flops/integral, 11.59 mops/integral)
}

// In: Gm[0 .. 11] (inclusive), Out: [r]^0, ordered as AngularComps(11)
// (moment 11, 78 entries)
void Mdrr11(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 10
   FReal2 r_100_10 = Rx * pIn[11];
   FReal2 r_010_10 = Ry * pIn[11];
   FReal2 r_001_10 = Rz * pIn[11];

   // m = 9
   FReal2 r_200_9 = Rx * r_100_10 - pIn[10];
   FReal2 r_020_9 = Ry * r_010_10 - pIn[10];
   FReal2 r_002_9 = Rz * r_001_10 - pIn[10];
   FReal2 r_100_9 = Rx * pIn[10];
   FReal2 r_010_9 = Ry * pIn[10];
   FReal2 r_001_9 = Rz * pIn[10];

   // m = 8
   FReal2 r_200_8 = Rx * r_100_9 - pIn[9];
   FReal2 r_020_8 = Ry * r_010_9 - pIn[9];
   FReal2 r_002_8 = Rz * r_001_9 - pIn[9];
   FReal2 r_300_8 = Rx * r_200_9 - 2 * r_100_9;
   FReal2 r_100_8 = Rx * pIn[9];
   FReal2 r_030_8 = Ry * r_020_9 - 2 * r_010_9;
   FReal2 r_010_8 = Ry * pIn[9];
   FReal2 r_003_8 = Rz * r_002_9 - 2 * r_001_9;
   FReal2 r_001_8 = Rz * pIn[9];

   // m = 7
   FReal2 r_400_7 = Rx * r_300_8 - 3 * r_200_8;
   FReal2 r_040_7 = Ry * r_030_8 - 3 * r_020_8;
   FReal2 r_004_7 = Rz * r_003_8 - 3 * r_002_8;
   FReal2 r_200_7 = Rx * r_100_8 - pIn[8];
   FReal2 r_020_7 = Ry * r_010_8 - pIn[8];
   FReal2 r_002_7 = Rz * r_001_8 - pIn[8];
   FReal2 r_300_7 = Rx * r_200_8 - 2 * r_100_8;
   FReal2 r_100_7 = Rx * pIn[8];
   FReal2 r_030_7 = Ry * r_020_8 - 2 * r_010_8;
   FReal2 r_010_7 = Ry * pIn[8];
   FReal2 r_003_7 = Rz * r_002_8 - 2 * r_001_8;
   FReal2 r_001_7 = Rz * pIn[8];

   // m = 6
   FReal2 r_400_6 = Rx * r_300_7 - 3 * r_200_7;
   FReal2 r_040_6 = Ry * r_030_7 - 3 * r_020_7;
   FReal2 r_004_6 = Rz * r_003_7 - 3 * r_002_7;
   FReal2 r_200_6 = Rx * r_100_7 - pIn[7];
   FReal2 r_020_6 = Ry * r_010_7 - pIn[7];
   FReal2 r_002_6 = Rz * r_001_7 - pIn[7];
   FReal2 r_500_6 = Rx * r_400_7 - 4 * r_300_7;
   FReal2 r_140_6 = Rx * r_040_7;
   FReal2 r_300_6 = Rx * r_200_7 - 2 * r_100_7;
   FReal2 r_100_6 = Rx * pIn[7];
   FReal2 r_050_6 = Ry * r_040_7 - 4 * r_030_7;
   FReal2 r_014_6 = Ry * r_004_7;
   FReal2 r_030_6 = Ry * r_020_7 - 2 * r_010_7;
   FReal2 r_010_6 = Ry * pIn[7];
   FReal2 r_401_6 = Rz * r_400_7;
   FReal2 r_005_6 = Rz * r_004_7 - 4 * r_003_7;
   FReal2 r_003_6 = Rz * r_002_7 - 2 * r_001_7;
   FReal2 r_001_6 = Rz * pIn[7];

   // m = 5
   FReal2 r_600_5 = Rx * r_500_6 - 5 * r_400_6;
   FReal2 r_402_5 = Rz * r_401_6 - r_400_6;
   FReal2 r_240_5 = Rx * r_140_6 - r_040_6;
   FReal2 r_060_5 = Ry * r_050_6 - 5 * r_040_6;
   FReal2 r_024_5 = Ry * r_014_6 - r_004_6;
   FReal2 r_006_5 = Rz * r_005_6 - 5 * r_004_6;
   FReal2 r_400_5 = Rx * r_300_6 - 3 * r_200_6;
   FReal2 r_040_5 = Ry * r_030_6 - 3 * r_020_6;
   FReal2 r_004_5 = Rz * r_003_6 - 3 * r_002_6;
   FReal2 r_200_5 = Rx * r_100_6 - pIn[6];
   FReal2 r_020_5 = Ry * r_010_6 - pIn[6];
   FReal2 r_002_5 = Rz * r_001_6 - pIn[6];
   FReal2 r_500_5 = Rx * r_400_6 - 4 * r_300_6;
   FReal2 r_140_5 = Rx * r_040_6;
   FReal2 r_300_5 = Rx * r_200_6 - 2 * r_100_6;
   FReal2 r_100_5 = Rx * pIn[6];
   FReal2 r_050_5 = Ry * r_040_6 - 4 * r_030_6;
   FReal2 r_014_5 = Ry * r_004_6;
   FReal2 r_030_5 = Ry * r_020_6 - 2 * r_010_6;
   FReal2 r_010_5 = Ry * pIn[6];
   FReal2 r_401_5 = Rz * r_400_6;
   FReal2 r_005_5 = Rz * r_004_6 - 4 * r_003_6;
   FReal2 r_003_5 = Rz * r_002_6 - 2 * r_001_6;
   FReal2 r_001_5 = Rz * pIn[6];
   FReal2 r_150_5 = Rx * r_050_6;
   FReal2 r_501_5 = Rz * r_500_6;
   FReal2 r_015_5 = Ry * r_005_6;

   // m = 4
   FReal2 r_600_4 = Rx * r_500_5 - 5 * r_400_5;
   FReal2 r_402_4 = Rz * r_401_5 - r_400_5;
   FReal2 r_240_4 = Rx * r_140_5 - r_040_5;
   FReal2 r_060_4 = Ry * r_050_5 - 5 * r_040_5;
   FReal2 r_024_4 = Ry * r_014_5 - r_004_5;
   FReal2 r_006_4 = Rz * r_005_5 - 5 * r_004_5;
   FReal2 r_400_4 = Rx * r_300_5 - 3 * r_200_5;
   FReal2 r_040_4 = Ry * r_030_5 - 3 * r_020_5;
   FReal2 r_004_4 = Rz * r_003_5 - 3 * r_002_5;
   FReal2 r_700_4 = Rx * r_600_5 - 6 * r_500_5;
   FReal2 r_502_4 = Rz * r_501_5 - r_500_5;
   FReal2 r_340_4 = Rx * r_240_5 - 2 * r_140_5;
   FReal2 r_106_4 = Rx * r_006_5;
   FReal2 r_500_4 = Rx * r_400_5 - 4 * r_300_5;
   FReal2 r_140_4 = Rx * r_040_5;
   FReal2 r_300_4 = Rx * r_200_5 - 2 * r_100_5;
   FReal2 r_610_4 = Ry * r_600_5;
   FReal2 r_250_4 = Rx * r_150_5 - r_050_5;
   FReal2 r_070_4 = Ry * r_060_5 - 6 * r_050_5;
   FReal2 r_034_4 = Ry * r_024_5 - 2 * r_014_5;
   FReal2 r_050_4 = Ry * r_040_5 - 4 * r_030_5;
   FReal2 r_014_4 = Ry * r_004_5;
   FReal2 r_030_4 = Ry * r_020_5 - 2 * r_010_5;
   FReal2 r_403_4 = Rz * r_402_5 - 2 * r_401_5;
   FReal2 r_061_4 = Rz * r_060_5;
   FReal2 r_025_4 = Ry * r_015_5 - r_005_5;
   FReal2 r_007_4 = Rz * r_006_5 - 6 * r_005_5;
   FReal2 r_401_4 = Rz * r_400_5;
   FReal2 r_005_4 = Rz * r_004_5 - 4 * r_003_5;
   FReal2 r_003_4 = Rz * r_002_5 - 2 * r_001_5;
   FReal2 r_510_4 = Ry * r_500_5;
   FReal2 r_150_4 = Rx * r_050_5;
   FReal2 r_501_4 = Rz * r_500_5;
   FReal2 r_105_4 = Rx * r_005_5;
   FReal2 r_051_4 = Rz * r_050_5;
   FReal2 r_015_4 = Ry * r_005_5;

   // m = 3
   FReal2 r_800_3 = Rx * r_700_4 - 7 * r_600_4;
   FReal2 r_620_3 = Ry * r_610_4 - r_600_4;
   FReal2 r_602_3 = Rx * r_502_4 - 5 * r_402_4;
   FReal2 r_440_3 = Rx * r_340_4 - 3 * r_240_4;
   FReal2 r_404_3 = Rz * r_403_4 - 3 * r_402_4;
   FReal2 r_260_3 = Ry * r_250_4 - 5 * r_240_4;
   FReal2 r_206_3 = Rx * r_106_4 - r_006_4;
   FReal2 r_080_3 = Ry * r_070_4 - 7 * r_060_4;
   FReal2 r_062_3 = Rz * r_061_4 - r_060_4;
   FReal2 r_044_3 = Ry * r_034_4 - 3 * r_024_4;
   FReal2 r_026_3 = Rz * r_025_4 - 5 * r_024_4;
   FReal2 r_008_3 = Rz * r_007_4 - 7 * r_006_4;
   FReal2 r_600_3 = Rx * r_500_4 - 5 * r_400_4;
   FReal2 r_402_3 = Rz * r_401_4 - r_400_4;
   FReal2 r_240_3 = Rx * r_140_4 - r_040_4;
   FReal2 r_060_3 = Ry * r_050_4 - 5 * r_040_4;
   FReal2 r_024_3 = Ry * r_014_4 - r_004_4;
   FReal2 r_006_3 = Rz * r_005_4 - 5 * r_004_4;
   FReal2 r_700_3 = Rx * r_600_4 - 6 * r_500_4;
   FReal2 r_520_3 = Ry * r_510_4 - r_500_4;
   FReal2 r_502_3 = Rz * r_501_4 - r_500_4;
   FReal2 r_340_3 = Rx * r_240_4 - 2 * r_140_4;
   FReal2 r_106_3 = Rx * r_006_4;
   FReal2 r_500_3 = Rx * r_400_4 - 4 * r_300_4;
   FReal2 r_610_3 = Ry * r_600_4;
   FReal2 r_250_3 = Rx * r_150_4 - r_050_4;
   FReal2 r_070_3 = Ry * r_060_4 - 6 * r_050_4;
   FReal2 r_052_3 = Rz * r_051_4 - r_050_4;
   FReal2 r_034_3 = Ry * r_024_4 - 2 * r_014_4;
   FReal2 r_050_3 = Ry * r_040_4 - 4 * r_030_4;
   FReal2 r_403_3 = Rz * r_402_4 - 2 * r_401_4;
   FReal2 r_205_3 = Rx * r_105_4 - r_005_4;
   FReal2 r_061_3 = Rz * r_060_4;
   FReal2 r_025_3 = Ry * r_015_4 - r_005_4;
   FReal2 r_007_3 = Rz * r_006_4 - 6 * r_005_4;
   FReal2 r_005_3 = Rz * r_004_4 - 4 * r_003_4;
   FReal2 r_350_3 = Rx * r_250_4 - 2 * r_150_4;
   FReal2 r_510_3 = Ry * r_500_4;
   FReal2 r_150_3 = Rx * r_050_4;
   FReal2 r_503_3 = Rz * r_502_4 - 2 * r_501_4;
   FReal2 r_501_3 = Rz * r_500_4;
   FReal2 r_105_3 = Rx * r_005_4;
   FReal2 r_035_3 = Ry * r_025_4 - 2 * r_015_4;
   FReal2 r_051_3 = Rz * r_050_4;
   FReal2 r_015_3 = Ry * r_005_4;

   // m = 2
   FReal2 r_800_2 = Rx * r_700_3 - 7 * r_600_3;
   FReal2 r_620_2 = Ry * r_610_3 - r_600_3;
   FReal2 r_602_2 = Rx * r_502_3 - 5 * r_402_3;
   FReal2 r_440_2 = Rx * r_340_3 - 3 * r_240_3;
   FReal2 r_404_2 = Rz * r_403_3 - 3 * r_402_3;
   FReal2 r_260_2 = Ry * r_250_3 - 5 * r_240_3;
   FReal2 r_206_2 = Rx * r_106_3 - r_006_3;
   FReal2 r_080_2 = Ry * r_070_3 - 7 * r_060_3;
   FReal2 r_062_2 = Rz * r_061_3 - r_060_3;
   FReal2 r_044_2 = Ry * r_034_3 - 3 * r_024_3;
   FReal2 r_026_2 = Rz * r_025_3 - 5 * r_024_3;
   FReal2 r_008_2 = Rz * r_007_3 - 7 * r_006_3;
   FReal2 r_900_2 = Rx * r_800_3 - 8 * r_700_3;
   FReal2 r_720_2 = Rx * r_620_3 - 6 * r_520_3;
   FReal2 r_702_2 = Rx * r_602_3 - 6 * r_502_3;
   FReal2 r_540_2 = Rx * r_440_3 - 4 * r_340_3;
   FReal2 r_504_2 = Rz * r_503_3 - 3 * r_502_3;
   FReal2 r_360_2 = Ry * r_350_3 - 5 * r_340_3;
   FReal2 r_306_2 = Rx * r_206_3 - 2 * r_106_3;
   FReal2 r_144_2 = Rx * r_044_3;
   FReal2 r_700_2 = Rx * r_600_3 - 6 * r_500_3;
   FReal2 r_520_2 = Ry * r_510_3 - r_500_3;
   FReal2 r_502_2 = Rz * r_501_3 - r_500_3;
   FReal2 r_160_2 = Rx * r_060_3;
   FReal2 r_106_2 = Rx * r_006_3;
   FReal2 r_630_2 = Ry * r_620_3 - 2 * r_610_3;
   FReal2 r_450_2 = Rx * r_350_3 - 3 * r_250_3;
   FReal2 r_414_2 = Ry * r_404_3;
   FReal2 r_270_2 = Ry * r_260_3 - 6 * r_250_3;
   FReal2 r_090_2 = Ry * r_080_3 - 8 * r_070_3;
   FReal2 r_072_2 = Ry * r_062_3 - 6 * r_052_3;
   FReal2 r_054_2 = Ry * r_044_3 - 4 * r_034_3;
   FReal2 r_036_2 = Rz * r_035_3 - 5 * r_034_3;
   FReal2 r_610_2 = Ry * r_600_3;
   FReal2 r_250_2 = Rx * r_150_3 - r_050_3;
   FReal2 r_070_2 = Ry * r_060_3 - 6 * r_050_3;
   FReal2 r_052_2 = Rz * r_051_3 - r_050_3;
   FReal2 r_016_2 = Ry * r_006_3;
   FReal2 r_603_2 = Rx * r_503_3 - 5 * r_403_3;
   FReal2 r_441_2 = Rz * r_440_3;
   FReal2 r_405_2 = Rz * r_404_3 - 4 * r_403_3;
   FReal2 r_207_2 = Rz * r_206_3 - 6 * r_205_3;
   FReal2 r_063_2 = Rz * r_062_3 - 2 * r_061_3;
   FReal2 r_045_2 = Ry * r_035_3 - 3 * r_025_3;
   FReal2 r_027_2 = Rz * r_026_3 - 6 * r_025_3;
   FReal2 r_009_2 = Rz * r_008_3 - 8 * r_007_3;
   FReal2 r_601_2 = Rz * r_600_3;
   FReal2 r_205_2 = Rx * r_105_3 - r_005_3;
   FReal2 r_061_2 = Rz * r_060_3;
   FReal2 r_025_2 = Ry * r_015_3 - r_005_3;
   FReal2 r_007_2 = Rz * r_006_3 - 6 * r_005_3;
   FReal2 r_530_2 = Ry * r_520_3 - 2 * r_510_3;
   FReal2 r_350_2 = Rx * r_250_3 - 2 * r_150_3;
   FReal2 r_503_2 = Rz * r_502_3 - 2 * r_501_3;
   FReal2 r_305_2 = Rx * r_205_3 - 2 * r_105_3;
   FReal2 r_053_2 = Rz * r_052_3 - 2 * r_051_3;
   FReal2 r_035_2 = Ry * r_025_3 - 2 * r_015_3;
   FReal2 r_513_2 = Ry * r_503_3;
   FReal2 r_351_2 = Rz * r_350_3;
   FReal2 r_135_2 = Rx * r_035_3;

   // m = 1
   FReal2 r_a00_1 = Rx * r_900_2 - 9 * r_800_2;
   FReal2 r_820_1 = Rx * r_720_2 - 7 * r_620_2;
   FReal2 r_802_1 = Rx * r_702_2 - 7 * r_602_2;
   FReal2 r_640_1 = Ry * r_630_2 - 3 * r_620_2;
   FReal2 r_604_1 = Rz * r_603_2 - 3 * r_602_2;
   FReal2 r_460_1 = Rx * r_360_2 - 3 * r_260_2;
   FReal2 r_442_1 = Rz * r_441_2 - r_440_2;
   FReal2 r_424_1 = Ry * r_414_2 - r_404_2;
   FReal2 r_406_1 = Rx * r_306_2 - 3 * r_206_2;
   FReal2 r_280_1 = Ry * r_270_2 - 7 * r_260_2;
   FReal2 r_244_1 = Rx * r_144_2 - r_044_2;
   FReal2 r_208_1 = Rz * r_207_2 - 7 * r_206_2;
   FReal2 r_0a0_1 = Ry * r_090_2 - 9 * r_080_2;
   FReal2 r_082_1 = Ry * r_072_2 - 7 * r_062_2;
   FReal2 r_064_1 = Rz * r_063_2 - 3 * r_062_2;
   FReal2 r_046_1 = Ry * r_036_2 - 3 * r_026_2;
   FReal2 r_028_1 = Rz * r_027_2 - 7 * r_026_2;
   FReal2 r_00a_1 = Rz * r_009_2 - 9 * r_008_2;
   FReal2 r_900_1 = Rx * r_800_2 - 8 * r_700_2;
   FReal2 r_702_1 = Rx * r_602_2 - 6 * r_502_2;
   FReal2 r_540_1 = Ry * r_530_2 - 3 * r_520_2;
   FReal2 r_504_1 = Rz * r_503_2 - 3 * r_502_2;
   FReal2 r_360_1 = Rx * r_260_2 - 2 * r_160_2;
   FReal2 r_306_1 = Rx * r_206_2 - 2 * r_106_2;
   FReal2 r_144_1 = Rx * r_044_2;
   FReal2 r_630_1 = Ry * r_620_2 - 2 * r_610_2;
   FReal2 r_450_1 = Rx * r_350_2 - 3 * r_250_2;
   FReal2 r_414_1 = Ry * r_404_2;
   FReal2 r_270_1 = Ry * r_260_2 - 6 * r_250_2;
   FReal2 r_090_1 = Ry * r_080_2 - 8 * r_070_2;
   FReal2 r_054_1 = Rz * r_053_2 - 3 * r_052_2;
   FReal2 r_036_1 = Ry * r_026_2 - 2 * r_016_2;
   FReal2 r_603_1 = Rz * r_602_2 - 2 * r_601_2;
   FReal2 r_441_1 = Rz * r_440_2;
   FReal2 r_405_1 = Rx * r_305_2 - 3 * r_205_2;
   FReal2 r_063_1 = Rz * r_062_2 - 2 * r_061_2;
   FReal2 r_045_1 = Ry * r_035_2 - 3 * r_025_2;
   FReal2 r_027_1 = Rz * r_026_2 - 6 * r_025_2;
   FReal2 r_009_1 = Rz * r_008_2 - 8 * r_007_2;
   FReal2 r_910_1 = Ry * r_900_2;
   FReal2 r_730_1 = Rx * r_630_2 - 6 * r_530_2;
   FReal2 r_712_1 = Ry * r_702_2;
   FReal2 r_514_1 = Ry * r_504_2;
   FReal2 r_370_1 = Ry * r_360_2 - 6 * r_350_2;
   FReal2 r_352_1 = Rz * r_351_2 - r_350_2;
   FReal2 r_316_1 = Ry * r_306_2;
   FReal2 r_190_1 = Rx * r_090_2;
   FReal2 r_154_1 = Rx * r_054_2;
   FReal2 r_136_1 = Rx * r_036_2;
   FReal2 r_901_1 = Rz * r_900_2;
   FReal2 r_703_1 = Rx * r_603_2 - 6 * r_503_2;
   FReal2 r_541_1 = Rz * r_540_2;
   FReal2 r_523_1 = Ry * r_513_2 - r_503_2;
   FReal2 r_361_1 = Rz * r_360_2;
   FReal2 r_307_1 = Rz * r_306_2 - 6 * r_305_2;
   FReal2 r_163_1 = Rx * r_063_2;
   FReal2 r_145_1 = Rx * r_045_2;
   FReal2 r_127_1 = Rx * r_027_2;
   FReal2 r_109_1 = Rx * r_009_2;
   FReal2 r_631_1 = Rz * r_630_2;
   FReal2 r_613_1 = Ry * r_603_2;
   FReal2 r_451_1 = Rz * r_450_2;
   FReal2 r_415_1 = Ry * r_405_2;
   FReal2 r_271_1 = Rz * r_270_2;
   FReal2 r_235_1 = Rx * r_135_2 - r_035_2;
   FReal2 r_091_1 = Rz * r_090_2;
   FReal2 r_073_1 = Ry * r_063_2 - 6 * r_053_2;
   FReal2 r_037_1 = Rz * r_036_2 - 6 * r_035_2;
   FReal2 r_019_1 = Ry * r_009_2;
   FReal2 r_513_1 = Ry * r_503_2;
   FReal2 r_351_1 = Rz * r_350_2;
   FReal2 r_135_1 = Rx * r_035_2;

   // m = 0
   pOut[0] = Rx * r_a00_1 - 10 * r_900_1;
   pOut[1] = Ry * r_910_1 - r_900_1;
   pOut[2] = Rx * r_802_1 - 8 * r_702_1;
   pOut[3] = Rx * r_640_1 - 6 * r_540_1;
   pOut[4] = Ry * r_712_1 - r_702_1;
   pOut[5] = Rx * r_604_1 - 6 * r_504_1;
   pOut[6] = Rx * r_460_1 - 4 * r_360_1;
   pOut[7] = Rz * r_541_1 - r_540_1;
   pOut[8] = Ry * r_514_1 - r_504_1;
   pOut[9] = Rx * r_406_1 - 4 * r_306_1;
   pOut[10] = Ry * r_370_1 - 7 * r_360_1;
   pOut[11] = Rz * r_361_1 - r_360_1;
   pOut[12] = Rx * r_244_1 - 2 * r_144_1;
   pOut[13] = Ry * r_316_1 - r_306_1;
   pOut[14] = Rz * r_307_1 - 7 * r_306_1;
   pOut[15] = Rx * r_0a0_1;
   pOut[16] = Rx * r_082_1;
   pOut[17] = Ry * r_154_1 - 5 * r_144_1;
   pOut[18] = Rz * r_145_1 - 5 * r_144_1;
   pOut[19] = Rx * r_028_1;
   pOut[20] = Rx * r_00a_1;
   pOut[21] = Ry * r_a00_1;
   pOut[22] = Rx * r_730_1 - 7 * r_630_1;
   pOut[23] = Ry * r_802_1;
   pOut[24] = Ry * r_640_1 - 4 * r_630_1;
   pOut[25] = Rz * r_631_1 - r_630_1;
   pOut[26] = Rx * r_514_1 - 5 * r_414_1;
   pOut[27] = Ry * r_460_1 - 6 * r_450_1;
   pOut[28] = Rz * r_451_1 - r_450_1;
   pOut[29] = Ry * r_424_1 - 2 * r_414_1;
   pOut[30] = Rz * r_415_1 - 5 * r_414_1;
   pOut[31] = Rx * r_190_1 - r_090_1;
   pOut[32] = Rz * r_271_1 - r_270_1;
   pOut[33] = Rx * r_154_1 - r_054_1;
   pOut[34] = Rx * r_136_1 - r_036_1;
   pOut[35] = Ry * r_208_1;
   pOut[36] = Ry * r_0a0_1 - 10 * r_090_1;
   pOut[37] = Rz * r_091_1 - r_090_1;
   pOut[38] = Ry * r_064_1 - 6 * r_054_1;
   pOut[39] = Ry * r_046_1 - 4 * r_036_1;
   pOut[40] = Rz * r_037_1 - 7 * r_036_1;
   pOut[41] = Ry * r_00a_1;
   pOut[42] = Rz * r_a00_1;
   pOut[43] = Rz * r_820_1;
   pOut[44] = Rx * r_703_1 - 7 * r_603_1;
   pOut[45] = Rx * r_541_1 - 5 * r_441_1;
   pOut[46] = Ry * r_613_1 - r_603_1;
   pOut[47] = Rz * r_604_1 - 4 * r_603_1;
   pOut[48] = Ry * r_451_1 - 5 * r_441_1;
   pOut[49] = Rz * r_442_1 - 2 * r_441_1;
   pOut[50] = Ry * r_415_1 - r_405_1;
   pOut[51] = Rz * r_406_1 - 6 * r_405_1;
   pOut[52] = Rz * r_280_1;
   pOut[53] = Rx * r_163_1 - r_063_1;
   pOut[54] = Rx * r_145_1 - r_045_1;
   pOut[55] = Rx * r_127_1 - r_027_1;
   pOut[56] = Rx * r_109_1 - r_009_1;
   pOut[57] = Rz * r_0a0_1;
   pOut[58] = Ry * r_073_1 - 7 * r_063_1;
   pOut[59] = Rz * r_064_1 - 4 * r_063_1;
   pOut[60] = Rz * r_046_1 - 6 * r_045_1;
   pOut[61] = Rz * r_028_1 - 8 * r_027_1;
   pOut[62] = Rz * r_00a_1 - 10 * r_009_1;
   pOut[63] = Ry * r_901_1;
   pOut[64] = Rz * r_730_1;
   pOut[65] = Rx * r_613_1 - 6 * r_513_1;
   pOut[66] = Rx * r_451_1 - 4 * r_351_1;
   pOut[67] = Ry * r_523_1 - 2 * r_513_1;
   pOut[68] = Rz * r_514_1 - 4 * r_513_1;
   pOut[69] = Ry * r_361_1 - 6 * r_351_1;
   pOut[70] = Rz * r_352_1 - 2 * r_351_1;
   pOut[71] = Rx * r_235_1 - 2 * r_135_1;
   pOut[72] = Ry * r_307_1;
   pOut[73] = Rz * r_190_1;
   pOut[74] = Rx * r_073_1;
   pOut[75] = Ry * r_145_1 - 4 * r_135_1;
   pOut[76] = Rz * r_136_1 - 6 * r_135_1;
   pOut[77] = Rx * r_019_1;
   // 557 flops, 996 mops   (7.14 flops/integral, 12.77 mops/integral)
}

// In: Gm[0 .. 12] (inclusive), Out: [r]^0, ordered as AngularComps(12)
// (moment 12, 91 entries)
void Mdrr12(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 11
   FReal2 r_100_11 = Rx * pIn[12];
   FReal2 r_010_11 = Ry * pIn[12];
   FReal2 r_001_11 = Rz * pIn[12];

   // m = 10
   FReal2 r_200_10 = Rx * r_100_11 - pIn[11];
   FReal2 r_020_10 = Ry * r_010_11 - pIn[11];
   FReal2 r_002_10 = Rz * r_001_11 - pIn[11];
   FReal2 r_100_10 = Rx * pIn[11];
   FReal2 r_010_10 = Ry * pIn[11];
   FReal2 r_001_10 = Rz * pIn[11];

   // m = 9
   FReal2 r_200_9 = Rx * r_100_10 - pIn[10];
   FReal2 r_020_9 = Ry * r_010_10 - pIn[10];
   FReal2 r_002_9 = Rz * r_001_10 - pIn[10];
   FReal2 r_300_9 = Rx * r_200_10 - 2 * r_100_10;
   FReal2 r_100_9 = Rx * pIn[10];
   FReal2 r_030_9 = Ry * r_020_10 - 2 * r_010_10;
   FReal2 r_010_9 = Ry * pIn[10];
   FReal2 r_003_9 = Rz * r_002_10 - 2 * r_001_10;
   FReal2 r_001_9 = Rz * pIn[10];

   // m = 8
   FReal2 r_400_8 = Rx * r_300_9 - 3 * r_200_9;
   FReal2 r_040_8 = Ry * r_030_9 - 3 * r_020_9;
   FReal2 r_004_8 = Rz * r_003_9 - 3 * r_002_9;
   FReal2 r_200_8 = Rx * r_100_9 - pIn[9];
   FReal2 r_020_8 = Ry * r_010_9 - pIn[9];
   FReal2 r_002_8 = Rz * r_001_9 - pIn[9];
   FReal2 r_300_8 = Rx * r_200_9 - 2 * r_100_9;
   FReal2 r_100_8 = Rx * pIn[9];
   FReal2 r_030_8 = Ry * r_020_9 - 2 * r_010_9;
   FReal2 r_010_8 = Ry * pIn[9];
   FReal2 r_003_8 = Rz * r_002_9 - 2 * r_001_9;
   FReal2 r_001_8 = Rz * pIn[9];

   // m = 7
   FReal2 r_400_7 = Rx * r_300_8 - 3 * r_200_8;
   FReal2 r_040_7 = Ry * r_030_8 - 3 * r_020_8;
   FReal2 r_004_7 = Rz * r_003_8 - 3 * r_002_8;
   FReal2 r_200_7 = Rx * r_100_8 - pIn[8];
   FReal2 r_020_7 = Ry * r_010_8 - pIn[8];
   FReal2 r_002_7 = Rz * r_001_8 - pIn[8];
   FReal2 r_500_7 = Rx * r_400_8 - 4 * r_300_8;
   FReal2 r_300_7 = Rx * r_200_8 - 2 * r_100_8;
   FReal2 r_100_7 = Rx * pIn[8];
   FReal2 r_050_7 = Ry * r_040_8 - 4 * r_030_8;
   FReal2 r_030_7 = Ry * r_020_8 - 2 * r_010_8;
   FReal2 r_010_7 = Ry * pIn[8];
   FReal2 r_401_7 = Rz * r_400_8;
   FReal2 r_005_7 = Rz * r_004_8 - 4 * r_003_8;
   FReal2 r_003_7 = Rz * r_002_8 - 2 * r_001_8;
   FReal2 r_001_7 = Rz * pIn[8];

   // m = 6
   FReal2 r_600_6 = Rx * r_500_7 - 5 * r_400_7;
   FReal2 r_402_6 = Rz * r_401_7 - r_400_7;
   FReal2 r_060_6 = Ry * r_050_7 - 5 * r_040_7;
   FReal2 r_006_6 = Rz * r_005_7 - 5 * r_004_7;
   FReal2 r_400_6 = Rx * r_300_7 - 3 * r_200_7;
   FReal2 r_040_6 = Ry * r_030_7 - 3 * r_020_7;
   FReal2 r_004_6 = Rz * r_003_7 - 3 * r_002_7;
   FReal2 r_200_6 = Rx * r_100_7 - pIn[7];
   FReal2 r_020_6 = Ry * r_010_7 - pIn[7];
   FReal2 r_002_6 = Rz * r_001_7 - pIn[7];
   FReal2 r_500_6 = Rx * r_400_7 - 4 * r_300_7;
   FReal2 r_300_6 = Rx * r_200_7 - 2 * r_100_7;
   FReal2 r_100_6 = Rx * pIn[7];
   FReal2 r_050_6 = Ry * r_040_7 - 4 * r_030_7;
   FReal2 r_030_6 = Ry * r_020_7 - 2 * r_010_7;
   FReal2 r_010_6 = Ry * pIn[7];
   FReal2 r_401_6 = Rz * r_400_7;
   FReal2 r_005_6 = Rz * r_004_7 - 4 * r_003_7;
   FReal2 r_003_6 = Rz * r_002_7 - 2 * r_001_7;
   FReal2 r_001_6 = Rz * pIn[7];
   FReal2 r_510_6 = Ry * r_500_7;
   FReal2 r_150_6 = Rx * r_050_7;
   FReal2 r_501_6 = Rz * r_500_7;
   FReal2 r_051_6 = Rz * r_050_7;
   FReal2 r_015_6 = Ry * r_005_7;

   // m = 5
   FReal2 r_600_5 = Rx * r_500_6 - 5 * r_400_6;
   FReal2 r_402_5 = Rz * r_401_6 - r_400_6;
   FReal2 r_060_5 = Ry * r_050_6 - 5 * r_040_6;
   FReal2 r_006_5 = Rz * r_005_6 - 5 * r_004_6;
   FReal2 r_400_5 = Rx * r_300_6 - 3 * r_200_6;
   FReal2 r_040_5 = Ry * r_030_6 - 3 * r_020_6;
   FReal2 r_004_5 = Rz * r_003_6 - 3 * r_002_6;
   FReal2 r_200_5 = Rx * r_100_6 - pIn[6];
   FReal2 r_020_5 = Ry * r_010_6 - pIn[6];
   FReal2 r_002_5 = Rz * r_001_6 - pIn[6];
   FReal2 r_700_5 = Rx * r_600_6 - 6 * r_500_6;
   FReal2 r_520_5 = Ry * r_510_6 - r_500_6;
   FReal2 r_502_5 = Rz * r_501_6 - r_500_6;
   FReal2 r_160_5 = Rx * r_060_6;
   FReal2 r_500_5 = Rx * r_400_6 - 4 * r_300_6;
   FReal2 r_300_5 = Rx * r_200_6 - 2 * r_100_6;
   FReal2 r_610_5 = Ry * r_600_6;
   FReal2 r_250_5 = Rx * r_150_6 - r_050_6;
   FReal2 r_070_5 = Ry * r_060_6 - 6 * r_050_6;
   FReal2 r_052_5 = Rz * r_051_6 - r_050_6;
   FReal2 r_016_5 = Ry * r_006_6;
   FReal2 r_050_5 = Ry * r_040_6 - 4 * r_030_6;
   FReal2 r_030_5 = Ry * r_020_6 - 2 * r_010_6;
   FReal2 r_403_5 = Rz * r_402_6 - 2 * r_401_6;
   FReal2 r_061_5 = Rz * r_060_6;
   FReal2 r_025_5 = Ry * r_015_6 - r_005_6;
   FReal2 r_007_5 = Rz * r_006_6 - 6 * r_005_6;
   FReal2 r_401_5 = Rz * r_400_6;
   FReal2 r_005_5 = Rz * r_004_6 - 4 * r_003_6;
   FReal2 r_003_5 = Rz * r_002_6 - 2 * r_001_6;
   FReal2 r_510_5 = Ry * r_500_6;
   FReal2 r_150_5 = Rx * r_050_6;
   FReal2 r_501_5 = Rz * r_500_6;
   FReal2 r_051_5 = Rz * r_050_6;
   FReal2 r_015_5 = Ry * r_005_6;

   // m = 4
   FReal2 r_800_4 = Rx * r_700_5 - 7 * r_600_5;
   FReal2 r_620_4 = Ry * r_610_5 - r_600_5;
   FReal2 r_404_4 = Rz * r_403_5 - 3 * r_402_5;
   FReal2 r_260_4 = Rx * r_160_5 - r_060_5;
   FReal2 r_080_4 = Ry * r_070_5 - 7 * r_060_5;
   FReal2 r_062_4 = Rz * r_061_5 - r_060_5;
   FReal2 r_026_4 = Ry * r_016_5 - r_006_5;
   FReal2 r_008_4 = Rz * r_007_5 - 7 * r_006_5;
   FReal2 r_600_4 = Rx * r_500_5 - 5 * r_400_5;
   FReal2 r_402_4 = Rz * r_401_5 - r_400_5;
   FReal2 r_060_4 = Ry * r_050_5 - 5 * r_040_5;
   FReal2 r_006_4 = Rz * r_005_5 - 5 * r_004_5;
   FReal2 r_400_4 = Rx * r_300_5 - 3 * r_200_5;
   FReal2 r_040_4 = Ry * r_030_5 - 3 * r_020_5;
   FReal2 r_004_4 = Rz * r_003_5 - 3 * r_002_5;
   FReal2 r_700_4 = Rx * r_600_5 - 6 * r_500_5;
   FReal2 r_520_4 = Ry * r_510_5 - r_500_5;
   FReal2 r_502_4 = Rz * r_501_5 - r_500_5;
   FReal2 r_160_4 = Rx * r_060_5;
   FReal2 r_500_4 = Rx * r_400_5 - 4 * r_300_5;
   FReal2 r_610_4 = Ry * r_600_5;
   FReal2 r_250_4 = Rx * r_150_5 - r_050_5;
   FReal2 r_070_4 = Ry * r_060_5 - 6 * r_050_5;
   FReal2 r_052_4 = Rz * r_051_5 - r_050_5;
   FReal2 r_016_4 = Ry * r_006_5;
   FReal2 r_050_4 = Ry * r_040_5 - 4 * r_030_5;
   FReal2 r_403_4 = Rz * r_402_5 - 2 * r_401_5;
   FReal2 r_061_4 = Rz * r_060_5;
   FReal2 r_025_4 = Ry * r_015_5 - r_005_5;
   FReal2 r_007_4 = Rz * r_006_5 - 6 * r_005_5;
   FReal2 r_401_4 = Rz * r_400_5;
   FReal2 r_005_4 = Rz * r_004_5 - 4 * r_003_5;
   FReal2 r_530_4 = Ry * r_520_5 - 2 * r_510_5;
   FReal2 r_350_4 = Rx * r_250_5 - 2 * r_150_5;
   FReal2 r_510_4 = Ry * r_500_5;
   FReal2 r_150_4 = Rx * r_050_5;
   FReal2 r_503_4 = Rz * r_502_5 - 2 * r_501_5;
   FReal2 r_107_4 = Rx * r_007_5;
   FReal2 r_501_4 = Rz * r_500_5;
   FReal2 r_053_4 = Rz * r_052_5 - 2 * r_051_5;
   FReal2 r_035_4 = Ry * r_025_5 - 2 * r_015_5;
   FReal2 r_051_4 = Rz * r_050_5;
   FReal2 r_015_4 = Ry * r_005_5;

   // m = 3
   FReal2 r_800_3 = Rx * r_700_4 - 7 * r_600_4;
   FReal2 r_620_3 = Ry * r_610_4 - r_600_4;
   FReal2 r_404_3 = Rz * r_403_4 - 3 * r_402_4;
   FReal2 r_260_3 = Rx * r_160_4 - r_060_4;
   FReal2 r_080_3 = Ry * r_070_4 - 7 * r_060_4;
   FReal2 r_062_3 = Rz * r_061_4 - r_060_4;
   FReal2 r_026_3 = Ry * r_016_4 - r_006_4;
   FReal2 r_008_3 = Rz * r_007_4 - 7 * r_006_4;
   FReal2 r_600_3 = Rx * r_500_4 - 5 * r_400_4;
   FReal2 r_402_3 = Rz * r_401_4 - r_400_4;
   FReal2 r_060_3 = Ry * r_050_4 - 5 * r_040_4;
   FReal2 r_006_3 = Rz * r_005_4 - 5 * r_004_4;
   FReal2 r_900_3 = Rx * r_800_4 - 8 * r_700_4;
   FReal2 r_540_3 = Ry * r_530_4 - 3 * r_520_4;
   FReal2 r_504_3 = Rz * r_503_4 - 3 * r_502_4;
   FReal2 r_360_3 = Rx * r_260_4 - 2 * r_160_4;
   FReal2 r_108_3 = Rx * r_008_4;
   FReal2 r_700_3 = Rx * r_600_4 - 6 * r_500_4;
   FReal2 r_520_3 = Ry * r_510_4 - r_500_4;
   FReal2 r_502_3 = Rz * r_501_4 - r_500_4;
   FReal2 r_160_3 = Rx * r_060_4;
   FReal2 r_630_3 = Ry * r_620_4 - 2 * r_610_4;
   FReal2 r_450_3 = Rx * r_350_4 - 3 * r_250_4;
   FReal2 r_414_3 = Ry * r_404_4;
   FReal2 r_270_3 = Ry * r_260_4 - 6 * r_250_4;
   FReal2 r_090_3 = Ry * r_080_4 - 8 * r_070_4;
   FReal2 r_072_3 = Ry * r_062_4 - 6 * r_052_4;
   FReal2 r_054_3 = Rz * r_053_4 - 3 * r_052_4;
   FReal2 r_036_3 = Ry * r_026_4 - 2 * r_016_4;
   FReal2 r_610_3 = Ry * r_600_4;
   FReal2 r_250_3 = Rx * r_150_4 - r_050_4;
   FReal2 r_070_3 = Ry * r_060_4 - 6 * r_050_4;
   FReal2 r_052_3 = Rz * r_051_4 - r_050_4;
   FReal2 r_016_3 = Ry * r_006_4;
   FReal2 r_801_3 = Rz * r_800_4;
   FReal2 r_603_3 = Rx * r_503_4 - 5 * r_403_4;
   FReal2 r_405_3 = Rz * r_404_4 - 4 * r_403_4;
   FReal2 r_207_3 = Rx * r_107_4 - r_007_4;
   FReal2 r_063_3 = Rz * r_062_4 - 2 * r_061_4;
   FReal2 r_045_3 = Ry * r_035_4 - 3 * r_025_4;
   FReal2 r_027_3 = Rz * r_026_4 - 6 * r_025_4;
   FReal2 r_009_3 = Rz * r_008_4 - 8 * r_007_4;
   FReal2 r_403_3 = Rz * r_402_4 - 2 * r_401_4;
   FReal2 r_061_3 = Rz * r_060_4;
   FReal2 r_025_3 = Ry * r_015_4 - r_005_4;
   FReal2 r_007_3 = Rz * r_006_4 - 6 * r_005_4;
   FReal2 r_530_3 = Ry * r_520_4 - 2 * r_510_4;
   FReal2 r_350_3 = Rx * r_250_4 - 2 * r_150_4;
   FReal2 r_150_3 = Rx * r_050_4;
   FReal2 r_503_3 = Rz * r_502_4 - 2 * r_501_4;
   FReal2 r_107_3 = Rx * r_007_4;
   FReal2 r_501_3 = Rz * r_500_4;
   FReal2 r_053_3 = Rz * r_052_4 - 2 * r_051_4;
   FReal2 r_035_3 = Ry * r_025_4 - 2 * r_015_4;
   FReal2 r_015_3 = Ry * r_005_4;

   // m = 2
   FReal2 r_a00_2 = Rx * r_900_3 - 9 * r_800_3;
   FReal2 r_802_2 = Rz * r_801_3 - r_800_3;
   FReal2 r_640_2 = Ry * r_630_3 - 3 * r_620_3;
   FReal2 r_604_2 = Rx * r_504_3 - 5 * r_404_3;
   FReal2 r_460_2 = Rx * r_360_3 - 3 * r_260_3;
   FReal2 r_424_2 = Ry * r_414_3 - r_404_3;
   FReal2 r_280_2 = Ry * r_270_3 - 7 * r_260_3;
   FReal2 r_208_2 = Rx * r_108_3 - r_008_3;
   FReal2 r_0a0_2 = Ry * r_090_3 - 9 * r_080_3;
   FReal2 r_082_2 = Ry * r_072_3 - 7 * r_062_3;
   FReal2 r_064_2 = Rz * r_063_3 - 3 * r_062_3;
   FReal2 r_046_2 = Ry * r_036_3 - 3 * r_026_3;
   FReal2 r_028_2 = Rz * r_027_3 - 7 * r_026_3;
   FReal2 r_00a_2 = Rz * r_009_3 - 9 * r_008_3;
   FReal2 r_800_2 = Rx * r_700_3 - 7 * r_600_3;
   FReal2 r_620_2 = Ry * r_610_3 - r_600_3;
   FReal2 r_404_2 = Rz * r_403_3 - 3 * r_402_3;
   FReal2 r_260_2 = Rx * r_160_3 - r_060_3;
   FReal2 r_080_2 = Ry * r_070_3 - 7 * r_060_3;
   FReal2 r_062_2 = Rz * r_061_3 - r_060_3;
   FReal2 r_026_2 = Ry * r_016_3 - r_006_3;
   FReal2 r_008_2 = Rz * r_007_3 - 7 * r_006_3;
   FReal2 r_900_2 = Rx * r_800_3 - 8 * r_700_3;
   FReal2 r_720_2 = Rx * r_620_3 - 6 * r_520_3;
   FReal2 r_540_2 = Ry * r_530_3 - 3 * r_520_3;
   FReal2 r_504_2 = Rz * r_503_3 - 3 * r_502_3;
   FReal2 r_360_2 = Rx * r_260_3 - 2 * r_160_3;
   FReal2 r_108_2 = Rx * r_008_3;
   FReal2 r_630_2 = Ry * r_620_3 - 2 * r_610_3;
   FReal2 r_450_2 = Rx * r_350_3 - 3 * r_250_3;
   FReal2 r_414_2 = Ry * r_404_3;
   FReal2 r_270_2 = Ry * r_260_3 - 6 * r_250_3;
   FReal2 r_090_2 = Ry * r_080_3 - 8 * r_070_3;
   FReal2 r_072_2 = Ry * r_062_3 - 6 * r_052_3;
   FReal2 r_054_2 = Rz * r_053_3 - 3 * r_052_3;
   FReal2 r_036_2 = Ry * r_026_3 - 2 * r_016_3;
   FReal2 r_801_2 = Rz * r_800_3;
   FReal2 r_603_2 = Rx * r_503_3 - 5 * r_403_3;
   FReal2 r_405_2 = Rz * r_404_3 - 4 * r_403_3;
   FReal2 r_207_2 = Rx * r_107_3 - r_007_3;
   FReal2 r_063_2 = Rz * r_062_3 - 2 * r_061_3;
   FReal2 r_045_2 = Ry * r_035_3 - 3 * r_025_3;
   FReal2 r_027_2 = Rz * r_026_3 - 6 * r_025_3;
   FReal2 r_009_2 = Rz * r_008_3 - 8 * r_007_3;
   FReal2 r_910_2 = Ry * r_900_3;
   FReal2 r_730_2 = Rx * r_630_3 - 6 * r_530_3;
   FReal2 r_550_2 = Rx * r_450_3 - 4 * r_350_3;
   FReal2 r_514_2 = Ry * r_504_3;
   FReal2 r_370_2 = Ry * r_360_3 - 6 * r_350_3;
   FReal2 r_154_2 = Rx * r_054_3;
   FReal2 r_136_2 = Rx * r_036_3;
   FReal2 r_710_2 = Ry * r_700_3;
   FReal2 r_350_2 = Rx * r_250_3 - 2 * r_150_3;
   FReal2 r_901_2 = Rz * r_900_3;
   FReal2 r_703_2 = Rx * r_603_3 - 6 * r_503_3;
   FReal2 r_541_2 = Rz * r_540_3;
   FReal2 r_505_2 = Rz * r_504_3 - 4 * r_503_3;
   FReal2 r_361_2 = Rz * r_360_3;
   FReal2 r_307_2 = Rx * r_207_3 - 2 * r_107_3;
   FReal2 r_145_2 = Rx * r_045_3;
   FReal2 r_503_2 = Rz * r_502_3 - 2 * r_501_3;
   FReal2 r_107_2 = Rx * r_007_3;
   FReal2 r_613_2 = Ry * r_603_3;
   FReal2 r_451_2 = Rz * r_450_3;
   FReal2 r_415_2 = Ry * r_405_3;
   FReal2 r_073_2 = Ry * r_063_3 - 6 * r_053_3;
   FReal2 r_055_2 = Ry * r_045_3 - 4 * r_035_3;
   FReal2 r_037_2 = Rz * r_036_3 - 6 * r_035_3;
   FReal2 r_071_2 = Rz * r_070_3;
   FReal2 r_035_2 = Ry * r_025_3 - 2 * r_015_3;

   // m = 1
   FReal2 r_a00_1 = Rx * r_900_2 - 9 * r_800_2;
   FReal2 r_802_1 = Rz * r_801_2 - r_800_2;
   FReal2 r_640_1 = Ry * r_630_2 - 3 * r_620_2;
   FReal2 r_604_1 = Rx * r_504_2 - 5 * r_404_2;
   FReal2 r_460_1 = Rx * r_360_2 - 3 * r_260_2;
   FReal2 r_424_1 = Ry * r_414_2 - r_404_2;
   FReal2 r_280_1 = Ry * r_270_2 - 7 * r_260_2;
   FReal2 r_208_1 = Rx * r_108_2 - r_008_2;
   FReal2 r_0a0_1 = Ry * r_090_2 - 9 * r_080_2;
   FReal2 r_064_1 = Rz * r_063_2 - 3 * r_062_2;
   FReal2 r_046_1 = Ry * r_036_2 - 3 * r_026_2;
   FReal2 r_028_1 = Rz * r_027_2 - 7 * r_026_2;
   FReal2 r_00a_1 = Rz * r_009_2 - 9 * r_008_2;
   FReal2 r_b00_1 = Rx * r_a00_2 - 10 * r_900_2;
   FReal2 r_920_1 = Ry * r_910_2 - r_900_2;
   FReal2 r_902_1 = Rz * r_901_2 - r_900_2;
   FReal2 r_740_1 = Ry * r_730_2 - 3 * r_720_2;
   FReal2 r_704_1 = Rx * r_604_2 - 6 * r_504_2;
   FReal2 r_560_1 = Rx * r_460_2 - 4 * r_360_2;
   FReal2 r_542_1 = Rz * r_541_2 - r_540_2;
   FReal2 r_524_1 = Ry * r_514_2 - r_504_2;
   FReal2 r_380_1 = Ry * r_370_2 - 7 * r_360_2;
   FReal2 r_362_1 = Rz * r_361_2 - r_360_2;
   FReal2 r_308_1 = Rx * r_208_2 - 2 * r_108_2;
   FReal2 r_1a0_1 = Rx * r_0a0_2;
   FReal2 r_164_1 = Rx * r_064_2;
   FReal2 r_146_1 = Rx * r_046_2;
   FReal2 r_128_1 = Rx * r_028_2;
   FReal2 r_a10_1 = Ry * r_a00_2;
   FReal2 r_830_1 = Rx * r_730_2 - 7 * r_630_2;
   FReal2 r_812_1 = Ry * r_802_2;
   FReal2 r_452_1 = Rz * r_451_2 - r_450_2;
   FReal2 r_434_1 = Ry * r_424_2 - 2 * r_414_2;
   FReal2 r_290_1 = Ry * r_280_2 - 8 * r_270_2;
   FReal2 r_254_1 = Rx * r_154_2 - r_054_2;
   FReal2 r_236_1 = Rx * r_136_2 - r_036_2;
   FReal2 r_0b0_1 = Ry * r_0a0_2 - 10 * r_090_2;
   FReal2 r_092_1 = Ry * r_082_2 - 8 * r_072_2;
   FReal2 r_074_1 = Rz * r_073_2 - 3 * r_072_2;
   FReal2 r_056_1 = Ry * r_046_2 - 4 * r_036_2;
   FReal2 r_038_1 = Rz * r_037_2 - 7 * r_036_2;
   FReal2 r_01a_1 = Ry * r_00a_2;
   FReal2 r_a01_1 = Rx * r_901_2 - 9 * r_801_2;
   FReal2 r_803_1 = Rz * r_802_2 - 2 * r_801_2;
   FReal2 r_641_1 = Rz * r_640_2;
   FReal2 r_623_1 = Ry * r_613_2 - r_603_2;
   FReal2 r_605_1 = Rx * r_505_2 - 5 * r_405_2;
   FReal2 r_461_1 = Rz * r_460_2;
   FReal2 r_425_1 = Ry * r_415_2 - r_405_2;
   FReal2 r_407_1 = Rx * r_307_2 - 3 * r_207_2;
   FReal2 r_281_1 = Rz * r_280_2;
   FReal2 r_245_1 = Rx * r_145_2 - r_045_2;
   FReal2 r_209_1 = Rz * r_208_2 - 8 * r_207_2;
   FReal2 r_0a1_1 = Rz * r_0a0_2;
   FReal2 r_083_1 = Ry * r_073_2 - 7 * r_063_2;
   FReal2 r_029_1 = Rz * r_028_2 - 8 * r_027_2;
   FReal2 r_00b_1 = Rz * r_00a_2 - 10 * r_009_2;
   FReal2 r_730_1 = Ry * r_720_2 - 2 * r_710_2;
   FReal2 r_550_1 = Rx * r_450_2 - 4 * r_350_2;
   FReal2 r_514_1 = Ry * r_504_2;
   FReal2 r_370_1 = Ry * r_360_2 - 6 * r_350_2;
   FReal2 r_154_1 = Rx * r_054_2;
   FReal2 r_136_1 = Rx * r_036_2;
   FReal2 r_703_1 = Rx * r_603_2 - 6 * r_503_2;
   FReal2 r_541_1 = Rz * r_540_2;
   FReal2 r_505_1 = Rz * r_504_2 - 4 * r_503_2;
   FReal2 r_361_1 = Rz * r_360_2;
   FReal2 r_307_1 = Rx * r_207_2 - 2 * r_107_2;
   FReal2 r_145_1 = Rx * r_045_2;
   FReal2 r_613_1 = Ry * r_603_2;
   FReal2 r_451_1 = Rz * r_450_2;
   FReal2 r_415_1 = Ry * r_405_2;
   FReal2 r_073_1 = Rz * r_072_2 - 2 * r_071_2;
   FReal2 r_055_1 = Ry * r_045_2 - 4 * r_035_2;
   FReal2 r_037_1 = Rz * r_036_2 - 6 * r_035_2;
   FReal2 r_731_1 = Rz * r_730_2;
   FReal2 r_713_1 = Ry * r_703_2;
   FReal2 r_551_1 = Rz * r_550_2;
   FReal2 r_515_1 = Ry * r_505_2;
   FReal2 r_371_1 = Rz * r_370_2;
   FReal2 r_317_1 = Ry * r_307_2;
   FReal2 r_173_1 = Rx * r_073_2;
   FReal2 r_155_1 = Rx * r_055_2;
   FReal2 r_137_1 = Rx * r_037_2;

   // m = 0
   pOut[0] = Rx * r_b00_1 - 11 * r_a00_1;
   pOut[1] = Ry * r_a10_1 - r_a00_1;
   pOut[2] = Rx * r_902_1 - 9 * r_802_1;
   pOut[3] = Rx * r_740_1 - 7 * r_640_1;
   pOut[4] = Ry * r_812_1 - r_802_1;
   pOut[5] = Rx * r_704_1 - 7 * r_604_1;
   pOut[6] = Rx * r_560_1 - 5 * r_460_1;
   pOut[7] = Rz * r_641_1 - r_640_1;
   pOut[8] = Rx * r_524_1 - 5 * r_424_1;
   pOut[9] = Rz * r_605_1 - 5 * r_604_1;
   pOut[10] = Rx * r_380_1 - 3 * r_280_1;
   pOut[11] = Rz * r_461_1 - r_460_1;
   pOut[12] = Ry * r_434_1 - 3 * r_424_1;
   pOut[13] = Rz * r_425_1 - 5 * r_424_1;
   pOut[14] = Rx * r_308_1 - 3 * r_208_1;
   pOut[15] = Rx * r_1a0_1 - r_0a0_1;
   pOut[16] = Rz * r_281_1 - r_280_1;
   pOut[17] = Rx * r_164_1 - r_064_1;
   pOut[18] = Rx * r_146_1 - r_046_1;
   pOut[19] = Rx * r_128_1 - r_028_1;
   pOut[20] = Rz * r_209_1 - 9 * r_208_1;
   pOut[21] = Ry * r_0b0_1 - 11 * r_0a0_1;
   pOut[22] = Rz * r_0a1_1 - r_0a0_1;
   pOut[23] = Ry * r_074_1 - 7 * r_064_1;
   pOut[24] = Ry * r_056_1 - 5 * r_046_1;
   pOut[25] = Ry * r_038_1 - 3 * r_028_1;
   pOut[26] = Rz * r_029_1 - 9 * r_028_1;
   pOut[27] = Rz * r_00b_1 - 11 * r_00a_1;
   pOut[28] = Ry * r_b00_1;
   pOut[29] = Rx * r_830_1 - 8 * r_730_1;
   pOut[30] = Ry * r_902_1;
   pOut[31] = Ry * r_740_1 - 4 * r_730_1;
   pOut[32] = Rz * r_731_1 - r_730_1;
   pOut[33] = Ry * r_704_1;
   pOut[34] = Ry * r_560_1 - 6 * r_550_1;
   pOut[35] = Rz * r_551_1 - r_550_1;
   pOut[36] = Ry * r_524_1 - 2 * r_514_1;
   pOut[37] = Rz * r_515_1 - 5 * r_514_1;
   pOut[38] = Ry * r_380_1 - 8 * r_370_1;
   pOut[39] = Rz * r_371_1 - r_370_1;
   pOut[40] = Rx * r_254_1 - 2 * r_154_1;
   pOut[41] = Rx * r_236_1 - 2 * r_136_1;
   pOut[42] = Ry * r_308_1;
   pOut[43] = Rx * r_0b0_1;
   pOut[44] = Rx * r_092_1;
   pOut[45] = Ry * r_164_1 - 6 * r_154_1;
   pOut[46] = Ry * r_146_1 - 4 * r_136_1;
   pOut[47] = Rz * r_137_1 - 7 * r_136_1;
   pOut[48] = Rx * r_01a_1;
   pOut[49] = Rz * r_b00_1;
   pOut[50] = Rz * r_920_1;
   pOut[51] = Rx * r_803_1 - 8 * r_703_1;
   pOut[52] = Rx * r_641_1 - 6 * r_541_1;
   pOut[53] = Ry * r_713_1 - r_703_1;
   pOut[54] = Rz * r_704_1 - 4 * r_703_1;
   pOut[55] = Rx * r_461_1 - 4 * r_361_1;
   pOut[56] = Rz * r_542_1 - 2 * r_541_1;
   pOut[57] = Ry * r_515_1 - r_505_1;
   pOut[58] = Rx * r_407_1 - 4 * r_307_1;
   pOut[59] = Ry * r_371_1 - 7 * r_361_1;
   pOut[60] = Rz * r_362_1 - 2 * r_361_1;
   pOut[61] = Rx * r_245_1 - 2 * r_145_1;
   pOut[62] = Ry * r_317_1 - r_307_1;
   pOut[63] = Rz * r_308_1 - 8 * r_307_1;
   pOut[64] = Rz * r_1a0_1;
   pOut[65] = Rx * r_083_1;
   pOut[66] = Ry * r_155_1 - 5 * r_145_1;
   pOut[67] = Rz * r_146_1 - 6 * r_145_1;
   pOut[68] = Rx * r_029_1;
   pOut[69] = Rx * r_00b_1;
   pOut[70] = Ry * r_a01_1;
   pOut[71] = Rz * r_830_1;
   pOut[72] = Rx * r_713_1 - 7 * r_613_1;
   pOut[73] = Rx * r_551_1 - 5 * r_451_1;
   pOut[74] = Ry * r_623_1 - 2 * r_613_1;
   pOut[75] = Rx * r_515_1 - 5 * r_415_1;
   pOut[76] = Ry * r_461_1 - 6 * r_451_1;
   pOut[77] = Rz * r_452_1 - 2 * r_451_1;
   pOut[78] = Ry * r_425_1 - 2 * r_415_1;
   pOut[79] = Ry * r_407_1;
   pOut[80] = Rz * r_290_1;
   pOut[81] = Rx * r_173_1 - r_073_1;
   pOut[82] = Rx * r_155_1 - r_055_1;
   pOut[83] = Rx * r_137_1 - r_037_1;
   pOut[84] = Ry * r_209_1;
   pOut[85] = Rz * r_0b0_1;
   pOut[86] = Ry * r_083_1 - 8 * r_073_1;
   pOut[87] = Rz * r_074_1 - 4 * r_073_1;
   pOut[88] = Rz * r_056_1 - 6 * r_055_1;
   pOut[89] = Rz * r_038_1 - 8 * r_037_1;
   pOut[90] = Ry * r_00b_1;
   // 688 flops, 1227 mops   (7.56 flops/integral, 13.48 mops/integral)
}

// In: Gm[0 .. 13] (inclusive), Out: [r]^0, ordered as AngularComps(13)
// (moment 13, 105 entries)
void Mdrr13(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 12
   FReal2 r_100_12 = Rx * pIn[13];
   FReal2 r_010_12 = Ry * pIn[13];
   FReal2 r_001_12 = Rz * pIn[13];

   // m = 11
   FReal2 r_200_11 = Rx * r_100_12 - pIn[12];
   FReal2 r_020_11 = Ry * r_010_12 - pIn[12];
   FReal2 r_002_11 = Rz * r_001_12 - pIn[12];
   FReal2 r_100_11 = Rx * pIn[12];
   FReal2 r_010_11 = Ry * pIn[12];
   FReal2 r_001_11 = Rz * pIn[12];

   // m = 10
   FReal2 r_200_10 = Rx * r_100_11 - pIn[11];
   FReal2 r_020_10 = Ry * r_010_11 - pIn[11];
   FReal2 r_002_10 = Rz * r_001_11 - pIn[11];
   FReal2 r_300_10 = Rx * r_200_11 - 2 * r_100_11;
   FReal2 r_100_10 = Rx * pIn[11];
   FReal2 r_030_10 = Ry * r_020_11 - 2 * r_010_11;
   FReal2 r_010_10 = Ry * pIn[11];
   FReal2 r_003_10 = Rz * r_002_11 - 2 * r_001_11;
   FReal2 r_001_10 = Rz * pIn[11];

   // m = 9
   FReal2 r_400_9 = Rx * r_300_10 - 3 * r_200_10;
   FReal2 r_040_9 = Ry * r_030_10 - 3 * r_020_10;
   FReal2 r_004_9 = Rz * r_003_10 - 3 * r_002_10;
   FReal2 r_200_9 = Rx * r_100_10 - pIn[10];
   FReal2 r_020_9 = Ry * r_010_10 - pIn[10];
   FReal2 r_002_9 = Rz * r_001_10 - pIn[10];
   FReal2 r_300_9 = Rx * r_200_10 - 2 * r_100_10;
   FReal2 r_100_9 = Rx * pIn[10];
   FReal2 r_030_9 = Ry * r_020_10 - 2 * r_010_10;
   FReal2 r_010_9 = Ry * pIn[10];
   FReal2 r_003_9 = Rz * r_002_10 - 2 * r_001_10;
   FReal2 r_001_9 = Rz * pIn[10];

   // m = 8
   FReal2 r_400_8 = Rx * r_300_9 - 3 * r_200_9;
   FReal2 r_040_8 = Ry * r_030_9 - 3 * r_020_9;
   FReal2 r_004_8 = Rz * r_003_9 - 3 * r_002_9;
   FReal2 r_200_8 = Rx * r_100_9 - pIn[9];
   FReal2 r_020_8 = Ry * r_010_9 - pIn[9];
   FReal2 r_002_8 = Rz * r_001_9 - pIn[9];
   FReal2 r_500_8 = Rx * r_400_9 - 4 * r_300_9;
   FReal2 r_300_8 = Rx * r_200_9 - 2 * r_100_9;
   FReal2 r_100_8 = Rx * pIn[9];
   FReal2 r_050_8 = Ry * r_040_9 - 4 * r_030_9;
   FReal2 r_030_8 = Ry * r_020_9 - 2 * r_010_9;
   FReal2 r_010_8 = Ry * pIn[9];
   FReal2 r_005_8 = Rz * r_004_9 - 4 * r_003_9;
   FReal2 r_003_8 = Rz * r_002_9 - 2 * r_001_9;
   FReal2 r_001_8 = Rz * pIn[9];

   // m = 7
   FReal2 r_600_7 = Rx * r_500_8 - 5 * r_400_8;
   FReal2 r_060_7 = Ry * r_050_8 - 5 * r_040_8;
   FReal2 r_006_7 = Rz * r_005_8 - 5 * r_004_8;
   FReal2 r_400_7 = Rx * r_300_8 - 3 * r_200_8;
   FReal2 r_040_7 = Ry * r_030_8 - 3 * r_020_8;
   FReal2 r_004_7 = Rz * r_003_8 - 3 * r_002_8;
   FReal2 r_200_7 = Rx * r_100_8 - pIn[8];
   FReal2 r_020_7 = Ry * r_010_8 - pIn[8];
   FReal2 r_002_7 = Rz * r_001_8 - pIn[8];
   FReal2 r_500_7 = Rx * r_400_8 - 4 * r_300_8;
   FReal2 r_300_7 = Rx * r_200_8 - 2 * r_100_8;
   FReal2 r_100_7 = Rx * pIn[8];
   FReal2 r_050_7 = Ry * r_040_8 - 4 * r_030_8;
   FReal2 r_030_7 = Ry * r_020_8 - 2 * r_010_8;
   FReal2 r_010_7 = Ry * pIn[8];
   FReal2 r_005_7 = Rz * r_004_8 - 4 * r_003_8;
   FReal2 r_003_7 = Rz * r_002_8 - 2 * r_001_8;
   FReal2 r_001_7 = Rz * pIn[8];
   FReal2 r_150_7 = Rx * r_050_8;
   FReal2 r_501_7 = Rz * r_500_8;
   FReal2 r_015_7 = Ry * r_005_8;

   // m = 6
   FReal2 r_600_6 = Rx * r_500_7 - 5 * r_400_7;
   FReal2 r_060_6 = Ry * r_050_7 - 5 * r_040_7;
   FReal2 r_006_6 = Rz * r_005_7 - 5 * r_004_7;
   FReal2 r_400_6 = Rx * r_300_7 - 3 * r_200_7;
   FReal2 r_040_6 = Ry * r_030_7 - 3 * r_020_7;
   FReal2 r_004_6 = Rz * r_003_7 - 3 * r_002_7;
   FReal2 r_200_6 = Rx * r_100_7 - pIn[7];
   FReal2 r_020_6 = Ry * r_010_7 - pIn[7];
   FReal2 r_002_6 = Rz * r_001_7 - pIn[7];
   FReal2 r_700_6 = Rx * r_600_7 - 6 * r_500_7;
   FReal2 r_502_6 = Rz * r_501_7 - r_500_7;
   FReal2 r_160_6 = Rx * r_060_7;
   FReal2 r_106_6 = Rx * r_006_7;
   FReal2 r_500_6 = Rx * r_400_7 - 4 * r_300_7;
   FReal2 r_300_6 = Rx * r_200_7 - 2 * r_100_7;
   FReal2 r_100_6 = Rx * pIn[7];
   FReal2 r_610_6 = Ry * r_600_7;
   FReal2 r_250_6 = Rx * r_150_7 - r_050_7;
   FReal2 r_070_6 = Ry * r_060_7 - 6 * r_050_7;
   FReal2 r_016_6 = Ry * r_006_7;
   FReal2 r_050_6 = Ry * r_040_7 - 4 * r_030_7;
   FReal2 r_030_6 = Ry * r_020_7 - 2 * r_010_7;
   FReal2 r_010_6 = Ry * pIn[7];
   FReal2 r_601_6 = Rz * r_600_7;
   FReal2 r_061_6 = Rz * r_060_7;
   FReal2 r_025_6 = Ry * r_015_7 - r_005_7;
   FReal2 r_007_6 = Rz * r_006_7 - 6 * r_005_7;
   FReal2 r_005_6 = Rz * r_004_7 - 4 * r_003_7;
   FReal2 r_003_6 = Rz * r_002_7 - 2 * r_001_7;
   FReal2 r_001_6 = Rz * pIn[7];
   FReal2 r_150_6 = Rx * r_050_7;
   FReal2 r_501_6 = Rz * r_500_7;
   FReal2 r_015_6 = Ry * r_005_7;

   // m = 5
   FReal2 r_800_5 = Rx * r_700_6 - 7 * r_600_6;
   FReal2 r_620_5 = Ry * r_610_6 - r_600_6;
   FReal2 r_602_5 = Rz * r_601_6 - r_600_6;
   FReal2 r_260_5 = Rx * r_160_6 - r_060_6;
   FReal2 r_206_5 = Rx * r_106_6 - r_006_6;
   FReal2 r_080_5 = Ry * r_070_6 - 7 * r_060_6;
   FReal2 r_062_5 = Rz * r_061_6 - r_060_6;
   FReal2 r_026_5 = Ry * r_016_6 - r_006_6;
   FReal2 r_008_5 = Rz * r_007_6 - 7 * r_006_6;
   FReal2 r_600_5 = Rx * r_500_6 - 5 * r_400_6;
   FReal2 r_060_5 = Ry * r_050_6 - 5 * r_040_6;
   FReal2 r_006_5 = Rz * r_005_6 - 5 * r_004_6;
   FReal2 r_400_5 = Rx * r_300_6 - 3 * r_200_6;
   FReal2 r_040_5 = Ry * r_030_6 - 3 * r_020_6;
   FReal2 r_004_5 = Rz * r_003_6 - 3 * r_002_6;
   FReal2 r_700_5 = Rx * r_600_6 - 6 * r_500_6;
   FReal2 r_502_5 = Rz * r_501_6 - r_500_6;
   FReal2 r_160_5 = Rx * r_060_6;
   FReal2 r_106_5 = Rx * r_006_6;
   FReal2 r_500_5 = Rx * r_400_6 - 4 * r_300_6;
   FReal2 r_300_5 = Rx * r_200_6 - 2 * r_100_6;
   FReal2 r_610_5 = Ry * r_600_6;
   FReal2 r_250_5 = Rx * r_150_6 - r_050_6;
   FReal2 r_070_5 = Ry * r_060_6 - 6 * r_050_6;
   FReal2 r_016_5 = Ry * r_006_6;
   FReal2 r_050_5 = Ry * r_040_6 - 4 * r_030_6;
   FReal2 r_030_5 = Ry * r_020_6 - 2 * r_010_6;
   FReal2 r_601_5 = Rz * r_600_6;
   FReal2 r_061_5 = Rz * r_060_6;
   FReal2 r_025_5 = Ry * r_015_6 - r_005_6;
   FReal2 r_007_5 = Rz * r_006_6 - 6 * r_005_6;
   FReal2 r_005_5 = Rz * r_004_6 - 4 * r_003_6;
   FReal2 r_003_5 = Rz * r_002_6 - 2 * r_001_6;
   FReal2 r_710_5 = Ry * r_700_6;
   FReal2 r_350_5 = Rx * r_250_6 - 2 * r_150_6;
   FReal2 r_150_5 = Rx * r_050_6;
   FReal2 r_503_5 = Rz * r_502_6 - 2 * r_501_6;
   FReal2 r_107_5 = Rx * r_007_6;
   FReal2 r_501_5 = Rz * r_500_6;
   FReal2 r_071_5 = Rz * r_070_6;
   FReal2 r_035_5 = Ry * r_025_6 - 2 * r_015_6;
   FReal2 r_015_5 = Ry * r_005_6;

   // m = 4
   FReal2 r_800_4 = Rx * r_700_5 - 7 * r_600_5;
   FReal2 r_620_4 = Ry * r_610_5 - r_600_5;
   FReal2 r_602_4 = Rz * r_601_5 - r_600_5;
   FReal2 r_260_4 = Rx * r_160_5 - r_060_5;
   FReal2 r_206_4 = Rx * r_106_5 - r_006_5;
   FReal2 r_080_4 = Ry * r_070_5 - 7 * r_060_5;
   FReal2 r_062_4 = Rz * r_061_5 - r_060_5;
   FReal2 r_026_4 = Ry * r_016_5 - r_006_5;
   FReal2 r_008_4 = Rz * r_007_5 - 7 * r_006_5;
   FReal2 r_600_4 = Rx * r_500_5 - 5 * r_400_5;
   FReal2 r_060_4 = Ry * r_050_5 - 5 * r_040_5;
   FReal2 r_006_4 = Rz * r_005_5 - 5 * r_004_5;
   FReal2 r_900_4 = Rx * r_800_5 - 8 * r_700_5;
   FReal2 r_720_4 = Ry * r_710_5 - r_700_5;
   FReal2 r_504_4 = Rz * r_503_5 - 3 * r_502_5;
   FReal2 r_360_4 = Rx * r_260_5 - 2 * r_160_5;
   FReal2 r_306_4 = Rx * r_206_5 - 2 * r_106_5;
   FReal2 r_700_4 = Rx * r_600_5 - 6 * r_500_5;
   FReal2 r_502_4 = Rz * r_501_5 - r_500_5;
   FReal2 r_160_4 = Rx * r_060_5;
   FReal2 r_106_4 = Rx * r_006_5;
   FReal2 r_500_4 = Rx * r_400_5 - 4 * r_300_5;
   FReal2 r_630_4 = Ry * r_620_5 - 2 * r_610_5;
   FReal2 r_450_4 = Rx * r_350_5 - 3 * r_250_5;
   FReal2 r_090_4 = Ry * r_080_5 - 8 * r_070_5;
   FReal2 r_072_4 = Rz * r_071_5 - r_070_5;
   FReal2 r_036_4 = Ry * r_026_5 - 2 * r_016_5;
   FReal2 r_610_4 = Ry * r_600_5;
   FReal2 r_250_4 = Rx * r_150_5 - r_050_5;
   FReal2 r_070_4 = Ry * r_060_5 - 6 * r_050_5;
   FReal2 r_016_4 = Ry * r_006_5;
   FReal2 r_050_4 = Ry * r_040_5 - 4 * r_030_5;
   FReal2 r_603_4 = Rz * r_602_5 - 2 * r_601_5;
   FReal2 r_207_4 = Rx * r_107_5 - r_007_5;
   FReal2 r_063_4 = Rz * r_062_5 - 2 * r_061_5;
   FReal2 r_045_4 = Ry * r_035_5 - 3 * r_025_5;
   FReal2 r_009_4 = Rz * r_008_5 - 8 * r_007_5;
   FReal2 r_601_4 = Rz * r_600_5;
   FReal2 r_061_4 = Rz * r_060_5;
   FReal2 r_025_4 = Ry * r_015_5 - r_005_5;
   FReal2 r_007_4 = Rz * r_006_5 - 6 * r_005_5;
   FReal2 r_005_4 = Rz * r_004_5 - 4 * r_003_5;
   FReal2 r_710_4 = Ry * r_700_5;
   FReal2 r_350_4 = Rx * r_250_5 - 2 * r_150_5;
   FReal2 r_150_4 = Rx * r_050_5;
   FReal2 r_503_4 = Rz * r_502_5 - 2 * r_501_5;
   FReal2 r_107_4 = Rx * r_007_5;
   FReal2 r_501_4 = Rz * r_500_5;
   FReal2 r_071_4 = Rz * r_070_5;
   FReal2 r_035_4 = Ry * r_025_5 - 2 * r_015_5;
   FReal2 r_015_4 = Ry * r_005_5;

   // m = 3
   FReal2 r_a00_3 = Rx * r_900_4 - 9 * r_800_4;
   FReal2 r_640_3 = Ry * r_630_4 - 3 * r_620_4;
   FReal2 r_604_3 = Rz * r_603_4 - 3 * r_602_4;
   FReal2 r_460_3 = Rx * r_360_4 - 3 * r_260_4;
   FReal2 r_406_3 = Rx * r_306_4 - 3 * r_206_4;
   FReal2 r_0a0_3 = Ry * r_090_4 - 9 * r_080_4;
   FReal2 r_064_3 = Rz * r_063_4 - 3 * r_062_4;
   FReal2 r_046_3 = Ry * r_036_4 - 3 * r_026_4;
   FReal2 r_00a_3 = Rz * r_009_4 - 9 * r_008_4;
   FReal2 r_800_3 = Rx * r_700_4 - 7 * r_600_4;
   FReal2 r_620_3 = Ry * r_610_4 - r_600_4;
   FReal2 r_602_3 = Rz * r_601_4 - r_600_4;
   FReal2 r_260_3 = Rx * r_160_4 - r_060_4;
   FReal2 r_206_3 = Rx * r_106_4 - r_006_4;
   FReal2 r_080_3 = Ry * r_070_4 - 7 * r_060_4;
   FReal2 r_062_3 = Rz * r_061_4 - r_060_4;
   FReal2 r_026_3 = Ry * r_016_4 - r_006_4;
   FReal2 r_008_3 = Rz * r_007_4 - 7 * r_006_4;
   FReal2 r_900_3 = Rx * r_800_4 - 8 * r_700_4;
   FReal2 r_720_3 = Ry * r_710_4 - r_700_4;
   FReal2 r_504_3 = Rz * r_503_4 - 3 * r_502_4;
   FReal2 r_360_3 = Rx * r_260_4 - 2 * r_160_4;
   FReal2 r_306_3 = Rx * r_206_4 - 2 * r_106_4;
   FReal2 r_700_3 = Rx * r_600_4 - 6 * r_500_4;
   FReal2 r_502_3 = Rz * r_501_4 - r_500_4;
   FReal2 r_160_3 = Rx * r_060_4;
   FReal2 r_630_3 = Ry * r_620_4 - 2 * r_610_4;
   FReal2 r_450_3 = Rx * r_350_4 - 3 * r_250_4;
   FReal2 r_090_3 = Ry * r_080_4 - 8 * r_070_4;
   FReal2 r_072_3 = Rz * r_071_4 - r_070_4;
   FReal2 r_036_3 = Ry * r_026_4 - 2 * r_016_4;
   FReal2 r_250_3 = Rx * r_150_4 - r_050_4;
   FReal2 r_070_3 = Ry * r_060_4 - 6 * r_050_4;
   FReal2 r_016_3 = Ry * r_006_4;
   FReal2 r_603_3 = Rz * r_602_4 - 2 * r_601_4;
   FReal2 r_207_3 = Rx * r_107_4 - r_007_4;
   FReal2 r_063_3 = Rz * r_062_4 - 2 * r_061_4;
   FReal2 r_045_3 = Ry * r_035_4 - 3 * r_025_4;
   FReal2 r_009_3 = Rz * r_008_4 - 8 * r_007_4;
   FReal2 r_601_3 = Rz * r_600_4;
   FReal2 r_025_3 = Ry * r_015_4 - r_005_4;
   FReal2 r_007_3 = Rz * r_006_4 - 6 * r_005_4;
   FReal2 r_730_3 = Ry * r_720_4 - 2 * r_710_4;
   FReal2 r_550_3 = Rx * r_450_4 - 4 * r_350_4;
   FReal2 r_514_3 = Ry * r_504_4;
   FReal2 r_370_3 = Ry * r_360_4 - 6 * r_350_4;
   FReal2 r_190_3 = Rx * r_090_4;
   FReal2 r_710_3 = Ry * r_700_4;
   FReal2 r_350_3 = Rx * r_250_4 - 2 * r_150_4;
   FReal2 r_901_3 = Rz * r_900_4;
   FReal2 r_703_3 = Rx * r_603_4 - 6 * r_503_4;
   FReal2 r_505_3 = Rz * r_504_4 - 4 * r_503_4;
   FReal2 r_307_3 = Rx * r_207_4 - 2 * r_107_4;
   FReal2 r_145_3 = Rx * r_045_4;
   FReal2 r_503_3 = Rz * r_502_4 - 2 * r_501_4;
   FReal2 r_107_3 = Rx * r_007_4;
   FReal2 r_451_3 = Rz * r_450_4;
   FReal2 r_091_3 = Rz * r_090_4;
   FReal2 r_073_3 = Rz * r_072_4 - 2 * r_071_4;
   FReal2 r_055_3 = Ry * r_045_4 - 4 * r_035_4;
   FReal2 r_037_3 = Rz * r_036_4 - 6 * r_035_4;
   FReal2 r_019_3 = Ry * r_009_4;
   FReal2 r_071_3 = Rz * r_070_4;
   FReal2 r_035_3 = Ry * r_025_4 - 2 * r_015_4;

   // m = 2
   FReal2 r_a00_2 = Rx * r_900_3 - 9 * r_800_3;
   FReal2 r_820_2 = Rx * r_720_3 - 7 * r_620_3;
   FReal2 r_640_2 = Ry * r_630_3 - 3 * r_620_3;
   FReal2 r_604_2 = Rz * r_603_3 - 3 * r_602_3;
   FReal2 r_460_2 = Rx * r_360_3 - 3 * r_260_3;
   FReal2 r_406_2 = Rx * r_306_3 - 3 * r_206_3;
   FReal2 r_208_2 = Rz * r_207_3 - 7 * r_206_3;
   FReal2 r_0a0_2 = Ry * r_090_3 - 9 * r_080_3;
   FReal2 r_082_2 = Ry * r_072_3 - 7 * r_062_3;
   FReal2 r_064_2 = Rz * r_063_3 - 3 * r_062_3;
   FReal2 r_046_2 = Ry * r_036_3 - 3 * r_026_3;
   FReal2 r_00a_2 = Rz * r_009_3 - 9 * r_008_3;
   FReal2 r_b00_2 = Rx * r_a00_3 - 10 * r_900_3;
   FReal2 r_902_2 = Rz * r_901_3 - r_900_3;
   FReal2 r_740_2 = Ry * r_730_3 - 3 * r_720_3;
   FReal2 r_704_2 = Rx * r_604_3 - 6 * r_504_3;
   FReal2 r_560_2 = Rx * r_460_3 - 4 * r_360_3;
   FReal2 r_524_2 = Ry * r_514_3 - r_504_3;
   FReal2 r_380_2 = Ry * r_370_3 - 7 * r_360_3;
   FReal2 r_308_2 = Rz * r_307_3 - 7 * r_306_3;
   FReal2 r_1a0_2 = Rx * r_0a0_3;
   FReal2 r_164_2 = Rx * r_064_3;
   FReal2 r_146_2 = Rx * r_046_3;
   FReal2 r_10a_2 = Rx * r_00a_3;
   FReal2 r_900_2 = Rx * r_800_3 - 8 * r_700_3;
   FReal2 r_720_2 = Ry * r_710_3 - r_700_3;
   FReal2 r_504_2 = Rz * r_503_3 - 3 * r_502_3;
   FReal2 r_360_2 = Rx * r_260_3 - 2 * r_160_3;
   FReal2 r_108_2 = Rx * r_008_3;
   FReal2 r_a10_2 = Ry * r_a00_3;
   FReal2 r_830_2 = Rx * r_730_3 - 7 * r_630_3;
   FReal2 r_614_2 = Ry * r_604_3;
   FReal2 r_470_2 = Ry * r_460_3 - 6 * r_450_3;
   FReal2 r_452_2 = Rz * r_451_3 - r_450_3;
   FReal2 r_416_2 = Ry * r_406_3;
   FReal2 r_290_2 = Rx * r_190_3 - r_090_3;
   FReal2 r_0b0_2 = Ry * r_0a0_3 - 10 * r_090_3;
   FReal2 r_092_2 = Rz * r_091_3 - r_090_3;
   FReal2 r_074_2 = Rz * r_073_3 - 3 * r_072_3;
   FReal2 r_056_2 = Ry * r_046_3 - 4 * r_036_3;
   FReal2 r_038_2 = Rz * r_037_3 - 7 * r_036_3;
   FReal2 r_01a_2 = Ry * r_00a_3;
   FReal2 r_810_2 = Ry * r_800_3;
   FReal2 r_450_2 = Rx * r_350_3 - 3 * r_250_3;
   FReal2 r_090_2 = Ry * r_080_3 - 8 * r_070_3;
   FReal2 r_072_2 = Rz * r_071_3 - r_070_3;
   FReal2 r_036_2 = Ry * r_026_3 - 2 * r_016_3;
   FReal2 r_a01_2 = Rz * r_a00_3;
   FReal2 r_803_2 = Rx * r_703_3 - 7 * r_603_3;
   FReal2 r_641_2 = Rz * r_640_3;
   FReal2 r_605_2 = Rz * r_604_3 - 4 * r_603_3;
   FReal2 r_461_2 = Rz * r_460_3;
   FReal2 r_407_2 = Rx * r_307_3 - 3 * r_207_3;
   FReal2 r_245_2 = Rx * r_145_3 - r_045_3;
   FReal2 r_083_2 = Ry * r_073_3 - 7 * r_063_3;
   FReal2 r_047_2 = Rz * r_046_3 - 6 * r_045_3;
   FReal2 r_029_2 = Ry * r_019_3 - r_009_3;
   FReal2 r_00b_2 = Rz * r_00a_3 - 10 * r_009_3;
   FReal2 r_603_2 = Rz * r_602_3 - 2 * r_601_3;
   FReal2 r_207_2 = Rx * r_107_3 - r_007_3;
   FReal2 r_081_2 = Rz * r_080_3;
   FReal2 r_045_2 = Ry * r_035_3 - 3 * r_025_3;
   FReal2 r_009_2 = Rz * r_008_3 - 8 * r_007_3;
   FReal2 r_730_2 = Ry * r_720_3 - 2 * r_710_3;
   FReal2 r_550_2 = Rx * r_450_3 - 4 * r_350_3;
   FReal2 r_514_2 = Ry * r_504_3;
   FReal2 r_370_2 = Ry * r_360_3 - 6 * r_350_3;
   FReal2 r_190_2 = Rx * r_090_3;
   FReal2 r_901_2 = Rz * r_900_3;
   FReal2 r_703_2 = Rx * r_603_3 - 6 * r_503_3;
   FReal2 r_505_2 = Rz * r_504_3 - 4 * r_503_3;
   FReal2 r_307_2 = Rx * r_207_3 - 2 * r_107_3;
   FReal2 r_145_2 = Rx * r_045_3;
   FReal2 r_451_2 = Rz * r_450_3;
   FReal2 r_073_2 = Rz * r_072_3 - 2 * r_071_3;
   FReal2 r_055_2 = Ry * r_045_3 - 4 * r_035_3;
   FReal2 r_037_2 = Rz * r_036_3 - 6 * r_035_3;
   FReal2 r_019_2 = Ry * r_009_3;
   FReal2 r_713_2 = Ry * r_703_3;
   FReal2 r_551_2 = Rz * r_550_3;
   FReal2 r_515_2 = Ry * r_505_3;
   FReal2 r_371_2 = Rz * r_370_3;
   FReal2 r_155_2 = Rx * r_055_3;
   FReal2 r_137_2 = Rx * r_037_3;

   // m = 1
   FReal2 r_c00_1 = Rx * r_b00_2 - 11 * r_a00_2;
   FReal2 r_a20_1 = Ry * r_a10_2 - r_a00_2;
   FReal2 r_a02_1 = Rz * r_a01_2 - r_a00_2;
   FReal2 r_840_1 = Ry * r_830_2 - 3 * r_820_2;
   FReal2 r_804_1 = Rx * r_704_2 - 7 * r_604_2;
   FReal2 r_642_1 = Rz * r_641_2 - r_640_2;
   FReal2 r_624_1 = Ry * r_614_2 - r_604_2;
   FReal2 r_480_1 = Ry * r_470_2 - 7 * r_460_2;
   FReal2 r_462_1 = Rz * r_461_2 - r_460_2;
   FReal2 r_426_1 = Ry * r_416_2 - r_406_2;
   FReal2 r_408_1 = Rx * r_308_2 - 3 * r_208_2;
   FReal2 r_2a0_1 = Rx * r_1a0_2 - r_0a0_2;
   FReal2 r_264_1 = Rx * r_164_2 - r_064_2;
   FReal2 r_246_1 = Rx * r_146_2 - r_046_2;
   FReal2 r_20a_1 = Rx * r_10a_2 - r_00a_2;
   FReal2 r_0c0_1 = Ry * r_0b0_2 - 11 * r_0a0_2;
   FReal2 r_0a2_1 = Ry * r_092_2 - 9 * r_082_2;
   FReal2 r_084_1 = Rz * r_083_2 - 3 * r_082_2;
   FReal2 r_048_1 = Rz * r_047_2 - 7 * r_046_2;
   FReal2 r_02a_1 = Ry * r_01a_2 - r_00a_2;
   FReal2 r_00c_1 = Rz * r_00b_2 - 11 * r_00a_2;
   FReal2 r_b00_1 = Rx * r_a00_2 - 10 * r_900_2;
   FReal2 r_902_1 = Rz * r_901_2 - r_900_2;
   FReal2 r_740_1 = Ry * r_730_2 - 3 * r_720_2;
   FReal2 r_704_1 = Rx * r_604_2 - 6 * r_504_2;
   FReal2 r_560_1 = Rx * r_460_2 - 4 * r_360_2;
   FReal2 r_524_1 = Ry * r_514_2 - r_504_2;
   FReal2 r_380_1 = Ry * r_370_2 - 7 * r_360_2;
   FReal2 r_308_1 = Rx * r_208_2 - 2 * r_108_2;
   FReal2 r_164_1 = Rx * r_064_2;
   FReal2 r_146_1 = Rx * r_046_2;
   FReal2 r_830_1 = Ry * r_820_2 - 2 * r_810_2;
   FReal2 r_614_1 = Ry * r_604_2;
   FReal2 r_470_1 = Ry * r_460_2 - 6 * r_450_2;
   FReal2 r_452_1 = Rz * r_451_2 - r_450_2;
   FReal2 r_416_1 = Ry * r_406_2;
   FReal2 r_290_1 = Rx * r_190_2 - r_090_2;
   FReal2 r_0b0_1 = Ry * r_0a0_2 - 10 * r_090_2;
   FReal2 r_074_1 = Rz * r_073_2 - 3 * r_072_2;
   FReal2 r_056_1 = Ry * r_046_2 - 4 * r_036_2;
   FReal2 r_038_1 = Rz * r_037_2 - 7 * r_036_2;
   FReal2 r_803_1 = Rx * r_703_2 - 7 * r_603_2;
   FReal2 r_641_1 = Rz * r_640_2;
   FReal2 r_605_1 = Rz * r_604_2 - 4 * r_603_2;
   FReal2 r_461_1 = Rz * r_460_2;
   FReal2 r_407_1 = Rx * r_307_2 - 3 * r_207_2;
   FReal2 r_245_1 = Rx * r_145_2 - r_045_2;
   FReal2 r_083_1 = Rz * r_082_2 - 2 * r_081_2;
   FReal2 r_047_1 = Rz * r_046_2 - 6 * r_045_2;
   FReal2 r_029_1 = Ry * r_019_2 - r_009_2;
   FReal2 r_00b_1 = Rz * r_00a_2 - 10 * r_009_2;
   FReal2 r_b10_1 = Ry * r_b00_2;
   FReal2 r_930_1 = Rx * r_830_2 - 8 * r_730_2;
   FReal2 r_912_1 = Ry * r_902_2;
   FReal2 r_750_1 = Ry * r_740_2 - 4 * r_730_2;
   FReal2 r_570_1 = Ry * r_560_2 - 6 * r_550_2;
   FReal2 r_552_1 = Rz * r_551_2 - r_550_2;
   FReal2 r_534_1 = Ry * r_524_2 - 2 * r_514_2;
   FReal2 r_390_1 = Rx * r_290_2 - 2 * r_190_2;
   FReal2 r_372_1 = Rz * r_371_2 - r_370_2;
   FReal2 r_318_1 = Ry * r_308_2;
   FReal2 r_1b0_1 = Rx * r_0b0_2;
   FReal2 r_174_1 = Rx * r_074_2;
   FReal2 r_156_1 = Rx * r_056_2;
   FReal2 r_138_1 = Rx * r_038_2;
   FReal2 r_b01_1 = Rx * r_a01_2 - 10 * r_901_2;
   FReal2 r_903_1 = Rz * r_902_2 - 2 * r_901_2;
   FReal2 r_741_1 = Rz * r_740_2;
   FReal2 r_723_1 = Ry * r_713_2 - r_703_2;
   FReal2 r_705_1 = Rx * r_605_2 - 6 * r_505_2;
   FReal2 r_561_1 = Rz * r_560_2;
   FReal2 r_525_1 = Ry * r_515_2 - r_505_2;
   FReal2 r_507_1 = Rx * r_407_2 - 4 * r_307_2;
   FReal2 r_381_1 = Rz * r_380_2;
   FReal2 r_345_1 = Rx * r_245_2 - 2 * r_145_2;
   FReal2 r_309_1 = Rz * r_308_2 - 8 * r_307_2;
   FReal2 r_183_1 = Rx * r_083_2;
   FReal2 r_129_1 = Rx * r_029_2;
   FReal2 r_10b_1 = Rx * r_00b_2;
   FReal2 r_831_1 = Rz * r_830_2;
   FReal2 r_813_1 = Ry * r_803_2;
   FReal2 r_615_1 = Ry * r_605_2;
   FReal2 r_453_1 = Rz * r_452_2 - 2 * r_451_2;
   FReal2 r_417_1 = Ry * r_407_2;
   FReal2 r_291_1 = Rz * r_290_2;
   FReal2 r_255_1 = Rx * r_155_2 - r_055_2;
   FReal2 r_237_1 = Rx * r_137_2 - r_037_2;
   FReal2 r_0b1_1 = Rz * r_0b0_2;
   FReal2 r_093_1 = Ry * r_083_2 - 8 * r_073_2;
   FReal2 r_075_1 = Rz * r_074_2 - 4 * r_073_2;
   FReal2 r_057_1 = Rz * r_056_2 - 6 * r_055_2;
   FReal2 r_039_1 = Ry * r_029_2 - 2 * r_019_2;
   FReal2 r_01b_1 = Rz * r_01a_2 - 10 * r_019_2;
   FReal2 r_713_1 = Ry * r_703_2;
   FReal2 r_551_1 = Rz * r_550_2;
   FReal2 r_515_1 = Ry * r_505_2;
   FReal2 r_371_1 = Rz * r_370_2;
   FReal2 r_155_1 = Rx * r_055_2;
   FReal2 r_137_1 = Rx * r_037_2;

   // m = 0
   pOut[0] = Rx * r_c00_1 - 12 * r_b00_1;
   pOut[1] = Ry * r_b10_1 - r_b00_1;
   pOut[2] = Rx * r_a02_1 - 10 * r_902_1;
   pOut[3] = Rx * r_840_1 - 8 * r_740_1;
   pOut[4] = Ry * r_912_1 - r_902_1;
   pOut[5] = Rx * r_804_1 - 8 * r_704_1;
   pOut[6] = Ry * r_750_1 - 5 * r_740_1;
   pOut[7] = Rz * r_741_1 - r_740_1;
   pOut[8] = Rx * r_624_1 - 6 * r_524_1;
   pOut[9] = Rz * r_705_1 - 5 * r_704_1;
   pOut[10] = Rx * r_480_1 - 4 * r_380_1;
   pOut[11] = Rz * r_561_1 - r_560_1;
   pOut[12] = Ry * r_534_1 - 3 * r_524_1;
   pOut[13] = Rz * r_525_1 - 5 * r_524_1;
   pOut[14] = Rx * r_408_1 - 4 * r_308_1;
   pOut[15] = Ry * r_390_1 - 9 * r_380_1;
   pOut[16] = Rz * r_381_1 - r_380_1;
   pOut[17] = Rx * r_264_1 - 2 * r_164_1;
   pOut[18] = Rx * r_246_1 - 2 * r_146_1;
   pOut[19] = Ry * r_318_1 - r_308_1;
   pOut[20] = Rz * r_309_1 - 9 * r_308_1;
   pOut[21] = Rx * r_0c0_1;
   pOut[22] = Rx * r_0a2_1;
   pOut[23] = Ry * r_174_1 - 7 * r_164_1;
   pOut[24] = Ry * r_156_1 - 5 * r_146_1;
   pOut[25] = Rx * r_048_1;
   pOut[26] = Rx * r_02a_1;
   pOut[27] = Rx * r_00c_1;
   pOut[28] = Ry * r_c00_1;
   pOut[29] = Rx * r_930_1 - 9 * r_830_1;
   pOut[30] = Ry * r_a02_1;
   pOut[31] = Ry * r_840_1 - 4 * r_830_1;
   pOut[32] = Rz * r_831_1 - r_830_1;
   pOut[33] = Ry * r_804_1;
   pOut[34] = Rx * r_570_1 - 5 * r_470_1;
   pOut[35] = Rx * r_552_1 - 5 * r_452_1;
   pOut[36] = Ry * r_624_1 - 2 * r_614_1;
   pOut[37] = Rz * r_615_1 - 5 * r_614_1;
   pOut[38] = Ry * r_480_1 - 8 * r_470_1;
   pOut[39] = Ry * r_462_1 - 6 * r_452_1;
   pOut[40] = Rz * r_453_1 - 3 * r_452_1;
   pOut[41] = Ry * r_426_1 - 2 * r_416_1;
   pOut[42] = Rz * r_417_1 - 7 * r_416_1;
   pOut[43] = Rx * r_1b0_1 - r_0b0_1;
   pOut[44] = Rz * r_291_1 - r_290_1;
   pOut[45] = Rx * r_174_1 - r_074_1;
   pOut[46] = Rx * r_156_1 - r_056_1;
   pOut[47] = Rx * r_138_1 - r_038_1;
   pOut[48] = Ry * r_20a_1;
   pOut[49] = Ry * r_0c0_1 - 12 * r_0b0_1;
   pOut[50] = Rz * r_0b1_1 - r_0b0_1;
   pOut[51] = Ry * r_084_1 - 8 * r_074_1;
   pOut[52] = Rz * r_075_1 - 5 * r_074_1;
   pOut[53] = Ry * r_048_1 - 4 * r_038_1;
   pOut[54] = Rz * r_039_1 - 9 * r_038_1;
   pOut[55] = Ry * r_00c_1;
   pOut[56] = Rz * r_c00_1;
   pOut[57] = Rz * r_a20_1;
   pOut[58] = Rx * r_903_1 - 9 * r_803_1;
   pOut[59] = Rx * r_741_1 - 7 * r_641_1;
   pOut[60] = Ry * r_813_1 - r_803_1;
   pOut[61] = Rz * r_804_1 - 4 * r_803_1;
   pOut[62] = Rx * r_561_1 - 5 * r_461_1;
   pOut[63] = Rz * r_642_1 - 2 * r_641_1;
   pOut[64] = Ry * r_615_1 - r_605_1;
   pOut[65] = Rx * r_507_1 - 5 * r_407_1;
   pOut[66] = Rz * r_480_1;
   pOut[67] = Rz * r_462_1 - 2 * r_461_1;
   pOut[68] = Rx * r_345_1 - 3 * r_245_1;
   pOut[69] = Ry * r_417_1 - r_407_1;
   pOut[70] = Rz * r_408_1 - 8 * r_407_1;
   pOut[71] = Rz * r_2a0_1;
   pOut[72] = Rx * r_183_1 - r_083_1;
   pOut[73] = Ry * r_255_1 - 5 * r_245_1;
   pOut[74] = Rz * r_246_1 - 6 * r_245_1;
   pOut[75] = Rx * r_129_1 - r_029_1;
   pOut[76] = Rx * r_10b_1 - r_00b_1;
   pOut[77] = Rz * r_0c0_1;
   pOut[78] = Ry * r_093_1 - 9 * r_083_1;
   pOut[79] = Rz * r_084_1 - 4 * r_083_1;
   pOut[80] = Ry * r_057_1 - 5 * r_047_1;
   pOut[81] = Rz * r_048_1 - 8 * r_047_1;
   pOut[82] = Rz * r_02a_1 - 10 * r_029_1;
   pOut[83] = Rz * r_00c_1 - 12 * r_00b_1;
   pOut[84] = Ry * r_b01_1;
   pOut[85] = Rz * r_930_1;
   pOut[86] = Rx * r_813_1 - 8 * r_713_1;
   pOut[87] = Rz * r_750_1;
   pOut[88] = Ry * r_723_1 - 2 * r_713_1;
   pOut[89] = Rx * r_615_1 - 6 * r_515_1;
   pOut[90] = Ry * r_561_1 - 6 * r_551_1;
   pOut[91] = Rz * r_552_1 - 2 * r_551_1;
   pOut[92] = Ry * r_525_1 - 2 * r_515_1;
   pOut[93] = Ry * r_507_1;
   pOut[94] = Ry * r_381_1 - 8 * r_371_1;
   pOut[95] = Rz * r_372_1 - 2 * r_371_1;
   pOut[96] = Rx * r_255_1 - 2 * r_155_1;
   pOut[97] = Rx * r_237_1 - 2 * r_137_1;
   pOut[98] = Ry * r_309_1;
   pOut[99] = Rz * r_1b0_1;
   pOut[100] = Rx * r_093_1;
   pOut[101] = Rx * r_075_1;
   pOut[102] = Rz * r_156_1 - 6 * r_155_1;
   pOut[103] = Rz * r_138_1 - 8 * r_137_1;
   pOut[104] = Rx * r_01b_1;
   // 831 flops, 1484 mops   (7.91 flops/integral, 14.13 mops/integral)
}

// In: Gm[0 .. 14] (inclusive), Out: [r]^0, ordered as AngularComps(14)
// (moment 14, 120 entries)
void Mdrr14(FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const R[3])
{
   FReal1
      Rx = R[0], Ry = R[1], Rz = R[2];
   // m = 13
   FReal2 r_100_13 = Rx * pIn[14];
   FReal2 r_010_13 = Ry * pIn[14];
   FReal2 r_001_13 = Rz * pIn[14];

   // m = 12
   FReal2 r_200_12 = Rx * r_100_13 - pIn[13];
   FReal2 r_020_12 = Ry * r_010_13 - pIn[13];
   FReal2 r_002_12 = Rz * r_001_13 - pIn[13];
   FReal2 r_100_12 = Rx * pIn[13];
   FReal2 r_010_12 = Ry * pIn[13];
   FReal2 r_001_12 = Rz * pIn[13];

   // m = 11
   FReal2 r_200_11 = Rx * r_100_12 - pIn[12];
   FReal2 r_020_11 = Ry * r_010_12 - pIn[12];
   FReal2 r_002_11 = Rz * r_001_12 - pIn[12];
   FReal2 r_300_11 = Rx * r_200_12 - 2 * r_100_12;
   FReal2 r_100_11 = Rx * pIn[12];
   FReal2 r_030_11 = Ry * r_020_12 - 2 * r_010_12;
   FReal2 r_010_11 = Ry * pIn[12];
   FReal2 r_003_11 = Rz * r_002_12 - 2 * r_001_12;
   FReal2 r_001_11 = Rz * pIn[12];

   // m = 10
   FReal2 r_400_10 = Rx * r_300_11 - 3 * r_200_11;
   FReal2 r_040_10 = Ry * r_030_11 - 3 * r_020_11;
   FReal2 r_004_10 = Rz * r_003_11 - 3 * r_002_11;
   FReal2 r_200_10 = Rx * r_100_11 - pIn[11];
   FReal2 r_020_10 = Ry * r_010_11 - pIn[11];
   FReal2 r_002_10 = Rz * r_001_11 - pIn[11];
   FReal2 r_300_10 = Rx * r_200_11 - 2 * r_100_11;
   FReal2 r_100_10 = Rx * pIn[11];
   FReal2 r_030_10 = Ry * r_020_11 - 2 * r_010_11;
   FReal2 r_010_10 = Ry * pIn[11];
   FReal2 r_003_10 = Rz * r_002_11 - 2 * r_001_11;
   FReal2 r_001_10 = Rz * pIn[11];

   // m = 9
   FReal2 r_400_9 = Rx * r_300_10 - 3 * r_200_10;
   FReal2 r_040_9 = Ry * r_030_10 - 3 * r_020_10;
   FReal2 r_004_9 = Rz * r_003_10 - 3 * r_002_10;
   FReal2 r_200_9 = Rx * r_100_10 - pIn[10];
   FReal2 r_020_9 = Ry * r_010_10 - pIn[10];
   FReal2 r_002_9 = Rz * r_001_10 - pIn[10];
   FReal2 r_500_9 = Rx * r_400_10 - 4 * r_300_10;
   FReal2 r_300_9 = Rx * r_200_10 - 2 * r_100_10;
   FReal2 r_100_9 = Rx * pIn[10];
   FReal2 r_050_9 = Ry * r_040_10 - 4 * r_030_10;
   FReal2 r_030_9 = Ry * r_020_10 - 2 * r_010_10;
   FReal2 r_010_9 = Ry * pIn[10];
   FReal2 r_005_9 = Rz * r_004_10 - 4 * r_003_10;
   FReal2 r_003_9 = Rz * r_002_10 - 2 * r_001_10;
   FReal2 r_001_9 = Rz * pIn[10];

   // m = 8
   FReal2 r_600_8 = Rx * r_500_9 - 5 * r_400_9;
   FReal2 r_060_8 = Ry * r_050_9 - 5 * r_040_9;
   FReal2 r_006_8 = Rz * r_005_9 - 5 * r_004_9;
   FReal2 r_400_8 = Rx * r_300_9 - 3 * r_200_9;
   FReal2 r_040_8 = Ry * r_030_9 - 3 * r_020_9;
   FReal2 r_004_8 = Rz * r_003_9 - 3 * r_002_9;
   FReal2 r_200_8 = Rx * r_100_9 - pIn[9];
   FReal2 r_020_8 = Ry * r_010_9 - pIn[9];
   FReal2 r_002_8 = Rz * r_001_9 - pIn[9];
   FReal2 r_500_8 = Rx * r_400_9 - 4 * r_300_9;
   FReal2 r_300_8 = Rx * r_200_9 - 2 * r_100_9;
   FReal2 r_100_8 = Rx * pIn[9];
   FReal2 r_050_8 = Ry * r_040_9 - 4 * r_030_9;
   FReal2 r_030_8 = Ry * r_020_9 - 2 * r_010_9;
   FReal2 r_010_8 = Ry * pIn[9];
   FReal2 r_005_8 = Rz * r_004_9 - 4 * r_003_9;
   FReal2 r_003_8 = Rz * r_002_9 - 2 * r_001_9;
   FReal2 r_001_8 = Rz * pIn[9];
   FReal2 r_150_8 = Rx * r_050_9;
   FReal2 r_501_8 = Rz * r_500_9;
   FReal2 r_015_8 = Ry * r_005_9;

   // m = 7
   FReal2 r_600_7 = Rx * r_500_8 - 5 * r_400_8;
   FReal2 r_060_7 = Ry * r_050_8 - 5 * r_040_8;
   FReal2 r_006_7 = Rz * r_005_8 - 5 * r_004_8;
   FReal2 r_400_7 = Rx * r_300_8 - 3 * r_200_8;
   FReal2 r_040_7 = Ry * r_030_8 - 3 * r_020_8;
   FReal2 r_004_7 = Rz * r_003_8 - 3 * r_002_8;
   FReal2 r_200_7 = Rx * r_100_8 - pIn[8];
   FReal2 r_020_7 = Ry * r_010_8 - pIn[8];
   FReal2 r_002_7 = Rz * r_001_8 - pIn[8];
   FReal2 r_700_7 = Rx * r_600_8 - 6 * r_500_8;
   FReal2 r_502_7 = Rz * r_501_8 - r_500_8;
   FReal2 r_160_7 = Rx * r_060_8;
   FReal2 r_500_7 = Rx * r_400_8 - 4 * r_300_8;
   FReal2 r_300_7 = Rx * r_200_8 - 2 * r_100_8;
   FReal2 r_100_7 = Rx * pIn[8];
   FReal2 r_250_7 = Rx * r_150_8 - r_050_8;
   FReal2 r_070_7 = Ry * r_060_8 - 6 * r_050_8;
   FReal2 r_016_7 = Ry * r_006_8;
   FReal2 r_050_7 = Ry * r_040_8 - 4 * r_030_8;
   FReal2 r_030_7 = Ry * r_020_8 - 2 * r_010_8;
   FReal2 r_010_7 = Ry * pIn[8];
   FReal2 r_601_7 = Rz * r_600_8;
   FReal2 r_025_7 = Ry * r_015_8 - r_005_8;
   FReal2 r_007_7 = Rz * r_006_8 - 6 * r_005_8;
   FReal2 r_005_7 = Rz * r_004_8 - 4 * r_003_8;
   FReal2 r_003_7 = Rz * r_002_8 - 2 * r_001_8;
   FReal2 r_001_7 = Rz * pIn[8];
   FReal2 r_150_7 = Rx * r_050_8;
   FReal2 r_501_7 = Rz * r_500_8;
   FReal2 r_015_7 = Ry * r_005_8;

   // m = 6
   FReal2 r_800_6 = Rx * r_700_7 - 7 * r_600_7;
   FReal2 r_602_6 = Rz * r_601_7 - r_600_7;
   FReal2 r_260_6 = Rx * r_160_7 - r_060_7;
   FReal2 r_080_6 = Ry * r_070_7 - 7 * r_060_7;
   FReal2 r_026_6 = Ry * r_016_7 - r_006_7;
   FReal2 r_008_6 = Rz * r_007_7 - 7 * r_006_7;
   FReal2 r_600_6 = Rx * r_500_7 - 5 * r_400_7;
   FReal2 r_060_6 = Ry * r_050_7 - 5 * r_040_7;
   FReal2 r_006_6 = Rz * r_005_7 - 5 * r_004_7;
   FReal2 r_400_6 = Rx * r_300_7 - 3 * r_200_7;
   FReal2 r_040_6 = Ry * r_030_7 - 3 * r_020_7;
   FReal2 r_004_6 = Rz * r_003_7 - 3 * r_002_7;
   FReal2 r_200_6 = Rx * r_100_7 - pIn[7];
   FReal2 r_020_6 = Ry * r_010_7 - pIn[7];
   FReal2 r_002_6 = Rz * r_001_7 - pIn[7];
   FReal2 r_700_6 = Rx * r_600_7 - 6 * r_500_7;
   FReal2 r_502_6 = Rz * r_501_7 - r_500_7;
   FReal2 r_160_6 = Rx * r_060_7;
   FReal2 r_500_6 = Rx * r_400_7 - 4 * r_300_7;
   FReal2 r_300_6 = Rx * r_200_7 - 2 * r_100_7;
   FReal2 r_250_6 = Rx * r_150_7 - r_050_7;
   FReal2 r_070_6 = Ry * r_060_7 - 6 * r_050_7;
   FReal2 r_016_6 = Ry * r_006_7;
   FReal2 r_050_6 = Ry * r_040_7 - 4 * r_030_7;
   FReal2 r_030_6 = Ry * r_020_7 - 2 * r_010_7;
   FReal2 r_601_6 = Rz * r_600_7;
   FReal2 r_025_6 = Ry * r_015_7 - r_005_7;
   FReal2 r_007_6 = Rz * r_006_7 - 6 * r_005_7;
   FReal2 r_005_6 = Rz * r_004_7 - 4 * r_003_7;
   FReal2 r_003_6 = Rz * r_002_7 - 2 * r_001_7;
   FReal2 r_710_6 = Ry * r_700_7;
   FReal2 r_350_6 = Rx * r_250_7 - 2 * r_150_7;
   FReal2 r_150_6 = Rx * r_050_7;
   FReal2 r_503_6 = Rz * r_502_7 - 2 * r_501_7;
   FReal2 r_107_6 = Rx * r_007_7;
   FReal2 r_501_6 = Rz * r_500_7;
   FReal2 r_071_6 = Rz * r_070_7;
   FReal2 r_035_6 = Ry * r_025_7 - 2 * r_015_7;
   FReal2 r_015_6 = Ry * r_005_7;

   // m = 5
   FReal2 r_800_5 = Rx * r_700_6 - 7 * r_600_6;
   FReal2 r_602_5 = Rz * r_601_6 - r_600_6;
   FReal2 r_260_5 = Rx * r_160_6 - r_060_6;
   FReal2 r_080_5 = Ry * r_070_6 - 7 * r_060_6;
   FReal2 r_026_5 = Ry * r_016_6 - r_006_6;
   FReal2 r_008_5 = Rz * r_007_6 - 7 * r_006_6;
   FReal2 r_600_5 = Rx * r_500_6 - 5 * r_400_6;
   FReal2 r_060_5 = Ry * r_050_6 - 5 * r_040_6;
   FReal2 r_006_5 = Rz * r_005_6 - 5 * r_004_6;
   FReal2 r_400_5 = Rx * r_300_6 - 3 * r_200_6;
   FReal2 r_040_5 = Ry * r_030_6 - 3 * r_020_6;
   FReal2 r_004_5 = Rz * r_003_6 - 3 * r_002_6;
   FReal2 r_900_5 = Rx * r_800_6 - 8 * r_700_6;
   FReal2 r_720_5 = Ry * r_710_6 - r_700_6;
   FReal2 r_504_5 = Rz * r_503_6 - 3 * r_502_6;
   FReal2 r_360_5 = Rx * r_260_6 - 2 * r_160_6;
   FReal2 r_108_5 = Rx * r_008_6;
   FReal2 r_700_5 = Rx * r_600_6 - 6 * r_500_6;
   FReal2 r_502_5 = Rz * r_501_6 - r_500_6;
   FReal2 r_160_5 = Rx * r_060_6;
   FReal2 r_500_5 = Rx * r_400_6 - 4 * r_300_6;
   FReal2 r_810_5 = Ry * r_800_6;
   FReal2 r_450_5 = Rx * r_350_6 - 3 * r_250_6;
   FReal2 r_090_5 = Ry * r_080_6 - 8 * r_070_6;
   FReal2 r_072_5 = Rz * r_071_6 - r_070_6;
   FReal2 r_036_5 = Ry * r_026_6 - 2 * r_016_6;
   FReal2 r_250_5 = Rx * r_150_6 - r_050_6;
   FReal2 r_070_5 = Ry * r_060_6 - 6 * r_050_6;
   FReal2 r_016_5 = Ry * r_006_6;
   FReal2 r_050_5 = Ry * r_040_6 - 4 * r_030_6;
   FReal2 r_603_5 = Rz * r_602_6 - 2 * r_601_6;
   FReal2 r_207_5 = Rx * r_107_6 - r_007_6;
   FReal2 r_081_5 = Rz * r_080_6;
   FReal2 r_045_5 = Ry * r_035_6 - 3 * r_025_6;
   FReal2 r_009_5 = Rz * r_008_6 - 8 * r_007_6;
   FReal2 r_601_5 = Rz * r_600_6;
   FReal2 r_025_5 = Ry * r_015_6 - r_005_6;
   FReal2 r_007_5 = Rz * r_006_6 - 6 * r_005_6;
   FReal2 r_005_5 = Rz * r_004_6 - 4 * r_003_6;
   FReal2 r_710_5 = Ry * r_700_6;
   FReal2 r_350_5 = Rx * r_250_6 - 2 * r_150_6;
   FReal2 r_150_5 = Rx * r_050_6;
   FReal2 r_503_5 = Rz * r_502_6 - 2 * r_501_6;
   FReal2 r_107_5 = Rx * r_007_6;
   FReal2 r_501_5 = Rz * r_500_6;
   FReal2 r_071_5 = Rz * r_070_6;
   FReal2 r_035_5 = Ry * r_025_6 - 2 * r_015_6;
   FReal2 r_015_5 = Ry * r_005_6;

   // m = 4
   FReal2 r_a00_4 = Rx * r_900_5 - 9 * r_800_5;
   FReal2 r_820_4 = Ry * r_810_5 - r_800_5;
   FReal2 r_604_4 = Rz * r_603_5 - 3 * r_602_5;
   FReal2 r_460_4 = Rx * r_360_5 - 3 * r_260_5;
   FReal2 r_208_4 = Rx * r_108_5 - r_008_5;
   FReal2 r_0a0_4 = Ry * r_090_5 - 9 * r_080_5;
   FReal2 r_082_4 = Rz * r_081_5 - r_080_5;
   FReal2 r_046_4 = Ry * r_036_5 - 3 * r_026_5;
   FReal2 r_00a_4 = Rz * r_009_5 - 9 * r_008_5;
   FReal2 r_800_4 = Rx * r_700_5 - 7 * r_600_5;
   FReal2 r_602_4 = Rz * r_601_5 - r_600_5;
   FReal2 r_260_4 = Rx * r_160_5 - r_060_5;
   FReal2 r_080_4 = Ry * r_070_5 - 7 * r_060_5;
   FReal2 r_026_4 = Ry * r_016_5 - r_006_5;
   FReal2 r_008_4 = Rz * r_007_5 - 7 * r_006_5;
   FReal2 r_600_4 = Rx * r_500_5 - 5 * r_400_5;
   FReal2 r_060_4 = Ry * r_050_5 - 5 * r_040_5;
   FReal2 r_006_4 = Rz * r_005_5 - 5 * r_004_5;
   FReal2 r_900_4 = Rx * r_800_5 - 8 * r_700_5;
   FReal2 r_720_4 = Ry * r_710_5 - r_700_5;
   FReal2 r_504_4 = Rz * r_503_5 - 3 * r_502_5;
   FReal2 r_360_4 = Rx * r_260_5 - 2 * r_160_5;
   FReal2 r_108_4 = Rx * r_008_5;
   FReal2 r_700_4 = Rx * r_600_5 - 6 * r_500_5;
   FReal2 r_502_4 = Rz * r_501_5 - r_500_5;
   FReal2 r_160_4 = Rx * r_060_5;
   FReal2 r_810_4 = Ry * r_800_5;
   FReal2 r_450_4 = Rx * r_350_5 - 3 * r_250_5;
   FReal2 r_090_4 = Ry * r_080_5 - 8 * r_070_5;
   FReal2 r_072_4 = Rz * r_071_5 - r_070_5;
   FReal2 r_036_4 = Ry * r_026_5 - 2 * r_016_5;
   FReal2 r_250_4 = Rx * r_150_5 - r_050_5;
   FReal2 r_070_4 = Ry * r_060_5 - 6 * r_050_5;
   FReal2 r_016_4 = Ry * r_006_5;
   FReal2 r_603_4 = Rz * r_602_5 - 2 * r_601_5;
   FReal2 r_207_4 = Rx * r_107_5 - r_007_5;
   FReal2 r_081_4 = Rz * r_080_5;
   FReal2 r_045_4 = Ry * r_035_5 - 3 * r_025_5;
   FReal2 r_009_4 = Rz * r_008_5 - 8 * r_007_5;
   FReal2 r_601_4 = Rz * r_600_5;
   FReal2 r_025_4 = Ry * r_015_5 - r_005_5;
   FReal2 r_007_4 = Rz * r_006_5 - 6 * r_005_5;
   FReal2 r_730_4 = Ry * r_720_5 - 2 * r_710_5;
   FReal2 r_550_4 = Rx * r_450_5 - 4 * r_350_5;
   FReal2 r_370_4 = Ry * r_360_5 - 6 * r_350_5;
   FReal2 r_710_4 = Ry * r_700_5;
   FReal2 r_350_4 = Rx * r_250_5 - 2 * r_150_5;
   FReal2 r_150_4 = Rx * r_050_5;
   FReal2 r_703_4 = Rx * r_603_5 - 6 * r_503_5;
   FReal2 r_505_4 = Rz * r_504_5 - 4 * r_503_5;
   FReal2 r_307_4 = Rx * r_207_5 - 2 * r_107_5;
   FReal2 r_503_4 = Rz * r_502_5 - 2 * r_501_5;
   FReal2 r_107_4 = Rx * r_007_5;
   FReal2 r_501_4 = Rz * r_500_5;
   FReal2 r_073_4 = Rz * r_072_5 - 2 * r_071_5;
   FReal2 r_055_4 = Ry * r_045_5 - 4 * r_035_5;
   FReal2 r_037_4 = Rz * r_036_5 - 6 * r_035_5;
   FReal2 r_071_4 = Rz * r_070_5;
   FReal2 r_035_4 = Ry * r_025_5 - 2 * r_015_5;
   FReal2 r_015_4 = Ry * r_005_5;

   // m = 3
   FReal2 r_a00_3 = Rx * r_900_4 - 9 * r_800_4;
   FReal2 r_820_3 = Ry * r_810_4 - r_800_4;
   FReal2 r_604_3 = Rz * r_603_4 - 3 * r_602_4;
   FReal2 r_460_3 = Rx * r_360_4 - 3 * r_260_4;
   FReal2 r_208_3 = Rx * r_108_4 - r_008_4;
   FReal2 r_0a0_3 = Ry * r_090_4 - 9 * r_080_4;
   FReal2 r_082_3 = Rz * r_081_4 - r_080_4;
   FReal2 r_046_3 = Ry * r_036_4 - 3 * r_026_4;
   FReal2 r_00a_3 = Rz * r_009_4 - 9 * r_008_4;
   FReal2 r_800_3 = Rx * r_700_4 - 7 * r_600_4;
   FReal2 r_602_3 = Rz * r_601_4 - r_600_4;
   FReal2 r_260_3 = Rx * r_160_4 - r_060_4;
   FReal2 r_080_3 = Ry * r_070_4 - 7 * r_060_4;
   FReal2 r_026_3 = Ry * r_016_4 - r_006_4;
   FReal2 r_008_3 = Rz * r_007_4 - 7 * r_006_4;
   FReal2 r_b00_3 = Rx * r_a00_4 - 10 * r_900_4;
   FReal2 r_740_3 = Ry * r_730_4 - 3 * r_720_4;
   FReal2 r_704_3 = Rx * r_604_4 - 6 * r_504_4;
   FReal2 r_560_3 = Rx * r_460_4 - 4 * r_360_4;
   FReal2 r_506_3 = Rz * r_505_4 - 5 * r_504_4;
   FReal2 r_380_3 = Ry * r_370_4 - 7 * r_360_4;
   FReal2 r_308_3 = Rx * r_208_4 - 2 * r_108_4;
   FReal2 r_1a0_3 = Rx * r_0a0_4;
   FReal2 r_146_3 = Rx * r_046_4;
   FReal2 r_10a_3 = Rx * r_00a_4;
   FReal2 r_900_3 = Rx * r_800_4 - 8 * r_700_4;
   FReal2 r_720_3 = Ry * r_710_4 - r_700_4;
   FReal2 r_504_3 = Rz * r_503_4 - 3 * r_502_4;
   FReal2 r_360_3 = Rx * r_260_4 - 2 * r_160_4;
   FReal2 r_108_3 = Rx * r_008_4;
   FReal2 r_a10_3 = Ry * r_a00_4;
   FReal2 r_830_3 = Ry * r_820_4 - 2 * r_810_4;
   FReal2 r_650_3 = Rx * r_550_4 - 5 * r_450_4;
   FReal2 r_614_3 = Ry * r_604_4;
   FReal2 r_470_3 = Ry * r_460_4 - 6 * r_450_4;
   FReal2 r_0b0_3 = Ry * r_0a0_4 - 10 * r_090_4;
   FReal2 r_074_3 = Rz * r_073_4 - 3 * r_072_4;
   FReal2 r_056_3 = Ry * r_046_4 - 4 * r_036_4;
   FReal2 r_038_3 = Rz * r_037_4 - 7 * r_036_4;
   FReal2 r_01a_3 = Ry * r_00a_4;
   FReal2 r_810_3 = Ry * r_800_4;
   FReal2 r_450_3 = Rx * r_350_4 - 3 * r_250_4;
   FReal2 r_090_3 = Ry * r_080_4 - 8 * r_070_4;
   FReal2 r_072_3 = Rz * r_071_4 - r_070_4;
   FReal2 r_036_3 = Ry * r_026_4 - 2 * r_016_4;
   FReal2 r_a01_3 = Rz * r_a00_4;
   FReal2 r_803_3 = Rx * r_703_4 - 7 * r_603_4;
   FReal2 r_605_3 = Rz * r_604_4 - 4 * r_603_4;
   FReal2 r_461_3 = Rz * r_460_4;
   FReal2 r_407_3 = Rx * r_307_4 - 3 * r_207_4;
   FReal2 r_0a1_3 = Rz * r_0a0_4;
   FReal2 r_083_3 = Rz * r_082_4 - 2 * r_081_4;
   FReal2 r_065_3 = Ry * r_055_4 - 5 * r_045_4;
   FReal2 r_047_3 = Rz * r_046_4 - 6 * r_045_4;
   FReal2 r_00b_3 = Rz * r_00a_4 - 10 * r_009_4;
   FReal2 r_603_3 = Rz * r_602_4 - 2 * r_601_4;
   FReal2 r_207_3 = Rx * r_107_4 - r_007_4;
   FReal2 r_081_3 = Rz * r_080_4;
   FReal2 r_045_3 = Ry * r_035_4 - 3 * r_025_4;
   FReal2 r_009_3 = Rz * r_008_4 - 8 * r_007_4;
   FReal2 r_730_3 = Ry * r_720_4 - 2 * r_710_4;
   FReal2 r_550_3 = Rx * r_450_4 - 4 * r_350_4;
   FReal2 r_370_3 = Ry * r_360_4 - 6 * r_350_4;
   FReal2 r_190_3 = Rx * r_090_4;
   FReal2 r_350_3 = Rx * r_250_4 - 2 * r_150_4;
   FReal2 r_901_3 = Rz * r_900_4;
   FReal2 r_703_3 = Rx * r_603_4 - 6 * r_503_4;
   FReal2 r_505_3 = Rz * r_504_4 - 4 * r_503_4;
   FReal2 r_307_3 = Rx * r_207_4 - 2 * r_107_4;
   FReal2 r_503_3 = Rz * r_502_4 - 2 * r_501_4;
   FReal2 r_073_3 = Rz * r_072_4 - 2 * r_071_4;
   FReal2 r_055_3 = Ry * r_045_4 - 4 * r_035_4;
   FReal2 r_037_3 = Rz * r_036_4 - 6 * r_035_4;
   FReal2 r_019_3 = Ry * r_009_4;
   FReal2 r_035_3 = Ry * r_025_4 - 2 * r_015_4;
   FReal2 r_551_3 = Rz * r_550_4;
   FReal2 r_515_3 = Ry * r_505_4;
   FReal2 r_155_3 = Rx * r_055_4;

   // m = 2
   FReal2 r_c00_2 = Rx * r_b00_3 - 11 * r_a00_3;
   FReal2 r_a20_2 = Ry * r_a10_3 - r_a00_3;
   FReal2 r_a02_2 = Rz * r_a01_3 - r_a00_3;
   FReal2 r_840_2 = Ry * r_830_3 - 3 * r_820_3;
   FReal2 r_624_2 = Ry * r_614_3 - r_604_3;
   FReal2 r_462_2 = Rz * r_461_3 - r_460_3;
   FReal2 r_408_2 = Rx * r_308_3 - 3 * r_208_3;
   FReal2 r_2a0_2 = Rx * r_1a0_3 - r_0a0_3;
   FReal2 r_246_2 = Rx * r_146_3 - r_046_3;
   FReal2 r_20a_2 = Rx * r_10a_3 - r_00a_3;
   FReal2 r_0c0_2 = Ry * r_0b0_3 - 11 * r_0a0_3;
   FReal2 r_0a2_2 = Rz * r_0a1_3 - r_0a0_3;
   FReal2 r_084_2 = Rz * r_083_3 - 3 * r_082_3;
   FReal2 r_02a_2 = Ry * r_01a_3 - r_00a_3;
   FReal2 r_00c_2 = Rz * r_00b_3 - 11 * r_00a_3;
   FReal2 r_a00_2 = Rx * r_900_3 - 9 * r_800_3;
   FReal2 r_820_2 = Ry * r_810_3 - r_800_3;
   FReal2 r_604_2 = Rz * r_603_3 - 3 * r_602_3;
   FReal2 r_460_2 = Rx * r_360_3 - 3 * r_260_3;
   FReal2 r_208_2 = Rx * r_108_3 - r_008_3;
   FReal2 r_0a0_2 = Ry * r_090_3 - 9 * r_080_3;
   FReal2 r_082_2 = Rz * r_081_3 - r_080_3;
   FReal2 r_046_2 = Ry * r_036_3 - 3 * r_026_3;
   FReal2 r_00a_2 = Rz * r_009_3 - 9 * r_008_3;
   FReal2 r_b00_2 = Rx * r_a00_3 - 10 * r_900_3;
   FReal2 r_920_2 = Rx * r_820_3 - 8 * r_720_3;
   FReal2 r_902_2 = Rz * r_901_3 - r_900_3;
   FReal2 r_740_2 = Ry * r_730_3 - 3 * r_720_3;
   FReal2 r_704_2 = Rx * r_604_3 - 6 * r_504_3;
   FReal2 r_560_2 = Rx * r_460_3 - 4 * r_360_3;
   FReal2 r_506_2 = Rz * r_505_3 - 5 * r_504_3;
   FReal2 r_380_2 = Ry * r_370_3 - 7 * r_360_3;
   FReal2 r_308_2 = Rx * r_208_3 - 2 * r_108_3;
   FReal2 r_1a0_2 = Rx * r_0a0_3;
   FReal2 r_146_2 = Rx * r_046_3;
   FReal2 r_830_2 = Ry * r_820_3 - 2 * r_810_3;
   FReal2 r_650_2 = Rx * r_550_3 - 5 * r_450_3;
   FReal2 r_614_2 = Ry * r_604_3;
   FReal2 r_470_2 = Ry * r_460_3 - 6 * r_450_3;
   FReal2 r_290_2 = Rx * r_190_3 - r_090_3;
   FReal2 r_0b0_2 = Ry * r_0a0_3 - 10 * r_090_3;
   FReal2 r_092_2 = Ry * r_082_3 - 8 * r_072_3;
   FReal2 r_074_2 = Rz * r_073_3 - 3 * r_072_3;
   FReal2 r_056_2 = Ry * r_046_3 - 4 * r_036_3;
   FReal2 r_038_2 = Rz * r_037_3 - 7 * r_036_3;
   FReal2 r_01a_2 = Ry * r_00a_3;
   FReal2 r_a01_2 = Rz * r_a00_3;
   FReal2 r_803_2 = Rx * r_703_3 - 7 * r_603_3;
   FReal2 r_605_2 = Rz * r_604_3 - 4 * r_603_3;
   FReal2 r_461_2 = Rz * r_460_3;
   FReal2 r_407_2 = Rx * r_307_3 - 3 * r_207_3;
   FReal2 r_209_2 = Rz * r_208_3 - 8 * r_207_3;
   FReal2 r_083_2 = Rz * r_082_3 - 2 * r_081_3;
   FReal2 r_065_2 = Ry * r_055_3 - 5 * r_045_3;
   FReal2 r_047_2 = Rz * r_046_3 - 6 * r_045_3;
   FReal2 r_029_2 = Ry * r_019_3 - r_009_3;
   FReal2 r_00b_2 = Rz * r_00a_3 - 10 * r_009_3;
   FReal2 r_930_2 = Rx * r_830_3 - 8 * r_730_3;
   FReal2 r_750_2 = Ry * r_740_3 - 4 * r_730_3;
   FReal2 r_714_2 = Ry * r_704_3;
   FReal2 r_570_2 = Rx * r_470_3 - 4 * r_370_3;
   FReal2 r_552_2 = Rz * r_551_3 - r_550_3;
   FReal2 r_516_2 = Ry * r_506_3;
   FReal2 r_390_2 = Ry * r_380_3 - 8 * r_370_3;
   FReal2 r_174_2 = Rx * r_074_3;
   FReal2 r_156_2 = Rx * r_056_3;
   FReal2 r_138_2 = Rx * r_038_3;
   FReal2 r_910_2 = Ry * r_900_3;
   FReal2 r_550_2 = Rx * r_450_3 - 4 * r_350_3;
   FReal2 r_190_2 = Rx * r_090_3;
   FReal2 r_903_2 = Rx * r_803_3 - 8 * r_703_3;
   FReal2 r_741_2 = Rz * r_740_3;
   FReal2 r_705_2 = Rz * r_704_3 - 4 * r_703_3;
   FReal2 r_561_2 = Rz * r_560_3;
   FReal2 r_525_2 = Ry * r_515_3 - r_505_3;
   FReal2 r_507_2 = Rx * r_407_3 - 4 * r_307_3;
   FReal2 r_381_2 = Rz * r_380_3;
   FReal2 r_309_2 = Rz * r_308_3 - 8 * r_307_3;
   FReal2 r_165_2 = Rx * r_065_3;
   FReal2 r_147_2 = Rx * r_047_3;
   FReal2 r_901_2 = Rz * r_900_3;
   FReal2 r_505_2 = Rz * r_504_3 - 4 * r_503_3;
   FReal2 r_109_2 = Rx * r_009_3;
   FReal2 r_813_2 = Ry * r_803_3;
   FReal2 r_651_2 = Rz * r_650_3;
   FReal2 r_615_2 = Ry * r_605_3;
   FReal2 r_471_2 = Rz * r_470_3;
   FReal2 r_417_2 = Ry * r_407_3;
   FReal2 r_255_2 = Rx * r_155_3 - r_055_3;
   FReal2 r_093_2 = Ry * r_083_3 - 8 * r_073_3;
   FReal2 r_075_2 = Rz * r_074_3 - 4 * r_073_3;
   FReal2 r_057_2 = Ry * r_047_3 - 4 * r_037_3;
   FReal2 r_039_2 = Rz * r_038_3 - 8 * r_037_3;
   FReal2 r_091_2 = Rz * r_090_3;
   FReal2 r_055_2 = Ry * r_045_3 - 4 * r_035_3;
   FReal2 r_019_2 = Ry * r_009_3;
   FReal2 r_551_2 = Rz * r_550_3;
   FReal2 r_515_2 = Ry * r_505_3;
   FReal2 r_155_2 = Rx * r_055_3;

   // m = 1
   FReal2 r_c00_1 = Rx * r_b00_2 - 11 * r_a00_2;
   FReal2 r_a02_1 = Rz * r_a01_2 - r_a00_2;
   FReal2 r_840_1 = Ry * r_830_2 - 3 * r_820_2;
   FReal2 r_804_1 = Rx * r_704_2 - 7 * r_604_2;
   FReal2 r_624_1 = Ry * r_614_2 - r_604_2;
   FReal2 r_480_1 = Ry * r_470_2 - 7 * r_460_2;
   FReal2 r_462_1 = Rz * r_461_2 - r_460_2;
   FReal2 r_408_1 = Rx * r_308_2 - 3 * r_208_2;
   FReal2 r_2a0_1 = Rx * r_1a0_2 - r_0a0_2;
   FReal2 r_246_1 = Rx * r_146_2 - r_046_2;
   FReal2 r_0c0_1 = Ry * r_0b0_2 - 11 * r_0a0_2;
   FReal2 r_084_1 = Rz * r_083_2 - 3 * r_082_2;
   FReal2 r_048_1 = Rz * r_047_2 - 7 * r_046_2;
   FReal2 r_02a_1 = Ry * r_01a_2 - r_00a_2;
   FReal2 r_00c_1 = Rz * r_00b_2 - 11 * r_00a_2;
   FReal2 r_d00_1 = Rx * r_c00_2 - 12 * r_b00_2;
   FReal2 r_b20_1 = Rx * r_a20_2 - 10 * r_920_2;
   FReal2 r_b02_1 = Rx * r_a02_2 - 10 * r_902_2;
   FReal2 r_940_1 = Ry * r_930_2 - 3 * r_920_2;
   FReal2 r_904_1 = Rz * r_903_2 - 3 * r_902_2;
   FReal2 r_760_1 = Ry * r_750_2 - 5 * r_740_2;
   FReal2 r_742_1 = Rz * r_741_2 - r_740_2;
   FReal2 r_724_1 = Ry * r_714_2 - r_704_2;
   FReal2 r_706_1 = Rz * r_705_2 - 5 * r_704_2;
   FReal2 r_580_1 = Ry * r_570_2 - 7 * r_560_2;
   FReal2 r_562_1 = Rz * r_561_2 - r_560_2;
   FReal2 r_526_1 = Ry * r_516_2 - r_506_2;
   FReal2 r_508_1 = Rx * r_408_2 - 4 * r_308_2;
   FReal2 r_3a0_1 = Rx * r_2a0_2 - 2 * r_1a0_2;
   FReal2 r_382_1 = Rz * r_381_2 - r_380_2;
   FReal2 r_346_1 = Rx * r_246_2 - 2 * r_146_2;
   FReal2 r_30a_1 = Rz * r_309_2 - 9 * r_308_2;
   FReal2 r_1c0_1 = Rx * r_0c0_2;
   FReal2 r_184_1 = Rx * r_084_2;
   FReal2 r_12a_1 = Rx * r_02a_2;
   FReal2 r_10c_1 = Rx * r_00c_2;
   FReal2 r_c10_1 = Ry * r_c00_2;
   FReal2 r_a30_1 = Rx * r_930_2 - 9 * r_830_2;
   FReal2 r_a12_1 = Ry * r_a02_2;
   FReal2 r_850_1 = Ry * r_840_2 - 4 * r_830_2;
   FReal2 r_670_1 = Rx * r_570_2 - 5 * r_470_2;
   FReal2 r_652_1 = Rz * r_651_2 - r_650_2;
   FReal2 r_634_1 = Ry * r_624_2 - 2 * r_614_2;
   FReal2 r_490_1 = Rx * r_390_2 - 3 * r_290_2;
   FReal2 r_472_1 = Rz * r_471_2 - r_470_2;
   FReal2 r_418_1 = Ry * r_408_2;
   FReal2 r_2b0_1 = Ry * r_2a0_2 - 10 * r_290_2;
   FReal2 r_274_1 = Rx * r_174_2 - r_074_2;
   FReal2 r_256_1 = Rx * r_156_2 - r_056_2;
   FReal2 r_238_1 = Rx * r_138_2 - r_038_2;
   FReal2 r_0d0_1 = Ry * r_0c0_2 - 12 * r_0b0_2;
   FReal2 r_0b2_1 = Ry * r_0a2_2 - 10 * r_092_2;
   FReal2 r_094_1 = Rz * r_093_2 - 3 * r_092_2;
   FReal2 r_076_1 = Rz * r_075_2 - 5 * r_074_2;
   FReal2 r_058_1 = Rz * r_057_2 - 7 * r_056_2;
   FReal2 r_03a_1 = Ry * r_02a_2 - 2 * r_01a_2;
   FReal2 r_01c_1 = Ry * r_00c_2;
   FReal2 r_c01_1 = Rz * r_c00_2;
   FReal2 r_a03_1 = Rz * r_a02_2 - 2 * r_a01_2;
   FReal2 r_841_1 = Rz * r_840_2;
   FReal2 r_823_1 = Ry * r_813_2 - r_803_2;
   FReal2 r_805_1 = Rx * r_705_2 - 7 * r_605_2;
   FReal2 r_625_1 = Ry * r_615_2 - r_605_2;
   FReal2 r_607_1 = Rx * r_507_2 - 5 * r_407_2;
   FReal2 r_463_1 = Rz * r_462_2 - 2 * r_461_2;
   FReal2 r_427_1 = Ry * r_417_2 - r_407_2;
   FReal2 r_409_1 = Rx * r_309_2 - 3 * r_209_2;
   FReal2 r_2a1_1 = Rz * r_2a0_2;
   FReal2 r_265_1 = Rx * r_165_2 - r_065_2;
   FReal2 r_247_1 = Rx * r_147_2 - r_047_2;
   FReal2 r_20b_1 = Rz * r_20a_2 - 10 * r_209_2;
   FReal2 r_0c1_1 = Rz * r_0c0_2;
   FReal2 r_0a3_1 = Ry * r_093_2 - 9 * r_083_2;
   FReal2 r_085_1 = Rz * r_084_2 - 4 * r_083_2;
   FReal2 r_067_1 = Ry * r_057_2 - 5 * r_047_2;
   FReal2 r_049_1 = Ry * r_039_2 - 3 * r_029_2;
   FReal2 r_02b_1 = Rz * r_02a_2 - 10 * r_029_2;
   FReal2 r_00d_1 = Rz * r_00c_2 - 12 * r_00b_2;
   FReal2 r_930_1 = Ry * r_920_2 - 2 * r_910_2;
   FReal2 r_714_1 = Ry * r_704_2;
   FReal2 r_570_1 = Ry * r_560_2 - 6 * r_550_2;
   FReal2 r_552_1 = Rz * r_551_2 - r_550_2;
   FReal2 r_516_1 = Ry * r_506_2;
   FReal2 r_390_1 = Rx * r_290_2 - 2 * r_190_2;
   FReal2 r_174_1 = Rx * r_074_2;
   FReal2 r_156_1 = Rx * r_056_2;
   FReal2 r_138_1 = Rx * r_038_2;
   FReal2 r_903_1 = Rz * r_902_2 - 2 * r_901_2;
   FReal2 r_741_1 = Rz * r_740_2;
   FReal2 r_705_1 = Rx * r_605_2 - 6 * r_505_2;
   FReal2 r_561_1 = Rz * r_560_2;
   FReal2 r_525_1 = Ry * r_515_2 - r_505_2;
   FReal2 r_381_1 = Rz * r_380_2;
   FReal2 r_309_1 = Rx * r_209_2 - 2 * r_109_2;
   FReal2 r_165_1 = Rx * r_065_2;
   FReal2 r_147_1 = Rx * r_047_2;
   FReal2 r_813_1 = Ry * r_803_2;
   FReal2 r_651_1 = Rz * r_650_2;
   FReal2 r_615_1 = Ry * r_605_2;
   FReal2 r_471_1 = Rz * r_470_2;
   FReal2 r_417_1 = Ry * r_407_2;
   FReal2 r_255_1 = Rx * r_155_2 - r_055_2;
   FReal2 r_093_1 = Rz * r_092_2 - 2 * r_091_2;
   FReal2 r_057_1 = Rz * r_056_2 - 6 * r_055_2;
   FReal2 r_039_1 = Ry * r_029_2 - 2 * r_019_2;
   FReal2 r_931_1 = Rz * r_930_2;
   FReal2 r_913_1 = Ry * r_903_2;
   FReal2 r_553_1 = Rz * r_552_2 - 2 * r_551_2;
   FReal2 r_535_1 = Ry * r_525_2 - 2 * r_515_2;
   FReal2 r_391_1 = Rz * r_390_2;
   FReal2 r_355_1 = Rx * r_255_2 - 2 * r_155_2;
   FReal2 r_319_1 = Ry * r_309_2;
   FReal2 r_193_1 = Rx * r_093_2;
   FReal2 r_139_1 = Rx * r_039_2;

   // m = 0
   pOut[0] = Rx * r_d00_1 - 13 * r_c00_1;
   pOut[1] = Ry * r_c10_1 - r_c00_1;
   pOut[2] = Rx * r_b02_1 - 11 * r_a02_1;
   pOut[3] = Rx * r_940_1 - 9 * r_840_1;
   pOut[4] = Ry * r_a12_1 - r_a02_1;
   pOut[5] = Rx * r_904_1 - 9 * r_804_1;
   pOut[6] = Ry * r_850_1 - 5 * r_840_1;
   pOut[7] = Rz * r_841_1 - r_840_1;
   pOut[8] = Rx * r_724_1 - 7 * r_624_1;
   pOut[9] = Rz * r_805_1 - 5 * r_804_1;
   pOut[10] = Rx * r_580_1 - 5 * r_480_1;
   pOut[11] = Rx * r_562_1 - 5 * r_462_1;
   pOut[12] = Ry * r_634_1 - 3 * r_624_1;
   pOut[13] = Rz * r_625_1 - 5 * r_624_1;
   pOut[14] = Rx * r_508_1 - 5 * r_408_1;
   pOut[15] = Ry * r_490_1 - 9 * r_480_1;
   pOut[16] = Ry * r_472_1 - 7 * r_462_1;
   pOut[17] = Rz * r_463_1 - 3 * r_462_1;
   pOut[18] = Rx * r_346_1 - 3 * r_246_1;
   pOut[19] = Ry * r_418_1 - r_408_1;
   pOut[20] = Rz * r_409_1 - 9 * r_408_1;
   pOut[21] = Rx * r_1c0_1 - r_0c0_1;
   pOut[22] = Rz * r_2a1_1 - r_2a0_1;
   pOut[23] = Rx * r_184_1 - r_084_1;
   pOut[24] = Ry * r_256_1 - 5 * r_246_1;
   pOut[25] = Rz * r_247_1 - 7 * r_246_1;
   pOut[26] = Rx * r_12a_1 - r_02a_1;
   pOut[27] = Rx * r_10c_1 - r_00c_1;
   pOut[28] = Ry * r_0d0_1 - 13 * r_0c0_1;
   pOut[29] = Rz * r_0c1_1 - r_0c0_1;
   pOut[30] = Ry * r_094_1 - 9 * r_084_1;
   pOut[31] = Rz * r_085_1 - 5 * r_084_1;
   pOut[32] = Ry * r_058_1 - 5 * r_048_1;
   pOut[33] = Rz * r_049_1 - 9 * r_048_1;
   pOut[34] = Rz * r_02b_1 - 11 * r_02a_1;
   pOut[35] = Rz * r_00d_1 - 13 * r_00c_1;
   pOut[36] = Ry * r_d00_1;
   pOut[37] = Rx * r_a30_1 - 10 * r_930_1;
   pOut[38] = Ry * r_b02_1;
   pOut[39] = Ry * r_940_1 - 4 * r_930_1;
   pOut[40] = Rz * r_931_1 - r_930_1;
   pOut[41] = Ry * r_904_1;
   pOut[42] = Rx * r_670_1 - 6 * r_570_1;
   pOut[43] = Rx * r_652_1 - 6 * r_552_1;
   pOut[44] = Ry * r_724_1 - 2 * r_714_1;
   pOut[45] = Ry * r_706_1;
   pOut[46] = Rx * r_490_1 - 4 * r_390_1;
   pOut[47] = Ry * r_562_1 - 6 * r_552_1;
   pOut[48] = Rz * r_553_1 - 3 * r_552_1;
   pOut[49] = Ry * r_526_1 - 2 * r_516_1;
   pOut[50] = Ry * r_508_1;
   pOut[51] = Ry * r_3a0_1 - 10 * r_390_1;
   pOut[52] = Rz * r_391_1 - r_390_1;
   pOut[53] = Rx * r_274_1 - 2 * r_174_1;
   pOut[54] = Rx * r_256_1 - 2 * r_156_1;
   pOut[55] = Rx * r_238_1 - 2 * r_138_1;
   pOut[56] = Ry * r_30a_1;
   pOut[57] = Rx * r_0d0_1;
   pOut[58] = Rx * r_0b2_1;
   pOut[59] = Ry * r_184_1 - 8 * r_174_1;
   pOut[60] = Rx * r_076_1;
   pOut[61] = Rx * r_058_1;
   pOut[62] = Rz * r_139_1 - 9 * r_138_1;
   pOut[63] = Rx * r_01c_1;
   pOut[64] = Rz * r_d00_1;
   pOut[65] = Rz * r_b20_1;
   pOut[66] = Rx * r_a03_1 - 10 * r_903_1;
   pOut[67] = Rx * r_841_1 - 8 * r_741_1;
   pOut[68] = Ry * r_913_1 - r_903_1;
   pOut[69] = Rz * r_904_1 - 4 * r_903_1;
   pOut[70] = Rz * r_760_1;
   pOut[71] = Rz * r_742_1 - 2 * r_741_1;
   pOut[72] = Rx * r_625_1 - 6 * r_525_1;
   pOut[73] = Rz * r_706_1 - 6 * r_705_1;
   pOut[74] = Rz * r_580_1;
   pOut[75] = Rz * r_562_1 - 2 * r_561_1;
   pOut[76] = Ry * r_535_1 - 3 * r_525_1;
   pOut[77] = Rz * r_526_1 - 6 * r_525_1;
   pOut[78] = Rx * r_409_1 - 4 * r_309_1;
   pOut[79] = Ry * r_391_1 - 9 * r_381_1;
   pOut[80] = Rz * r_382_1 - 2 * r_381_1;
   pOut[81] = Rx * r_265_1 - 2 * r_165_1;
   pOut[82] = Rx * r_247_1 - 2 * r_147_1;
   pOut[83] = Ry * r_319_1 - r_309_1;
   pOut[84] = Rz * r_30a_1 - 10 * r_309_1;
   pOut[85] = Rz * r_1c0_1;
   pOut[86] = Rx * r_0a3_1;
   pOut[87] = Rx * r_085_1;
   pOut[88] = Rx * r_067_1;
   pOut[89] = Rx * r_049_1;
   pOut[90] = Rx * r_02b_1;
   pOut[91] = Rx * r_00d_1;
   pOut[92] = Ry * r_c01_1;
   pOut[93] = Rz * r_a30_1;
   pOut[94] = Rx * r_913_1 - 9 * r_813_1;
   pOut[95] = Rz * r_850_1;
   pOut[96] = Ry * r_823_1 - 2 * r_813_1;
   pOut[97] = Ry * r_805_1;
   pOut[98] = Rz * r_670_1;
   pOut[99] = Rz * r_652_1 - 2 * r_651_1;
   pOut[100] = Ry * r_625_1 - 2 * r_615_1;
   pOut[101] = Ry * r_607_1;
   pOut[102] = Rz * r_490_1;
   pOut[103] = Rz * r_472_1 - 2 * r_471_1;
   pOut[104] = Rx * r_355_1 - 3 * r_255_1;
   pOut[105] = Ry * r_427_1 - 2 * r_417_1;
   pOut[106] = Rz * r_418_1 - 8 * r_417_1;
   pOut[107] = Rz * r_2b0_1;
   pOut[108] = Rx * r_193_1 - r_093_1;
   pOut[109] = Ry * r_265_1 - 6 * r_255_1;
   pOut[110] = Rz * r_256_1 - 6 * r_255_1;
   pOut[111] = Rx * r_139_1 - r_039_1;
   pOut[112] = Ry * r_20b_1;
   pOut[113] = Rz * r_0d0_1;
   pOut[114] = Ry * r_0a3_1 - 10 * r_093_1;
   pOut[115] = Rz * r_094_1 - 4 * r_093_1;
   pOut[116] = Ry * r_067_1 - 6 * r_057_1;
   pOut[117] = Ry * r_049_1 - 4 * r_039_1;
   pOut[118] = Rz * r_03a_1 - 10 * r_039_1;
   pOut[119] = Ry * r_00d_1;
   // 1013 flops, 1788 mops   (8.44 flops/integral, 14.90 mops/integral)
}


FMdrrFn const
   ShellMdrr[MaxL_MDRR+1] = { &Mdrr0, &Mdrr1, &Mdrr2, &Mdrr3, &Mdrr4, &Mdrr5, &Mdrr6, &Mdrr7, &Mdrr8, &Mdrr9, &Mdrr10, &Mdrr11, &Mdrr12, &Mdrr13, &Mdrr14 };


FScatter2e2cFn const
   Scatter2e2c[2] = {Scatter2e2c_O, Scatter2e2c_A};

// indices of [r0] x [r0] in [R0]
static unsigned short iCartP00[1][1] = {
   {0},
};
// indices of [r0] x [r1] in [R1]
static unsigned short iCartP01[3][1] = {
   {0},{1},{2},
};
// indices of [r0] x [r2] in [R2]
static unsigned short iCartP02[6][1] = {
   {0},{1},{2},{3},{4},{5},
};
// indices of [r0] x [r3] in [R3]
static unsigned short iCartP03[10][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},
};
// indices of [r0] x [r4] in [R4]
static unsigned short iCartP04[15][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},
};
// indices of [r0] x [r5] in [R5]
static unsigned short iCartP05[21][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},
};
// indices of [r0] x [r6] in [R6]
static unsigned short iCartP06[28][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},
};
// indices of [r0] x [r7] in [R7]
static unsigned short iCartP07[36][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},
};
// indices of [r0] x [r8] in [R8]
static unsigned short iCartP08[45][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},
};
// indices of [r0] x [r9] in [R9]
static unsigned short iCartP09[55][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},
};
// indices of [r0] x [r10] in [R10]
static unsigned short iCartP0a[66][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60},{61},{62},{63},{64},{65},
};
// indices of [r0] x [r11] in [R11]
static unsigned short iCartP0b[78][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60},{61},{62},{63},{64},{65},{66},{67},{68},{69},{70},{71},{72},{73},{74},{75},{76},{77},
};
// indices of [r0] x [r12] in [R12]
static unsigned short iCartP0c[91][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60},{61},{62},{63},{64},{65},{66},{67},{68},{69},{70},{71},{72},{73},{74},{75},{76},{77},{78},{79},{80},{81},{82},{83},{84},{85},{86},{87},{88},{89},{90},
};
// indices of [r0] x [r13] in [R13]
static unsigned short iCartP0d[105][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60},{61},{62},{63},{64},{65},{66},{67},{68},{69},{70},{71},{72},{73},{74},{75},{76},{77},{78},{79},{80},{81},{82},{83},{84},{85},{86},{87},{88},{89},{90},{91},{92},{93},{94},{95},{96},{97},{98},{99},{100},{101},{102},{103},{104},
};
// indices of [r0] x [r14] in [R14]
static unsigned short iCartP0e[120][1] = {
   {0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60},{61},{62},{63},{64},{65},{66},{67},{68},{69},{70},{71},{72},{73},{74},{75},{76},{77},{78},{79},{80},{81},{82},{83},{84},{85},{86},{87},{88},{89},{90},{91},{92},{93},{94},{95},{96},{97},{98},{99},{100},{101},{102},{103},{104},{105},{106},{107},{108},{109},{110},{111},{112},{113},{114},{115},{116},{117},{118},{119}
};
// indices of [r1] x [r0] in [R1]
static unsigned short iCartP10[1][3] = {
   {0,1,2},
};
// indices of [r1] x [r1] in [R2]
static unsigned short iCartP11[3][3] = {
   {0,3,4},
   {3,1,5},
   {4,5,2},
};
// indices of [r1] x [r2] in [R3]
static unsigned short iCartP12[6][3] = {
   {0,3,6},
   {1,4,7},
   {2,5,8},
   {3,1,9},
   {6,9,2},
   {9,7,5},
};
// indices of [r1] x [r3] in [R4]
static unsigned short iCartP13[10][3] = {
   {0,6,9},
   {1,7,10},
   {2,8,11},
   {6,1,12},
   {7,3,13},
   {8,4,14},
   {9,12,2},
   {10,13,4},
   {11,14,5},
   {12,10,8},
};
// indices of [r1] x [r4] in [R5]
static unsigned short iCartP14[15][3] = {
   {0,6,12},
   {1,7,13},
   {2,8,14},
   {3,9,15},
   {4,10,16},
   {5,11,17},
   {6,1,18},
   {7,3,19},
   {8,4,20},
   {12,18,2},
   {13,19,4},
   {14,20,5},
   {18,13,8},
   {19,15,10},
   {20,16,11},
};
// indices of [r1] x [r5] in [R6]
static unsigned short iCartP15[21][3] = {
   {0,10,16},
   {1,11,17},
   {2,12,18},
   {3,13,19},
   {4,14,20},
   {5,15,21},
   {10,1,22},
   {11,3,23},
   {12,4,24},
   {13,6,25},
   {14,7,26},
   {15,8,27},
   {16,22,2},
   {17,23,4},
   {18,24,5},
   {19,25,7},
   {20,26,8},
   {21,27,9},
   {22,17,12},
   {23,19,14},
   {24,20,15},
};
// indices of [r1] x [r6] in [R7]
static unsigned short iCartP16[28][3] = {
   {0,10,20},
   {1,11,21},
   {2,12,22},
   {3,13,23},
   {4,14,24},
   {5,15,25},
   {6,16,26},
   {7,17,27},
   {8,18,28},
   {9,19,29},
   {10,1,30},
   {11,3,31},
   {12,4,32},
   {13,6,33},
   {14,7,34},
   {15,8,35},
   {20,30,2},
   {21,31,4},
   {22,32,5},
   {23,33,7},
   {24,34,8},
   {25,35,9},
   {30,21,12},
   {31,23,14},
   {32,24,15},
   {33,26,17},
   {34,27,18},
   {35,28,19},
};
// indices of [r1] x [r7] in [R8]
static unsigned short iCartP17[36][3] = {
   {0,15,25},
   {1,16,26},
   {2,17,27},
   {3,18,28},
   {4,19,29},
   {5,20,30},
   {6,21,31},
   {7,22,32},
   {8,23,33},
   {9,24,34},
   {15,1,35},
   {16,3,36},
   {17,4,37},
   {18,6,38},
   {19,7,39},
   {20,8,40},
   {21,10,41},
   {22,11,42},
   {23,12,43},
   {24,13,44},
   {25,35,2},
   {26,36,4},
   {27,37,5},
   {28,38,7},
   {29,39,8},
   {30,40,9},
   {31,41,11},
   {32,42,12},
   {33,43,13},
   {34,44,14},
   {35,26,17},
   {36,28,19},
   {37,29,20},
   {38,31,22},
   {39,32,23},
   {40,33,24},
};
// indices of [r1] x [r8] in [R9]
static unsigned short iCartP18[45][3] = {
   {0,15,30},
   {1,16,31},
   {2,17,32},
   {3,18,33},
   {4,19,34},
   {5,20,35},
   {6,21,36},
   {7,22,37},
   {8,23,38},
   {9,24,39},
   {10,25,40},
   {11,26,41},
   {12,27,42},
   {13,28,43},
   {14,29,44},
   {15,1,45},
   {16,3,46},
   {17,4,47},
   {18,6,48},
   {19,7,49},
   {20,8,50},
   {21,10,51},
   {22,11,52},
   {23,12,53},
   {24,13,54},
   {30,45,2},
   {31,46,4},
   {32,47,5},
   {33,48,7},
   {34,49,8},
   {35,50,9},
   {36,51,11},
   {37,52,12},
   {38,53,13},
   {39,54,14},
   {45,31,17},
   {46,33,19},
   {47,34,20},
   {48,36,22},
   {49,37,23},
   {50,38,24},
   {51,40,26},
   {52,41,27},
   {53,42,28},
   {54,43,29},
};
// indices of [r1] x [r9] in [R10]
static unsigned short iCartP19[55][3] = {
   {0,21,36},
   {1,22,37},
   {2,23,38},
   {3,24,39},
   {4,25,40},
   {5,26,41},
   {6,27,42},
   {7,28,43},
   {8,29,44},
   {9,30,45},
   {10,31,46},
   {11,32,47},
   {12,33,48},
   {13,34,49},
   {14,35,50},
   {21,1,51},
   {22,3,52},
   {23,4,53},
   {24,6,54},
   {25,7,55},
   {26,8,56},
   {27,10,57},
   {28,11,58},
   {29,12,59},
   {30,13,60},
   {31,15,61},
   {32,16,62},
   {33,17,63},
   {34,18,64},
   {35,19,65},
   {36,51,2},
   {37,52,4},
   {38,53,5},
   {39,54,7},
   {40,55,8},
   {41,56,9},
   {42,57,11},
   {43,58,12},
   {44,59,13},
   {45,60,14},
   {46,61,16},
   {47,62,17},
   {48,63,18},
   {49,64,19},
   {50,65,20},
   {51,37,23},
   {52,39,25},
   {53,40,26},
   {54,42,28},
   {55,43,29},
   {56,44,30},
   {57,46,32},
   {58,47,33},
   {59,48,34},
   {60,49,35},
};
// indices of [r1] x [r10] in [R11]
static unsigned short iCartP1a[66][3] = {
   {0,21,42},
   {1,22,43},
   {2,23,44},
   {3,24,45},
   {4,25,46},
   {5,26,47},
   {6,27,48},
   {7,28,49},
   {8,29,50},
   {9,30,51},
   {10,31,52},
   {11,32,53},
   {12,33,54},
   {13,34,55},
   {14,35,56},
   {15,36,57},
   {16,37,58},
   {17,38,59},
   {18,39,60},
   {19,40,61},
   {20,41,62},
   {21,1,63},
   {22,3,64},
   {23,4,65},
   {24,6,66},
   {25,7,67},
   {26,8,68},
   {27,10,69},
   {28,11,70},
   {29,12,71},
   {30,13,72},
   {31,15,73},
   {32,16,74},
   {33,17,75},
   {34,18,76},
   {35,19,77},
   {42,63,2},
   {43,64,4},
   {44,65,5},
   {45,66,7},
   {46,67,8},
   {47,68,9},
   {48,69,11},
   {49,70,12},
   {50,71,13},
   {51,72,14},
   {52,73,16},
   {53,74,17},
   {54,75,18},
   {55,76,19},
   {56,77,20},
   {63,43,23},
   {64,45,25},
   {65,46,26},
   {66,48,28},
   {67,49,29},
   {68,50,30},
   {69,52,32},
   {70,53,33},
   {71,54,34},
   {72,55,35},
   {73,57,37},
   {74,58,38},
   {75,59,39},
   {76,60,40},
   {77,61,41},
};
// indices of [r1] x [r11] in [R12]
static unsigned short iCartP1b[78][3] = {
   {0,28,49},
   {1,29,50},
   {2,30,51},
   {3,31,52},
   {4,32,53},
   {5,33,54},
   {6,34,55},
   {7,35,56},
   {8,36,57},
   {9,37,58},
   {10,38,59},
   {11,39,60},
   {12,40,61},
   {13,41,62},
   {14,42,63},
   {15,43,64},
   {16,44,65},
   {17,45,66},
   {18,46,67},
   {19,47,68},
   {20,48,69},
   {28,1,70},
   {29,3,71},
   {30,4,72},
   {31,6,73},
   {32,7,74},
   {33,8,75},
   {34,10,76},
   {35,11,77},
   {36,12,78},
   {37,13,79},
   {38,15,80},
   {39,16,81},
   {40,17,82},
   {41,18,83},
   {42,19,84},
   {43,21,85},
   {44,22,86},
   {45,23,87},
   {46,24,88},
   {47,25,89},
   {48,26,90},
   {49,70,2},
   {50,71,4},
   {51,72,5},
   {52,73,7},
   {53,74,8},
   {54,75,9},
   {55,76,11},
   {56,77,12},
   {57,78,13},
   {58,79,14},
   {59,80,16},
   {60,81,17},
   {61,82,18},
   {62,83,19},
   {63,84,20},
   {64,85,22},
   {65,86,23},
   {66,87,24},
   {67,88,25},
   {68,89,26},
   {69,90,27},
   {70,50,30},
   {71,52,32},
   {72,53,33},
   {73,55,35},
   {74,56,36},
   {75,57,37},
   {76,59,39},
   {77,60,40},
   {78,61,41},
   {79,62,42},
   {80,64,44},
   {81,65,45},
   {82,66,46},
   {83,67,47},
   {84,68,48},
};
// indices of [r1] x [r12] in [R13]
static unsigned short iCartP1c[91][3] = {
   {0,28,56},
   {1,29,57},
   {2,30,58},
   {3,31,59},
   {4,32,60},
   {5,33,61},
   {6,34,62},
   {7,35,63},
   {8,36,64},
   {9,37,65},
   {10,38,66},
   {11,39,67},
   {12,40,68},
   {13,41,69},
   {14,42,70},
   {15,43,71},
   {16,44,72},
   {17,45,73},
   {18,46,74},
   {19,47,75},
   {20,48,76},
   {21,49,77},
   {22,50,78},
   {23,51,79},
   {24,52,80},
   {25,53,81},
   {26,54,82},
   {27,55,83},
   {28,1,84},
   {29,3,85},
   {30,4,86},
   {31,6,87},
   {32,7,88},
   {33,8,89},
   {34,10,90},
   {35,11,91},
   {36,12,92},
   {37,13,93},
   {38,15,94},
   {39,16,95},
   {40,17,96},
   {41,18,97},
   {42,19,98},
   {43,21,99},
   {44,22,100},
   {45,23,101},
   {46,24,102},
   {47,25,103},
   {48,26,104},
   {56,84,2},
   {57,85,4},
   {58,86,5},
   {59,87,7},
   {60,88,8},
   {61,89,9},
   {62,90,11},
   {63,91,12},
   {64,92,13},
   {65,93,14},
   {66,94,16},
   {67,95,17},
   {68,96,18},
   {69,97,19},
   {70,98,20},
   {71,99,22},
   {72,100,23},
   {73,101,24},
   {74,102,25},
   {75,103,26},
   {76,104,27},
   {84,57,30},
   {85,59,32},
   {86,60,33},
   {87,62,35},
   {88,63,36},
   {89,64,37},
   {90,66,39},
   {91,67,40},
   {92,68,41},
   {93,69,42},
   {94,71,44},
   {95,72,45},
   {96,73,46},
   {97,74,47},
   {98,75,48},
   {99,77,50},
   {100,78,51},
   {101,79,52},
   {102,80,53},
   {103,81,54},
   {104,82,55},
};
// indices of [r1] x [r13] in [R14]
static unsigned short iCartP1d[105][3] = {
   {0,36,64},
   {1,37,65},
   {2,38,66},
   {3,39,67},
   {4,40,68},
   {5,41,69},
   {6,42,70},
   {7,43,71},
   {8,44,72},
   {9,45,73},
   {10,46,74},
   {11,47,75},
   {12,48,76},
   {13,49,77},
   {14,50,78},
   {15,51,79},
   {16,52,80},
   {17,53,81},
   {18,54,82},
   {19,55,83},
   {20,56,84},
   {21,57,85},
   {22,58,86},
   {23,59,87},
   {24,60,88},
   {25,61,89},
   {26,62,90},
   {27,63,91},
   {36,1,92},
   {37,3,93},
   {38,4,94},
   {39,6,95},
   {40,7,96},
   {41,8,97},
   {42,10,98},
   {43,11,99},
   {44,12,100},
   {45,13,101},
   {46,15,102},
   {47,16,103},
   {48,17,104},
   {49,18,105},
   {50,19,106},
   {51,21,107},
   {52,22,108},
   {53,23,109},
   {54,24,110},
   {55,25,111},
   {56,26,112},
   {57,28,113},
   {58,29,114},
   {59,30,115},
   {60,31,116},
   {61,32,117},
   {62,33,118},
   {63,34,119},
   {64,92,2},
   {65,93,4},
   {66,94,5},
   {67,95,7},
   {68,96,8},
   {69,97,9},
   {70,98,11},
   {71,99,12},
   {72,100,13},
   {73,101,14},
   {74,102,16},
   {75,103,17},
   {76,104,18},
   {77,105,19},
   {78,106,20},
   {79,107,22},
   {80,108,23},
   {81,109,24},
   {82,110,25},
   {83,111,26},
   {84,112,27},
   {85,113,29},
   {86,114,30},
   {87,115,31},
   {88,116,32},
   {89,117,33},
   {90,118,34},
   {91,119,35},
   {92,65,38},
   {93,67,40},
   {94,68,41},
   {95,70,43},
   {96,71,44},
   {97,72,45},
   {98,74,47},
   {99,75,48},
   {100,76,49},
   {101,77,50},
   {102,79,52},
   {103,80,53},
   {104,81,54},
   {105,82,55},
   {106,83,56},
   {107,85,58},
   {108,86,59},
   {109,87,60},
   {110,88,61},
   {111,89,62},
   {112,90,63}
};
// indices of [r2] x [r0] in [R2]
static unsigned short iCartP20[1][6] = {
   {0,1,2,3,4,5},
};
// indices of [r2] x [r1] in [R3]
static unsigned short iCartP21[3][6] = {
   {0,1,2,3,6,9},
   {3,4,5,1,9,7},
   {6,7,8,9,2,5},
};
// indices of [r2] x [r2] in [R4]
static unsigned short iCartP22[6][6] = {
   {0,1,2,6,9,12},
   {1,3,4,7,10,13},
   {2,4,5,8,11,14},
   {6,7,8,1,12,10},
   {9,10,11,12,2,8},
   {12,13,14,10,8,4},
};
// indices of [r2] x [r3] in [R5]
static unsigned short iCartP23[10][6] = {
   {0,1,2,6,12,18},
   {1,3,4,7,13,19},
   {2,4,5,8,14,20},
   {6,7,8,1,18,13},
   {7,9,10,3,19,15},
   {8,10,11,4,20,16},
   {12,13,14,18,2,8},
   {13,15,16,19,4,10},
   {14,16,17,20,5,11},
   {18,19,20,13,8,4},
};
// indices of [r2] x [r4] in [R6]
static unsigned short iCartP24[15][6] = {
   {0,1,2,10,16,22},
   {1,3,4,11,17,23},
   {2,4,5,12,18,24},
   {3,6,7,13,19,25},
   {4,7,8,14,20,26},
   {5,8,9,15,21,27},
   {10,11,12,1,22,17},
   {11,13,14,3,23,19},
   {12,14,15,4,24,20},
   {16,17,18,22,2,12},
   {17,19,20,23,4,14},
   {18,20,21,24,5,15},
   {22,23,24,17,12,4},
   {23,25,26,19,14,7},
   {24,26,27,20,15,8},
};
// indices of [r2] x [r5] in [R7]
static unsigned short iCartP25[21][6] = {
   {0,1,2,10,20,30},
   {1,3,4,11,21,31},
   {2,4,5,12,22,32},
   {3,6,7,13,23,33},
   {4,7,8,14,24,34},
   {5,8,9,15,25,35},
   {10,11,12,1,30,21},
   {11,13,14,3,31,23},
   {12,14,15,4,32,24},
   {13,16,17,6,33,26},
   {14,17,18,7,34,27},
   {15,18,19,8,35,28},
   {20,21,22,30,2,12},
   {21,23,24,31,4,14},
   {22,24,25,32,5,15},
   {23,26,27,33,7,17},
   {24,27,28,34,8,18},
   {25,28,29,35,9,19},
   {30,31,32,21,12,4},
   {31,33,34,23,14,7},
   {32,34,35,24,15,8},
};
// indices of [r2] x [r6] in [R8]
static unsigned short iCartP26[28][6] = {
   {0,1,2,15,25,35},
   {1,3,4,16,26,36},
   {2,4,5,17,27,37},
   {3,6,7,18,28,38},
   {4,7,8,19,29,39},
   {5,8,9,20,30,40},
   {6,10,11,21,31,41},
   {7,11,12,22,32,42},
   {8,12,13,23,33,43},
   {9,13,14,24,34,44},
   {15,16,17,1,35,26},
   {16,18,19,3,36,28},
   {17,19,20,4,37,29},
   {18,21,22,6,38,31},
   {19,22,23,7,39,32},
   {20,23,24,8,40,33},
   {25,26,27,35,2,17},
   {26,28,29,36,4,19},
   {27,29,30,37,5,20},
   {28,31,32,38,7,22},
   {29,32,33,39,8,23},
   {30,33,34,40,9,24},
   {35,36,37,26,17,4},
   {36,38,39,28,19,7},
   {37,39,40,29,20,8},
   {38,41,42,31,22,11},
   {39,42,43,32,23,12},
   {40,43,44,33,24,13},
};
// indices of [r2] x [r7] in [R9]
static unsigned short iCartP27[36][6] = {
   {0,1,2,15,30,45},
   {1,3,4,16,31,46},
   {2,4,5,17,32,47},
   {3,6,7,18,33,48},
   {4,7,8,19,34,49},
   {5,8,9,20,35,50},
   {6,10,11,21,36,51},
   {7,11,12,22,37,52},
   {8,12,13,23,38,53},
   {9,13,14,24,39,54},
   {15,16,17,1,45,31},
   {16,18,19,3,46,33},
   {17,19,20,4,47,34},
   {18,21,22,6,48,36},
   {19,22,23,7,49,37},
   {20,23,24,8,50,38},
   {21,25,26,10,51,40},
   {22,26,27,11,52,41},
   {23,27,28,12,53,42},
   {24,28,29,13,54,43},
   {30,31,32,45,2,17},
   {31,33,34,46,4,19},
   {32,34,35,47,5,20},
   {33,36,37,48,7,22},
   {34,37,38,49,8,23},
   {35,38,39,50,9,24},
   {36,40,41,51,11,26},
   {37,41,42,52,12,27},
   {38,42,43,53,13,28},
   {39,43,44,54,14,29},
   {45,46,47,31,17,4},
   {46,48,49,33,19,7},
   {47,49,50,34,20,8},
   {48,51,52,36,22,11},
   {49,52,53,37,23,12},
   {50,53,54,38,24,13},
};
// indices of [r2] x [r8] in [R10]
static unsigned short iCartP28[45][6] = {
   {0,1,2,21,36,51},
   {1,3,4,22,37,52},
   {2,4,5,23,38,53},
   {3,6,7,24,39,54},
   {4,7,8,25,40,55},
   {5,8,9,26,41,56},
   {6,10,11,27,42,57},
   {7,11,12,28,43,58},
   {8,12,13,29,44,59},
   {9,13,14,30,45,60},
   {10,15,16,31,46,61},
   {11,16,17,32,47,62},
   {12,17,18,33,48,63},
   {13,18,19,34,49,64},
   {14,19,20,35,50,65},
   {21,22,23,1,51,37},
   {22,24,25,3,52,39},
   {23,25,26,4,53,40},
   {24,27,28,6,54,42},
   {25,28,29,7,55,43},
   {26,29,30,8,56,44},
   {27,31,32,10,57,46},
   {28,32,33,11,58,47},
   {29,33,34,12,59,48},
   {30,34,35,13,60,49},
   {36,37,38,51,2,23},
   {37,39,40,52,4,25},
   {38,40,41,53,5,26},
   {39,42,43,54,7,28},
   {40,43,44,55,8,29},
   {41,44,45,56,9,30},
   {42,46,47,57,11,32},
   {43,47,48,58,12,33},
   {44,48,49,59,13,34},
   {45,49,50,60,14,35},
   {51,52,53,37,23,4},
   {52,54,55,39,25,7},
   {53,55,56,40,26,8},
   {54,57,58,42,28,11},
   {55,58,59,43,29,12},
   {56,59,60,44,30,13},
   {57,61,62,46,32,16},
   {58,62,63,47,33,17},
   {59,63,64,48,34,18},
   {60,64,65,49,35,19},
};
// indices of [r2] x [r9] in [R11]
static unsigned short iCartP29[55][6] = {
   {0,1,2,21,42,63},
   {1,3,4,22,43,64},
   {2,4,5,23,44,65},
   {3,6,7,24,45,66},
   {4,7,8,25,46,67},
   {5,8,9,26,47,68},
   {6,10,11,27,48,69},
   {7,11,12,28,49,70},
   {8,12,13,29,50,71},
   {9,13,14,30,51,72},
   {10,15,16,31,52,73},
   {11,16,17,32,53,74},
   {12,17,18,33,54,75},
   {13,18,19,34,55,76},
   {14,19,20,35,56,77},
   {21,22,23,1,63,43},
   {22,24,25,3,64,45},
   {23,25,26,4,65,46},
   {24,27,28,6,66,48},
   {25,28,29,7,67,49},
   {26,29,30,8,68,50},
   {27,31,32,10,69,52},
   {28,32,33,11,70,53},
   {29,33,34,12,71,54},
   {30,34,35,13,72,55},
   {31,36,37,15,73,57},
   {32,37,38,16,74,58},
   {33,38,39,17,75,59},
   {34,39,40,18,76,60},
   {35,40,41,19,77,61},
   {42,43,44,63,2,23},
   {43,45,46,64,4,25},
   {44,46,47,65,5,26},
   {45,48,49,66,7,28},
   {46,49,50,67,8,29},
   {47,50,51,68,9,30},
   {48,52,53,69,11,32},
   {49,53,54,70,12,33},
   {50,54,55,71,13,34},
   {51,55,56,72,14,35},
   {52,57,58,73,16,37},
   {53,58,59,74,17,38},
   {54,59,60,75,18,39},
   {55,60,61,76,19,40},
   {56,61,62,77,20,41},
   {63,64,65,43,23,4},
   {64,66,67,45,25,7},
   {65,67,68,46,26,8},
   {66,69,70,48,28,11},
   {67,70,71,49,29,12},
   {68,71,72,50,30,13},
   {69,73,74,52,32,16},
   {70,74,75,53,33,17},
   {71,75,76,54,34,18},
   {72,76,77,55,35,19},
};
// indices of [r2] x [r10] in [R12]
static unsigned short iCartP2a[66][6] = {
   {0,1,2,28,49,70},
   {1,3,4,29,50,71},
   {2,4,5,30,51,72},
   {3,6,7,31,52,73},
   {4,7,8,32,53,74},
   {5,8,9,33,54,75},
   {6,10,11,34,55,76},
   {7,11,12,35,56,77},
   {8,12,13,36,57,78},
   {9,13,14,37,58,79},
   {10,15,16,38,59,80},
   {11,16,17,39,60,81},
   {12,17,18,40,61,82},
   {13,18,19,41,62,83},
   {14,19,20,42,63,84},
   {15,21,22,43,64,85},
   {16,22,23,44,65,86},
   {17,23,24,45,66,87},
   {18,24,25,46,67,88},
   {19,25,26,47,68,89},
   {20,26,27,48,69,90},
   {28,29,30,1,70,50},
   {29,31,32,3,71,52},
   {30,32,33,4,72,53},
   {31,34,35,6,73,55},
   {32,35,36,7,74,56},
   {33,36,37,8,75,57},
   {34,38,39,10,76,59},
   {35,39,40,11,77,60},
   {36,40,41,12,78,61},
   {37,41,42,13,79,62},
   {38,43,44,15,80,64},
   {39,44,45,16,81,65},
   {40,45,46,17,82,66},
   {41,46,47,18,83,67},
   {42,47,48,19,84,68},
   {49,50,51,70,2,30},
   {50,52,53,71,4,32},
   {51,53,54,72,5,33},
   {52,55,56,73,7,35},
   {53,56,57,74,8,36},
   {54,57,58,75,9,37},
   {55,59,60,76,11,39},
   {56,60,61,77,12,40},
   {57,61,62,78,13,41},
   {58,62,63,79,14,42},
   {59,64,65,80,16,44},
   {60,65,66,81,17,45},
   {61,66,67,82,18,46},
   {62,67,68,83,19,47},
   {63,68,69,84,20,48},
   {70,71,72,50,30,4},
   {71,73,74,52,32,7},
   {72,74,75,53,33,8},
   {73,76,77,55,35,11},
   {74,77,78,56,36,12},
   {75,78,79,57,37,13},
   {76,80,81,59,39,16},
   {77,81,82,60,40,17},
   {78,82,83,61,41,18},
   {79,83,84,62,42,19},
   {80,85,86,64,44,22},
   {81,86,87,65,45,23},
   {82,87,88,66,46,24},
   {83,88,89,67,47,25},
   {84,89,90,68,48,26},
};
// indices of [r2] x [r11] in [R13]
static unsigned short iCartP2b[78][6] = {
   {0,1,2,28,56,84},
   {1,3,4,29,57,85},
   {2,4,5,30,58,86},
   {3,6,7,31,59,87},
   {4,7,8,32,60,88},
   {5,8,9,33,61,89},
   {6,10,11,34,62,90},
   {7,11,12,35,63,91},
   {8,12,13,36,64,92},
   {9,13,14,37,65,93},
   {10,15,16,38,66,94},
   {11,16,17,39,67,95},
   {12,17,18,40,68,96},
   {13,18,19,41,69,97},
   {14,19,20,42,70,98},
   {15,21,22,43,71,99},
   {16,22,23,44,72,100},
   {17,23,24,45,73,101},
   {18,24,25,46,74,102},
   {19,25,26,47,75,103},
   {20,26,27,48,76,104},
   {28,29,30,1,84,57},
   {29,31,32,3,85,59},
   {30,32,33,4,86,60},
   {31,34,35,6,87,62},
   {32,35,36,7,88,63},
   {33,36,37,8,89,64},
   {34,38,39,10,90,66},
   {35,39,40,11,91,67},
   {36,40,41,12,92,68},
   {37,41,42,13,93,69},
   {38,43,44,15,94,71},
   {39,44,45,16,95,72},
   {40,45,46,17,96,73},
   {41,46,47,18,97,74},
   {42,47,48,19,98,75},
   {43,49,50,21,99,77},
   {44,50,51,22,100,78},
   {45,51,52,23,101,79},
   {46,52,53,24,102,80},
   {47,53,54,25,103,81},
   {48,54,55,26,104,82},
   {56,57,58,84,2,30},
   {57,59,60,85,4,32},
   {58,60,61,86,5,33},
   {59,62,63,87,7,35},
   {60,63,64,88,8,36},
   {61,64,65,89,9,37},
   {62,66,67,90,11,39},
   {63,67,68,91,12,40},
   {64,68,69,92,13,41},
   {65,69,70,93,14,42},
   {66,71,72,94,16,44},
   {67,72,73,95,17,45},
   {68,73,74,96,18,46},
   {69,74,75,97,19,47},
   {70,75,76,98,20,48},
   {71,77,78,99,22,50},
   {72,78,79,100,23,51},
   {73,79,80,101,24,52},
   {74,80,81,102,25,53},
   {75,81,82,103,26,54},
   {76,82,83,104,27,55},
   {84,85,86,57,30,4},
   {85,87,88,59,32,7},
   {86,88,89,60,33,8},
   {87,90,91,62,35,11},
   {88,91,92,63,36,12},
   {89,92,93,64,37,13},
   {90,94,95,66,39,16},
   {91,95,96,67,40,17},
   {92,96,97,68,41,18},
   {93,97,98,69,42,19},
   {94,99,100,71,44,22},
   {95,100,101,72,45,23},
   {96,101,102,73,46,24},
   {97,102,103,74,47,25},
   {98,103,104,75,48,26},
};
// indices of [r2] x [r12] in [R14]
static unsigned short iCartP2c[91][6] = {
   {0,1,2,36,64,92},
   {1,3,4,37,65,93},
   {2,4,5,38,66,94},
   {3,6,7,39,67,95},
   {4,7,8,40,68,96},
   {5,8,9,41,69,97},
   {6,10,11,42,70,98},
   {7,11,12,43,71,99},
   {8,12,13,44,72,100},
   {9,13,14,45,73,101},
   {10,15,16,46,74,102},
   {11,16,17,47,75,103},
   {12,17,18,48,76,104},
   {13,18,19,49,77,105},
   {14,19,20,50,78,106},
   {15,21,22,51,79,107},
   {16,22,23,52,80,108},
   {17,23,24,53,81,109},
   {18,24,25,54,82,110},
   {19,25,26,55,83,111},
   {20,26,27,56,84,112},
   {21,28,29,57,85,113},
   {22,29,30,58,86,114},
   {23,30,31,59,87,115},
   {24,31,32,60,88,116},
   {25,32,33,61,89,117},
   {26,33,34,62,90,118},
   {27,34,35,63,91,119},
   {36,37,38,1,92,65},
   {37,39,40,3,93,67},
   {38,40,41,4,94,68},
   {39,42,43,6,95,70},
   {40,43,44,7,96,71},
   {41,44,45,8,97,72},
   {42,46,47,10,98,74},
   {43,47,48,11,99,75},
   {44,48,49,12,100,76},
   {45,49,50,13,101,77},
   {46,51,52,15,102,79},
   {47,52,53,16,103,80},
   {48,53,54,17,104,81},
   {49,54,55,18,105,82},
   {50,55,56,19,106,83},
   {51,57,58,21,107,85},
   {52,58,59,22,108,86},
   {53,59,60,23,109,87},
   {54,60,61,24,110,88},
   {55,61,62,25,111,89},
   {56,62,63,26,112,90},
   {64,65,66,92,2,38},
   {65,67,68,93,4,40},
   {66,68,69,94,5,41},
   {67,70,71,95,7,43},
   {68,71,72,96,8,44},
   {69,72,73,97,9,45},
   {70,74,75,98,11,47},
   {71,75,76,99,12,48},
   {72,76,77,100,13,49},
   {73,77,78,101,14,50},
   {74,79,80,102,16,52},
   {75,80,81,103,17,53},
   {76,81,82,104,18,54},
   {77,82,83,105,19,55},
   {78,83,84,106,20,56},
   {79,85,86,107,22,58},
   {80,86,87,108,23,59},
   {81,87,88,109,24,60},
   {82,88,89,110,25,61},
   {83,89,90,111,26,62},
   {84,90,91,112,27,63},
   {92,93,94,65,38,4},
   {93,95,96,67,40,7},
   {94,96,97,68,41,8},
   {95,98,99,70,43,11},
   {96,99,100,71,44,12},
   {97,100,101,72,45,13},
   {98,102,103,74,47,16},
   {99,103,104,75,48,17},
   {100,104,105,76,49,18},
   {101,105,106,77,50,19},
   {102,107,108,79,52,22},
   {103,108,109,80,53,23},
   {104,109,110,81,54,24},
   {105,110,111,82,55,25},
   {106,111,112,83,56,26},
   {107,113,114,85,58,29},
   {108,114,115,86,59,30},
   {109,115,116,87,60,31},
   {110,116,117,88,61,32},
   {111,117,118,89,62,33},
   {112,118,119,90,63,34}
};
// indices of [r3] x [r0] in [R3]
static unsigned short iCartP30[1][10] = {
   {0,1,2,3,4,5,6,7,8,9},
};
// indices of [r3] x [r1] in [R4]
static unsigned short iCartP31[3][10] = {
   {0,1,2,6,7,8,9,10,11,12},
   {6,7,8,1,3,4,12,13,14,10},
   {9,10,11,12,13,14,2,4,5,8},
};
// indices of [r3] x [r2] in [R5]
static unsigned short iCartP32[6][10] = {
   {0,1,2,6,7,8,12,13,14,18},
   {1,3,4,7,9,10,13,15,16,19},
   {2,4,5,8,10,11,14,16,17,20},
   {6,7,8,1,3,4,18,19,20,13},
   {12,13,14,18,19,20,2,4,5,8},
   {18,19,20,13,15,16,8,10,11,4},
};
// indices of [r3] x [r3] in [R6]
static unsigned short iCartP33[10][10] = {
   {0,1,2,10,11,12,16,17,18,22},
   {1,3,4,11,13,14,17,19,20,23},
   {2,4,5,12,14,15,18,20,21,24},
   {10,11,12,1,3,4,22,23,24,17},
   {11,13,14,3,6,7,23,25,26,19},
   {12,14,15,4,7,8,24,26,27,20},
   {16,17,18,22,23,24,2,4,5,12},
   {17,19,20,23,25,26,4,7,8,14},
   {18,20,21,24,26,27,5,8,9,15},
   {22,23,24,17,19,20,12,14,15,4},
};
// indices of [r3] x [r4] in [R7]
static unsigned short iCartP34[15][10] = {
   {0,1,2,10,11,12,20,21,22,30},
   {1,3,4,11,13,14,21,23,24,31},
   {2,4,5,12,14,15,22,24,25,32},
   {3,6,7,13,16,17,23,26,27,33},
   {4,7,8,14,17,18,24,27,28,34},
   {5,8,9,15,18,19,25,28,29,35},
   {10,11,12,1,3,4,30,31,32,21},
   {11,13,14,3,6,7,31,33,34,23},
   {12,14,15,4,7,8,32,34,35,24},
   {20,21,22,30,31,32,2,4,5,12},
   {21,23,24,31,33,34,4,7,8,14},
   {22,24,25,32,34,35,5,8,9,15},
   {30,31,32,21,23,24,12,14,15,4},
   {31,33,34,23,26,27,14,17,18,7},
   {32,34,35,24,27,28,15,18,19,8},
};
// indices of [r3] x [r5] in [R8]
static unsigned short iCartP35[21][10] = {
   {0,1,2,15,16,17,25,26,27,35},
   {1,3,4,16,18,19,26,28,29,36},
   {2,4,5,17,19,20,27,29,30,37},
   {3,6,7,18,21,22,28,31,32,38},
   {4,7,8,19,22,23,29,32,33,39},
   {5,8,9,20,23,24,30,33,34,40},
   {15,16,17,1,3,4,35,36,37,26},
   {16,18,19,3,6,7,36,38,39,28},
   {17,19,20,4,7,8,37,39,40,29},
   {18,21,22,6,10,11,38,41,42,31},
   {19,22,23,7,11,12,39,42,43,32},
   {20,23,24,8,12,13,40,43,44,33},
   {25,26,27,35,36,37,2,4,5,17},
   {26,28,29,36,38,39,4,7,8,19},
   {27,29,30,37,39,40,5,8,9,20},
   {28,31,32,38,41,42,7,11,12,22},
   {29,32,33,39,42,43,8,12,13,23},
   {30,33,34,40,43,44,9,13,14,24},
   {35,36,37,26,28,29,17,19,20,4},
   {36,38,39,28,31,32,19,22,23,7},
   {37,39,40,29,32,33,20,23,24,8},
};
// indices of [r3] x [r6] in [R9]
static unsigned short iCartP36[28][10] = {
   {0,1,2,15,16,17,30,31,32,45},
   {1,3,4,16,18,19,31,33,34,46},
   {2,4,5,17,19,20,32,34,35,47},
   {3,6,7,18,21,22,33,36,37,48},
   {4,7,8,19,22,23,34,37,38,49},
   {5,8,9,20,23,24,35,38,39,50},
   {6,10,11,21,25,26,36,40,41,51},
   {7,11,12,22,26,27,37,41,42,52},
   {8,12,13,23,27,28,38,42,43,53},
   {9,13,14,24,28,29,39,43,44,54},
   {15,16,17,1,3,4,45,46,47,31},
   {16,18,19,3,6,7,46,48,49,33},
   {17,19,20,4,7,8,47,49,50,34},
   {18,21,22,6,10,11,48,51,52,36},
   {19,22,23,7,11,12,49,52,53,37},
   {20,23,24,8,12,13,50,53,54,38},
   {30,31,32,45,46,47,2,4,5,17},
   {31,33,34,46,48,49,4,7,8,19},
   {32,34,35,47,49,50,5,8,9,20},
   {33,36,37,48,51,52,7,11,12,22},
   {34,37,38,49,52,53,8,12,13,23},
   {35,38,39,50,53,54,9,13,14,24},
   {45,46,47,31,33,34,17,19,20,4},
   {46,48,49,33,36,37,19,22,23,7},
   {47,49,50,34,37,38,20,23,24,8},
   {48,51,52,36,40,41,22,26,27,11},
   {49,52,53,37,41,42,23,27,28,12},
   {50,53,54,38,42,43,24,28,29,13},
};
// indices of [r3] x [r7] in [R10]
static unsigned short iCartP37[36][10] = {
   {0,1,2,21,22,23,36,37,38,51},
   {1,3,4,22,24,25,37,39,40,52},
   {2,4,5,23,25,26,38,40,41,53},
   {3,6,7,24,27,28,39,42,43,54},
   {4,7,8,25,28,29,40,43,44,55},
   {5,8,9,26,29,30,41,44,45,56},
   {6,10,11,27,31,32,42,46,47,57},
   {7,11,12,28,32,33,43,47,48,58},
   {8,12,13,29,33,34,44,48,49,59},
   {9,13,14,30,34,35,45,49,50,60},
   {21,22,23,1,3,4,51,52,53,37},
   {22,24,25,3,6,7,52,54,55,39},
   {23,25,26,4,7,8,53,55,56,40},
   {24,27,28,6,10,11,54,57,58,42},
   {25,28,29,7,11,12,55,58,59,43},
   {26,29,30,8,12,13,56,59,60,44},
   {27,31,32,10,15,16,57,61,62,46},
   {28,32,33,11,16,17,58,62,63,47},
   {29,33,34,12,17,18,59,63,64,48},
   {30,34,35,13,18,19,60,64,65,49},
   {36,37,38,51,52,53,2,4,5,23},
   {37,39,40,52,54,55,4,7,8,25},
   {38,40,41,53,55,56,5,8,9,26},
   {39,42,43,54,57,58,7,11,12,28},
   {40,43,44,55,58,59,8,12,13,29},
   {41,44,45,56,59,60,9,13,14,30},
   {42,46,47,57,61,62,11,16,17,32},
   {43,47,48,58,62,63,12,17,18,33},
   {44,48,49,59,63,64,13,18,19,34},
   {45,49,50,60,64,65,14,19,20,35},
   {51,52,53,37,39,40,23,25,26,4},
   {52,54,55,39,42,43,25,28,29,7},
   {53,55,56,40,43,44,26,29,30,8},
   {54,57,58,42,46,47,28,32,33,11},
   {55,58,59,43,47,48,29,33,34,12},
   {56,59,60,44,48,49,30,34,35,13},
};
// indices of [r3] x [r8] in [R11]
static unsigned short iCartP38[45][10] = {
   {0,1,2,21,22,23,42,43,44,63},
   {1,3,4,22,24,25,43,45,46,64},
   {2,4,5,23,25,26,44,46,47,65},
   {3,6,7,24,27,28,45,48,49,66},
   {4,7,8,25,28,29,46,49,50,67},
   {5,8,9,26,29,30,47,50,51,68},
   {6,10,11,27,31,32,48,52,53,69},
   {7,11,12,28,32,33,49,53,54,70},
   {8,12,13,29,33,34,50,54,55,71},
   {9,13,14,30,34,35,51,55,56,72},
   {10,15,16,31,36,37,52,57,58,73},
   {11,16,17,32,37,38,53,58,59,74},
   {12,17,18,33,38,39,54,59,60,75},
   {13,18,19,34,39,40,55,60,61,76},
   {14,19,20,35,40,41,56,61,62,77},
   {21,22,23,1,3,4,63,64,65,43},
   {22,24,25,3,6,7,64,66,67,45},
   {23,25,26,4,7,8,65,67,68,46},
   {24,27,28,6,10,11,66,69,70,48},
   {25,28,29,7,11,12,67,70,71,49},
   {26,29,30,8,12,13,68,71,72,50},
   {27,31,32,10,15,16,69,73,74,52},
   {28,32,33,11,16,17,70,74,75,53},
   {29,33,34,12,17,18,71,75,76,54},
   {30,34,35,13,18,19,72,76,77,55},
   {42,43,44,63,64,65,2,4,5,23},
   {43,45,46,64,66,67,4,7,8,25},
   {44,46,47,65,67,68,5,8,9,26},
   {45,48,49,66,69,70,7,11,12,28},
   {46,49,50,67,70,71,8,12,13,29},
   {47,50,51,68,71,72,9,13,14,30},
   {48,52,53,69,73,74,11,16,17,32},
   {49,53,54,70,74,75,12,17,18,33},
   {50,54,55,71,75,76,13,18,19,34},
   {51,55,56,72,76,77,14,19,20,35},
   {63,64,65,43,45,46,23,25,26,4},
   {64,66,67,45,48,49,25,28,29,7},
   {65,67,68,46,49,50,26,29,30,8},
   {66,69,70,48,52,53,28,32,33,11},
   {67,70,71,49,53,54,29,33,34,12},
   {68,71,72,50,54,55,30,34,35,13},
   {69,73,74,52,57,58,32,37,38,16},
   {70,74,75,53,58,59,33,38,39,17},
   {71,75,76,54,59,60,34,39,40,18},
   {72,76,77,55,60,61,35,40,41,19},
};
// indices of [r3] x [r9] in [R12]
static unsigned short iCartP39[55][10] = {
   {0,1,2,28,29,30,49,50,51,70},
   {1,3,4,29,31,32,50,52,53,71},
   {2,4,5,30,32,33,51,53,54,72},
   {3,6,7,31,34,35,52,55,56,73},
   {4,7,8,32,35,36,53,56,57,74},
   {5,8,9,33,36,37,54,57,58,75},
   {6,10,11,34,38,39,55,59,60,76},
   {7,11,12,35,39,40,56,60,61,77},
   {8,12,13,36,40,41,57,61,62,78},
   {9,13,14,37,41,42,58,62,63,79},
   {10,15,16,38,43,44,59,64,65,80},
   {11,16,17,39,44,45,60,65,66,81},
   {12,17,18,40,45,46,61,66,67,82},
   {13,18,19,41,46,47,62,67,68,83},
   {14,19,20,42,47,48,63,68,69,84},
   {28,29,30,1,3,4,70,71,72,50},
   {29,31,32,3,6,7,71,73,74,52},
   {30,32,33,4,7,8,72,74,75,53},
   {31,34,35,6,10,11,73,76,77,55},
   {32,35,36,7,11,12,74,77,78,56},
   {33,36,37,8,12,13,75,78,79,57},
   {34,38,39,10,15,16,76,80,81,59},
   {35,39,40,11,16,17,77,81,82,60},
   {36,40,41,12,17,18,78,82,83,61},
   {37,41,42,13,18,19,79,83,84,62},
   {38,43,44,15,21,22,80,85,86,64},
   {39,44,45,16,22,23,81,86,87,65},
   {40,45,46,17,23,24,82,87,88,66},
   {41,46,47,18,24,25,83,88,89,67},
   {42,47,48,19,25,26,84,89,90,68},
   {49,50,51,70,71,72,2,4,5,30},
   {50,52,53,71,73,74,4,7,8,32},
   {51,53,54,72,74,75,5,8,9,33},
   {52,55,56,73,76,77,7,11,12,35},
   {53,56,57,74,77,78,8,12,13,36},
   {54,57,58,75,78,79,9,13,14,37},
   {55,59,60,76,80,81,11,16,17,39},
   {56,60,61,77,81,82,12,17,18,40},
   {57,61,62,78,82,83,13,18,19,41},
   {58,62,63,79,83,84,14,19,20,42},
   {59,64,65,80,85,86,16,22,23,44},
   {60,65,66,81,86,87,17,23,24,45},
   {61,66,67,82,87,88,18,24,25,46},
   {62,67,68,83,88,89,19,25,26,47},
   {63,68,69,84,89,90,20,26,27,48},
   {70,71,72,50,52,53,30,32,33,4},
   {71,73,74,52,55,56,32,35,36,7},
   {72,74,75,53,56,57,33,36,37,8},
   {73,76,77,55,59,60,35,39,40,11},
   {74,77,78,56,60,61,36,40,41,12},
   {75,78,79,57,61,62,37,41,42,13},
   {76,80,81,59,64,65,39,44,45,16},
   {77,81,82,60,65,66,40,45,46,17},
   {78,82,83,61,66,67,41,46,47,18},
   {79,83,84,62,67,68,42,47,48,19},
};
// indices of [r3] x [r10] in [R13]
static unsigned short iCartP3a[66][10] = {
   {0,1,2,28,29,30,56,57,58,84},
   {1,3,4,29,31,32,57,59,60,85},
   {2,4,5,30,32,33,58,60,61,86},
   {3,6,7,31,34,35,59,62,63,87},
   {4,7,8,32,35,36,60,63,64,88},
   {5,8,9,33,36,37,61,64,65,89},
   {6,10,11,34,38,39,62,66,67,90},
   {7,11,12,35,39,40,63,67,68,91},
   {8,12,13,36,40,41,64,68,69,92},
   {9,13,14,37,41,42,65,69,70,93},
   {10,15,16,38,43,44,66,71,72,94},
   {11,16,17,39,44,45,67,72,73,95},
   {12,17,18,40,45,46,68,73,74,96},
   {13,18,19,41,46,47,69,74,75,97},
   {14,19,20,42,47,48,70,75,76,98},
   {15,21,22,43,49,50,71,77,78,99},
   {16,22,23,44,50,51,72,78,79,100},
   {17,23,24,45,51,52,73,79,80,101},
   {18,24,25,46,52,53,74,80,81,102},
   {19,25,26,47,53,54,75,81,82,103},
   {20,26,27,48,54,55,76,82,83,104},
   {28,29,30,1,3,4,84,85,86,57},
   {29,31,32,3,6,7,85,87,88,59},
   {30,32,33,4,7,8,86,88,89,60},
   {31,34,35,6,10,11,87,90,91,62},
   {32,35,36,7,11,12,88,91,92,63},
   {33,36,37,8,12,13,89,92,93,64},
   {34,38,39,10,15,16,90,94,95,66},
   {35,39,40,11,16,17,91,95,96,67},
   {36,40,41,12,17,18,92,96,97,68},
   {37,41,42,13,18,19,93,97,98,69},
   {38,43,44,15,21,22,94,99,100,71},
   {39,44,45,16,22,23,95,100,101,72},
   {40,45,46,17,23,24,96,101,102,73},
   {41,46,47,18,24,25,97,102,103,74},
   {42,47,48,19,25,26,98,103,104,75},
   {56,57,58,84,85,86,2,4,5,30},
   {57,59,60,85,87,88,4,7,8,32},
   {58,60,61,86,88,89,5,8,9,33},
   {59,62,63,87,90,91,7,11,12,35},
   {60,63,64,88,91,92,8,12,13,36},
   {61,64,65,89,92,93,9,13,14,37},
   {62,66,67,90,94,95,11,16,17,39},
   {63,67,68,91,95,96,12,17,18,40},
   {64,68,69,92,96,97,13,18,19,41},
   {65,69,70,93,97,98,14,19,20,42},
   {66,71,72,94,99,100,16,22,23,44},
   {67,72,73,95,100,101,17,23,24,45},
   {68,73,74,96,101,102,18,24,25,46},
   {69,74,75,97,102,103,19,25,26,47},
   {70,75,76,98,103,104,20,26,27,48},
   {84,85,86,57,59,60,30,32,33,4},
   {85,87,88,59,62,63,32,35,36,7},
   {86,88,89,60,63,64,33,36,37,8},
   {87,90,91,62,66,67,35,39,40,11},
   {88,91,92,63,67,68,36,40,41,12},
   {89,92,93,64,68,69,37,41,42,13},
   {90,94,95,66,71,72,39,44,45,16},
   {91,95,96,67,72,73,40,45,46,17},
   {92,96,97,68,73,74,41,46,47,18},
   {93,97,98,69,74,75,42,47,48,19},
   {94,99,100,71,77,78,44,50,51,22},
   {95,100,101,72,78,79,45,51,52,23},
   {96,101,102,73,79,80,46,52,53,24},
   {97,102,103,74,80,81,47,53,54,25},
   {98,103,104,75,81,82,48,54,55,26},
};
// indices of [r3] x [r11] in [R14]
static unsigned short iCartP3b[78][10] = {
   {0,1,2,36,37,38,64,65,66,92},
   {1,3,4,37,39,40,65,67,68,93},
   {2,4,5,38,40,41,66,68,69,94},
   {3,6,7,39,42,43,67,70,71,95},
   {4,7,8,40,43,44,68,71,72,96},
   {5,8,9,41,44,45,69,72,73,97},
   {6,10,11,42,46,47,70,74,75,98},
   {7,11,12,43,47,48,71,75,76,99},
   {8,12,13,44,48,49,72,76,77,100},
   {9,13,14,45,49,50,73,77,78,101},
   {10,15,16,46,51,52,74,79,80,102},
   {11,16,17,47,52,53,75,80,81,103},
   {12,17,18,48,53,54,76,81,82,104},
   {13,18,19,49,54,55,77,82,83,105},
   {14,19,20,50,55,56,78,83,84,106},
   {15,21,22,51,57,58,79,85,86,107},
   {16,22,23,52,58,59,80,86,87,108},
   {17,23,24,53,59,60,81,87,88,109},
   {18,24,25,54,60,61,82,88,89,110},
   {19,25,26,55,61,62,83,89,90,111},
   {20,26,27,56,62,63,84,90,91,112},
   {36,37,38,1,3,4,92,93,94,65},
   {37,39,40,3,6,7,93,95,96,67},
   {38,40,41,4,7,8,94,96,97,68},
   {39,42,43,6,10,11,95,98,99,70},
   {40,43,44,7,11,12,96,99,100,71},
   {41,44,45,8,12,13,97,100,101,72},
   {42,46,47,10,15,16,98,102,103,74},
   {43,47,48,11,16,17,99,103,104,75},
   {44,48,49,12,17,18,100,104,105,76},
   {45,49,50,13,18,19,101,105,106,77},
   {46,51,52,15,21,22,102,107,108,79},
   {47,52,53,16,22,23,103,108,109,80},
   {48,53,54,17,23,24,104,109,110,81},
   {49,54,55,18,24,25,105,110,111,82},
   {50,55,56,19,25,26,106,111,112,83},
   {51,57,58,21,28,29,107,113,114,85},
   {52,58,59,22,29,30,108,114,115,86},
   {53,59,60,23,30,31,109,115,116,87},
   {54,60,61,24,31,32,110,116,117,88},
   {55,61,62,25,32,33,111,117,118,89},
   {56,62,63,26,33,34,112,118,119,90},
   {64,65,66,92,93,94,2,4,5,38},
   {65,67,68,93,95,96,4,7,8,40},
   {66,68,69,94,96,97,5,8,9,41},
   {67,70,71,95,98,99,7,11,12,43},
   {68,71,72,96,99,100,8,12,13,44},
   {69,72,73,97,100,101,9,13,14,45},
   {70,74,75,98,102,103,11,16,17,47},
   {71,75,76,99,103,104,12,17,18,48},
   {72,76,77,100,104,105,13,18,19,49},
   {73,77,78,101,105,106,14,19,20,50},
   {74,79,80,102,107,108,16,22,23,52},
   {75,80,81,103,108,109,17,23,24,53},
   {76,81,82,104,109,110,18,24,25,54},
   {77,82,83,105,110,111,19,25,26,55},
   {78,83,84,106,111,112,20,26,27,56},
   {79,85,86,107,113,114,22,29,30,58},
   {80,86,87,108,114,115,23,30,31,59},
   {81,87,88,109,115,116,24,31,32,60},
   {82,88,89,110,116,117,25,32,33,61},
   {83,89,90,111,117,118,26,33,34,62},
   {84,90,91,112,118,119,27,34,35,63},
   {92,93,94,65,67,68,38,40,41,4},
   {93,95,96,67,70,71,40,43,44,7},
   {94,96,97,68,71,72,41,44,45,8},
   {95,98,99,70,74,75,43,47,48,11},
   {96,99,100,71,75,76,44,48,49,12},
   {97,100,101,72,76,77,45,49,50,13},
   {98,102,103,74,79,80,47,52,53,16},
   {99,103,104,75,80,81,48,53,54,17},
   {100,104,105,76,81,82,49,54,55,18},
   {101,105,106,77,82,83,50,55,56,19},
   {102,107,108,79,85,86,52,58,59,22},
   {103,108,109,80,86,87,53,59,60,23},
   {104,109,110,81,87,88,54,60,61,24},
   {105,110,111,82,88,89,55,61,62,25},
   {106,111,112,83,89,90,56,62,63,26}
};
// indices of [r4] x [r0] in [R4]
static unsigned short iCartP40[1][15] = {
   {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14},
};
// indices of [r4] x [r1] in [R5]
static unsigned short iCartP41[3][15] = {
   {0,1,2,3,4,5,6,7,8,12,13,14,18,19,20},
   {6,7,8,9,10,11,1,3,4,18,19,20,13,15,16},
   {12,13,14,15,16,17,18,19,20,2,4,5,8,10,11},
};
// indices of [r4] x [r2] in [R6]
static unsigned short iCartP42[6][15] = {
   {0,1,2,3,4,5,10,11,12,16,17,18,22,23,24},
   {1,3,4,6,7,8,11,13,14,17,19,20,23,25,26},
   {2,4,5,7,8,9,12,14,15,18,20,21,24,26,27},
   {10,11,12,13,14,15,1,3,4,22,23,24,17,19,20},
   {16,17,18,19,20,21,22,23,24,2,4,5,12,14,15},
   {22,23,24,25,26,27,17,19,20,12,14,15,4,7,8},
};
// indices of [r4] x [r3] in [R7]
static unsigned short iCartP43[10][15] = {
   {0,1,2,3,4,5,10,11,12,20,21,22,30,31,32},
   {1,3,4,6,7,8,11,13,14,21,23,24,31,33,34},
   {2,4,5,7,8,9,12,14,15,22,24,25,32,34,35},
   {10,11,12,13,14,15,1,3,4,30,31,32,21,23,24},
   {11,13,14,16,17,18,3,6,7,31,33,34,23,26,27},
   {12,14,15,17,18,19,4,7,8,32,34,35,24,27,28},
   {20,21,22,23,24,25,30,31,32,2,4,5,12,14,15},
   {21,23,24,26,27,28,31,33,34,4,7,8,14,17,18},
   {22,24,25,27,28,29,32,34,35,5,8,9,15,18,19},
   {30,31,32,33,34,35,21,23,24,12,14,15,4,7,8},
};
// indices of [r4] x [r4] in [R8]
static unsigned short iCartP44[15][15] = {
   {0,1,2,3,4,5,15,16,17,25,26,27,35,36,37},
   {1,3,4,6,7,8,16,18,19,26,28,29,36,38,39},
   {2,4,5,7,8,9,17,19,20,27,29,30,37,39,40},
   {3,6,7,10,11,12,18,21,22,28,31,32,38,41,42},
   {4,7,8,11,12,13,19,22,23,29,32,33,39,42,43},
   {5,8,9,12,13,14,20,23,24,30,33,34,40,43,44},
   {15,16,17,18,19,20,1,3,4,35,36,37,26,28,29},
   {16,18,19,21,22,23,3,6,7,36,38,39,28,31,32},
   {17,19,20,22,23,24,4,7,8,37,39,40,29,32,33},
   {25,26,27,28,29,30,35,36,37,2,4,5,17,19,20},
   {26,28,29,31,32,33,36,38,39,4,7,8,19,22,23},
   {27,29,30,32,33,34,37,39,40,5,8,9,20,23,24},
   {35,36,37,38,39,40,26,28,29,17,19,20,4,7,8},
   {36,38,39,41,42,43,28,31,32,19,22,23,7,11,12},
   {37,39,40,42,43,44,29,32,33,20,23,24,8,12,13},
};
// indices of [r4] x [r5] in [R9]
static unsigned short iCartP45[21][15] = {
   {0,1,2,3,4,5,15,16,17,30,31,32,45,46,47},
   {1,3,4,6,7,8,16,18,19,31,33,34,46,48,49},
   {2,4,5,7,8,9,17,19,20,32,34,35,47,49,50},
   {3,6,7,10,11,12,18,21,22,33,36,37,48,51,52},
   {4,7,8,11,12,13,19,22,23,34,37,38,49,52,53},
   {5,8,9,12,13,14,20,23,24,35,38,39,50,53,54},
   {15,16,17,18,19,20,1,3,4,45,46,47,31,33,34},
   {16,18,19,21,22,23,3,6,7,46,48,49,33,36,37},
   {17,19,20,22,23,24,4,7,8,47,49,50,34,37,38},
   {18,21,22,25,26,27,6,10,11,48,51,52,36,40,41},
   {19,22,23,26,27,28,7,11,12,49,52,53,37,41,42},
   {20,23,24,27,28,29,8,12,13,50,53,54,38,42,43},
   {30,31,32,33,34,35,45,46,47,2,4,5,17,19,20},
   {31,33,34,36,37,38,46,48,49,4,7,8,19,22,23},
   {32,34,35,37,38,39,47,49,50,5,8,9,20,23,24},
   {33,36,37,40,41,42,48,51,52,7,11,12,22,26,27},
   {34,37,38,41,42,43,49,52,53,8,12,13,23,27,28},
   {35,38,39,42,43,44,50,53,54,9,13,14,24,28,29},
   {45,46,47,48,49,50,31,33,34,17,19,20,4,7,8},
   {46,48,49,51,52,53,33,36,37,19,22,23,7,11,12},
   {47,49,50,52,53,54,34,37,38,20,23,24,8,12,13},
};
// indices of [r4] x [r6] in [R10]
static unsigned short iCartP46[28][15] = {
   {0,1,2,3,4,5,21,22,23,36,37,38,51,52,53},
   {1,3,4,6,7,8,22,24,25,37,39,40,52,54,55},
   {2,4,5,7,8,9,23,25,26,38,40,41,53,55,56},
   {3,6,7,10,11,12,24,27,28,39,42,43,54,57,58},
   {4,7,8,11,12,13,25,28,29,40,43,44,55,58,59},
   {5,8,9,12,13,14,26,29,30,41,44,45,56,59,60},
   {6,10,11,15,16,17,27,31,32,42,46,47,57,61,62},
   {7,11,12,16,17,18,28,32,33,43,47,48,58,62,63},
   {8,12,13,17,18,19,29,33,34,44,48,49,59,63,64},
   {9,13,14,18,19,20,30,34,35,45,49,50,60,64,65},
   {21,22,23,24,25,26,1,3,4,51,52,53,37,39,40},
   {22,24,25,27,28,29,3,6,7,52,54,55,39,42,43},
   {23,25,26,28,29,30,4,7,8,53,55,56,40,43,44},
   {24,27,28,31,32,33,6,10,11,54,57,58,42,46,47},
   {25,28,29,32,33,34,7,11,12,55,58,59,43,47,48},
   {26,29,30,33,34,35,8,12,13,56,59,60,44,48,49},
   {36,37,38,39,40,41,51,52,53,2,4,5,23,25,26},
   {37,39,40,42,43,44,52,54,55,4,7,8,25,28,29},
   {38,40,41,43,44,45,53,55,56,5,8,9,26,29,30},
   {39,42,43,46,47,48,54,57,58,7,11,12,28,32,33},
   {40,43,44,47,48,49,55,58,59,8,12,13,29,33,34},
   {41,44,45,48,49,50,56,59,60,9,13,14,30,34,35},
   {51,52,53,54,55,56,37,39,40,23,25,26,4,7,8},
   {52,54,55,57,58,59,39,42,43,25,28,29,7,11,12},
   {53,55,56,58,59,60,40,43,44,26,29,30,8,12,13},
   {54,57,58,61,62,63,42,46,47,28,32,33,11,16,17},
   {55,58,59,62,63,64,43,47,48,29,33,34,12,17,18},
   {56,59,60,63,64,65,44,48,49,30,34,35,13,18,19},
};
// indices of [r4] x [r7] in [R11]
static unsigned short iCartP47[36][15] = {
   {0,1,2,3,4,5,21,22,23,42,43,44,63,64,65},
   {1,3,4,6,7,8,22,24,25,43,45,46,64,66,67},
   {2,4,5,7,8,9,23,25,26,44,46,47,65,67,68},
   {3,6,7,10,11,12,24,27,28,45,48,49,66,69,70},
   {4,7,8,11,12,13,25,28,29,46,49,50,67,70,71},
   {5,8,9,12,13,14,26,29,30,47,50,51,68,71,72},
   {6,10,11,15,16,17,27,31,32,48,52,53,69,73,74},
   {7,11,12,16,17,18,28,32,33,49,53,54,70,74,75},
   {8,12,13,17,18,19,29,33,34,50,54,55,71,75,76},
   {9,13,14,18,19,20,30,34,35,51,55,56,72,76,77},
   {21,22,23,24,25,26,1,3,4,63,64,65,43,45,46},
   {22,24,25,27,28,29,3,6,7,64,66,67,45,48,49},
   {23,25,26,28,29,30,4,7,8,65,67,68,46,49,50},
   {24,27,28,31,32,33,6,10,11,66,69,70,48,52,53},
   {25,28,29,32,33,34,7,11,12,67,70,71,49,53,54},
   {26,29,30,33,34,35,8,12,13,68,71,72,50,54,55},
   {27,31,32,36,37,38,10,15,16,69,73,74,52,57,58},
   {28,32,33,37,38,39,11,16,17,70,74,75,53,58,59},
   {29,33,34,38,39,40,12,17,18,71,75,76,54,59,60},
   {30,34,35,39,40,41,13,18,19,72,76,77,55,60,61},
   {42,43,44,45,46,47,63,64,65,2,4,5,23,25,26},
   {43,45,46,48,49,50,64,66,67,4,7,8,25,28,29},
   {44,46,47,49,50,51,65,67,68,5,8,9,26,29,30},
   {45,48,49,52,53,54,66,69,70,7,11,12,28,32,33},
   {46,49,50,53,54,55,67,70,71,8,12,13,29,33,34},
   {47,50,51,54,55,56,68,71,72,9,13,14,30,34,35},
   {48,52,53,57,58,59,69,73,74,11,16,17,32,37,38},
   {49,53,54,58,59,60,70,74,75,12,17,18,33,38,39},
   {50,54,55,59,60,61,71,75,76,13,18,19,34,39,40},
   {51,55,56,60,61,62,72,76,77,14,19,20,35,40,41},
   {63,64,65,66,67,68,43,45,46,23,25,26,4,7,8},
   {64,66,67,69,70,71,45,48,49,25,28,29,7,11,12},
   {65,67,68,70,71,72,46,49,50,26,29,30,8,12,13},
   {66,69,70,73,74,75,48,52,53,28,32,33,11,16,17},
   {67,70,71,74,75,76,49,53,54,29,33,34,12,17,18},
   {68,71,72,75,76,77,50,54,55,30,34,35,13,18,19},
};
// indices of [r4] x [r8] in [R12]
static unsigned short iCartP48[45][15] = {
   {0,1,2,3,4,5,28,29,30,49,50,51,70,71,72},
   {1,3,4,6,7,8,29,31,32,50,52,53,71,73,74},
   {2,4,5,7,8,9,30,32,33,51,53,54,72,74,75},
   {3,6,7,10,11,12,31,34,35,52,55,56,73,76,77},
   {4,7,8,11,12,13,32,35,36,53,56,57,74,77,78},
   {5,8,9,12,13,14,33,36,37,54,57,58,75,78,79},
   {6,10,11,15,16,17,34,38,39,55,59,60,76,80,81},
   {7,11,12,16,17,18,35,39,40,56,60,61,77,81,82},
   {8,12,13,17,18,19,36,40,41,57,61,62,78,82,83},
   {9,13,14,18,19,20,37,41,42,58,62,63,79,83,84},
   {10,15,16,21,22,23,38,43,44,59,64,65,80,85,86},
   {11,16,17,22,23,24,39,44,45,60,65,66,81,86,87},
   {12,17,18,23,24,25,40,45,46,61,66,67,82,87,88},
   {13,18,19,24,25,26,41,46,47,62,67,68,83,88,89},
   {14,19,20,25,26,27,42,47,48,63,68,69,84,89,90},
   {28,29,30,31,32,33,1,3,4,70,71,72,50,52,53},
   {29,31,32,34,35,36,3,6,7,71,73,74,52,55,56},
   {30,32,33,35,36,37,4,7,8,72,74,75,53,56,57},
   {31,34,35,38,39,40,6,10,11,73,76,77,55,59,60},
   {32,35,36,39,40,41,7,11,12,74,77,78,56,60,61},
   {33,36,37,40,41,42,8,12,13,75,78,79,57,61,62},
   {34,38,39,43,44,45,10,15,16,76,80,81,59,64,65},
   {35,39,40,44,45,46,11,16,17,77,81,82,60,65,66},
   {36,40,41,45,46,47,12,17,18,78,82,83,61,66,67},
   {37,41,42,46,47,48,13,18,19,79,83,84,62,67,68},
   {49,50,51,52,53,54,70,71,72,2,4,5,30,32,33},
   {50,52,53,55,56,57,71,73,74,4,7,8,32,35,36},
   {51,53,54,56,57,58,72,74,75,5,8,9,33,36,37},
   {52,55,56,59,60,61,73,76,77,7,11,12,35,39,40},
   {53,56,57,60,61,62,74,77,78,8,12,13,36,40,41},
   {54,57,58,61,62,63,75,78,79,9,13,14,37,41,42},
   {55,59,60,64,65,66,76,80,81,11,16,17,39,44,45},
   {56,60,61,65,66,67,77,81,82,12,17,18,40,45,46},
   {57,61,62,66,67,68,78,82,83,13,18,19,41,46,47},
   {58,62,63,67,68,69,79,83,84,14,19,20,42,47,48},
   {70,71,72,73,74,75,50,52,53,30,32,33,4,7,8},
   {71,73,74,76,77,78,52,55,56,32,35,36,7,11,12},
   {72,74,75,77,78,79,53,56,57,33,36,37,8,12,13},
   {73,76,77,80,81,82,55,59,60,35,39,40,11,16,17},
   {74,77,78,81,82,83,56,60,61,36,40,41,12,17,18},
   {75,78,79,82,83,84,57,61,62,37,41,42,13,18,19},
   {76,80,81,85,86,87,59,64,65,39,44,45,16,22,23},
   {77,81,82,86,87,88,60,65,66,40,45,46,17,23,24},
   {78,82,83,87,88,89,61,66,67,41,46,47,18,24,25},
   {79,83,84,88,89,90,62,67,68,42,47,48,19,25,26},
};
// indices of [r4] x [r9] in [R13]
static unsigned short iCartP49[55][15] = {
   {0,1,2,3,4,5,28,29,30,56,57,58,84,85,86},
   {1,3,4,6,7,8,29,31,32,57,59,60,85,87,88},
   {2,4,5,7,8,9,30,32,33,58,60,61,86,88,89},
   {3,6,7,10,11,12,31,34,35,59,62,63,87,90,91},
   {4,7,8,11,12,13,32,35,36,60,63,64,88,91,92},
   {5,8,9,12,13,14,33,36,37,61,64,65,89,92,93},
   {6,10,11,15,16,17,34,38,39,62,66,67,90,94,95},
   {7,11,12,16,17,18,35,39,40,63,67,68,91,95,96},
   {8,12,13,17,18,19,36,40,41,64,68,69,92,96,97},
   {9,13,14,18,19,20,37,41,42,65,69,70,93,97,98},
   {10,15,16,21,22,23,38,43,44,66,71,72,94,99,100},
   {11,16,17,22,23,24,39,44,45,67,72,73,95,100,101},
   {12,17,18,23,24,25,40,45,46,68,73,74,96,101,102},
   {13,18,19,24,25,26,41,46,47,69,74,75,97,102,103},
   {14,19,20,25,26,27,42,47,48,70,75,76,98,103,104},
   {28,29,30,31,32,33,1,3,4,84,85,86,57,59,60},
   {29,31,32,34,35,36,3,6,7,85,87,88,59,62,63},
   {30,32,33,35,36,37,4,7,8,86,88,89,60,63,64},
   {31,34,35,38,39,40,6,10,11,87,90,91,62,66,67},
   {32,35,36,39,40,41,7,11,12,88,91,92,63,67,68},
   {33,36,37,40,41,42,8,12,13,89,92,93,64,68,69},
   {34,38,39,43,44,45,10,15,16,90,94,95,66,71,72},
   {35,39,40,44,45,46,11,16,17,91,95,96,67,72,73},
   {36,40,41,45,46,47,12,17,18,92,96,97,68,73,74},
   {37,41,42,46,47,48,13,18,19,93,97,98,69,74,75},
   {38,43,44,49,50,51,15,21,22,94,99,100,71,77,78},
   {39,44,45,50,51,52,16,22,23,95,100,101,72,78,79},
   {40,45,46,51,52,53,17,23,24,96,101,102,73,79,80},
   {41,46,47,52,53,54,18,24,25,97,102,103,74,80,81},
   {42,47,48,53,54,55,19,25,26,98,103,104,75,81,82},
   {56,57,58,59,60,61,84,85,86,2,4,5,30,32,33},
   {57,59,60,62,63,64,85,87,88,4,7,8,32,35,36},
   {58,60,61,63,64,65,86,88,89,5,8,9,33,36,37},
   {59,62,63,66,67,68,87,90,91,7,11,12,35,39,40},
   {60,63,64,67,68,69,88,91,92,8,12,13,36,40,41},
   {61,64,65,68,69,70,89,92,93,9,13,14,37,41,42},
   {62,66,67,71,72,73,90,94,95,11,16,17,39,44,45},
   {63,67,68,72,73,74,91,95,96,12,17,18,40,45,46},
   {64,68,69,73,74,75,92,96,97,13,18,19,41,46,47},
   {65,69,70,74,75,76,93,97,98,14,19,20,42,47,48},
   {66,71,72,77,78,79,94,99,100,16,22,23,44,50,51},
   {67,72,73,78,79,80,95,100,101,17,23,24,45,51,52},
   {68,73,74,79,80,81,96,101,102,18,24,25,46,52,53},
   {69,74,75,80,81,82,97,102,103,19,25,26,47,53,54},
   {70,75,76,81,82,83,98,103,104,20,26,27,48,54,55},
   {84,85,86,87,88,89,57,59,60,30,32,33,4,7,8},
   {85,87,88,90,91,92,59,62,63,32,35,36,7,11,12},
   {86,88,89,91,92,93,60,63,64,33,36,37,8,12,13},
   {87,90,91,94,95,96,62,66,67,35,39,40,11,16,17},
   {88,91,92,95,96,97,63,67,68,36,40,41,12,17,18},
   {89,92,93,96,97,98,64,68,69,37,41,42,13,18,19},
   {90,94,95,99,100,101,66,71,72,39,44,45,16,22,23},
   {91,95,96,100,101,102,67,72,73,40,45,46,17,23,24},
   {92,96,97,101,102,103,68,73,74,41,46,47,18,24,25},
   {93,97,98,102,103,104,69,74,75,42,47,48,19,25,26},
};
// indices of [r4] x [r10] in [R14]
static unsigned short iCartP4a[66][15] = {
   {0,1,2,3,4,5,36,37,38,64,65,66,92,93,94},
   {1,3,4,6,7,8,37,39,40,65,67,68,93,95,96},
   {2,4,5,7,8,9,38,40,41,66,68,69,94,96,97},
   {3,6,7,10,11,12,39,42,43,67,70,71,95,98,99},
   {4,7,8,11,12,13,40,43,44,68,71,72,96,99,100},
   {5,8,9,12,13,14,41,44,45,69,72,73,97,100,101},
   {6,10,11,15,16,17,42,46,47,70,74,75,98,102,103},
   {7,11,12,16,17,18,43,47,48,71,75,76,99,103,104},
   {8,12,13,17,18,19,44,48,49,72,76,77,100,104,105},
   {9,13,14,18,19,20,45,49,50,73,77,78,101,105,106},
   {10,15,16,21,22,23,46,51,52,74,79,80,102,107,108},
   {11,16,17,22,23,24,47,52,53,75,80,81,103,108,109},
   {12,17,18,23,24,25,48,53,54,76,81,82,104,109,110},
   {13,18,19,24,25,26,49,54,55,77,82,83,105,110,111},
   {14,19,20,25,26,27,50,55,56,78,83,84,106,111,112},
   {15,21,22,28,29,30,51,57,58,79,85,86,107,113,114},
   {16,22,23,29,30,31,52,58,59,80,86,87,108,114,115},
   {17,23,24,30,31,32,53,59,60,81,87,88,109,115,116},
   {18,24,25,31,32,33,54,60,61,82,88,89,110,116,117},
   {19,25,26,32,33,34,55,61,62,83,89,90,111,117,118},
   {20,26,27,33,34,35,56,62,63,84,90,91,112,118,119},
   {36,37,38,39,40,41,1,3,4,92,93,94,65,67,68},
   {37,39,40,42,43,44,3,6,7,93,95,96,67,70,71},
   {38,40,41,43,44,45,4,7,8,94,96,97,68,71,72},
   {39,42,43,46,47,48,6,10,11,95,98,99,70,74,75},
   {40,43,44,47,48,49,7,11,12,96,99,100,71,75,76},
   {41,44,45,48,49,50,8,12,13,97,100,101,72,76,77},
   {42,46,47,51,52,53,10,15,16,98,102,103,74,79,80},
   {43,47,48,52,53,54,11,16,17,99,103,104,75,80,81},
   {44,48,49,53,54,55,12,17,18,100,104,105,76,81,82},
   {45,49,50,54,55,56,13,18,19,101,105,106,77,82,83},
   {46,51,52,57,58,59,15,21,22,102,107,108,79,85,86},
   {47,52,53,58,59,60,16,22,23,103,108,109,80,86,87},
   {48,53,54,59,60,61,17,23,24,104,109,110,81,87,88},
   {49,54,55,60,61,62,18,24,25,105,110,111,82,88,89},
   {50,55,56,61,62,63,19,25,26,106,111,112,83,89,90},
   {64,65,66,67,68,69,92,93,94,2,4,5,38,40,41},
   {65,67,68,70,71,72,93,95,96,4,7,8,40,43,44},
   {66,68,69,71,72,73,94,96,97,5,8,9,41,44,45},
   {67,70,71,74,75,76,95,98,99,7,11,12,43,47,48},
   {68,71,72,75,76,77,96,99,100,8,12,13,44,48,49},
   {69,72,73,76,77,78,97,100,101,9,13,14,45,49,50},
   {70,74,75,79,80,81,98,102,103,11,16,17,47,52,53},
   {71,75,76,80,81,82,99,103,104,12,17,18,48,53,54},
   {72,76,77,81,82,83,100,104,105,13,18,19,49,54,55},
   {73,77,78,82,83,84,101,105,106,14,19,20,50,55,56},
   {74,79,80,85,86,87,102,107,108,16,22,23,52,58,59},
   {75,80,81,86,87,88,103,108,109,17,23,24,53,59,60},
   {76,81,82,87,88,89,104,109,110,18,24,25,54,60,61},
   {77,82,83,88,89,90,105,110,111,19,25,26,55,61,62},
   {78,83,84,89,90,91,106,111,112,20,26,27,56,62,63},
   {92,93,94,95,96,97,65,67,68,38,40,41,4,7,8},
   {93,95,96,98,99,100,67,70,71,40,43,44,7,11,12},
   {94,96,97,99,100,101,68,71,72,41,44,45,8,12,13},
   {95,98,99,102,103,104,70,74,75,43,47,48,11,16,17},
   {96,99,100,103,104,105,71,75,76,44,48,49,12,17,18},
   {97,100,101,104,105,106,72,76,77,45,49,50,13,18,19},
   {98,102,103,107,108,109,74,79,80,47,52,53,16,22,23},
   {99,103,104,108,109,110,75,80,81,48,53,54,17,23,24},
   {100,104,105,109,110,111,76,81,82,49,54,55,18,24,25},
   {101,105,106,110,111,112,77,82,83,50,55,56,19,25,26},
   {102,107,108,113,114,115,79,85,86,52,58,59,22,29,30},
   {103,108,109,114,115,116,80,86,87,53,59,60,23,30,31},
   {104,109,110,115,116,117,81,87,88,54,60,61,24,31,32},
   {105,110,111,116,117,118,82,88,89,55,61,62,25,32,33},
   {106,111,112,117,118,119,83,89,90,56,62,63,26,33,34}
};
// indices of [r5] x [r0] in [R5]
static unsigned short iCartP50[1][21] = {
   {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20},
};
// indices of [r5] x [r1] in [R6]
static unsigned short iCartP51[3][21] = {
   {0,1,2,3,4,5,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24},
   {10,11,12,13,14,15,1,3,4,6,7,8,22,23,24,25,26,27,17,19,20},
   {16,17,18,19,20,21,22,23,24,25,26,27,2,4,5,7,8,9,12,14,15},
};
// indices of [r5] x [r2] in [R7]
static unsigned short iCartP52[6][21] = {
   {0,1,2,3,4,5,10,11,12,13,14,15,20,21,22,23,24,25,30,31,32},
   {1,3,4,6,7,8,11,13,14,16,17,18,21,23,24,26,27,28,31,33,34},
   {2,4,5,7,8,9,12,14,15,17,18,19,22,24,25,27,28,29,32,34,35},
   {10,11,12,13,14,15,1,3,4,6,7,8,30,31,32,33,34,35,21,23,24},
   {20,21,22,23,24,25,30,31,32,33,34,35,2,4,5,7,8,9,12,14,15},
   {30,31,32,33,34,35,21,23,24,26,27,28,12,14,15,17,18,19,4,7,8},
};
// indices of [r5] x [r3] in [R8]
static unsigned short iCartP53[10][21] = {
   {0,1,2,3,4,5,15,16,17,18,19,20,25,26,27,28,29,30,35,36,37},
   {1,3,4,6,7,8,16,18,19,21,22,23,26,28,29,31,32,33,36,38,39},
   {2,4,5,7,8,9,17,19,20,22,23,24,27,29,30,32,33,34,37,39,40},
   {15,16,17,18,19,20,1,3,4,6,7,8,35,36,37,38,39,40,26,28,29},
   {16,18,19,21,22,23,3,6,7,10,11,12,36,38,39,41,42,43,28,31,32},
   {17,19,20,22,23,24,4,7,8,11,12,13,37,39,40,42,43,44,29,32,33},
   {25,26,27,28,29,30,35,36,37,38,39,40,2,4,5,7,8,9,17,19,20},
   {26,28,29,31,32,33,36,38,39,41,42,43,4,7,8,11,12,13,19,22,23},
   {27,29,30,32,33,34,37,39,40,42,43,44,5,8,9,12,13,14,20,23,24},
   {35,36,37,38,39,40,26,28,29,31,32,33,17,19,20,22,23,24,4,7,8},
};
// indices of [r5] x [r4] in [R9]
static unsigned short iCartP54[15][21] = {
   {0,1,2,3,4,5,15,16,17,18,19,20,30,31,32,33,34,35,45,46,47},
   {1,3,4,6,7,8,16,18,19,21,22,23,31,33,34,36,37,38,46,48,49},
   {2,4,5,7,8,9,17,19,20,22,23,24,32,34,35,37,38,39,47,49,50},
   {3,6,7,10,11,12,18,21,22,25,26,27,33,36,37,40,41,42,48,51,52},
   {4,7,8,11,12,13,19,22,23,26,27,28,34,37,38,41,42,43,49,52,53},
   {5,8,9,12,13,14,20,23,24,27,28,29,35,38,39,42,43,44,50,53,54},
   {15,16,17,18,19,20,1,3,4,6,7,8,45,46,47,48,49,50,31,33,34},
   {16,18,19,21,22,23,3,6,7,10,11,12,46,48,49,51,52,53,33,36,37},
   {17,19,20,22,23,24,4,7,8,11,12,13,47,49,50,52,53,54,34,37,38},
   {30,31,32,33,34,35,45,46,47,48,49,50,2,4,5,7,8,9,17,19,20},
   {31,33,34,36,37,38,46,48,49,51,52,53,4,7,8,11,12,13,19,22,23},
   {32,34,35,37,38,39,47,49,50,52,53,54,5,8,9,12,13,14,20,23,24},
   {45,46,47,48,49,50,31,33,34,36,37,38,17,19,20,22,23,24,4,7,8},
   {46,48,49,51,52,53,33,36,37,40,41,42,19,22,23,26,27,28,7,11,12},
   {47,49,50,52,53,54,34,37,38,41,42,43,20,23,24,27,28,29,8,12,13},
};
// indices of [r5] x [r5] in [R10]
static unsigned short iCartP55[21][21] = {
   {0,1,2,3,4,5,21,22,23,24,25,26,36,37,38,39,40,41,51,52,53},
   {1,3,4,6,7,8,22,24,25,27,28,29,37,39,40,42,43,44,52,54,55},
   {2,4,5,7,8,9,23,25,26,28,29,30,38,40,41,43,44,45,53,55,56},
   {3,6,7,10,11,12,24,27,28,31,32,33,39,42,43,46,47,48,54,57,58},
   {4,7,8,11,12,13,25,28,29,32,33,34,40,43,44,47,48,49,55,58,59},
   {5,8,9,12,13,14,26,29,30,33,34,35,41,44,45,48,49,50,56,59,60},
   {21,22,23,24,25,26,1,3,4,6,7,8,51,52,53,54,55,56,37,39,40},
   {22,24,25,27,28,29,3,6,7,10,11,12,52,54,55,57,58,59,39,42,43},
   {23,25,26,28,29,30,4,7,8,11,12,13,53,55,56,58,59,60,40,43,44},
   {24,27,28,31,32,33,6,10,11,15,16,17,54,57,58,61,62,63,42,46,47},
   {25,28,29,32,33,34,7,11,12,16,17,18,55,58,59,62,63,64,43,47,48},
   {26,29,30,33,34,35,8,12,13,17,18,19,56,59,60,63,64,65,44,48,49},
   {36,37,38,39,40,41,51,52,53,54,55,56,2,4,5,7,8,9,23,25,26},
   {37,39,40,42,43,44,52,54,55,57,58,59,4,7,8,11,12,13,25,28,29},
   {38,40,41,43,44,45,53,55,56,58,59,60,5,8,9,12,13,14,26,29,30},
   {39,42,43,46,47,48,54,57,58,61,62,63,7,11,12,16,17,18,28,32,33},
   {40,43,44,47,48,49,55,58,59,62,63,64,8,12,13,17,18,19,29,33,34},
   {41,44,45,48,49,50,56,59,60,63,64,65,9,13,14,18,19,20,30,34,35},
   {51,52,53,54,55,56,37,39,40,42,43,44,23,25,26,28,29,30,4,7,8},
   {52,54,55,57,58,59,39,42,43,46,47,48,25,28,29,32,33,34,7,11,12},
   {53,55,56,58,59,60,40,43,44,47,48,49,26,29,30,33,34,35,8,12,13},
};
// indices of [r5] x [r6] in [R11]
static unsigned short iCartP56[28][21] = {
   {0,1,2,3,4,5,21,22,23,24,25,26,42,43,44,45,46,47,63,64,65},
   {1,3,4,6,7,8,22,24,25,27,28,29,43,45,46,48,49,50,64,66,67},
   {2,4,5,7,8,9,23,25,26,28,29,30,44,46,47,49,50,51,65,67,68},
   {3,6,7,10,11,12,24,27,28,31,32,33,45,48,49,52,53,54,66,69,70},
   {4,7,8,11,12,13,25,28,29,32,33,34,46,49,50,53,54,55,67,70,71},
   {5,8,9,12,13,14,26,29,30,33,34,35,47,50,51,54,55,56,68,71,72},
   {6,10,11,15,16,17,27,31,32,36,37,38,48,52,53,57,58,59,69,73,74},
   {7,11,12,16,17,18,28,32,33,37,38,39,49,53,54,58,59,60,70,74,75},
   {8,12,13,17,18,19,29,33,34,38,39,40,50,54,55,59,60,61,71,75,76},
   {9,13,14,18,19,20,30,34,35,39,40,41,51,55,56,60,61,62,72,76,77},
   {21,22,23,24,25,26,1,3,4,6,7,8,63,64,65,66,67,68,43,45,46},
   {22,24,25,27,28,29,3,6,7,10,11,12,64,66,67,69,70,71,45,48,49},
   {23,25,26,28,29,30,4,7,8,11,12,13,65,67,68,70,71,72,46,49,50},
   {24,27,28,31,32,33,6,10,11,15,16,17,66,69,70,73,74,75,48,52,53},
   {25,28,29,32,33,34,7,11,12,16,17,18,67,70,71,74,75,76,49,53,54},
   {26,29,30,33,34,35,8,12,13,17,18,19,68,71,72,75,76,77,50,54,55},
   {42,43,44,45,46,47,63,64,65,66,67,68,2,4,5,7,8,9,23,25,26},
   {43,45,46,48,49,50,64,66,67,69,70,71,4,7,8,11,12,13,25,28,29},
   {44,46,47,49,50,51,65,67,68,70,71,72,5,8,9,12,13,14,26,29,30},
   {45,48,49,52,53,54,66,69,70,73,74,75,7,11,12,16,17,18,28,32,33},
   {46,49,50,53,54,55,67,70,71,74,75,76,8,12,13,17,18,19,29,33,34},
   {47,50,51,54,55,56,68,71,72,75,76,77,9,13,14,18,19,20,30,34,35},
   {63,64,65,66,67,68,43,45,46,48,49,50,23,25,26,28,29,30,4,7,8},
   {64,66,67,69,70,71,45,48,49,52,53,54,25,28,29,32,33,34,7,11,12},
   {65,67,68,70,71,72,46,49,50,53,54,55,26,29,30,33,34,35,8,12,13},
   {66,69,70,73,74,75,48,52,53,57,58,59,28,32,33,37,38,39,11,16,17},
   {67,70,71,74,75,76,49,53,54,58,59,60,29,33,34,38,39,40,12,17,18},
   {68,71,72,75,76,77,50,54,55,59,60,61,30,34,35,39,40,41,13,18,19},
};
// indices of [r5] x [r7] in [R12]
static unsigned short iCartP57[36][21] = {
   {0,1,2,3,4,5,28,29,30,31,32,33,49,50,51,52,53,54,70,71,72},
   {1,3,4,6,7,8,29,31,32,34,35,36,50,52,53,55,56,57,71,73,74},
   {2,4,5,7,8,9,30,32,33,35,36,37,51,53,54,56,57,58,72,74,75},
   {3,6,7,10,11,12,31,34,35,38,39,40,52,55,56,59,60,61,73,76,77},
   {4,7,8,11,12,13,32,35,36,39,40,41,53,56,57,60,61,62,74,77,78},
   {5,8,9,12,13,14,33,36,37,40,41,42,54,57,58,61,62,63,75,78,79},
   {6,10,11,15,16,17,34,38,39,43,44,45,55,59,60,64,65,66,76,80,81},
   {7,11,12,16,17,18,35,39,40,44,45,46,56,60,61,65,66,67,77,81,82},
   {8,12,13,17,18,19,36,40,41,45,46,47,57,61,62,66,67,68,78,82,83},
   {9,13,14,18,19,20,37,41,42,46,47,48,58,62,63,67,68,69,79,83,84},
   {28,29,30,31,32,33,1,3,4,6,7,8,70,71,72,73,74,75,50,52,53},
   {29,31,32,34,35,36,3,6,7,10,11,12,71,73,74,76,77,78,52,55,56},
   {30,32,33,35,36,37,4,7,8,11,12,13,72,74,75,77,78,79,53,56,57},
   {31,34,35,38,39,40,6,10,11,15,16,17,73,76,77,80,81,82,55,59,60},
   {32,35,36,39,40,41,7,11,12,16,17,18,74,77,78,81,82,83,56,60,61},
   {33,36,37,40,41,42,8,12,13,17,18,19,75,78,79,82,83,84,57,61,62},
   {34,38,39,43,44,45,10,15,16,21,22,23,76,80,81,85,86,87,59,64,65},
   {35,39,40,44,45,46,11,16,17,22,23,24,77,81,82,86,87,88,60,65,66},
   {36,40,41,45,46,47,12,17,18,23,24,25,78,82,83,87,88,89,61,66,67},
   {37,41,42,46,47,48,13,18,19,24,25,26,79,83,84,88,89,90,62,67,68},
   {49,50,51,52,53,54,70,71,72,73,74,75,2,4,5,7,8,9,30,32,33},
   {50,52,53,55,56,57,71,73,74,76,77,78,4,7,8,11,12,13,32,35,36},
   {51,53,54,56,57,58,72,74,75,77,78,79,5,8,9,12,13,14,33,36,37},
   {52,55,56,59,60,61,73,76,77,80,81,82,7,11,12,16,17,18,35,39,40},
   {53,56,57,60,61,62,74,77,78,81,82,83,8,12,13,17,18,19,36,40,41},
   {54,57,58,61,62,63,75,78,79,82,83,84,9,13,14,18,19,20,37,41,42},
   {55,59,60,64,65,66,76,80,81,85,86,87,11,16,17,22,23,24,39,44,45},
   {56,60,61,65,66,67,77,81,82,86,87,88,12,17,18,23,24,25,40,45,46},
   {57,61,62,66,67,68,78,82,83,87,88,89,13,18,19,24,25,26,41,46,47},
   {58,62,63,67,68,69,79,83,84,88,89,90,14,19,20,25,26,27,42,47,48},
   {70,71,72,73,74,75,50,52,53,55,56,57,30,32,33,35,36,37,4,7,8},
   {71,73,74,76,77,78,52,55,56,59,60,61,32,35,36,39,40,41,7,11,12},
   {72,74,75,77,78,79,53,56,57,60,61,62,33,36,37,40,41,42,8,12,13},
   {73,76,77,80,81,82,55,59,60,64,65,66,35,39,40,44,45,46,11,16,17},
   {74,77,78,81,82,83,56,60,61,65,66,67,36,40,41,45,46,47,12,17,18},
   {75,78,79,82,83,84,57,61,62,66,67,68,37,41,42,46,47,48,13,18,19},
};
// indices of [r5] x [r8] in [R13]
static unsigned short iCartP58[45][21] = {
   {0,1,2,3,4,5,28,29,30,31,32,33,56,57,58,59,60,61,84,85,86},
   {1,3,4,6,7,8,29,31,32,34,35,36,57,59,60,62,63,64,85,87,88},
   {2,4,5,7,8,9,30,32,33,35,36,37,58,60,61,63,64,65,86,88,89},
   {3,6,7,10,11,12,31,34,35,38,39,40,59,62,63,66,67,68,87,90,91},
   {4,7,8,11,12,13,32,35,36,39,40,41,60,63,64,67,68,69,88,91,92},
   {5,8,9,12,13,14,33,36,37,40,41,42,61,64,65,68,69,70,89,92,93},
   {6,10,11,15,16,17,34,38,39,43,44,45,62,66,67,71,72,73,90,94,95},
   {7,11,12,16,17,18,35,39,40,44,45,46,63,67,68,72,73,74,91,95,96},
   {8,12,13,17,18,19,36,40,41,45,46,47,64,68,69,73,74,75,92,96,97},
   {9,13,14,18,19,20,37,41,42,46,47,48,65,69,70,74,75,76,93,97,98},
   {10,15,16,21,22,23,38,43,44,49,50,51,66,71,72,77,78,79,94,99,100},
   {11,16,17,22,23,24,39,44,45,50,51,52,67,72,73,78,79,80,95,100,101},
   {12,17,18,23,24,25,40,45,46,51,52,53,68,73,74,79,80,81,96,101,102},
   {13,18,19,24,25,26,41,46,47,52,53,54,69,74,75,80,81,82,97,102,103},
   {14,19,20,25,26,27,42,47,48,53,54,55,70,75,76,81,82,83,98,103,104},
   {28,29,30,31,32,33,1,3,4,6,7,8,84,85,86,87,88,89,57,59,60},
   {29,31,32,34,35,36,3,6,7,10,11,12,85,87,88,90,91,92,59,62,63},
   {30,32,33,35,36,37,4,7,8,11,12,13,86,88,89,91,92,93,60,63,64},
   {31,34,35,38,39,40,6,10,11,15,16,17,87,90,91,94,95,96,62,66,67},
   {32,35,36,39,40,41,7,11,12,16,17,18,88,91,92,95,96,97,63,67,68},
   {33,36,37,40,41,42,8,12,13,17,18,19,89,92,93,96,97,98,64,68,69},
   {34,38,39,43,44,45,10,15,16,21,22,23,90,94,95,99,100,101,66,71,72},
   {35,39,40,44,45,46,11,16,17,22,23,24,91,95,96,100,101,102,67,72,73},
   {36,40,41,45,46,47,12,17,18,23,24,25,92,96,97,101,102,103,68,73,74},
   {37,41,42,46,47,48,13,18,19,24,25,26,93,97,98,102,103,104,69,74,75},
   {56,57,58,59,60,61,84,85,86,87,88,89,2,4,5,7,8,9,30,32,33},
   {57,59,60,62,63,64,85,87,88,90,91,92,4,7,8,11,12,13,32,35,36},
   {58,60,61,63,64,65,86,88,89,91,92,93,5,8,9,12,13,14,33,36,37},
   {59,62,63,66,67,68,87,90,91,94,95,96,7,11,12,16,17,18,35,39,40},
   {60,63,64,67,68,69,88,91,92,95,96,97,8,12,13,17,18,19,36,40,41},
   {61,64,65,68,69,70,89,92,93,96,97,98,9,13,14,18,19,20,37,41,42},
   {62,66,67,71,72,73,90,94,95,99,100,101,11,16,17,22,23,24,39,44,45},
   {63,67,68,72,73,74,91,95,96,100,101,102,12,17,18,23,24,25,40,45,46},
   {64,68,69,73,74,75,92,96,97,101,102,103,13,18,19,24,25,26,41,46,47},
   {65,69,70,74,75,76,93,97,98,102,103,104,14,19,20,25,26,27,42,47,48},
   {84,85,86,87,88,89,57,59,60,62,63,64,30,32,33,35,36,37,4,7,8},
   {85,87,88,90,91,92,59,62,63,66,67,68,32,35,36,39,40,41,7,11,12},
   {86,88,89,91,92,93,60,63,64,67,68,69,33,36,37,40,41,42,8,12,13},
   {87,90,91,94,95,96,62,66,67,71,72,73,35,39,40,44,45,46,11,16,17},
   {88,91,92,95,96,97,63,67,68,72,73,74,36,40,41,45,46,47,12,17,18},
   {89,92,93,96,97,98,64,68,69,73,74,75,37,41,42,46,47,48,13,18,19},
   {90,94,95,99,100,101,66,71,72,77,78,79,39,44,45,50,51,52,16,22,23},
   {91,95,96,100,101,102,67,72,73,78,79,80,40,45,46,51,52,53,17,23,24},
   {92,96,97,101,102,103,68,73,74,79,80,81,41,46,47,52,53,54,18,24,25},
   {93,97,98,102,103,104,69,74,75,80,81,82,42,47,48,53,54,55,19,25,26},
};
// indices of [r5] x [r9] in [R14]
static unsigned short iCartP59[55][21] = {
   {0,1,2,3,4,5,36,37,38,39,40,41,64,65,66,67,68,69,92,93,94},
   {1,3,4,6,7,8,37,39,40,42,43,44,65,67,68,70,71,72,93,95,96},
   {2,4,5,7,8,9,38,40,41,43,44,45,66,68,69,71,72,73,94,96,97},
   {3,6,7,10,11,12,39,42,43,46,47,48,67,70,71,74,75,76,95,98,99},
   {4,7,8,11,12,13,40,43,44,47,48,49,68,71,72,75,76,77,96,99,100},
   {5,8,9,12,13,14,41,44,45,48,49,50,69,72,73,76,77,78,97,100,101},
   {6,10,11,15,16,17,42,46,47,51,52,53,70,74,75,79,80,81,98,102,103},
   {7,11,12,16,17,18,43,47,48,52,53,54,71,75,76,80,81,82,99,103,104},
   {8,12,13,17,18,19,44,48,49,53,54,55,72,76,77,81,82,83,100,104,105},
   {9,13,14,18,19,20,45,49,50,54,55,56,73,77,78,82,83,84,101,105,106},
   {10,15,16,21,22,23,46,51,52,57,58,59,74,79,80,85,86,87,102,107,108},
   {11,16,17,22,23,24,47,52,53,58,59,60,75,80,81,86,87,88,103,108,109},
   {12,17,18,23,24,25,48,53,54,59,60,61,76,81,82,87,88,89,104,109,110},
   {13,18,19,24,25,26,49,54,55,60,61,62,77,82,83,88,89,90,105,110,111},
   {14,19,20,25,26,27,50,55,56,61,62,63,78,83,84,89,90,91,106,111,112},
   {36,37,38,39,40,41,1,3,4,6,7,8,92,93,94,95,96,97,65,67,68},
   {37,39,40,42,43,44,3,6,7,10,11,12,93,95,96,98,99,100,67,70,71},
   {38,40,41,43,44,45,4,7,8,11,12,13,94,96,97,99,100,101,68,71,72},
   {39,42,43,46,47,48,6,10,11,15,16,17,95,98,99,102,103,104,70,74,75},
   {40,43,44,47,48,49,7,11,12,16,17,18,96,99,100,103,104,105,71,75,76},
   {41,44,45,48,49,50,8,12,13,17,18,19,97,100,101,104,105,106,72,76,77},
   {42,46,47,51,52,53,10,15,16,21,22,23,98,102,103,107,108,109,74,79,80},
   {43,47,48,52,53,54,11,16,17,22,23,24,99,103,104,108,109,110,75,80,81},
   {44,48,49,53,54,55,12,17,18,23,24,25,100,104,105,109,110,111,76,81,82},
   {45,49,50,54,55,56,13,18,19,24,25,26,101,105,106,110,111,112,77,82,83},
   {46,51,52,57,58,59,15,21,22,28,29,30,102,107,108,113,114,115,79,85,86},
   {47,52,53,58,59,60,16,22,23,29,30,31,103,108,109,114,115,116,80,86,87},
   {48,53,54,59,60,61,17,23,24,30,31,32,104,109,110,115,116,117,81,87,88},
   {49,54,55,60,61,62,18,24,25,31,32,33,105,110,111,116,117,118,82,88,89},
   {50,55,56,61,62,63,19,25,26,32,33,34,106,111,112,117,118,119,83,89,90},
   {64,65,66,67,68,69,92,93,94,95,96,97,2,4,5,7,8,9,38,40,41},
   {65,67,68,70,71,72,93,95,96,98,99,100,4,7,8,11,12,13,40,43,44},
   {66,68,69,71,72,73,94,96,97,99,100,101,5,8,9,12,13,14,41,44,45},
   {67,70,71,74,75,76,95,98,99,102,103,104,7,11,12,16,17,18,43,47,48},
   {68,71,72,75,76,77,96,99,100,103,104,105,8,12,13,17,18,19,44,48,49},
   {69,72,73,76,77,78,97,100,101,104,105,106,9,13,14,18,19,20,45,49,50},
   {70,74,75,79,80,81,98,102,103,107,108,109,11,16,17,22,23,24,47,52,53},
   {71,75,76,80,81,82,99,103,104,108,109,110,12,17,18,23,24,25,48,53,54},
   {72,76,77,81,82,83,100,104,105,109,110,111,13,18,19,24,25,26,49,54,55},
   {73,77,78,82,83,84,101,105,106,110,111,112,14,19,20,25,26,27,50,55,56},
   {74,79,80,85,86,87,102,107,108,113,114,115,16,22,23,29,30,31,52,58,59},
   {75,80,81,86,87,88,103,108,109,114,115,116,17,23,24,30,31,32,53,59,60},
   {76,81,82,87,88,89,104,109,110,115,116,117,18,24,25,31,32,33,54,60,61},
   {77,82,83,88,89,90,105,110,111,116,117,118,19,25,26,32,33,34,55,61,62},
   {78,83,84,89,90,91,106,111,112,117,118,119,20,26,27,33,34,35,56,62,63},
   {92,93,94,95,96,97,65,67,68,70,71,72,38,40,41,43,44,45,4,7,8},
   {93,95,96,98,99,100,67,70,71,74,75,76,40,43,44,47,48,49,7,11,12},
   {94,96,97,99,100,101,68,71,72,75,76,77,41,44,45,48,49,50,8,12,13},
   {95,98,99,102,103,104,70,74,75,79,80,81,43,47,48,52,53,54,11,16,17},
   {96,99,100,103,104,105,71,75,76,80,81,82,44,48,49,53,54,55,12,17,18},
   {97,100,101,104,105,106,72,76,77,81,82,83,45,49,50,54,55,56,13,18,19},
   {98,102,103,107,108,109,74,79,80,85,86,87,47,52,53,58,59,60,16,22,23},
   {99,103,104,108,109,110,75,80,81,86,87,88,48,53,54,59,60,61,17,23,24},
   {100,104,105,109,110,111,76,81,82,87,88,89,49,54,55,60,61,62,18,24,25},
   {101,105,106,110,111,112,77,82,83,88,89,90,50,55,56,61,62,63,19,25,26}
};
// indices of [r6] x [r0] in [R6]
static unsigned short iCartP60[1][28] = {
   {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27},
};
// indices of [r6] x [r1] in [R7]
static unsigned short iCartP61[3][28] = {
   {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,21,22,23,24,25,30,31,32,33,34,35},
   {10,11,12,13,14,15,16,17,18,19,1,3,4,6,7,8,30,31,32,33,34,35,21,23,24,26,27,28},
   {20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,2,4,5,7,8,9,12,14,15,17,18,19},
};
// indices of [r6] x [r2] in [R8]
static unsigned short iCartP62[6][28] = {
   {0,1,2,3,4,5,6,7,8,9,15,16,17,18,19,20,25,26,27,28,29,30,35,36,37,38,39,40},
   {1,3,4,6,7,8,10,11,12,13,16,18,19,21,22,23,26,28,29,31,32,33,36,38,39,41,42,43},
   {2,4,5,7,8,9,11,12,13,14,17,19,20,22,23,24,27,29,30,32,33,34,37,39,40,42,43,44},
   {15,16,17,18,19,20,21,22,23,24,1,3,4,6,7,8,35,36,37,38,39,40,26,28,29,31,32,33},
   {25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,2,4,5,7,8,9,17,19,20,22,23,24},
   {35,36,37,38,39,40,41,42,43,44,26,28,29,31,32,33,17,19,20,22,23,24,4,7,8,11,12,13},
};
// indices of [r6] x [r3] in [R9]
static unsigned short iCartP63[10][28] = {
   {0,1,2,3,4,5,6,7,8,9,15,16,17,18,19,20,30,31,32,33,34,35,45,46,47,48,49,50},
   {1,3,4,6,7,8,10,11,12,13,16,18,19,21,22,23,31,33,34,36,37,38,46,48,49,51,52,53},
   {2,4,5,7,8,9,11,12,13,14,17,19,20,22,23,24,32,34,35,37,38,39,47,49,50,52,53,54},
   {15,16,17,18,19,20,21,22,23,24,1,3,4,6,7,8,45,46,47,48,49,50,31,33,34,36,37,38},
   {16,18,19,21,22,23,25,26,27,28,3,6,7,10,11,12,46,48,49,51,52,53,33,36,37,40,41,42},
   {17,19,20,22,23,24,26,27,28,29,4,7,8,11,12,13,47,49,50,52,53,54,34,37,38,41,42,43},
   {30,31,32,33,34,35,36,37,38,39,45,46,47,48,49,50,2,4,5,7,8,9,17,19,20,22,23,24},
   {31,33,34,36,37,38,40,41,42,43,46,48,49,51,52,53,4,7,8,11,12,13,19,22,23,26,27,28},
   {32,34,35,37,38,39,41,42,43,44,47,49,50,52,53,54,5,8,9,12,13,14,20,23,24,27,28,29},
   {45,46,47,48,49,50,51,52,53,54,31,33,34,36,37,38,17,19,20,22,23,24,4,7,8,11,12,13},
};
// indices of [r6] x [r4] in [R10]
static unsigned short iCartP64[15][28] = {
   {0,1,2,3,4,5,6,7,8,9,21,22,23,24,25,26,36,37,38,39,40,41,51,52,53,54,55,56},
   {1,3,4,6,7,8,10,11,12,13,22,24,25,27,28,29,37,39,40,42,43,44,52,54,55,57,58,59},
   {2,4,5,7,8,9,11,12,13,14,23,25,26,28,29,30,38,40,41,43,44,45,53,55,56,58,59,60},
   {3,6,7,10,11,12,15,16,17,18,24,27,28,31,32,33,39,42,43,46,47,48,54,57,58,61,62,63},
   {4,7,8,11,12,13,16,17,18,19,25,28,29,32,33,34,40,43,44,47,48,49,55,58,59,62,63,64},
   {5,8,9,12,13,14,17,18,19,20,26,29,30,33,34,35,41,44,45,48,49,50,56,59,60,63,64,65},
   {21,22,23,24,25,26,27,28,29,30,1,3,4,6,7,8,51,52,53,54,55,56,37,39,40,42,43,44},
   {22,24,25,27,28,29,31,32,33,34,3,6,7,10,11,12,52,54,55,57,58,59,39,42,43,46,47,48},
   {23,25,26,28,29,30,32,33,34,35,4,7,8,11,12,13,53,55,56,58,59,60,40,43,44,47,48,49},
   {36,37,38,39,40,41,42,43,44,45,51,52,53,54,55,56,2,4,5,7,8,9,23,25,26,28,29,30},
   {37,39,40,42,43,44,46,47,48,49,52,54,55,57,58,59,4,7,8,11,12,13,25,28,29,32,33,34},
   {38,40,41,43,44,45,47,48,49,50,53,55,56,58,59,60,5,8,9,12,13,14,26,29,30,33,34,35},
   {51,52,53,54,55,56,57,58,59,60,37,39,40,42,43,44,23,25,26,28,29,30,4,7,8,11,12,13},
   {52,54,55,57,58,59,61,62,63,64,39,42,43,46,47,48,25,28,29,32,33,34,7,11,12,16,17,18},
   {53,55,56,58,59,60,62,63,64,65,40,43,44,47,48,49,26,29,30,33,34,35,8,12,13,17,18,19},
};
// indices of [r6] x [r5] in [R11]
static unsigned short iCartP65[21][28] = {
   {0,1,2,3,4,5,6,7,8,9,21,22,23,24,25,26,42,43,44,45,46,47,63,64,65,66,67,68},
   {1,3,4,6,7,8,10,11,12,13,22,24,25,27,28,29,43,45,46,48,49,50,64,66,67,69,70,71},
   {2,4,5,7,8,9,11,12,13,14,23,25,26,28,29,30,44,46,47,49,50,51,65,67,68,70,71,72},
   {3,6,7,10,11,12,15,16,17,18,24,27,28,31,32,33,45,48,49,52,53,54,66,69,70,73,74,75},
   {4,7,8,11,12,13,16,17,18,19,25,28,29,32,33,34,46,49,50,53,54,55,67,70,71,74,75,76},
   {5,8,9,12,13,14,17,18,19,20,26,29,30,33,34,35,47,50,51,54,55,56,68,71,72,75,76,77},
   {21,22,23,24,25,26,27,28,29,30,1,3,4,6,7,8,63,64,65,66,67,68,43,45,46,48,49,50},
   {22,24,25,27,28,29,31,32,33,34,3,6,7,10,11,12,64,66,67,69,70,71,45,48,49,52,53,54},
   {23,25,26,28,29,30,32,33,34,35,4,7,8,11,12,13,65,67,68,70,71,72,46,49,50,53,54,55},
   {24,27,28,31,32,33,36,37,38,39,6,10,11,15,16,17,66,69,70,73,74,75,48,52,53,57,58,59},
   {25,28,29,32,33,34,37,38,39,40,7,11,12,16,17,18,67,70,71,74,75,76,49,53,54,58,59,60},
   {26,29,30,33,34,35,38,39,40,41,8,12,13,17,18,19,68,71,72,75,76,77,50,54,55,59,60,61},
   {42,43,44,45,46,47,48,49,50,51,63,64,65,66,67,68,2,4,5,7,8,9,23,25,26,28,29,30},
   {43,45,46,48,49,50,52,53,54,55,64,66,67,69,70,71,4,7,8,11,12,13,25,28,29,32,33,34},
   {44,46,47,49,50,51,53,54,55,56,65,67,68,70,71,72,5,8,9,12,13,14,26,29,30,33,34,35},
   {45,48,49,52,53,54,57,58,59,60,66,69,70,73,74,75,7,11,12,16,17,18,28,32,33,37,38,39},
   {46,49,50,53,54,55,58,59,60,61,67,70,71,74,75,76,8,12,13,17,18,19,29,33,34,38,39,40},
   {47,50,51,54,55,56,59,60,61,62,68,71,72,75,76,77,9,13,14,18,19,20,30,34,35,39,40,41},
   {63,64,65,66,67,68,69,70,71,72,43,45,46,48,49,50,23,25,26,28,29,30,4,7,8,11,12,13},
   {64,66,67,69,70,71,73,74,75,76,45,48,49,52,53,54,25,28,29,32,33,34,7,11,12,16,17,18},
   {65,67,68,70,71,72,74,75,76,77,46,49,50,53,54,55,26,29,30,33,34,35,8,12,13,17,18,19},
};
// indices of [r6] x [r6] in [R12]
static unsigned short iCartP66[28][28] = {
   {0,1,2,3,4,5,6,7,8,9,28,29,30,31,32,33,49,50,51,52,53,54,70,71,72,73,74,75},
   {1,3,4,6,7,8,10,11,12,13,29,31,32,34,35,36,50,52,53,55,56,57,71,73,74,76,77,78},
   {2,4,5,7,8,9,11,12,13,14,30,32,33,35,36,37,51,53,54,56,57,58,72,74,75,77,78,79},
   {3,6,7,10,11,12,15,16,17,18,31,34,35,38,39,40,52,55,56,59,60,61,73,76,77,80,81,82},
   {4,7,8,11,12,13,16,17,18,19,32,35,36,39,40,41,53,56,57,60,61,62,74,77,78,81,82,83},
   {5,8,9,12,13,14,17,18,19,20,33,36,37,40,41,42,54,57,58,61,62,63,75,78,79,82,83,84},
   {6,10,11,15,16,17,21,22,23,24,34,38,39,43,44,45,55,59,60,64,65,66,76,80,81,85,86,87},
   {7,11,12,16,17,18,22,23,24,25,35,39,40,44,45,46,56,60,61,65,66,67,77,81,82,86,87,88},
   {8,12,13,17,18,19,23,24,25,26,36,40,41,45,46,47,57,61,62,66,67,68,78,82,83,87,88,89},
   {9,13,14,18,19,20,24,25,26,27,37,41,42,46,47,48,58,62,63,67,68,69,79,83,84,88,89,90},
   {28,29,30,31,32,33,34,35,36,37,1,3,4,6,7,8,70,71,72,73,74,75,50,52,53,55,56,57},
   {29,31,32,34,35,36,38,39,40,41,3,6,7,10,11,12,71,73,74,76,77,78,52,55,56,59,60,61},
   {30,32,33,35,36,37,39,40,41,42,4,7,8,11,12,13,72,74,75,77,78,79,53,56,57,60,61,62},
   {31,34,35,38,39,40,43,44,45,46,6,10,11,15,16,17,73,76,77,80,81,82,55,59,60,64,65,66},
   {32,35,36,39,40,41,44,45,46,47,7,11,12,16,17,18,74,77,78,81,82,83,56,60,61,65,66,67},
   {33,36,37,40,41,42,45,46,47,48,8,12,13,17,18,19,75,78,79,82,83,84,57,61,62,66,67,68},
   {49,50,51,52,53,54,55,56,57,58,70,71,72,73,74,75,2,4,5,7,8,9,30,32,33,35,36,37},
   {50,52,53,55,56,57,59,60,61,62,71,73,74,76,77,78,4,7,8,11,12,13,32,35,36,39,40,41},
   {51,53,54,56,57,58,60,61,62,63,72,74,75,77,78,79,5,8,9,12,13,14,33,36,37,40,41,42},
   {52,55,56,59,60,61,64,65,66,67,73,76,77,80,81,82,7,11,12,16,17,18,35,39,40,44,45,46},
   {53,56,57,60,61,62,65,66,67,68,74,77,78,81,82,83,8,12,13,17,18,19,36,40,41,45,46,47},
   {54,57,58,61,62,63,66,67,68,69,75,78,79,82,83,84,9,13,14,18,19,20,37,41,42,46,47,48},
   {70,71,72,73,74,75,76,77,78,79,50,52,53,55,56,57,30,32,33,35,36,37,4,7,8,11,12,13},
   {71,73,74,76,77,78,80,81,82,83,52,55,56,59,60,61,32,35,36,39,40,41,7,11,12,16,17,18},
   {72,74,75,77,78,79,81,82,83,84,53,56,57,60,61,62,33,36,37,40,41,42,8,12,13,17,18,19},
   {73,76,77,80,81,82,85,86,87,88,55,59,60,64,65,66,35,39,40,44,45,46,11,16,17,22,23,24},
   {74,77,78,81,82,83,86,87,88,89,56,60,61,65,66,67,36,40,41,45,46,47,12,17,18,23,24,25},
   {75,78,79,82,83,84,87,88,89,90,57,61,62,66,67,68,37,41,42,46,47,48,13,18,19,24,25,26},
};
// indices of [r6] x [r7] in [R13]
static unsigned short iCartP67[36][28] = {
   {0,1,2,3,4,5,6,7,8,9,28,29,30,31,32,33,56,57,58,59,60,61,84,85,86,87,88,89},
   {1,3,4,6,7,8,10,11,12,13,29,31,32,34,35,36,57,59,60,62,63,64,85,87,88,90,91,92},
   {2,4,5,7,8,9,11,12,13,14,30,32,33,35,36,37,58,60,61,63,64,65,86,88,89,91,92,93},
   {3,6,7,10,11,12,15,16,17,18,31,34,35,38,39,40,59,62,63,66,67,68,87,90,91,94,95,96},
   {4,7,8,11,12,13,16,17,18,19,32,35,36,39,40,41,60,63,64,67,68,69,88,91,92,95,96,97},
   {5,8,9,12,13,14,17,18,19,20,33,36,37,40,41,42,61,64,65,68,69,70,89,92,93,96,97,98},
   {6,10,11,15,16,17,21,22,23,24,34,38,39,43,44,45,62,66,67,71,72,73,90,94,95,99,100,101},
   {7,11,12,16,17,18,22,23,24,25,35,39,40,44,45,46,63,67,68,72,73,74,91,95,96,100,101,102},
   {8,12,13,17,18,19,23,24,25,26,36,40,41,45,46,47,64,68,69,73,74,75,92,96,97,101,102,103},
   {9,13,14,18,19,20,24,25,26,27,37,41,42,46,47,48,65,69,70,74,75,76,93,97,98,102,103,104},
   {28,29,30,31,32,33,34,35,36,37,1,3,4,6,7,8,84,85,86,87,88,89,57,59,60,62,63,64},
   {29,31,32,34,35,36,38,39,40,41,3,6,7,10,11,12,85,87,88,90,91,92,59,62,63,66,67,68},
   {30,32,33,35,36,37,39,40,41,42,4,7,8,11,12,13,86,88,89,91,92,93,60,63,64,67,68,69},
   {31,34,35,38,39,40,43,44,45,46,6,10,11,15,16,17,87,90,91,94,95,96,62,66,67,71,72,73},
   {32,35,36,39,40,41,44,45,46,47,7,11,12,16,17,18,88,91,92,95,96,97,63,67,68,72,73,74},
   {33,36,37,40,41,42,45,46,47,48,8,12,13,17,18,19,89,92,93,96,97,98,64,68,69,73,74,75},
   {34,38,39,43,44,45,49,50,51,52,10,15,16,21,22,23,90,94,95,99,100,101,66,71,72,77,78,79},
   {35,39,40,44,45,46,50,51,52,53,11,16,17,22,23,24,91,95,96,100,101,102,67,72,73,78,79,80},
   {36,40,41,45,46,47,51,52,53,54,12,17,18,23,24,25,92,96,97,101,102,103,68,73,74,79,80,81},
   {37,41,42,46,47,48,52,53,54,55,13,18,19,24,25,26,93,97,98,102,103,104,69,74,75,80,81,82},
   {56,57,58,59,60,61,62,63,64,65,84,85,86,87,88,89,2,4,5,7,8,9,30,32,33,35,36,37},
   {57,59,60,62,63,64,66,67,68,69,85,87,88,90,91,92,4,7,8,11,12,13,32,35,36,39,40,41},
   {58,60,61,63,64,65,67,68,69,70,86,88,89,91,92,93,5,8,9,12,13,14,33,36,37,40,41,42},
   {59,62,63,66,67,68,71,72,73,74,87,90,91,94,95,96,7,11,12,16,17,18,35,39,40,44,45,46},
   {60,63,64,67,68,69,72,73,74,75,88,91,92,95,96,97,8,12,13,17,18,19,36,40,41,45,46,47},
   {61,64,65,68,69,70,73,74,75,76,89,92,93,96,97,98,9,13,14,18,19,20,37,41,42,46,47,48},
   {62,66,67,71,72,73,77,78,79,80,90,94,95,99,100,101,11,16,17,22,23,24,39,44,45,50,51,52},
   {63,67,68,72,73,74,78,79,80,81,91,95,96,100,101,102,12,17,18,23,24,25,40,45,46,51,52,53},
   {64,68,69,73,74,75,79,80,81,82,92,96,97,101,102,103,13,18,19,24,25,26,41,46,47,52,53,54},
   {65,69,70,74,75,76,80,81,82,83,93,97,98,102,103,104,14,19,20,25,26,27,42,47,48,53,54,55},
   {84,85,86,87,88,89,90,91,92,93,57,59,60,62,63,64,30,32,33,35,36,37,4,7,8,11,12,13},
   {85,87,88,90,91,92,94,95,96,97,59,62,63,66,67,68,32,35,36,39,40,41,7,11,12,16,17,18},
   {86,88,89,91,92,93,95,96,97,98,60,63,64,67,68,69,33,36,37,40,41,42,8,12,13,17,18,19},
   {87,90,91,94,95,96,99,100,101,102,62,66,67,71,72,73,35,39,40,44,45,46,11,16,17,22,23,24},
   {88,91,92,95,96,97,100,101,102,103,63,67,68,72,73,74,36,40,41,45,46,47,12,17,18,23,24,25},
   {89,92,93,96,97,98,101,102,103,104,64,68,69,73,74,75,37,41,42,46,47,48,13,18,19,24,25,26},
};
// indices of [r6] x [r8] in [R14]
static unsigned short iCartP68[45][28] = {
   {0,1,2,3,4,5,6,7,8,9,36,37,38,39,40,41,64,65,66,67,68,69,92,93,94,95,96,97},
   {1,3,4,6,7,8,10,11,12,13,37,39,40,42,43,44,65,67,68,70,71,72,93,95,96,98,99,100},
   {2,4,5,7,8,9,11,12,13,14,38,40,41,43,44,45,66,68,69,71,72,73,94,96,97,99,100,101},
   {3,6,7,10,11,12,15,16,17,18,39,42,43,46,47,48,67,70,71,74,75,76,95,98,99,102,103,104},
   {4,7,8,11,12,13,16,17,18,19,40,43,44,47,48,49,68,71,72,75,76,77,96,99,100,103,104,105},
   {5,8,9,12,13,14,17,18,19,20,41,44,45,48,49,50,69,72,73,76,77,78,97,100,101,104,105,106},
   {6,10,11,15,16,17,21,22,23,24,42,46,47,51,52,53,70,74,75,79,80,81,98,102,103,107,108,109},
   {7,11,12,16,17,18,22,23,24,25,43,47,48,52,53,54,71,75,76,80,81,82,99,103,104,108,109,110},
   {8,12,13,17,18,19,23,24,25,26,44,48,49,53,54,55,72,76,77,81,82,83,100,104,105,109,110,111},
   {9,13,14,18,19,20,24,25,26,27,45,49,50,54,55,56,73,77,78,82,83,84,101,105,106,110,111,112},
   {10,15,16,21,22,23,28,29,30,31,46,51,52,57,58,59,74,79,80,85,86,87,102,107,108,113,114,115},
   {11,16,17,22,23,24,29,30,31,32,47,52,53,58,59,60,75,80,81,86,87,88,103,108,109,114,115,116},
   {12,17,18,23,24,25,30,31,32,33,48,53,54,59,60,61,76,81,82,87,88,89,104,109,110,115,116,117},
   {13,18,19,24,25,26,31,32,33,34,49,54,55,60,61,62,77,82,83,88,89,90,105,110,111,116,117,118},
   {14,19,20,25,26,27,32,33,34,35,50,55,56,61,62,63,78,83,84,89,90,91,106,111,112,117,118,119},
   {36,37,38,39,40,41,42,43,44,45,1,3,4,6,7,8,92,93,94,95,96,97,65,67,68,70,71,72},
   {37,39,40,42,43,44,46,47,48,49,3,6,7,10,11,12,93,95,96,98,99,100,67,70,71,74,75,76},
   {38,40,41,43,44,45,47,48,49,50,4,7,8,11,12,13,94,96,97,99,100,101,68,71,72,75,76,77},
   {39,42,43,46,47,48,51,52,53,54,6,10,11,15,16,17,95,98,99,102,103,104,70,74,75,79,80,81},
   {40,43,44,47,48,49,52,53,54,55,7,11,12,16,17,18,96,99,100,103,104,105,71,75,76,80,81,82},
   {41,44,45,48,49,50,53,54,55,56,8,12,13,17,18,19,97,100,101,104,105,106,72,76,77,81,82,83},
   {42,46,47,51,52,53,57,58,59,60,10,15,16,21,22,23,98,102,103,107,108,109,74,79,80,85,86,87},
   {43,47,48,52,53,54,58,59,60,61,11,16,17,22,23,24,99,103,104,108,109,110,75,80,81,86,87,88},
   {44,48,49,53,54,55,59,60,61,62,12,17,18,23,24,25,100,104,105,109,110,111,76,81,82,87,88,89},
   {45,49,50,54,55,56,60,61,62,63,13,18,19,24,25,26,101,105,106,110,111,112,77,82,83,88,89,90},
   {64,65,66,67,68,69,70,71,72,73,92,93,94,95,96,97,2,4,5,7,8,9,38,40,41,43,44,45},
   {65,67,68,70,71,72,74,75,76,77,93,95,96,98,99,100,4,7,8,11,12,13,40,43,44,47,48,49},
   {66,68,69,71,72,73,75,76,77,78,94,96,97,99,100,101,5,8,9,12,13,14,41,44,45,48,49,50},
   {67,70,71,74,75,76,79,80,81,82,95,98,99,102,103,104,7,11,12,16,17,18,43,47,48,52,53,54},
   {68,71,72,75,76,77,80,81,82,83,96,99,100,103,104,105,8,12,13,17,18,19,44,48,49,53,54,55},
   {69,72,73,76,77,78,81,82,83,84,97,100,101,104,105,106,9,13,14,18,19,20,45,49,50,54,55,56},
   {70,74,75,79,80,81,85,86,87,88,98,102,103,107,108,109,11,16,17,22,23,24,47,52,53,58,59,60},
   {71,75,76,80,81,82,86,87,88,89,99,103,104,108,109,110,12,17,18,23,24,25,48,53,54,59,60,61},
   {72,76,77,81,82,83,87,88,89,90,100,104,105,109,110,111,13,18,19,24,25,26,49,54,55,60,61,62},
   {73,77,78,82,83,84,88,89,90,91,101,105,106,110,111,112,14,19,20,25,26,27,50,55,56,61,62,63},
   {92,93,94,95,96,97,98,99,100,101,65,67,68,70,71,72,38,40,41,43,44,45,4,7,8,11,12,13},
   {93,95,96,98,99,100,102,103,104,105,67,70,71,74,75,76,40,43,44,47,48,49,7,11,12,16,17,18},
   {94,96,97,99,100,101,103,104,105,106,68,71,72,75,76,77,41,44,45,48,49,50,8,12,13,17,18,19},
   {95,98,99,102,103,104,107,108,109,110,70,74,75,79,80,81,43,47,48,52,53,54,11,16,17,22,23,24},
   {96,99,100,103,104,105,108,109,110,111,71,75,76,80,81,82,44,48,49,53,54,55,12,17,18,23,24,25},
   {97,100,101,104,105,106,109,110,111,112,72,76,77,81,82,83,45,49,50,54,55,56,13,18,19,24,25,26},
   {98,102,103,107,108,109,113,114,115,116,74,79,80,85,86,87,47,52,53,58,59,60,16,22,23,29,30,31},
   {99,103,104,108,109,110,114,115,116,117,75,80,81,86,87,88,48,53,54,59,60,61,17,23,24,30,31,32},
   {100,104,105,109,110,111,115,116,117,118,76,81,82,87,88,89,49,54,55,60,61,62,18,24,25,31,32,33},
   {101,105,106,110,111,112,116,117,118,119,77,82,83,88,89,90,50,55,56,61,62,63,19,25,26,32,33,34}
};
unsigned short *iCartPxx[15][7] = {
   {&iCartP00[0][0], &iCartP10[0][0], &iCartP20[0][0], &iCartP30[0][0], &iCartP40[0][0], &iCartP50[0][0], &iCartP60[0][0]},
   {&iCartP01[0][0], &iCartP11[0][0], &iCartP21[0][0], &iCartP31[0][0], &iCartP41[0][0], &iCartP51[0][0], &iCartP61[0][0]},
   {&iCartP02[0][0], &iCartP12[0][0], &iCartP22[0][0], &iCartP32[0][0], &iCartP42[0][0], &iCartP52[0][0], &iCartP62[0][0]},
   {&iCartP03[0][0], &iCartP13[0][0], &iCartP23[0][0], &iCartP33[0][0], &iCartP43[0][0], &iCartP53[0][0], &iCartP63[0][0]},
   {&iCartP04[0][0], &iCartP14[0][0], &iCartP24[0][0], &iCartP34[0][0], &iCartP44[0][0], &iCartP54[0][0], &iCartP64[0][0]},
   {&iCartP05[0][0], &iCartP15[0][0], &iCartP25[0][0], &iCartP35[0][0], &iCartP45[0][0], &iCartP55[0][0], &iCartP65[0][0]},
   {&iCartP06[0][0], &iCartP16[0][0], &iCartP26[0][0], &iCartP36[0][0], &iCartP46[0][0], &iCartP56[0][0], &iCartP66[0][0]},
   {&iCartP07[0][0], &iCartP17[0][0], &iCartP27[0][0], &iCartP37[0][0], &iCartP47[0][0], &iCartP57[0][0], &iCartP67[0][0]},
   {&iCartP08[0][0], &iCartP18[0][0], &iCartP28[0][0], &iCartP38[0][0], &iCartP48[0][0], &iCartP58[0][0], &iCartP68[0][0]},
   {&iCartP09[0][0], &iCartP19[0][0], &iCartP29[0][0], &iCartP39[0][0], &iCartP49[0][0], &iCartP59[0][0], 0},
   {&iCartP0a[0][0], &iCartP1a[0][0], &iCartP2a[0][0], &iCartP3a[0][0], &iCartP4a[0][0], 0, 0},
   {&iCartP0b[0][0], &iCartP1b[0][0], &iCartP2b[0][0], &iCartP3b[0][0], 0, 0, 0},
   {&iCartP0c[0][0], &iCartP1c[0][0], &iCartP2c[0][0], 0, 0, 0, 0},
   {&iCartP0d[0][0], &iCartP1d[0][0], 0, 0, 0, 0, 0},
   {&iCartP0e[0][0], 0, 0, 0, 0, 0, 0}
};
} // namespace aic

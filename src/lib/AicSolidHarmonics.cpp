/* Copyright (c) 2012  Gerald Knizia
 * 
 * This file is part of the bfint program
 * (See http://www.theochem.uni-stuttgart.de/~knizia)
 * 
 * bfint is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 * 
 * bfint is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with bfint (LICENSE). If not, see http://www.gnu.org/licenses/
 */

/* AicSolidHarmonics.cpp v20121011 EST [charge, Gerald Knizia] */

#include "AicSolidHarmonics.h"

namespace aic {
unsigned char SlmSymSig[49] = {
   0, 1, 2, 4, 0, 3, 5, 0, 6, 1, 2, 4, 1, 7, 2, 4, 0, 3, 5, 0, 6, 0, 3, 5, 6, 1, 2, 4, 1, 7, 2, 4,
   2, 4, 1, 7, 0, 3, 5, 0, 6, 0, 3, 5, 3, 0, 6, 6, 5
};

// constants for direct expansion of Slms in terms of cartesians
static double const sd0 = 5.e-01;
static double const sd1 = 1.7320508075688772;
static double const sd2 = 8.660254037844386e-01;
static double const sd3 = 6.1237243569579447e-01;
static double const sd4 = 2.4494897427831779;
static double const sd5 = 1.5;
static double const sd6 = 2.3717082451262841;
static double const sd7 = 7.9056941504209477e-01;
static double const sd8 = 3.8729833462074166;
static double const sd9 = 1.9364916731037083;
static double const sd10 = 3.;
static double const sd11 = 3.75e-01;
static double const sd12 = 7.5e-01;
static double const sd13 = 1.1180339887498947;
static double const sd14 = 6.7082039324993685;
static double const sd15 = 3.1622776601683791;
static double const sd16 = 7.3950997288745191e-01;
static double const sd17 = 4.4370598373247114;
static double const sd18 = 3.3541019662496843;
static double const sd19 = 5.5901699437494734e-01;
static double const sd20 = 2.9580398915498076;
static double const sd21 = 6.2749501990055663;
static double const sd22 = 2.0916500663351885;
static double const sd23 = 5.8094750193111251;
static double const sd24 = 9.6824583655185414e-01;
static double const sd25 = 4.8412291827592707e-01;
static double const sd26 = 2.5617376914898995;
static double const sd27 = 5.123475382979799;
static double const sd28 = 1.2549900398011133e+01;
static double const sd29 = 4.1833001326703769;
static double const sd30 = 1.0458250331675942;
static double const sd31 = 1.5687375497513916;
static double const sd32 = 5.2291251658379712e-01;
static double const sd33 = 8.8741196746494229;
static double const sd34 = 2.2185299186623557;
static double const sd35 = 1.3311179511974135e+01;
static double const sd36 = 7.0156076002011392e-01;
static double const sd37 = 3.5078038001005698;
static double const sd38 = 7.0156076002011396;
static double const sd39 = 1.875;
static double const sd40 = 5.;
static double const sd41 = 3.75;
static double const sd42 = 1.0246950765959598e+01;
static double const sd43 = 6.7169328938139605e-01;
static double const sd44 = 1.0075399340720942e+01;
static double const sd45 = 1.4491376746189438e+01;
static double const sd46 = 9.0571104663683988e-01;
static double const sd47 = 1.8114220932736798;
static double const sd48 = 1.1634069043116428e+01;
static double const sd49 = 2.3268138086232852;
static double const sd50 = 2.3268138086232856e+01;
static double const sd51 = 2.9764702249476642e+01;
static double const sd52 = 4.9607837082461073;
static double const sd53 = 2.4803918541230536;
static double const sd54 = 4.960783708246107e-01;
static double const sd55 = 7.245688373094719;
static double const sd56 = 4.5285552331841994e-01;
static double const sd57 = 1.3433865787627923e+01;
static double const sd58 = 4.0301597362883763;
static double const sd59 = 8.1513994197315576;
static double const sd60 = 5.4342662798210393;
static double const sd61 = 2.7171331399105196;
static double const sd62 = 2.1737065119284157e+01;
static double const sd63 = 1.9843134832984428;
static double const sd64 = 1.9843134832984429e+01;
static double const sd65 = 1.125e+01;
static double const sd66 = 5.625;
static double const sd67 = 7.5;
static double const sd68 = 9.375e-01;
static double const sd69 = 3.125e-01;
static double const sd70 = 1.1456439237389599e+01;
static double const sd71 = 2.8641098093473998;
static double const sd72 = 4.5825756949558398;
static double const sd73 = 5.7282196186947996;

// constants for incremental expansion of Slms in terms of lower Slms.
static double const si0 = 0.;
static double const si1 = 8.660254037844386e-01;
static double const si2 = 1.5;
static double const si3 = 5.e-01;
static double const si4 = 1.7320508075688774;
static double const si5 = 9.128709291752769e-01;
static double const si6 = 1.7677669529663687;
static double const si7 = 6.1237243569579447e-01;
static double const si8 = 1.6666666666666667;
static double const si9 = 6.6666666666666663e-01;
static double const si10 = 2.2360679774997898;
static double const si11 = 9.3541434669348533e-01;
static double const si12 = 1.75;
static double const si13 = 7.5e-01;
static double const si14 = 2.0207259421636903;
static double const si15 = 6.454972243679028e-01;
static double const si16 = 1.8073922282301278;
static double const si17 = 7.3029674334022143e-01;
static double const si18 = 2.6457513110645903;
static double const si19 = 9.4868329805051377e-01;
static double const si20 = 1.8371173070873836;
static double const si21 = 7.9056941504209488e-01;
static double const si22 = 1.9639610121239315;
static double const si23 = 7.559289460184544e-01;
static double const si24 = 2.25;
static double const si25 = 6.6143782776614768e-01;
static double const si26 = 3.;
static double const si27 = 1.8;
static double const si28 = 8.0000000000000004e-01;
static double const si29 = 9.574271077563381e-01;
static double const si30 = 1.9445436482630056;
static double const si31 = 8.1009258730098255e-01;
static double const si32 = 3.3166247903553998;
static double const si33 = 2.4596747752497685;
static double const si34 = 6.7082039324993692e-01;
static double const si35 = 2.1169509870286278;
static double const si36 = 7.6980035891950105e-01;
static double const si37 = 1.8333333333333333;
static double const si38 = 8.3333333333333337e-01;
static double const si39 = 1.8593393604027364;
static double const si40 = 8.2807867121082501e-01;

// evaluate all S^l_m(R) for l = 0 .. L. Results stored in Out[iSlmA(l,m)].
void EvalSlmX( double *AIC_RP Out, double const *AIC_RP R, int L )
{
   // S(l=0, m=-0..0)
   Out[0] = 1;
   if ( 0 == L ) return;

   double
       x = R[0], y = R[1], z = R[2];
   // S(l=1, m=-1..1)
   Out[1] = x;
   Out[2] = y;
   Out[3] = z;
   if ( 1 == L ) return;

   double
       rsq = x*x + y*y + z*z;
   // S(l=2, m=-2..2)
   Out[4] = si2 * z * z - si3 * rsq;
   Out[5] = si1 * (y * x + x * y);
   Out[6] = si4 * z * x;
   Out[7] = si1 * (x * x - y * y);
   Out[8] = si4 * z * y;
   if ( 2 == L ) return;

   // S(l=3, m=-3..3)
   Out[9] = si6 * z * Out[6] - si7 * rsq * x;
   Out[10] = si6 * z * Out[8] - si7 * rsq * y;
   Out[11] = si8 * z * Out[4] - si9 * rsq * z;
   Out[12] = si5 * (x * Out[7] - y * Out[5]);
   Out[13] = si10 * z * Out[5];
   Out[14] = si5 * (y * Out[7] + x * Out[5]);
   Out[15] = si10 * z * Out[7];
   if ( 3 == L ) return;

   // S(l=4, m=-4..4)
   Out[16] = si12 * z * Out[11] - si13 * rsq * Out[4];
   Out[17] = si14 * z * Out[13] - si15 * rsq * Out[5];
   Out[18] = si16 * z * Out[9] - si17 * rsq * Out[6];
   Out[19] = si11 * (x * Out[12] - y * Out[14]);
   Out[20] = si16 * z * Out[10] - si17 * rsq * Out[8];
   Out[21] = si14 * z * Out[15] - si15 * rsq * Out[7];
   Out[22] = si11 * (y * Out[12] + x * Out[14]);
   Out[23] = si18 * z * Out[12];
   Out[24] = si18 * z * Out[14];
   if ( 4 == L ) return;

   // S(l=5, m=-5..5)
   Out[25] = si20 * z * Out[18] - si21 * rsq * Out[9];
   Out[26] = si20 * z * Out[20] - si21 * rsq * Out[10];
   Out[27] = si22 * z * Out[21] - si23 * rsq * Out[15];
   Out[28] = si24 * z * Out[23] - si25 * rsq * Out[12];
   Out[29] = si26 * z * Out[22];
   Out[30] = si24 * z * Out[24] - si25 * rsq * Out[14];
   Out[31] = si26 * z * Out[19];
   Out[32] = si19 * (y * Out[19] + x * Out[22]);
   Out[33] = si27 * z * Out[16] - si28 * rsq * Out[11];
   Out[34] = si19 * (x * Out[19] - y * Out[22]);
   Out[35] = si22 * z * Out[17] - si23 * rsq * Out[13];
   if ( 5 == L ) return;

   // S(l=6, m=-6..6)
   Out[36] = si29 * (x * Out[34] - y * Out[32]);
   Out[37] = si30 * z * Out[35] - si31 * rsq * Out[17];
   Out[38] = si32 * z * Out[34];
   Out[39] = si33 * z * Out[31] - si34 * rsq * Out[19];
   Out[40] = si32 * z * Out[32];
   Out[41] = si30 * z * Out[27] - si31 * rsq * Out[21];
   Out[42] = si29 * (y * Out[34] + x * Out[32]);
   Out[43] = si35 * z * Out[28] - si36 * rsq * Out[23];
   Out[44] = si33 * z * Out[29] - si34 * rsq * Out[22];
   Out[45] = si37 * z * Out[33] - si38 * rsq * Out[16];
   Out[46] = si35 * z * Out[30] - si36 * rsq * Out[24];
   Out[47] = si39 * z * Out[26] - si40 * rsq * Out[20];
   Out[48] = si39 * z * Out[25] - si40 * rsq * Out[18];
   if ( 6 == L ) return;

   AicAssertFail("L too large.", __FILE__, __LINE__);
}
// evaluate all S^l_m(R) and d/dx_i S^l_m(R) for l = 0 .. L.
// Results stored in Out[4*iSlmA(l,m) + iComp]. iComp = 0: value, =1,2,3: d/dx_i.
void EvalSlmX_Grd1( double *AIC_RP Out, double const *AIC_RP R, int L )
{
   // [1, d/dx, d/dy, d/dz] S(l=0, m=-0..0)
   Out[0] = 1.;
   Out[1] = 0.;
   Out[2] = 0.;
   Out[3] = 0.;
   if ( 0 == L ) return;

   double
       x = R[0], y = R[1], z = R[2];
   // [1, d/dx, d/dy, d/dz] S(l=1, m=-1..1)
   Out[4] = x;
   Out[5] = 1.;
   Out[6] = 0.;
   Out[7] = 0.;
   Out[8] = y;
   Out[9] = 0.;
   Out[10] = 1.;
   Out[11] = 0.;
   Out[12] = z;
   Out[13] = 0.;
   Out[14] = 0.;
   Out[15] = 1.;
   if ( 1 == L ) return;

   double
       rsq = x*x + y*y + z*z;
   // [1, d/dx, d/dy, d/dz] S(l=2, m=-2..2)
   Out[16] = si2 * z * z - si3 * rsq;
   Out[17] = -si3 * 2*x;
   Out[18] = -si3 * 2*y;
   Out[19] = si2 * 2*z - si3 * 2*z;
   Out[20] = si1 * (y * x + x * y);
   Out[21] = si1 * (y + y);
   Out[22] = si1 * (x + x);
   Out[23] = 0.;
   Out[24] = si4 * z * x;
   Out[25] = si4 * z;
   Out[26] = 0.;
   Out[27] = si4 * x;
   Out[28] = si1 * (x * x - y * y);
   Out[29] = si1 * (2*x);
   Out[30] = si1 * (-2*y);
   Out[31] = 0.;
   Out[32] = si4 * z * y;
   Out[33] = 0.;
   Out[34] = si4 * z;
   Out[35] = si4 * y;
   if ( 2 == L ) return;

   // [1, d/dx, d/dy, d/dz] S(l=3, m=-3..3)
   Out[36] = si6 * z * Out[24] - si7 * rsq * x;
   Out[37] = si6 * z * Out[25] - (si7 * (rsq + 2*x * x));
   Out[38] = -si7 * 2*y * x;
   Out[39] = (si6 * (z * Out[27] + Out[24])) - si7 * 2*z * x;
   Out[40] = si6 * z * Out[32] - si7 * rsq * y;
   Out[41] = -si7 * 2*x * y;
   Out[42] = si6 * z * Out[34] - (si7 * (rsq + 2*y * y));
   Out[43] = (si6 * (z * Out[35] + Out[32])) - si7 * 2*z * y;
   Out[44] = si8 * z * Out[16] - si9 * rsq * z;
   Out[45] = si8 * z * Out[17] - si9 * 2*x * z;
   Out[46] = si8 * z * Out[18] - si9 * 2*y * z;
   Out[47] = (si8 * (z * Out[19] + Out[16])) - (si9 * (rsq + 2*z * z));
   Out[48] = si5 * (x * Out[28] - y * Out[20]);
   Out[49] = si5 * ((x * Out[29] + Out[28]) - y * Out[21]);
   Out[50] = si5 * (x * Out[30] - (y * Out[22] + Out[20]));
   Out[51] = 0.;
   Out[52] = si10 * z * Out[20];
   Out[53] = si10 * z * Out[21];
   Out[54] = si10 * z * Out[22];
   Out[55] = si10 * Out[20];
   Out[56] = si5 * (y * Out[28] + x * Out[20]);
   Out[57] = si5 * (y * Out[29] + (x * Out[21] + Out[20]));
   Out[58] = si5 * ((y * Out[30] + Out[28]) + x * Out[22]);
   Out[59] = 0.;
   Out[60] = si10 * z * Out[28];
   Out[61] = si10 * z * Out[29];
   Out[62] = si10 * z * Out[30];
   Out[63] = si10 * Out[28];
   if ( 3 == L ) return;

   // [1, d/dx, d/dy, d/dz] S(l=4, m=-4..4)
   Out[64] = si12 * z * Out[44] - si13 * rsq * Out[16];
   Out[65] = si12 * z * Out[45] - (si13 * (rsq * Out[17] + 2*x * Out[16]));
   Out[66] = si12 * z * Out[46] - (si13 * (rsq * Out[18] + 2*y * Out[16]));
   Out[67] = (si12 * (z * Out[47] + Out[44])) - (si13 * (rsq * Out[19] + 2*z * Out[16]));
   Out[68] = si14 * z * Out[52] - si15 * rsq * Out[20];
   Out[69] = si14 * z * Out[53] - (si15 * (rsq * Out[21] + 2*x * Out[20]));
   Out[70] = si14 * z * Out[54] - (si15 * (rsq * Out[22] + 2*y * Out[20]));
   Out[71] = (si14 * (z * Out[55] + Out[52])) - si15 * 2*z * Out[20];
   Out[72] = si16 * z * Out[36] - si17 * rsq * Out[24];
   Out[73] = si16 * z * Out[37] - (si17 * (rsq * Out[25] + 2*x * Out[24]));
   Out[74] = si16 * z * Out[38] - si17 * 2*y * Out[24];
   Out[75] = (si16 * (z * Out[39] + Out[36])) - (si17 * (rsq * Out[27] + 2*z * Out[24]));
   Out[76] = si11 * (x * Out[48] - y * Out[56]);
   Out[77] = si11 * ((x * Out[49] + Out[48]) - y * Out[57]);
   Out[78] = si11 * (x * Out[50] - (y * Out[58] + Out[56]));
   Out[79] = 0.;
   Out[80] = si16 * z * Out[40] - si17 * rsq * Out[32];
   Out[81] = si16 * z * Out[41] - si17 * 2*x * Out[32];
   Out[82] = si16 * z * Out[42] - (si17 * (rsq * Out[34] + 2*y * Out[32]));
   Out[83] = (si16 * (z * Out[43] + Out[40])) - (si17 * (rsq * Out[35] + 2*z * Out[32]));
   Out[84] = si14 * z * Out[60] - si15 * rsq * Out[28];
   Out[85] = si14 * z * Out[61] - (si15 * (rsq * Out[29] + 2*x * Out[28]));
   Out[86] = si14 * z * Out[62] - (si15 * (rsq * Out[30] + 2*y * Out[28]));
   Out[87] = (si14 * (z * Out[63] + Out[60])) - si15 * 2*z * Out[28];
   Out[88] = si11 * (y * Out[48] + x * Out[56]);
   Out[89] = si11 * (y * Out[49] + (x * Out[57] + Out[56]));
   Out[90] = si11 * ((y * Out[50] + Out[48]) + x * Out[58]);
   Out[91] = 0.;
   Out[92] = si18 * z * Out[48];
   Out[93] = si18 * z * Out[49];
   Out[94] = si18 * z * Out[50];
   Out[95] = si18 * Out[48];
   Out[96] = si18 * z * Out[56];
   Out[97] = si18 * z * Out[57];
   Out[98] = si18 * z * Out[58];
   Out[99] = si18 * Out[56];
   if ( 4 == L ) return;

   // [1, d/dx, d/dy, d/dz] S(l=5, m=-5..5)
   Out[100] = si20 * z * Out[72] - si21 * rsq * Out[36];
   Out[101] = si20 * z * Out[73] - (si21 * (rsq * Out[37] + 2*x * Out[36]));
   Out[102] = si20 * z * Out[74] - (si21 * (rsq * Out[38] + 2*y * Out[36]));
   Out[103] = (si20 * (z * Out[75] + Out[72])) - (si21 * (rsq * Out[39] + 2*z * Out[36]));
   Out[104] = si20 * z * Out[80] - si21 * rsq * Out[40];
   Out[105] = si20 * z * Out[81] - (si21 * (rsq * Out[41] + 2*x * Out[40]));
   Out[106] = si20 * z * Out[82] - (si21 * (rsq * Out[42] + 2*y * Out[40]));
   Out[107] = (si20 * (z * Out[83] + Out[80])) - (si21 * (rsq * Out[43] + 2*z * Out[40]));
   Out[108] = si22 * z * Out[84] - si23 * rsq * Out[60];
   Out[109] = si22 * z * Out[85] - (si23 * (rsq * Out[61] + 2*x * Out[60]));
   Out[110] = si22 * z * Out[86] - (si23 * (rsq * Out[62] + 2*y * Out[60]));
   Out[111] = (si22 * (z * Out[87] + Out[84])) - (si23 * (rsq * Out[63] + 2*z * Out[60]));
   Out[112] = si24 * z * Out[92] - si25 * rsq * Out[48];
   Out[113] = si24 * z * Out[93] - (si25 * (rsq * Out[49] + 2*x * Out[48]));
   Out[114] = si24 * z * Out[94] - (si25 * (rsq * Out[50] + 2*y * Out[48]));
   Out[115] = (si24 * (z * Out[95] + Out[92])) - si25 * 2*z * Out[48];
   Out[116] = si26 * z * Out[88];
   Out[117] = si26 * z * Out[89];
   Out[118] = si26 * z * Out[90];
   Out[119] = si26 * Out[88];
   Out[120] = si24 * z * Out[96] - si25 * rsq * Out[56];
   Out[121] = si24 * z * Out[97] - (si25 * (rsq * Out[57] + 2*x * Out[56]));
   Out[122] = si24 * z * Out[98] - (si25 * (rsq * Out[58] + 2*y * Out[56]));
   Out[123] = (si24 * (z * Out[99] + Out[96])) - si25 * 2*z * Out[56];
   Out[124] = si26 * z * Out[76];
   Out[125] = si26 * z * Out[77];
   Out[126] = si26 * z * Out[78];
   Out[127] = si26 * Out[76];
   Out[128] = si19 * (y * Out[76] + x * Out[88]);
   Out[129] = si19 * (y * Out[77] + (x * Out[89] + Out[88]));
   Out[130] = si19 * ((y * Out[78] + Out[76]) + x * Out[90]);
   Out[131] = 0.;
   Out[132] = si27 * z * Out[64] - si28 * rsq * Out[44];
   Out[133] = si27 * z * Out[65] - (si28 * (rsq * Out[45] + 2*x * Out[44]));
   Out[134] = si27 * z * Out[66] - (si28 * (rsq * Out[46] + 2*y * Out[44]));
   Out[135] = (si27 * (z * Out[67] + Out[64])) - (si28 * (rsq * Out[47] + 2*z * Out[44]));
   Out[136] = si19 * (x * Out[76] - y * Out[88]);
   Out[137] = si19 * ((x * Out[77] + Out[76]) - y * Out[89]);
   Out[138] = si19 * (x * Out[78] - (y * Out[90] + Out[88]));
   Out[139] = 0.;
   Out[140] = si22 * z * Out[68] - si23 * rsq * Out[52];
   Out[141] = si22 * z * Out[69] - (si23 * (rsq * Out[53] + 2*x * Out[52]));
   Out[142] = si22 * z * Out[70] - (si23 * (rsq * Out[54] + 2*y * Out[52]));
   Out[143] = (si22 * (z * Out[71] + Out[68])) - (si23 * (rsq * Out[55] + 2*z * Out[52]));
   if ( 5 == L ) return;

   // [1, d/dx, d/dy, d/dz] S(l=6, m=-6..6)
   Out[144] = si29 * (x * Out[136] - y * Out[128]);
   Out[145] = si29 * ((x * Out[137] + Out[136]) - y * Out[129]);
   Out[146] = si29 * (x * Out[138] - (y * Out[130] + Out[128]));
   Out[147] = 0.;
   Out[148] = si30 * z * Out[140] - si31 * rsq * Out[68];
   Out[149] = si30 * z * Out[141] - (si31 * (rsq * Out[69] + 2*x * Out[68]));
   Out[150] = si30 * z * Out[142] - (si31 * (rsq * Out[70] + 2*y * Out[68]));
   Out[151] = (si30 * (z * Out[143] + Out[140])) - (si31 * (rsq * Out[71] + 2*z * Out[68]));
   Out[152] = si32 * z * Out[136];
   Out[153] = si32 * z * Out[137];
   Out[154] = si32 * z * Out[138];
   Out[155] = si32 * Out[136];
   Out[156] = si33 * z * Out[124] - si34 * rsq * Out[76];
   Out[157] = si33 * z * Out[125] - (si34 * (rsq * Out[77] + 2*x * Out[76]));
   Out[158] = si33 * z * Out[126] - (si34 * (rsq * Out[78] + 2*y * Out[76]));
   Out[159] = (si33 * (z * Out[127] + Out[124])) - si34 * 2*z * Out[76];
   Out[160] = si32 * z * Out[128];
   Out[161] = si32 * z * Out[129];
   Out[162] = si32 * z * Out[130];
   Out[163] = si32 * Out[128];
   Out[164] = si30 * z * Out[108] - si31 * rsq * Out[84];
   Out[165] = si30 * z * Out[109] - (si31 * (rsq * Out[85] + 2*x * Out[84]));
   Out[166] = si30 * z * Out[110] - (si31 * (rsq * Out[86] + 2*y * Out[84]));
   Out[167] = (si30 * (z * Out[111] + Out[108])) - (si31 * (rsq * Out[87] + 2*z * Out[84]));
   Out[168] = si29 * (y * Out[136] + x * Out[128]);
   Out[169] = si29 * (y * Out[137] + (x * Out[129] + Out[128]));
   Out[170] = si29 * ((y * Out[138] + Out[136]) + x * Out[130]);
   Out[171] = 0.;
   Out[172] = si35 * z * Out[112] - si36 * rsq * Out[92];
   Out[173] = si35 * z * Out[113] - (si36 * (rsq * Out[93] + 2*x * Out[92]));
   Out[174] = si35 * z * Out[114] - (si36 * (rsq * Out[94] + 2*y * Out[92]));
   Out[175] = (si35 * (z * Out[115] + Out[112])) - (si36 * (rsq * Out[95] + 2*z * Out[92]));
   Out[176] = si33 * z * Out[116] - si34 * rsq * Out[88];
   Out[177] = si33 * z * Out[117] - (si34 * (rsq * Out[89] + 2*x * Out[88]));
   Out[178] = si33 * z * Out[118] - (si34 * (rsq * Out[90] + 2*y * Out[88]));
   Out[179] = (si33 * (z * Out[119] + Out[116])) - si34 * 2*z * Out[88];
   Out[180] = si37 * z * Out[132] - si38 * rsq * Out[64];
   Out[181] = si37 * z * Out[133] - (si38 * (rsq * Out[65] + 2*x * Out[64]));
   Out[182] = si37 * z * Out[134] - (si38 * (rsq * Out[66] + 2*y * Out[64]));
   Out[183] = (si37 * (z * Out[135] + Out[132])) - (si38 * (rsq * Out[67] + 2*z * Out[64]));
   Out[184] = si35 * z * Out[120] - si36 * rsq * Out[96];
   Out[185] = si35 * z * Out[121] - (si36 * (rsq * Out[97] + 2*x * Out[96]));
   Out[186] = si35 * z * Out[122] - (si36 * (rsq * Out[98] + 2*y * Out[96]));
   Out[187] = (si35 * (z * Out[123] + Out[120])) - (si36 * (rsq * Out[99] + 2*z * Out[96]));
   Out[188] = si39 * z * Out[104] - si40 * rsq * Out[80];
   Out[189] = si39 * z * Out[105] - (si40 * (rsq * Out[81] + 2*x * Out[80]));
   Out[190] = si39 * z * Out[106] - (si40 * (rsq * Out[82] + 2*y * Out[80]));
   Out[191] = (si39 * (z * Out[107] + Out[104])) - (si40 * (rsq * Out[83] + 2*z * Out[80]));
   Out[192] = si39 * z * Out[100] - si40 * rsq * Out[72];
   Out[193] = si39 * z * Out[101] - (si40 * (rsq * Out[73] + 2*x * Out[72]));
   Out[194] = si39 * z * Out[102] - (si40 * (rsq * Out[74] + 2*y * Out[72]));
   Out[195] = (si39 * (z * Out[103] + Out[100])) - (si40 * (rsq * Out[75] + 2*z * Out[72]));
   if ( 6 == L ) return;

   AicAssertFail("L too large.", __FILE__, __LINE__);
}
// evaluate all S^l_m(R), d/dx_i S^l_m(R), and d/dx_i d/dx_j S^l_m(R) for l = 0 .. L.
// Results stored in Out[10*iSlmA(l,m) + iComp]. iComp = 0: value, =1,2,3: d/dx_i; 4..10: d/dx_i d/dx_j
void EvalSlmX_Grd2( double *AIC_RP Out, double const *AIC_RP R, int L )
{
   // [1, d/dx, d/dy, d/dz, d^2/d(x x), d^2/d(y y), d^2/d(z z), d^2/d(x y), d^2/d(x z), d^2/d(y z)], S(l=0, m=-0..0)
   Out[0] = 1.;
   Out[1] = 0.;
   Out[2] = 0.;
   Out[3] = 0.;
   Out[4] = 0.;
   Out[5] = 0.;
   Out[6] = 0.;
   Out[7] = 0.;
   Out[8] = 0.;
   Out[9] = 0.;
   if ( 0 == L ) return;

   double
       x = R[0], y = R[1], z = R[2];
   // [1, d/dx, d/dy, d/dz, d^2/d(x x), d^2/d(y y), d^2/d(z z), d^2/d(x y), d^2/d(x z), d^2/d(y z)], S(l=1, m=-1..1)
   Out[10] = x;
   Out[11] = 1.;
   Out[12] = 0.;
   Out[13] = 0.;
   Out[14] = 0.;
   Out[15] = 0.;
   Out[16] = 0.;
   Out[17] = 0.;
   Out[18] = 0.;
   Out[19] = 0.;
   Out[20] = y;
   Out[21] = 0.;
   Out[22] = 1.;
   Out[23] = 0.;
   Out[24] = 0.;
   Out[25] = 0.;
   Out[26] = 0.;
   Out[27] = 0.;
   Out[28] = 0.;
   Out[29] = 0.;
   Out[30] = z;
   Out[31] = 0.;
   Out[32] = 0.;
   Out[33] = 1.;
   Out[34] = 0.;
   Out[35] = 0.;
   Out[36] = 0.;
   Out[37] = 0.;
   Out[38] = 0.;
   Out[39] = 0.;
   if ( 1 == L ) return;

   double
       rsq = x*x + y*y + z*z;
   // [1, d/dx, d/dy, d/dz, d^2/d(x x), d^2/d(y y), d^2/d(z z), d^2/d(x y), d^2/d(x z), d^2/d(y z)], S(l=2, m=-2..2)
   Out[40] = si2 * z * z - si3 * rsq;
   Out[41] = -si3 * 2*x;
   Out[42] = -si3 * 2*y;
   Out[43] = si2 * 2*z - si3 * 2*z;
   Out[44] = -si3 * 2;
   Out[45] = -si3 * 2;
   Out[46] = si2 * 2 - si3 * 2;
   Out[47] = 0.;
   Out[48] = 0.;
   Out[49] = 0.;
   Out[50] = si1 * (y * x + x * y);
   Out[51] = si1 * (y + y);
   Out[52] = si1 * (x + x);
   Out[53] = 0.;
   Out[54] = 0.;
   Out[55] = 0.;
   Out[56] = 0.;
   Out[57] = si1 * (1. + 1.);
   Out[58] = 0.;
   Out[59] = 0.;
   Out[60] = si4 * z * x;
   Out[61] = si4 * z;
   Out[62] = 0.;
   Out[63] = si4 * x;
   Out[64] = 0.;
   Out[65] = 0.;
   Out[66] = 0.;
   Out[67] = 0.;
   Out[68] = si4;
   Out[69] = 0.;
   Out[70] = si1 * (x * x - y * y);
   Out[71] = si1 * (2*x);
   Out[72] = si1 * (-2*y);
   Out[73] = 0.;
   Out[74] = si1 * (2);
   Out[75] = si1 * (-2);
   Out[76] = 0.;
   Out[77] = 0.;
   Out[78] = 0.;
   Out[79] = 0.;
   Out[80] = si4 * z * y;
   Out[81] = 0.;
   Out[82] = si4 * z;
   Out[83] = si4 * y;
   Out[84] = 0.;
   Out[85] = 0.;
   Out[86] = 0.;
   Out[87] = 0.;
   Out[88] = 0.;
   Out[89] = si4;
   if ( 2 == L ) return;

   // [1, d/dx, d/dy, d/dz, d^2/d(x x), d^2/d(y y), d^2/d(z z), d^2/d(x y), d^2/d(x z), d^2/d(y z)], S(l=3, m=-3..3)
   Out[90] = si6 * z * Out[60] - si7 * rsq * x;
   Out[91] = si6 * z * Out[61] - (si7 * (rsq + 2*x * x));
   Out[92] = -si7 * 2*y * x;
   Out[93] = (si6 * (z * Out[63] + Out[60])) - si7 * 2*z * x;
   Out[94] = -(si7 * (4*x + 2 * x));
   Out[95] = -si7 * 2 * x;
   Out[96] = si6 * 2 * Out[63] - si7 * 2 * x;
   Out[97] = -si7 * 2*y;
   Out[98] = (si6 * (z * Out[68] + Out[61])) - si7 * 2*z;
   Out[99] = 0.;
   Out[100] = si6 * z * Out[80] - si7 * rsq * y;
   Out[101] = -si7 * 2*x * y;
   Out[102] = si6 * z * Out[82] - (si7 * (rsq + 2*y * y));
   Out[103] = (si6 * (z * Out[83] + Out[80])) - si7 * 2*z * y;
   Out[104] = -si7 * 2 * y;
   Out[105] = -(si7 * (4*y + 2 * y));
   Out[106] = si6 * 2 * Out[83] - si7 * 2 * y;
   Out[107] = -si7 * 2*x;
   Out[108] = 0.;
   Out[109] = (si6 * (z * Out[89] + Out[82])) - si7 * 2*z;
   Out[110] = si8 * z * Out[40] - si9 * rsq * z;
   Out[111] = si8 * z * Out[41] - si9 * 2*x * z;
   Out[112] = si8 * z * Out[42] - si9 * 2*y * z;
   Out[113] = (si8 * (z * Out[43] + Out[40])) - (si9 * (rsq + 2*z * z));
   Out[114] = si8 * z * Out[44] - si9 * 2 * z;
   Out[115] = si8 * z * Out[45] - si9 * 2 * z;
   Out[116] = (si8 * (z * Out[46] + 2 * Out[43])) - (si9 * (4*z + 2 * z));
   Out[117] = 0.;
   Out[118] = si8 * Out[41] - si9 * 2*x;
   Out[119] = si8 * Out[42] - si9 * 2*y;
   Out[120] = si5 * (x * Out[70] - y * Out[50]);
   Out[121] = si5 * ((x * Out[71] + Out[70]) - y * Out[51]);
   Out[122] = si5 * (x * Out[72] - (y * Out[52] + Out[50]));
   Out[123] = 0.;
   Out[124] = si5 * ((x * Out[74] + 2 * Out[71]));
   Out[125] = si5 * (x * Out[75] - 2 * Out[52]);
   Out[126] = 0.;
   Out[127] = si5 * (Out[72] - (y * Out[57] + Out[51]));
   Out[128] = 0.;
   Out[129] = 0.;
   Out[130] = si10 * z * Out[50];
   Out[131] = si10 * z * Out[51];
   Out[132] = si10 * z * Out[52];
   Out[133] = si10 * Out[50];
   Out[134] = 0.;
   Out[135] = 0.;
   Out[136] = 0.;
   Out[137] = si10 * z * Out[57];
   Out[138] = si10 * Out[51];
   Out[139] = si10 * Out[52];
   Out[140] = si5 * (y * Out[70] + x * Out[50]);
   Out[141] = si5 * (y * Out[71] + (x * Out[51] + Out[50]));
   Out[142] = si5 * ((y * Out[72] + Out[70]) + x * Out[52]);
   Out[143] = 0.;
   Out[144] = si5 * (y * Out[74] + 2 * Out[51]);
   Out[145] = si5 * ((y * Out[75] + 2 * Out[72]));
   Out[146] = 0.;
   Out[147] = si5 * (Out[71] + (x * Out[57] + Out[52]));
   Out[148] = 0.;
   Out[149] = 0.;
   Out[150] = si10 * z * Out[70];
   Out[151] = si10 * z * Out[71];
   Out[152] = si10 * z * Out[72];
   Out[153] = si10 * Out[70];
   Out[154] = si10 * z * Out[74];
   Out[155] = si10 * z * Out[75];
   Out[156] = 0.;
   Out[157] = 0.;
   Out[158] = si10 * Out[71];
   Out[159] = si10 * Out[72];
   if ( 3 == L ) return;

   // [1, d/dx, d/dy, d/dz, d^2/d(x x), d^2/d(y y), d^2/d(z z), d^2/d(x y), d^2/d(x z), d^2/d(y z)], S(l=4, m=-4..4)
   Out[160] = si12 * z * Out[110] - si13 * rsq * Out[40];
   Out[161] = si12 * z * Out[111] - (si13 * (rsq * Out[41] + 2*x * Out[40]));
   Out[162] = si12 * z * Out[112] - (si13 * (rsq * Out[42] + 2*y * Out[40]));
   Out[163] = (si12 * (z * Out[113] + Out[110])) - (si13 * (rsq * Out[43] + 2*z * Out[40]));
   Out[164] = si12 * z * Out[114] - (si13 * (rsq * Out[44] + 4*x * Out[41] + 2 * Out[40]));
   Out[165] = si12 * z * Out[115] - (si13 * (rsq * Out[45] + 4*y * Out[42] + 2 * Out[40]));
   Out[166] = (si12 * (z * Out[116] + 2 * Out[113])) - (si13 * (rsq * Out[46] + 4*z * Out[43] + 2 * Out[40]));
   Out[167] = -(si13 * (2*x * Out[42] + 2*y * Out[41]));
   Out[168] = (si12 * (z * Out[118] + Out[111])) - (si13 * (2*x * Out[43] + 2*z * Out[41]));
   Out[169] = (si12 * (z * Out[119] + Out[112])) - (si13 * (2*y * Out[43] + 2*z * Out[42]));
   Out[170] = si14 * z * Out[130] - si15 * rsq * Out[50];
   Out[171] = si14 * z * Out[131] - (si15 * (rsq * Out[51] + 2*x * Out[50]));
   Out[172] = si14 * z * Out[132] - (si15 * (rsq * Out[52] + 2*y * Out[50]));
   Out[173] = (si14 * (z * Out[133] + Out[130])) - si15 * 2*z * Out[50];
   Out[174] = -(si15 * (4*x * Out[51] + 2 * Out[50]));
   Out[175] = -(si15 * (4*y * Out[52] + 2 * Out[50]));
   Out[176] = si14 * 2 * Out[133] - si15 * 2 * Out[50];
   Out[177] = si14 * z * Out[137] - (si15 * (rsq * Out[57] + 2*x * Out[52] + 2*y * Out[51]));
   Out[178] = (si14 * (z * Out[138] + Out[131])) - si15 * 2*z * Out[51];
   Out[179] = (si14 * (z * Out[139] + Out[132])) - si15 * 2*z * Out[52];
   Out[180] = si16 * z * Out[90] - si17 * rsq * Out[60];
   Out[181] = si16 * z * Out[91] - (si17 * (rsq * Out[61] + 2*x * Out[60]));
   Out[182] = si16 * z * Out[92] - si17 * 2*y * Out[60];
   Out[183] = (si16 * (z * Out[93] + Out[90])) - (si17 * (rsq * Out[63] + 2*z * Out[60]));
   Out[184] = si16 * z * Out[94] - (si17 * (4*x * Out[61] + 2 * Out[60]));
   Out[185] = si16 * z * Out[95] - si17 * 2 * Out[60];
   Out[186] = (si16 * (z * Out[96] + 2 * Out[93])) - (si17 * (4*z * Out[63] + 2 * Out[60]));
   Out[187] = si16 * z * Out[97] - si17 * 2*y * Out[61];
   Out[188] = (si16 * (z * Out[98] + Out[91])) - (si17 * (rsq * Out[68] + 2*x * Out[63] + 2*z * Out[61]));
   Out[189] = si16 * Out[92] - si17 * 2*y * Out[63];
   Out[190] = si11 * (x * Out[120] - y * Out[140]);
   Out[191] = si11 * ((x * Out[121] + Out[120]) - y * Out[141]);
   Out[192] = si11 * (x * Out[122] - (y * Out[142] + Out[140]));
   Out[193] = 0.;
   Out[194] = si11 * ((x * Out[124] + 2 * Out[121]) - y * Out[144]);
   Out[195] = si11 * (x * Out[125] - (y * Out[145] + 2 * Out[142]));
   Out[196] = 0.;
   Out[197] = si11 * ((x * Out[127] + Out[122]) - (y * Out[147] + Out[141]));
   Out[198] = 0.;
   Out[199] = 0.;
   Out[200] = si16 * z * Out[100] - si17 * rsq * Out[80];
   Out[201] = si16 * z * Out[101] - si17 * 2*x * Out[80];
   Out[202] = si16 * z * Out[102] - (si17 * (rsq * Out[82] + 2*y * Out[80]));
   Out[203] = (si16 * (z * Out[103] + Out[100])) - (si17 * (rsq * Out[83] + 2*z * Out[80]));
   Out[204] = si16 * z * Out[104] - si17 * 2 * Out[80];
   Out[205] = si16 * z * Out[105] - (si17 * (4*y * Out[82] + 2 * Out[80]));
   Out[206] = (si16 * (z * Out[106] + 2 * Out[103])) - (si17 * (4*z * Out[83] + 2 * Out[80]));
   Out[207] = si16 * z * Out[107] - si17 * 2*x * Out[82];
   Out[208] = si16 * Out[101] - si17 * 2*x * Out[83];
   Out[209] = (si16 * (z * Out[109] + Out[102])) - (si17 * (rsq * Out[89] + 2*y * Out[83] + 2*z * Out[82]));
   Out[210] = si14 * z * Out[150] - si15 * rsq * Out[70];
   Out[211] = si14 * z * Out[151] - (si15 * (rsq * Out[71] + 2*x * Out[70]));
   Out[212] = si14 * z * Out[152] - (si15 * (rsq * Out[72] + 2*y * Out[70]));
   Out[213] = (si14 * (z * Out[153] + Out[150])) - si15 * 2*z * Out[70];
   Out[214] = si14 * z * Out[154] - (si15 * (rsq * Out[74] + 4*x * Out[71] + 2 * Out[70]));
   Out[215] = si14 * z * Out[155] - (si15 * (rsq * Out[75] + 4*y * Out[72] + 2 * Out[70]));
   Out[216] = si14 * 2 * Out[153] - si15 * 2 * Out[70];
   Out[217] = -(si15 * (2*x * Out[72] + 2*y * Out[71]));
   Out[218] = (si14 * (z * Out[158] + Out[151])) - si15 * 2*z * Out[71];
   Out[219] = (si14 * (z * Out[159] + Out[152])) - si15 * 2*z * Out[72];
   Out[220] = si11 * (y * Out[120] + x * Out[140]);
   Out[221] = si11 * (y * Out[121] + (x * Out[141] + Out[140]));
   Out[222] = si11 * ((y * Out[122] + Out[120]) + x * Out[142]);
   Out[223] = 0.;
   Out[224] = si11 * (y * Out[124] + (x * Out[144] + 2 * Out[141]));
   Out[225] = si11 * ((y * Out[125] + 2 * Out[122]) + x * Out[145]);
   Out[226] = 0.;
   Out[227] = si11 * ((y * Out[127] + Out[121]) + (x * Out[147] + Out[142]));
   Out[228] = 0.;
   Out[229] = 0.;
   Out[230] = si18 * z * Out[120];
   Out[231] = si18 * z * Out[121];
   Out[232] = si18 * z * Out[122];
   Out[233] = si18 * Out[120];
   Out[234] = si18 * z * Out[124];
   Out[235] = si18 * z * Out[125];
   Out[236] = 0.;
   Out[237] = si18 * z * Out[127];
   Out[238] = si18 * Out[121];
   Out[239] = si18 * Out[122];
   Out[240] = si18 * z * Out[140];
   Out[241] = si18 * z * Out[141];
   Out[242] = si18 * z * Out[142];
   Out[243] = si18 * Out[140];
   Out[244] = si18 * z * Out[144];
   Out[245] = si18 * z * Out[145];
   Out[246] = 0.;
   Out[247] = si18 * z * Out[147];
   Out[248] = si18 * Out[141];
   Out[249] = si18 * Out[142];
   if ( 4 == L ) return;

   // [1, d/dx, d/dy, d/dz, d^2/d(x x), d^2/d(y y), d^2/d(z z), d^2/d(x y), d^2/d(x z), d^2/d(y z)], S(l=5, m=-5..5)
   Out[250] = si20 * z * Out[180] - si21 * rsq * Out[90];
   Out[251] = si20 * z * Out[181] - (si21 * (rsq * Out[91] + 2*x * Out[90]));
   Out[252] = si20 * z * Out[182] - (si21 * (rsq * Out[92] + 2*y * Out[90]));
   Out[253] = (si20 * (z * Out[183] + Out[180])) - (si21 * (rsq * Out[93] + 2*z * Out[90]));
   Out[254] = si20 * z * Out[184] - (si21 * (rsq * Out[94] + 4*x * Out[91] + 2 * Out[90]));
   Out[255] = si20 * z * Out[185] - (si21 * (rsq * Out[95] + 4*y * Out[92] + 2 * Out[90]));
   Out[256] = (si20 * (z * Out[186] + 2 * Out[183])) - (si21 * (rsq * Out[96] + 4*z * Out[93] + 2 * Out[90]));
   Out[257] = si20 * z * Out[187] - (si21 * (rsq * Out[97] + 2*x * Out[92] + 2*y * Out[91]));
   Out[258] = (si20 * (z * Out[188] + Out[181])) - (si21 * (rsq * Out[98] + 2*x * Out[93] + 2*z * Out[91]));
   Out[259] = (si20 * (z * Out[189] + Out[182])) - (si21 * (2*y * Out[93] + 2*z * Out[92]));
   Out[260] = si20 * z * Out[200] - si21 * rsq * Out[100];
   Out[261] = si20 * z * Out[201] - (si21 * (rsq * Out[101] + 2*x * Out[100]));
   Out[262] = si20 * z * Out[202] - (si21 * (rsq * Out[102] + 2*y * Out[100]));
   Out[263] = (si20 * (z * Out[203] + Out[200])) - (si21 * (rsq * Out[103] + 2*z * Out[100]));
   Out[264] = si20 * z * Out[204] - (si21 * (rsq * Out[104] + 4*x * Out[101] + 2 * Out[100]));
   Out[265] = si20 * z * Out[205] - (si21 * (rsq * Out[105] + 4*y * Out[102] + 2 * Out[100]));
   Out[266] = (si20 * (z * Out[206] + 2 * Out[203])) - (si21 * (rsq * Out[106] + 4*z * Out[103] + 2 * Out[100]));
   Out[267] = si20 * z * Out[207] - (si21 * (rsq * Out[107] + 2*x * Out[102] + 2*y * Out[101]));
   Out[268] = (si20 * (z * Out[208] + Out[201])) - (si21 * (2*x * Out[103] + 2*z * Out[101]));
   Out[269] = (si20 * (z * Out[209] + Out[202])) - (si21 * (rsq * Out[109] + 2*y * Out[103] + 2*z * Out[102]));
   Out[270] = si22 * z * Out[210] - si23 * rsq * Out[150];
   Out[271] = si22 * z * Out[211] - (si23 * (rsq * Out[151] + 2*x * Out[150]));
   Out[272] = si22 * z * Out[212] - (si23 * (rsq * Out[152] + 2*y * Out[150]));
   Out[273] = (si22 * (z * Out[213] + Out[210])) - (si23 * (rsq * Out[153] + 2*z * Out[150]));
   Out[274] = si22 * z * Out[214] - (si23 * (rsq * Out[154] + 4*x * Out[151] + 2 * Out[150]));
   Out[275] = si22 * z * Out[215] - (si23 * (rsq * Out[155] + 4*y * Out[152] + 2 * Out[150]));
   Out[276] = (si22 * (z * Out[216] + 2 * Out[213])) - (si23 * (4*z * Out[153] + 2 * Out[150]));
   Out[277] = si22 * z * Out[217] - (si23 * (2*x * Out[152] + 2*y * Out[151]));
   Out[278] = (si22 * (z * Out[218] + Out[211])) - (si23 * (rsq * Out[158] + 2*x * Out[153] + 2*z * Out[151]));
   Out[279] = (si22 * (z * Out[219] + Out[212])) - (si23 * (rsq * Out[159] + 2*y * Out[153] + 2*z * Out[152]));
   Out[280] = si24 * z * Out[230] - si25 * rsq * Out[120];
   Out[281] = si24 * z * Out[231] - (si25 * (rsq * Out[121] + 2*x * Out[120]));
   Out[282] = si24 * z * Out[232] - (si25 * (rsq * Out[122] + 2*y * Out[120]));
   Out[283] = (si24 * (z * Out[233] + Out[230])) - si25 * 2*z * Out[120];
   Out[284] = si24 * z * Out[234] - (si25 * (rsq * Out[124] + 4*x * Out[121] + 2 * Out[120]));
   Out[285] = si24 * z * Out[235] - (si25 * (rsq * Out[125] + 4*y * Out[122] + 2 * Out[120]));
   Out[286] = si24 * 2 * Out[233] - si25 * 2 * Out[120];
   Out[287] = si24 * z * Out[237] - (si25 * (rsq * Out[127] + 2*x * Out[122] + 2*y * Out[121]));
   Out[288] = (si24 * (z * Out[238] + Out[231])) - si25 * 2*z * Out[121];
   Out[289] = (si24 * (z * Out[239] + Out[232])) - si25 * 2*z * Out[122];
   Out[290] = si26 * z * Out[220];
   Out[291] = si26 * z * Out[221];
   Out[292] = si26 * z * Out[222];
   Out[293] = si26 * Out[220];
   Out[294] = si26 * z * Out[224];
   Out[295] = si26 * z * Out[225];
   Out[296] = 0.;
   Out[297] = si26 * z * Out[227];
   Out[298] = si26 * Out[221];
   Out[299] = si26 * Out[222];
   Out[300] = si24 * z * Out[240] - si25 * rsq * Out[140];
   Out[301] = si24 * z * Out[241] - (si25 * (rsq * Out[141] + 2*x * Out[140]));
   Out[302] = si24 * z * Out[242] - (si25 * (rsq * Out[142] + 2*y * Out[140]));
   Out[303] = (si24 * (z * Out[243] + Out[240])) - si25 * 2*z * Out[140];
   Out[304] = si24 * z * Out[244] - (si25 * (rsq * Out[144] + 4*x * Out[141] + 2 * Out[140]));
   Out[305] = si24 * z * Out[245] - (si25 * (rsq * Out[145] + 4*y * Out[142] + 2 * Out[140]));
   Out[306] = si24 * 2 * Out[243] - si25 * 2 * Out[140];
   Out[307] = si24 * z * Out[247] - (si25 * (rsq * Out[147] + 2*x * Out[142] + 2*y * Out[141]));
   Out[308] = (si24 * (z * Out[248] + Out[241])) - si25 * 2*z * Out[141];
   Out[309] = (si24 * (z * Out[249] + Out[242])) - si25 * 2*z * Out[142];
   Out[310] = si26 * z * Out[190];
   Out[311] = si26 * z * Out[191];
   Out[312] = si26 * z * Out[192];
   Out[313] = si26 * Out[190];
   Out[314] = si26 * z * Out[194];
   Out[315] = si26 * z * Out[195];
   Out[316] = 0.;
   Out[317] = si26 * z * Out[197];
   Out[318] = si26 * Out[191];
   Out[319] = si26 * Out[192];
   Out[320] = si19 * (y * Out[190] + x * Out[220]);
   Out[321] = si19 * (y * Out[191] + (x * Out[221] + Out[220]));
   Out[322] = si19 * ((y * Out[192] + Out[190]) + x * Out[222]);
   Out[323] = 0.;
   Out[324] = si19 * (y * Out[194] + (x * Out[224] + 2 * Out[221]));
   Out[325] = si19 * ((y * Out[195] + 2 * Out[192]) + x * Out[225]);
   Out[326] = 0.;
   Out[327] = si19 * ((y * Out[197] + Out[191]) + (x * Out[227] + Out[222]));
   Out[328] = 0.;
   Out[329] = 0.;
   Out[330] = si27 * z * Out[160] - si28 * rsq * Out[110];
   Out[331] = si27 * z * Out[161] - (si28 * (rsq * Out[111] + 2*x * Out[110]));
   Out[332] = si27 * z * Out[162] - (si28 * (rsq * Out[112] + 2*y * Out[110]));
   Out[333] = (si27 * (z * Out[163] + Out[160])) - (si28 * (rsq * Out[113] + 2*z * Out[110]));
   Out[334] = si27 * z * Out[164] - (si28 * (rsq * Out[114] + 4*x * Out[111] + 2 * Out[110]));
   Out[335] = si27 * z * Out[165] - (si28 * (rsq * Out[115] + 4*y * Out[112] + 2 * Out[110]));
   Out[336] = (si27 * (z * Out[166] + 2 * Out[163])) - (si28 * (rsq * Out[116] + 4*z * Out[113] + 2 * Out[110]));
   Out[337] = si27 * z * Out[167] - (si28 * (2*x * Out[112] + 2*y * Out[111]));
   Out[338] = (si27 * (z * Out[168] + Out[161])) - (si28 * (rsq * Out[118] + 2*x * Out[113] + 2*z * Out[111]));
   Out[339] = (si27 * (z * Out[169] + Out[162])) - (si28 * (rsq * Out[119] + 2*y * Out[113] + 2*z * Out[112]));
   Out[340] = si19 * (x * Out[190] - y * Out[220]);
   Out[341] = si19 * ((x * Out[191] + Out[190]) - y * Out[221]);
   Out[342] = si19 * (x * Out[192] - (y * Out[222] + Out[220]));
   Out[343] = 0.;
   Out[344] = si19 * ((x * Out[194] + 2 * Out[191]) - y * Out[224]);
   Out[345] = si19 * (x * Out[195] - (y * Out[225] + 2 * Out[222]));
   Out[346] = 0.;
   Out[347] = si19 * ((x * Out[197] + Out[192]) - (y * Out[227] + Out[221]));
   Out[348] = 0.;
   Out[349] = 0.;
   Out[350] = si22 * z * Out[170] - si23 * rsq * Out[130];
   Out[351] = si22 * z * Out[171] - (si23 * (rsq * Out[131] + 2*x * Out[130]));
   Out[352] = si22 * z * Out[172] - (si23 * (rsq * Out[132] + 2*y * Out[130]));
   Out[353] = (si22 * (z * Out[173] + Out[170])) - (si23 * (rsq * Out[133] + 2*z * Out[130]));
   Out[354] = si22 * z * Out[174] - (si23 * (4*x * Out[131] + 2 * Out[130]));
   Out[355] = si22 * z * Out[175] - (si23 * (4*y * Out[132] + 2 * Out[130]));
   Out[356] = (si22 * (z * Out[176] + 2 * Out[173])) - (si23 * (4*z * Out[133] + 2 * Out[130]));
   Out[357] = si22 * z * Out[177] - (si23 * (rsq * Out[137] + 2*x * Out[132] + 2*y * Out[131]));
   Out[358] = (si22 * (z * Out[178] + Out[171])) - (si23 * (rsq * Out[138] + 2*x * Out[133] + 2*z * Out[131]));
   Out[359] = (si22 * (z * Out[179] + Out[172])) - (si23 * (rsq * Out[139] + 2*y * Out[133] + 2*z * Out[132]));
   if ( 5 == L ) return;

   // [1, d/dx, d/dy, d/dz, d^2/d(x x), d^2/d(y y), d^2/d(z z), d^2/d(x y), d^2/d(x z), d^2/d(y z)], S(l=6, m=-6..6)
   Out[360] = si29 * (x * Out[340] - y * Out[320]);
   Out[361] = si29 * ((x * Out[341] + Out[340]) - y * Out[321]);
   Out[362] = si29 * (x * Out[342] - (y * Out[322] + Out[320]));
   Out[363] = 0.;
   Out[364] = si29 * ((x * Out[344] + 2 * Out[341]) - y * Out[324]);
   Out[365] = si29 * (x * Out[345] - (y * Out[325] + 2 * Out[322]));
   Out[366] = 0.;
   Out[367] = si29 * ((x * Out[347] + Out[342]) - (y * Out[327] + Out[321]));
   Out[368] = 0.;
   Out[369] = 0.;
   Out[370] = si30 * z * Out[350] - si31 * rsq * Out[170];
   Out[371] = si30 * z * Out[351] - (si31 * (rsq * Out[171] + 2*x * Out[170]));
   Out[372] = si30 * z * Out[352] - (si31 * (rsq * Out[172] + 2*y * Out[170]));
   Out[373] = (si30 * (z * Out[353] + Out[350])) - (si31 * (rsq * Out[173] + 2*z * Out[170]));
   Out[374] = si30 * z * Out[354] - (si31 * (rsq * Out[174] + 4*x * Out[171] + 2 * Out[170]));
   Out[375] = si30 * z * Out[355] - (si31 * (rsq * Out[175] + 4*y * Out[172] + 2 * Out[170]));
   Out[376] = (si30 * (z * Out[356] + 2 * Out[353])) - (si31 * (rsq * Out[176] + 4*z * Out[173] + 2 * Out[170]));
   Out[377] = si30 * z * Out[357] - (si31 * (rsq * Out[177] + 2*x * Out[172] + 2*y * Out[171]));
   Out[378] = (si30 * (z * Out[358] + Out[351])) - (si31 * (rsq * Out[178] + 2*x * Out[173] + 2*z * Out[171]));
   Out[379] = (si30 * (z * Out[359] + Out[352])) - (si31 * (rsq * Out[179] + 2*y * Out[173] + 2*z * Out[172]));
   Out[380] = si32 * z * Out[340];
   Out[381] = si32 * z * Out[341];
   Out[382] = si32 * z * Out[342];
   Out[383] = si32 * Out[340];
   Out[384] = si32 * z * Out[344];
   Out[385] = si32 * z * Out[345];
   Out[386] = 0.;
   Out[387] = si32 * z * Out[347];
   Out[388] = si32 * Out[341];
   Out[389] = si32 * Out[342];
   Out[390] = si33 * z * Out[310] - si34 * rsq * Out[190];
   Out[391] = si33 * z * Out[311] - (si34 * (rsq * Out[191] + 2*x * Out[190]));
   Out[392] = si33 * z * Out[312] - (si34 * (rsq * Out[192] + 2*y * Out[190]));
   Out[393] = (si33 * (z * Out[313] + Out[310])) - si34 * 2*z * Out[190];
   Out[394] = si33 * z * Out[314] - (si34 * (rsq * Out[194] + 4*x * Out[191] + 2 * Out[190]));
   Out[395] = si33 * z * Out[315] - (si34 * (rsq * Out[195] + 4*y * Out[192] + 2 * Out[190]));
   Out[396] = si33 * 2 * Out[313] - si34 * 2 * Out[190];
   Out[397] = si33 * z * Out[317] - (si34 * (rsq * Out[197] + 2*x * Out[192] + 2*y * Out[191]));
   Out[398] = (si33 * (z * Out[318] + Out[311])) - si34 * 2*z * Out[191];
   Out[399] = (si33 * (z * Out[319] + Out[312])) - si34 * 2*z * Out[192];
   Out[400] = si32 * z * Out[320];
   Out[401] = si32 * z * Out[321];
   Out[402] = si32 * z * Out[322];
   Out[403] = si32 * Out[320];
   Out[404] = si32 * z * Out[324];
   Out[405] = si32 * z * Out[325];
   Out[406] = 0.;
   Out[407] = si32 * z * Out[327];
   Out[408] = si32 * Out[321];
   Out[409] = si32 * Out[322];
   Out[410] = si30 * z * Out[270] - si31 * rsq * Out[210];
   Out[411] = si30 * z * Out[271] - (si31 * (rsq * Out[211] + 2*x * Out[210]));
   Out[412] = si30 * z * Out[272] - (si31 * (rsq * Out[212] + 2*y * Out[210]));
   Out[413] = (si30 * (z * Out[273] + Out[270])) - (si31 * (rsq * Out[213] + 2*z * Out[210]));
   Out[414] = si30 * z * Out[274] - (si31 * (rsq * Out[214] + 4*x * Out[211] + 2 * Out[210]));
   Out[415] = si30 * z * Out[275] - (si31 * (rsq * Out[215] + 4*y * Out[212] + 2 * Out[210]));
   Out[416] = (si30 * (z * Out[276] + 2 * Out[273])) - (si31 * (rsq * Out[216] + 4*z * Out[213] + 2 * Out[210]));
   Out[417] = si30 * z * Out[277] - (si31 * (rsq * Out[217] + 2*x * Out[212] + 2*y * Out[211]));
   Out[418] = (si30 * (z * Out[278] + Out[271])) - (si31 * (rsq * Out[218] + 2*x * Out[213] + 2*z * Out[211]));
   Out[419] = (si30 * (z * Out[279] + Out[272])) - (si31 * (rsq * Out[219] + 2*y * Out[213] + 2*z * Out[212]));
   Out[420] = si29 * (y * Out[340] + x * Out[320]);
   Out[421] = si29 * (y * Out[341] + (x * Out[321] + Out[320]));
   Out[422] = si29 * ((y * Out[342] + Out[340]) + x * Out[322]);
   Out[423] = 0.;
   Out[424] = si29 * (y * Out[344] + (x * Out[324] + 2 * Out[321]));
   Out[425] = si29 * ((y * Out[345] + 2 * Out[342]) + x * Out[325]);
   Out[426] = 0.;
   Out[427] = si29 * ((y * Out[347] + Out[341]) + (x * Out[327] + Out[322]));
   Out[428] = 0.;
   Out[429] = 0.;
   Out[430] = si35 * z * Out[280] - si36 * rsq * Out[230];
   Out[431] = si35 * z * Out[281] - (si36 * (rsq * Out[231] + 2*x * Out[230]));
   Out[432] = si35 * z * Out[282] - (si36 * (rsq * Out[232] + 2*y * Out[230]));
   Out[433] = (si35 * (z * Out[283] + Out[280])) - (si36 * (rsq * Out[233] + 2*z * Out[230]));
   Out[434] = si35 * z * Out[284] - (si36 * (rsq * Out[234] + 4*x * Out[231] + 2 * Out[230]));
   Out[435] = si35 * z * Out[285] - (si36 * (rsq * Out[235] + 4*y * Out[232] + 2 * Out[230]));
   Out[436] = (si35 * (z * Out[286] + 2 * Out[283])) - (si36 * (4*z * Out[233] + 2 * Out[230]));
   Out[437] = si35 * z * Out[287] - (si36 * (rsq * Out[237] + 2*x * Out[232] + 2*y * Out[231]));
   Out[438] = (si35 * (z * Out[288] + Out[281])) - (si36 * (rsq * Out[238] + 2*x * Out[233] + 2*z * Out[231]));
   Out[439] = (si35 * (z * Out[289] + Out[282])) - (si36 * (rsq * Out[239] + 2*y * Out[233] + 2*z * Out[232]));
   Out[440] = si33 * z * Out[290] - si34 * rsq * Out[220];
   Out[441] = si33 * z * Out[291] - (si34 * (rsq * Out[221] + 2*x * Out[220]));
   Out[442] = si33 * z * Out[292] - (si34 * (rsq * Out[222] + 2*y * Out[220]));
   Out[443] = (si33 * (z * Out[293] + Out[290])) - si34 * 2*z * Out[220];
   Out[444] = si33 * z * Out[294] - (si34 * (rsq * Out[224] + 4*x * Out[221] + 2 * Out[220]));
   Out[445] = si33 * z * Out[295] - (si34 * (rsq * Out[225] + 4*y * Out[222] + 2 * Out[220]));
   Out[446] = si33 * 2 * Out[293] - si34 * 2 * Out[220];
   Out[447] = si33 * z * Out[297] - (si34 * (rsq * Out[227] + 2*x * Out[222] + 2*y * Out[221]));
   Out[448] = (si33 * (z * Out[298] + Out[291])) - si34 * 2*z * Out[221];
   Out[449] = (si33 * (z * Out[299] + Out[292])) - si34 * 2*z * Out[222];
   Out[450] = si37 * z * Out[330] - si38 * rsq * Out[160];
   Out[451] = si37 * z * Out[331] - (si38 * (rsq * Out[161] + 2*x * Out[160]));
   Out[452] = si37 * z * Out[332] - (si38 * (rsq * Out[162] + 2*y * Out[160]));
   Out[453] = (si37 * (z * Out[333] + Out[330])) - (si38 * (rsq * Out[163] + 2*z * Out[160]));
   Out[454] = si37 * z * Out[334] - (si38 * (rsq * Out[164] + 4*x * Out[161] + 2 * Out[160]));
   Out[455] = si37 * z * Out[335] - (si38 * (rsq * Out[165] + 4*y * Out[162] + 2 * Out[160]));
   Out[456] = (si37 * (z * Out[336] + 2 * Out[333])) - (si38 * (rsq * Out[166] + 4*z * Out[163] + 2 * Out[160]));
   Out[457] = si37 * z * Out[337] - (si38 * (rsq * Out[167] + 2*x * Out[162] + 2*y * Out[161]));
   Out[458] = (si37 * (z * Out[338] + Out[331])) - (si38 * (rsq * Out[168] + 2*x * Out[163] + 2*z * Out[161]));
   Out[459] = (si37 * (z * Out[339] + Out[332])) - (si38 * (rsq * Out[169] + 2*y * Out[163] + 2*z * Out[162]));
   Out[460] = si35 * z * Out[300] - si36 * rsq * Out[240];
   Out[461] = si35 * z * Out[301] - (si36 * (rsq * Out[241] + 2*x * Out[240]));
   Out[462] = si35 * z * Out[302] - (si36 * (rsq * Out[242] + 2*y * Out[240]));
   Out[463] = (si35 * (z * Out[303] + Out[300])) - (si36 * (rsq * Out[243] + 2*z * Out[240]));
   Out[464] = si35 * z * Out[304] - (si36 * (rsq * Out[244] + 4*x * Out[241] + 2 * Out[240]));
   Out[465] = si35 * z * Out[305] - (si36 * (rsq * Out[245] + 4*y * Out[242] + 2 * Out[240]));
   Out[466] = (si35 * (z * Out[306] + 2 * Out[303])) - (si36 * (4*z * Out[243] + 2 * Out[240]));
   Out[467] = si35 * z * Out[307] - (si36 * (rsq * Out[247] + 2*x * Out[242] + 2*y * Out[241]));
   Out[468] = (si35 * (z * Out[308] + Out[301])) - (si36 * (rsq * Out[248] + 2*x * Out[243] + 2*z * Out[241]));
   Out[469] = (si35 * (z * Out[309] + Out[302])) - (si36 * (rsq * Out[249] + 2*y * Out[243] + 2*z * Out[242]));
   Out[470] = si39 * z * Out[260] - si40 * rsq * Out[200];
   Out[471] = si39 * z * Out[261] - (si40 * (rsq * Out[201] + 2*x * Out[200]));
   Out[472] = si39 * z * Out[262] - (si40 * (rsq * Out[202] + 2*y * Out[200]));
   Out[473] = (si39 * (z * Out[263] + Out[260])) - (si40 * (rsq * Out[203] + 2*z * Out[200]));
   Out[474] = si39 * z * Out[264] - (si40 * (rsq * Out[204] + 4*x * Out[201] + 2 * Out[200]));
   Out[475] = si39 * z * Out[265] - (si40 * (rsq * Out[205] + 4*y * Out[202] + 2 * Out[200]));
   Out[476] = (si39 * (z * Out[266] + 2 * Out[263])) - (si40 * (rsq * Out[206] + 4*z * Out[203] + 2 * Out[200]));
   Out[477] = si39 * z * Out[267] - (si40 * (rsq * Out[207] + 2*x * Out[202] + 2*y * Out[201]));
   Out[478] = (si39 * (z * Out[268] + Out[261])) - (si40 * (rsq * Out[208] + 2*x * Out[203] + 2*z * Out[201]));
   Out[479] = (si39 * (z * Out[269] + Out[262])) - (si40 * (rsq * Out[209] + 2*y * Out[203] + 2*z * Out[202]));
   Out[480] = si39 * z * Out[250] - si40 * rsq * Out[180];
   Out[481] = si39 * z * Out[251] - (si40 * (rsq * Out[181] + 2*x * Out[180]));
   Out[482] = si39 * z * Out[252] - (si40 * (rsq * Out[182] + 2*y * Out[180]));
   Out[483] = (si39 * (z * Out[253] + Out[250])) - (si40 * (rsq * Out[183] + 2*z * Out[180]));
   Out[484] = si39 * z * Out[254] - (si40 * (rsq * Out[184] + 4*x * Out[181] + 2 * Out[180]));
   Out[485] = si39 * z * Out[255] - (si40 * (rsq * Out[185] + 4*y * Out[182] + 2 * Out[180]));
   Out[486] = (si39 * (z * Out[256] + 2 * Out[253])) - (si40 * (rsq * Out[186] + 4*z * Out[183] + 2 * Out[180]));
   Out[487] = si39 * z * Out[257] - (si40 * (rsq * Out[187] + 2*x * Out[182] + 2*y * Out[181]));
   Out[488] = (si39 * (z * Out[258] + Out[251])) - (si40 * (rsq * Out[188] + 2*x * Out[183] + 2*z * Out[181]));
   Out[489] = (si39 * (z * Out[259] + Out[252])) - (si40 * (rsq * Out[189] + 2*y * Out[183] + 2*z * Out[182]));
   if ( 6 == L ) return;

   AicAssertFail("L too large.", __FILE__, __LINE__);
}

// transform one cartesian coordinate to solid harmonics.
void ShTrN( double *AIC_RP pOut, unsigned ns, double const *AIC_RP pIn, unsigned dummy, unsigned L, unsigned Count )
{
   assert(ns == dummy);
   switch(L) {
      case 0: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            pOut[0*ns+i] = pIn[0*ns+i];
         }
         break;
      }
      case 1: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            pOut[0*ns+i] = pIn[0*ns+i];
            pOut[1*ns+i] = pIn[1*ns+i];
            pOut[2*ns+i] = pIn[2*ns+i];
         }
         break;
      }
      case 2: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[0*ns+i]; double d1 = pIn[1*ns+i];
            pOut[0*ns+i] = -sd0*d0 - sd0*d1 + pIn[2*ns+i];
            pOut[1*ns+i] = sd1*pIn[3*ns+i];
            pOut[2*ns+i] = sd1*pIn[4*ns+i];
            pOut[3*ns+i] = sd2*d0 - sd2*d1;
            pOut[4*ns+i] = sd1*pIn[5*ns+i];
         }
         break;
      }
      case 3: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[0*ns+i]; double d1 = pIn[1*ns+i]; double d3 = pIn[3*ns+i]; double d4 = pIn[4*ns+i]; double d6 = pIn[6*ns+i]; double d7 = pIn[7*ns+i];
            pOut[0*ns+i] = -sd3*d0 - sd3*d1 + sd4*pIn[2*ns+i];
            pOut[1*ns+i] = -sd3*d3 - sd3*d4 + sd4*pIn[5*ns+i];
            pOut[2*ns+i] = -sd5*d6 - sd5*d7 + pIn[8*ns+i];
            pOut[3*ns+i] = sd7*d0 - sd6*d1;
            pOut[4*ns+i] = sd8*pIn[9*ns+i];
            pOut[5*ns+i] = sd6*d3 - sd7*d4;
            pOut[6*ns+i] = sd9*d6 - sd9*d7;
         }
         break;
      }
      case 4: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[0*ns+i]; double d1 = pIn[1*ns+i]; double d2 = pIn[2*ns+i]; double d3 = pIn[3*ns+i]; double d4 = pIn[4*ns+i]; double d6 = pIn[6*ns+i]; double d7 = pIn[7*ns+i]; double d9 = pIn[9*ns+i]; double d10 = pIn[10*ns+i]; double d12 = pIn[12*ns+i]; double d13 = pIn[13*ns+i];
            pOut[0*ns+i] = sd11*d0 + sd12*d1 - sd10*d2 + sd11*d3 - sd10*d4 + pIn[5*ns+i];
            pOut[1*ns+i] = -sd13*d6 - sd13*d7 + sd14*pIn[8*ns+i];
            pOut[2*ns+i] = -sd6*d9 - sd6*d10 + sd15*pIn[11*ns+i];
            pOut[3*ns+i] = sd16*d0 - sd17*d1 + sd16*d3;
            pOut[4*ns+i] = -sd6*d12 - sd6*d13 + sd15*pIn[14*ns+i];
            pOut[5*ns+i] = -sd19*d0 + sd18*d2 + sd19*d3 - sd18*d4;
            pOut[6*ns+i] = sd20*d6 - sd20*d7;
            pOut[7*ns+i] = sd22*d9 - sd21*d10;
            pOut[8*ns+i] = sd21*d12 - sd22*d13;
         }
         break;
      }
      case 5: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[0*ns+i]; double d1 = pIn[1*ns+i]; double d2 = pIn[2*ns+i]; double d3 = pIn[3*ns+i]; double d4 = pIn[4*ns+i]; double d6 = pIn[6*ns+i]; double d7 = pIn[7*ns+i]; double d8 = pIn[8*ns+i]; double d9 = pIn[9*ns+i]; double d10 = pIn[10*ns+i]; double d12 = pIn[12*ns+i]; double d13 = pIn[13*ns+i]; double d14 = pIn[14*ns+i]; double d15 = pIn[15*ns+i]; double d16 = pIn[16*ns+i]; double d18 = pIn[18*ns+i]; double d19 = pIn[19*ns+i];
            pOut[0*ns+i] = sd25*d0 + sd24*d1 - sd23*d2 + sd25*d3 - sd23*d4 + sd8*pIn[5*ns+i];
            pOut[1*ns+i] = sd25*d6 + sd24*d7 - sd23*d8 + sd25*d9 - sd23*d10 + sd8*pIn[11*ns+i];
            pOut[2*ns+i] = -sd26*d12 + sd27*d14 + sd26*d15 - sd27*d16;
            pOut[3*ns+i] = -sd32*d0 + sd30*d1 + sd29*d2 + sd31*d3 - sd28*d4;
            pOut[4*ns+i] = sd33*d18 - sd33*d19;
            pOut[5*ns+i] = -sd31*d6 - sd30*d7 + sd28*d8 + sd32*d9 - sd29*d10;
            pOut[6*ns+i] = sd34*d12 - sd35*d13 + sd34*d15;
            pOut[7*ns+i] = sd37*d6 - sd38*d7 + sd36*d9;
            pOut[8*ns+i] = sd39*d12 + sd41*d13 - sd40*d14 + sd39*d15 - sd40*d16 + pIn[17*ns+i];
            pOut[9*ns+i] = sd36*d0 - sd38*d1 + sd37*d3;
            pOut[10*ns+i] = -sd27*d18 - sd27*d19 + sd42*pIn[20*ns+i];
         }
         break;
      }
      case 6: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[0*ns+i]; double d1 = pIn[1*ns+i]; double d2 = pIn[2*ns+i]; double d3 = pIn[3*ns+i]; double d4 = pIn[4*ns+i]; double d5 = pIn[5*ns+i]; double d6 = pIn[6*ns+i]; double d7 = pIn[7*ns+i]; double d8 = pIn[8*ns+i]; double d10 = pIn[10*ns+i]; double d11 = pIn[11*ns+i]; double d12 = pIn[12*ns+i]; double d13 = pIn[13*ns+i]; double d14 = pIn[14*ns+i]; double d16 = pIn[16*ns+i]; double d17 = pIn[17*ns+i]; double d18 = pIn[18*ns+i]; double d19 = pIn[19*ns+i]; double d20 = pIn[20*ns+i]; double d22 = pIn[22*ns+i]; double d23 = pIn[23*ns+i]; double d24 = pIn[24*ns+i]; double d25 = pIn[25*ns+i]; double d26 = pIn[26*ns+i];
            pOut[0*ns+i] = sd43*d0 - sd44*d1 + sd44*d3 - sd43*d6;
            pOut[1*ns+i] = sd46*d10 + sd47*d11 - sd45*d12 + sd46*d13 - sd45*d14 + sd45*pIn[15*ns+i];
            pOut[2*ns+i] = sd49*d16 - sd50*d17 + sd48*d19;
            pOut[3*ns+i] = -sd54*d0 + sd53*d1 + sd52*d2 + sd53*d3 - sd51*d4 - sd54*d6 + sd52*d7;
            pOut[4*ns+i] = sd48*d22 - sd50*d23 + sd49*d25;
            pOut[5*ns+i] = sd56*d0 + sd56*d1 - sd55*d2 - sd56*d3 + sd55*d5 - sd56*d6 + sd55*d7 - sd55*d8;
            pOut[6*ns+i] = sd58*d10 - sd57*d11 + sd58*d13;
            pOut[7*ns+i] = -sd61*d16 + sd60*d17 + sd55*d18 + sd59*d19 - sd62*d20;
            pOut[8*ns+i] = -sd63*d10 + sd64*d12 + sd63*d13 - sd64*d14;
            pOut[9*ns+i] = -sd69*d0 - sd68*d1 + sd66*d2 - sd68*d3 + sd65*d4 - sd67*d5 - sd69*d6 + sd66*d7 - sd67*d8 + pIn[9*ns+i];
            pOut[10*ns+i] = -sd59*d22 - sd60*d23 + sd62*d24 + sd61*d25 - sd55*d26;
            pOut[11*ns+i] = sd71*d22 + sd73*d23 - sd70*d24 + sd71*d25 - sd70*d26 + sd72*pIn[27*ns+i];
            pOut[12*ns+i] = sd71*d16 + sd73*d17 - sd70*d18 + sd71*d19 - sd70*d20 + sd72*pIn[21*ns+i];
         }
         break;
      }
   }
}

// transform one cartesian coordinate to solid harmonics. input indices explicitly given.
void ShTrN_Indirect( double *AIC_RP pOut, unsigned ns, double const *AIC_RP pIn, unsigned short const *AIC_RP ii, unsigned L, unsigned Count )
{
   switch(L) {
      case 0: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[ii[1*i+0]];
            pOut[0*ns+i] = d0;
         }
         break;
      }
      case 1: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[ii[3*i+0]]; double d1 = pIn[ii[3*i+1]]; double d2 = pIn[ii[3*i+2]];
            pOut[0*ns+i] = d0;
            pOut[1*ns+i] = d1;
            pOut[2*ns+i] = d2;
         }
         break;
      }
      case 2: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[ii[6*i+0]]; double d1 = pIn[ii[6*i+1]]; double d2 = pIn[ii[6*i+2]]; double d3 = pIn[ii[6*i+3]]; double d4 = pIn[ii[6*i+4]]; double d5 = pIn[ii[6*i+5]];
            pOut[0*ns+i] = -sd0*d0 - sd0*d1 + d2;
            pOut[1*ns+i] = sd1*d3;
            pOut[2*ns+i] = sd1*d4;
            pOut[3*ns+i] = sd2*d0 - sd2*d1;
            pOut[4*ns+i] = sd1*d5;
         }
         break;
      }
      case 3: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[ii[10*i+0]]; double d1 = pIn[ii[10*i+1]]; double d2 = pIn[ii[10*i+2]]; double d3 = pIn[ii[10*i+3]]; double d4 = pIn[ii[10*i+4]]; double d5 = pIn[ii[10*i+5]]; double d6 = pIn[ii[10*i+6]]; double d7 = pIn[ii[10*i+7]]; double d8 = pIn[ii[10*i+8]]; double d9 = pIn[ii[10*i+9]];
            pOut[0*ns+i] = -sd3*d0 - sd3*d1 + sd4*d2;
            pOut[1*ns+i] = -sd3*d3 - sd3*d4 + sd4*d5;
            pOut[2*ns+i] = -sd5*d6 - sd5*d7 + d8;
            pOut[3*ns+i] = sd7*d0 - sd6*d1;
            pOut[4*ns+i] = sd8*d9;
            pOut[5*ns+i] = sd6*d3 - sd7*d4;
            pOut[6*ns+i] = sd9*d6 - sd9*d7;
         }
         break;
      }
      case 4: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[ii[15*i+0]]; double d1 = pIn[ii[15*i+1]]; double d2 = pIn[ii[15*i+2]]; double d3 = pIn[ii[15*i+3]]; double d4 = pIn[ii[15*i+4]]; double d5 = pIn[ii[15*i+5]]; double d6 = pIn[ii[15*i+6]]; double d7 = pIn[ii[15*i+7]]; double d8 = pIn[ii[15*i+8]]; double d9 = pIn[ii[15*i+9]]; double d10 = pIn[ii[15*i+10]]; double d11 = pIn[ii[15*i+11]]; double d12 = pIn[ii[15*i+12]]; double d13 = pIn[ii[15*i+13]]; double d14 = pIn[ii[15*i+14]];
            pOut[0*ns+i] = sd11*d0 + sd12*d1 - sd10*d2 + sd11*d3 - sd10*d4 + d5;
            pOut[1*ns+i] = -sd13*d6 - sd13*d7 + sd14*d8;
            pOut[2*ns+i] = -sd6*d9 - sd6*d10 + sd15*d11;
            pOut[3*ns+i] = sd16*d0 - sd17*d1 + sd16*d3;
            pOut[4*ns+i] = -sd6*d12 - sd6*d13 + sd15*d14;
            pOut[5*ns+i] = -sd19*d0 + sd18*d2 + sd19*d3 - sd18*d4;
            pOut[6*ns+i] = sd20*d6 - sd20*d7;
            pOut[7*ns+i] = sd22*d9 - sd21*d10;
            pOut[8*ns+i] = sd21*d12 - sd22*d13;
         }
         break;
      }
      case 5: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[ii[21*i+0]]; double d1 = pIn[ii[21*i+1]]; double d2 = pIn[ii[21*i+2]]; double d3 = pIn[ii[21*i+3]]; double d4 = pIn[ii[21*i+4]]; double d5 = pIn[ii[21*i+5]]; double d6 = pIn[ii[21*i+6]]; double d7 = pIn[ii[21*i+7]]; double d8 = pIn[ii[21*i+8]]; double d9 = pIn[ii[21*i+9]]; double d10 = pIn[ii[21*i+10]]; double d11 = pIn[ii[21*i+11]]; double d12 = pIn[ii[21*i+12]]; double d13 = pIn[ii[21*i+13]]; double d14 = pIn[ii[21*i+14]]; double d15 = pIn[ii[21*i+15]]; double d16 = pIn[ii[21*i+16]]; double d17 = pIn[ii[21*i+17]]; double d18 = pIn[ii[21*i+18]]; double d19 = pIn[ii[21*i+19]]; double d20 = pIn[ii[21*i+20]];
            pOut[0*ns+i] = sd25*d0 + sd24*d1 - sd23*d2 + sd25*d3 - sd23*d4 + sd8*d5;
            pOut[1*ns+i] = sd25*d6 + sd24*d7 - sd23*d8 + sd25*d9 - sd23*d10 + sd8*d11;
            pOut[2*ns+i] = -sd26*d12 + sd27*d14 + sd26*d15 - sd27*d16;
            pOut[3*ns+i] = -sd32*d0 + sd30*d1 + sd29*d2 + sd31*d3 - sd28*d4;
            pOut[4*ns+i] = sd33*d18 - sd33*d19;
            pOut[5*ns+i] = -sd31*d6 - sd30*d7 + sd28*d8 + sd32*d9 - sd29*d10;
            pOut[6*ns+i] = sd34*d12 - sd35*d13 + sd34*d15;
            pOut[7*ns+i] = sd37*d6 - sd38*d7 + sd36*d9;
            pOut[8*ns+i] = sd39*d12 + sd41*d13 - sd40*d14 + sd39*d15 - sd40*d16 + d17;
            pOut[9*ns+i] = sd36*d0 - sd38*d1 + sd37*d3;
            pOut[10*ns+i] = -sd27*d18 - sd27*d19 + sd42*d20;
         }
         break;
      }
      case 6: {
         for ( unsigned i = 0; i != Count; ++ i ) {
            double d0 = pIn[ii[28*i+0]]; double d1 = pIn[ii[28*i+1]]; double d2 = pIn[ii[28*i+2]]; double d3 = pIn[ii[28*i+3]]; double d4 = pIn[ii[28*i+4]]; double d5 = pIn[ii[28*i+5]]; double d6 = pIn[ii[28*i+6]]; double d7 = pIn[ii[28*i+7]]; double d8 = pIn[ii[28*i+8]]; double d9 = pIn[ii[28*i+9]]; double d10 = pIn[ii[28*i+10]]; double d11 = pIn[ii[28*i+11]]; double d12 = pIn[ii[28*i+12]]; double d13 = pIn[ii[28*i+13]]; double d14 = pIn[ii[28*i+14]]; double d15 = pIn[ii[28*i+15]]; double d16 = pIn[ii[28*i+16]]; double d17 = pIn[ii[28*i+17]]; double d18 = pIn[ii[28*i+18]]; double d19 = pIn[ii[28*i+19]]; double d20 = pIn[ii[28*i+20]]; double d21 = pIn[ii[28*i+21]]; double d22 = pIn[ii[28*i+22]]; double d23 = pIn[ii[28*i+23]]; double d24 = pIn[ii[28*i+24]]; double d25 = pIn[ii[28*i+25]]; double d26 = pIn[ii[28*i+26]]; double d27 = pIn[ii[28*i+27]];
            pOut[0*ns+i] = sd43*d0 - sd44*d1 + sd44*d3 - sd43*d6;
            pOut[1*ns+i] = sd46*d10 + sd47*d11 - sd45*d12 + sd46*d13 - sd45*d14 + sd45*d15;
            pOut[2*ns+i] = sd49*d16 - sd50*d17 + sd48*d19;
            pOut[3*ns+i] = -sd54*d0 + sd53*d1 + sd52*d2 + sd53*d3 - sd51*d4 - sd54*d6 + sd52*d7;
            pOut[4*ns+i] = sd48*d22 - sd50*d23 + sd49*d25;
            pOut[5*ns+i] = sd56*d0 + sd56*d1 - sd55*d2 - sd56*d3 + sd55*d5 - sd56*d6 + sd55*d7 - sd55*d8;
            pOut[6*ns+i] = sd58*d10 - sd57*d11 + sd58*d13;
            pOut[7*ns+i] = -sd61*d16 + sd60*d17 + sd55*d18 + sd59*d19 - sd62*d20;
            pOut[8*ns+i] = -sd63*d10 + sd64*d12 + sd63*d13 - sd64*d14;
            pOut[9*ns+i] = -sd69*d0 - sd68*d1 + sd66*d2 - sd68*d3 + sd65*d4 - sd67*d5 - sd69*d6 + sd66*d7 - sd67*d8 + d9;
            pOut[10*ns+i] = -sd59*d22 - sd60*d23 + sd62*d24 + sd61*d25 - sd55*d26;
            pOut[11*ns+i] = sd71*d22 + sd73*d23 - sd70*d24 + sd71*d25 - sd70*d26 + sd72*d27;
            pOut[12*ns+i] = sd71*d16 + sd73*d17 - sd70*d18 + sd71*d19 - sd70*d20 + sd72*d21;
         }
         break;
      }
   }
}

// transform one solid harmonic coordinate to cartesians.
void CaTrN1( double *AIC_RP pOut, double const *AIC_RP pIn, unsigned si, unsigned L )
{
   switch(L) {
      case 0: {
         pOut[0] = pIn[0*si];
         break;
      }
      case 1: {
         pOut[0] = pIn[0*si];
         pOut[1] = pIn[1*si];
         pOut[2] = pIn[2*si];
         break;
      }
      case 2: {
         pOut[0] = -sd0*pIn[0*si] + sd2*pIn[3*si];
         pOut[1] = -sd0*pIn[0*si] - sd2*pIn[3*si];
         pOut[2] = pIn[0*si];
         pOut[3] = sd1*pIn[1*si];
         pOut[4] = sd1*pIn[2*si];
         pOut[5] = sd1*pIn[4*si];
         break;
      }
      case 3: {
         pOut[0] = -sd3*pIn[0*si] + sd7*pIn[3*si];
         pOut[1] = -sd3*pIn[0*si] - sd6*pIn[3*si];
         pOut[2] = sd4*pIn[0*si];
         pOut[3] = -sd3*pIn[1*si] + sd6*pIn[5*si];
         pOut[4] = -sd3*pIn[1*si] - sd7*pIn[5*si];
         pOut[5] = sd4*pIn[1*si];
         pOut[6] = -sd5*pIn[2*si] + sd9*pIn[6*si];
         pOut[7] = -sd5*pIn[2*si] - sd9*pIn[6*si];
         pOut[8] = pIn[2*si];
         pOut[9] = sd8*pIn[4*si];
         break;
      }
      case 4: {
         pOut[0] = sd11*pIn[0*si] + sd16*pIn[3*si] - sd19*pIn[5*si];
         pOut[1] = sd12*pIn[0*si] - sd17*pIn[3*si];
         pOut[2] = -sd10*pIn[0*si] + sd18*pIn[5*si];
         pOut[3] = sd11*pIn[0*si] + sd16*pIn[3*si] + sd19*pIn[5*si];
         pOut[4] = -sd10*pIn[0*si] - sd18*pIn[5*si];
         pOut[5] = pIn[0*si];
         pOut[6] = -sd13*pIn[1*si] + sd20*pIn[6*si];
         pOut[7] = -sd13*pIn[1*si] - sd20*pIn[6*si];
         pOut[8] = sd14*pIn[1*si];
         pOut[9] = -sd6*pIn[2*si] + sd22*pIn[7*si];
         pOut[10] = -sd6*pIn[2*si] - sd21*pIn[7*si];
         pOut[11] = sd15*pIn[2*si];
         pOut[12] = -sd6*pIn[4*si] + sd21*pIn[8*si];
         pOut[13] = -sd6*pIn[4*si] - sd22*pIn[8*si];
         pOut[14] = sd15*pIn[4*si];
         break;
      }
      case 5: {
         pOut[0] = sd25*pIn[0*si] - sd32*pIn[3*si] + sd36*pIn[9*si];
         pOut[1] = sd24*pIn[0*si] + sd30*pIn[3*si] - sd38*pIn[9*si];
         pOut[2] = -sd23*pIn[0*si] + sd29*pIn[3*si];
         pOut[3] = sd25*pIn[0*si] + sd31*pIn[3*si] + sd37*pIn[9*si];
         pOut[4] = -sd23*pIn[0*si] - sd28*pIn[3*si];
         pOut[5] = sd8*pIn[0*si];
         pOut[6] = sd25*pIn[1*si] - sd31*pIn[5*si] + sd37*pIn[7*si];
         pOut[7] = sd24*pIn[1*si] - sd30*pIn[5*si] - sd38*pIn[7*si];
         pOut[8] = -sd23*pIn[1*si] + sd28*pIn[5*si];
         pOut[9] = sd25*pIn[1*si] + sd32*pIn[5*si] + sd36*pIn[7*si];
         pOut[10] = -sd23*pIn[1*si] - sd29*pIn[5*si];
         pOut[11] = sd8*pIn[1*si];
         pOut[12] = -sd26*pIn[2*si] + sd34*pIn[6*si] + sd39*pIn[8*si];
         pOut[13] = -sd35*pIn[6*si] + sd41*pIn[8*si];
         pOut[14] = sd27*pIn[2*si] - sd40*pIn[8*si];
         pOut[15] = sd26*pIn[2*si] + sd34*pIn[6*si] + sd39*pIn[8*si];
         pOut[16] = -sd27*pIn[2*si] - sd40*pIn[8*si];
         pOut[17] = pIn[8*si];
         pOut[18] = sd33*pIn[4*si] - sd27*pIn[10*si];
         pOut[19] = -sd33*pIn[4*si] - sd27*pIn[10*si];
         pOut[20] = sd42*pIn[10*si];
         break;
      }
      case 6: {
         pOut[0] = sd43*pIn[0*si] - sd54*pIn[3*si] + sd56*pIn[5*si] - sd69*pIn[9*si];
         pOut[1] = -sd44*pIn[0*si] + sd53*pIn[3*si] + sd56*pIn[5*si] - sd68*pIn[9*si];
         pOut[2] = sd52*pIn[3*si] - sd55*pIn[5*si] + sd66*pIn[9*si];
         pOut[3] = sd44*pIn[0*si] + sd53*pIn[3*si] - sd56*pIn[5*si] - sd68*pIn[9*si];
         pOut[4] = -sd51*pIn[3*si] + sd65*pIn[9*si];
         pOut[5] = sd55*pIn[5*si] - sd67*pIn[9*si];
         pOut[6] = -sd43*pIn[0*si] - sd54*pIn[3*si] - sd56*pIn[5*si] - sd69*pIn[9*si];
         pOut[7] = sd52*pIn[3*si] + sd55*pIn[5*si] + sd66*pIn[9*si];
         pOut[8] = -sd55*pIn[5*si] - sd67*pIn[9*si];
         pOut[9] = pIn[9*si];
         pOut[10] = sd46*pIn[1*si] + sd58*pIn[6*si] - sd63*pIn[8*si];
         pOut[11] = sd47*pIn[1*si] - sd57*pIn[6*si];
         pOut[12] = -sd45*pIn[1*si] + sd64*pIn[8*si];
         pOut[13] = sd46*pIn[1*si] + sd58*pIn[6*si] + sd63*pIn[8*si];
         pOut[14] = -sd45*pIn[1*si] - sd64*pIn[8*si];
         pOut[15] = sd45*pIn[1*si];
         pOut[16] = sd49*pIn[2*si] - sd61*pIn[7*si] + sd71*pIn[12*si];
         pOut[17] = -sd50*pIn[2*si] + sd60*pIn[7*si] + sd73*pIn[12*si];
         pOut[18] = sd55*pIn[7*si] - sd70*pIn[12*si];
         pOut[19] = sd48*pIn[2*si] + sd59*pIn[7*si] + sd71*pIn[12*si];
         pOut[20] = -sd62*pIn[7*si] - sd70*pIn[12*si];
         pOut[21] = sd72*pIn[12*si];
         pOut[22] = sd48*pIn[4*si] - sd59*pIn[10*si] + sd71*pIn[11*si];
         pOut[23] = -sd50*pIn[4*si] - sd60*pIn[10*si] + sd73*pIn[11*si];
         pOut[24] = sd62*pIn[10*si] - sd70*pIn[11*si];
         pOut[25] = sd49*pIn[4*si] + sd61*pIn[10*si] + sd71*pIn[11*si];
         pOut[26] = -sd55*pIn[10*si] - sd70*pIn[11*si];
         pOut[27] = sd72*pIn[11*si];
         break;
      }
   }
}
} // namespace aic

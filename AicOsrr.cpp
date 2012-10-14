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

/* AicOsrr.cpp v20121011 EST [charge, Gerald Knizia] */
#include "AicOsrr.h"
#ifdef _DEBUG
   #include <assert.h>
#endif // _DEBUG
namespace aic{

void CaTrN1( double *AIC_RP pOut, double const *AIC_RP pIn, unsigned StrideIn, unsigned L );
void ShTrN_Indirect( double *AIC_RP pOut, unsigned ns, double const *AIC_RP pIn, unsigned short const *AIC_RP ii, unsigned L, unsigned Count );

void OsrrKrnB_3cen( FReal2 *AIC_RP pOut, FReal2 const *AIC_RP pIn, unsigned iDir, FOsrrParamsB const &P )
{
   double
      PmQf = P.PmQf[iDir],
      InvEtaABC2 = P.InvEtaABC2;
   switch( iDir * 10 + P.AngMomAB_Min )
   {
      case (0*10 + 0):
         pOut[0] = PmQf * pIn[0];
         if ( P.AngMomAB_Max == 0 ) break;
      case (0*10 + 1):
         pOut[1] = PmQf * pIn[1] + InvEtaABC2 * pIn[0];
         pOut[2] = PmQf * pIn[2];
         pOut[3] = PmQf * pIn[3];
         if ( P.AngMomAB_Max == 1 ) break;
      case (0*10 + 2):
         pOut[4] = PmQf * pIn[4] + 2 * InvEtaABC2 * pIn[1];
         pOut[5] = PmQf * pIn[5];
         pOut[6] = PmQf * pIn[6];
         pOut[7] = PmQf * pIn[7] + InvEtaABC2 * pIn[2];
         pOut[8] = PmQf * pIn[8] + InvEtaABC2 * pIn[3];
         pOut[9] = PmQf * pIn[9];
         if ( P.AngMomAB_Max == 2 ) break;
      case (0*10 + 3):
         pOut[10] = PmQf * pIn[10] + 3 * InvEtaABC2 * pIn[4];
         pOut[11] = PmQf * pIn[11] + InvEtaABC2 * pIn[5];
         pOut[12] = PmQf * pIn[12] + InvEtaABC2 * pIn[6];
         pOut[13] = PmQf * pIn[13] + 2 * InvEtaABC2 * pIn[7];
         pOut[14] = PmQf * pIn[14];
         pOut[15] = PmQf * pIn[15];
         pOut[16] = PmQf * pIn[16] + 2 * InvEtaABC2 * pIn[8];
         pOut[17] = PmQf * pIn[17];
         pOut[18] = PmQf * pIn[18];
         pOut[19] = PmQf * pIn[19] + InvEtaABC2 * pIn[9];
         if ( P.AngMomAB_Max == 3 ) break;
      case (0*10 + 4):
         pOut[20] = PmQf * pIn[20] + 4 * InvEtaABC2 * pIn[10];
         pOut[21] = PmQf * pIn[21] + 2 * InvEtaABC2 * pIn[11];
         pOut[22] = PmQf * pIn[22] + 2 * InvEtaABC2 * pIn[12];
         pOut[23] = PmQf * pIn[23];
         pOut[24] = PmQf * pIn[24];
         pOut[25] = PmQf * pIn[25];
         pOut[26] = PmQf * pIn[26] + 3 * InvEtaABC2 * pIn[13];
         pOut[27] = PmQf * pIn[27] + InvEtaABC2 * pIn[14];
         pOut[28] = PmQf * pIn[28] + InvEtaABC2 * pIn[15];
         pOut[29] = PmQf * pIn[29] + 3 * InvEtaABC2 * pIn[16];
         pOut[30] = PmQf * pIn[30] + InvEtaABC2 * pIn[17];
         pOut[31] = PmQf * pIn[31] + InvEtaABC2 * pIn[18];
         pOut[32] = PmQf * pIn[32] + 2 * InvEtaABC2 * pIn[19];
         pOut[33] = PmQf * pIn[33];
         pOut[34] = PmQf * pIn[34];
         if ( P.AngMomAB_Max == 4 ) break;
      case (0*10 + 5):
         pOut[35] = PmQf * pIn[35] + 5 * InvEtaABC2 * pIn[20];
         pOut[36] = PmQf * pIn[36] + 3 * InvEtaABC2 * pIn[21];
         pOut[37] = PmQf * pIn[37] + 3 * InvEtaABC2 * pIn[22];
         pOut[38] = PmQf * pIn[38] + InvEtaABC2 * pIn[23];
         pOut[39] = PmQf * pIn[39] + InvEtaABC2 * pIn[24];
         pOut[40] = PmQf * pIn[40] + InvEtaABC2 * pIn[25];
         pOut[41] = PmQf * pIn[41] + 4 * InvEtaABC2 * pIn[26];
         pOut[42] = PmQf * pIn[42] + 2 * InvEtaABC2 * pIn[27];
         pOut[43] = PmQf * pIn[43] + 2 * InvEtaABC2 * pIn[28];
         pOut[44] = PmQf * pIn[44];
         pOut[45] = PmQf * pIn[45];
         pOut[46] = PmQf * pIn[46];
         pOut[47] = PmQf * pIn[47] + 4 * InvEtaABC2 * pIn[29];
         pOut[48] = PmQf * pIn[48] + 2 * InvEtaABC2 * pIn[30];
         pOut[49] = PmQf * pIn[49] + 2 * InvEtaABC2 * pIn[31];
         pOut[50] = PmQf * pIn[50];
         pOut[51] = PmQf * pIn[51];
         pOut[52] = PmQf * pIn[52];
         pOut[53] = PmQf * pIn[53] + 3 * InvEtaABC2 * pIn[32];
         pOut[54] = PmQf * pIn[54] + InvEtaABC2 * pIn[33];
         pOut[55] = PmQf * pIn[55] + InvEtaABC2 * pIn[34];
         if ( P.AngMomAB_Max == 5 ) break;
      case (0*10 + 6):
         pOut[56] = PmQf * pIn[56] + 6 * InvEtaABC2 * pIn[35];
         pOut[57] = PmQf * pIn[57] + 4 * InvEtaABC2 * pIn[36];
         pOut[58] = PmQf * pIn[58] + 4 * InvEtaABC2 * pIn[37];
         pOut[59] = PmQf * pIn[59] + 2 * InvEtaABC2 * pIn[38];
         pOut[60] = PmQf * pIn[60] + 2 * InvEtaABC2 * pIn[39];
         pOut[61] = PmQf * pIn[61] + 2 * InvEtaABC2 * pIn[40];
         pOut[62] = PmQf * pIn[62];
         pOut[63] = PmQf * pIn[63];
         pOut[64] = PmQf * pIn[64];
         pOut[65] = PmQf * pIn[65];
         pOut[66] = PmQf * pIn[66] + 5 * InvEtaABC2 * pIn[41];
         pOut[67] = PmQf * pIn[67] + 3 * InvEtaABC2 * pIn[42];
         pOut[68] = PmQf * pIn[68] + 3 * InvEtaABC2 * pIn[43];
         pOut[69] = PmQf * pIn[69] + InvEtaABC2 * pIn[44];
         pOut[70] = PmQf * pIn[70] + InvEtaABC2 * pIn[45];
         pOut[71] = PmQf * pIn[71] + InvEtaABC2 * pIn[46];
         pOut[72] = PmQf * pIn[72] + 5 * InvEtaABC2 * pIn[47];
         pOut[73] = PmQf * pIn[73] + 3 * InvEtaABC2 * pIn[48];
         pOut[74] = PmQf * pIn[74] + 3 * InvEtaABC2 * pIn[49];
         pOut[75] = PmQf * pIn[75] + InvEtaABC2 * pIn[50];
         pOut[76] = PmQf * pIn[76] + InvEtaABC2 * pIn[51];
         pOut[77] = PmQf * pIn[77] + InvEtaABC2 * pIn[52];
         pOut[78] = PmQf * pIn[78] + 4 * InvEtaABC2 * pIn[53];
         pOut[79] = PmQf * pIn[79] + 2 * InvEtaABC2 * pIn[54];
         pOut[80] = PmQf * pIn[80] + 2 * InvEtaABC2 * pIn[55];
         pOut[81] = PmQf * pIn[81];
         pOut[82] = PmQf * pIn[82];
         pOut[83] = PmQf * pIn[83];
         if ( P.AngMomAB_Max == 6 ) break;
      case (0*10 + 7):
         pOut[84] = PmQf * pIn[84] + 7 * InvEtaABC2 * pIn[56];
         pOut[85] = PmQf * pIn[85] + 5 * InvEtaABC2 * pIn[57];
         pOut[86] = PmQf * pIn[86] + 5 * InvEtaABC2 * pIn[58];
         pOut[87] = PmQf * pIn[87] + 3 * InvEtaABC2 * pIn[59];
         pOut[88] = PmQf * pIn[88] + 3 * InvEtaABC2 * pIn[60];
         pOut[89] = PmQf * pIn[89] + 3 * InvEtaABC2 * pIn[61];
         pOut[90] = PmQf * pIn[90] + InvEtaABC2 * pIn[62];
         pOut[91] = PmQf * pIn[91] + InvEtaABC2 * pIn[63];
         pOut[92] = PmQf * pIn[92] + InvEtaABC2 * pIn[64];
         pOut[93] = PmQf * pIn[93] + InvEtaABC2 * pIn[65];
         pOut[94] = PmQf * pIn[94] + 6 * InvEtaABC2 * pIn[66];
         pOut[95] = PmQf * pIn[95] + 4 * InvEtaABC2 * pIn[67];
         pOut[96] = PmQf * pIn[96] + 4 * InvEtaABC2 * pIn[68];
         pOut[97] = PmQf * pIn[97] + 2 * InvEtaABC2 * pIn[69];
         pOut[98] = PmQf * pIn[98] + 2 * InvEtaABC2 * pIn[70];
         pOut[99] = PmQf * pIn[99] + 2 * InvEtaABC2 * pIn[71];
         pOut[100] = PmQf * pIn[100];
         pOut[101] = PmQf * pIn[101];
         pOut[102] = PmQf * pIn[102];
         pOut[103] = PmQf * pIn[103];
         pOut[104] = PmQf * pIn[104] + 6 * InvEtaABC2 * pIn[72];
         pOut[105] = PmQf * pIn[105] + 4 * InvEtaABC2 * pIn[73];
         pOut[106] = PmQf * pIn[106] + 4 * InvEtaABC2 * pIn[74];
         pOut[107] = PmQf * pIn[107] + 2 * InvEtaABC2 * pIn[75];
         pOut[108] = PmQf * pIn[108] + 2 * InvEtaABC2 * pIn[76];
         pOut[109] = PmQf * pIn[109] + 2 * InvEtaABC2 * pIn[77];
         pOut[110] = PmQf * pIn[110];
         pOut[111] = PmQf * pIn[111];
         pOut[112] = PmQf * pIn[112];
         pOut[113] = PmQf * pIn[113];
         pOut[114] = PmQf * pIn[114] + 5 * InvEtaABC2 * pIn[78];
         pOut[115] = PmQf * pIn[115] + 3 * InvEtaABC2 * pIn[79];
         pOut[116] = PmQf * pIn[116] + 3 * InvEtaABC2 * pIn[80];
         pOut[117] = PmQf * pIn[117] + InvEtaABC2 * pIn[81];
         pOut[118] = PmQf * pIn[118] + InvEtaABC2 * pIn[82];
         pOut[119] = PmQf * pIn[119] + InvEtaABC2 * pIn[83];
         if ( P.AngMomAB_Max == 7 ) break;
      case (0*10 + 8):
         pOut[120] = PmQf * pIn[120] + 8 * InvEtaABC2 * pIn[84];
         pOut[121] = PmQf * pIn[121] + 6 * InvEtaABC2 * pIn[85];
         pOut[122] = PmQf * pIn[122] + 6 * InvEtaABC2 * pIn[86];
         pOut[123] = PmQf * pIn[123] + 4 * InvEtaABC2 * pIn[87];
         pOut[124] = PmQf * pIn[124] + 4 * InvEtaABC2 * pIn[88];
         pOut[125] = PmQf * pIn[125] + 4 * InvEtaABC2 * pIn[89];
         pOut[126] = PmQf * pIn[126] + 2 * InvEtaABC2 * pIn[90];
         pOut[127] = PmQf * pIn[127] + 2 * InvEtaABC2 * pIn[91];
         pOut[128] = PmQf * pIn[128] + 2 * InvEtaABC2 * pIn[92];
         pOut[129] = PmQf * pIn[129] + 2 * InvEtaABC2 * pIn[93];
         pOut[130] = PmQf * pIn[130];
         pOut[131] = PmQf * pIn[131];
         pOut[132] = PmQf * pIn[132];
         pOut[133] = PmQf * pIn[133];
         pOut[134] = PmQf * pIn[134];
         pOut[135] = PmQf * pIn[135] + 7 * InvEtaABC2 * pIn[94];
         pOut[136] = PmQf * pIn[136] + 5 * InvEtaABC2 * pIn[95];
         pOut[137] = PmQf * pIn[137] + 5 * InvEtaABC2 * pIn[96];
         pOut[138] = PmQf * pIn[138] + 3 * InvEtaABC2 * pIn[97];
         pOut[139] = PmQf * pIn[139] + 3 * InvEtaABC2 * pIn[98];
         pOut[140] = PmQf * pIn[140] + 3 * InvEtaABC2 * pIn[99];
         pOut[141] = PmQf * pIn[141] + InvEtaABC2 * pIn[100];
         pOut[142] = PmQf * pIn[142] + InvEtaABC2 * pIn[101];
         pOut[143] = PmQf * pIn[143] + InvEtaABC2 * pIn[102];
         pOut[144] = PmQf * pIn[144] + InvEtaABC2 * pIn[103];
         pOut[145] = PmQf * pIn[145] + 7 * InvEtaABC2 * pIn[104];
         pOut[146] = PmQf * pIn[146] + 5 * InvEtaABC2 * pIn[105];
         pOut[147] = PmQf * pIn[147] + 5 * InvEtaABC2 * pIn[106];
         pOut[148] = PmQf * pIn[148] + 3 * InvEtaABC2 * pIn[107];
         pOut[149] = PmQf * pIn[149] + 3 * InvEtaABC2 * pIn[108];
         pOut[150] = PmQf * pIn[150] + 3 * InvEtaABC2 * pIn[109];
         pOut[151] = PmQf * pIn[151] + InvEtaABC2 * pIn[110];
         pOut[152] = PmQf * pIn[152] + InvEtaABC2 * pIn[111];
         pOut[153] = PmQf * pIn[153] + InvEtaABC2 * pIn[112];
         pOut[154] = PmQf * pIn[154] + InvEtaABC2 * pIn[113];
         pOut[155] = PmQf * pIn[155] + 6 * InvEtaABC2 * pIn[114];
         pOut[156] = PmQf * pIn[156] + 4 * InvEtaABC2 * pIn[115];
         pOut[157] = PmQf * pIn[157] + 4 * InvEtaABC2 * pIn[116];
         pOut[158] = PmQf * pIn[158] + 2 * InvEtaABC2 * pIn[117];
         pOut[159] = PmQf * pIn[159] + 2 * InvEtaABC2 * pIn[118];
         pOut[160] = PmQf * pIn[160] + 2 * InvEtaABC2 * pIn[119];
         pOut[161] = PmQf * pIn[161];
         pOut[162] = PmQf * pIn[162];
         pOut[163] = PmQf * pIn[163];
         pOut[164] = PmQf * pIn[164];
         if ( P.AngMomAB_Max == 8 ) break;
      case (0*10 + 9):
         pOut[165] = PmQf * pIn[165] + 9 * InvEtaABC2 * pIn[120];
         pOut[166] = PmQf * pIn[166] + 7 * InvEtaABC2 * pIn[121];
         pOut[167] = PmQf * pIn[167] + 7 * InvEtaABC2 * pIn[122];
         pOut[168] = PmQf * pIn[168] + 5 * InvEtaABC2 * pIn[123];
         pOut[169] = PmQf * pIn[169] + 5 * InvEtaABC2 * pIn[124];
         pOut[170] = PmQf * pIn[170] + 5 * InvEtaABC2 * pIn[125];
         pOut[171] = PmQf * pIn[171] + 3 * InvEtaABC2 * pIn[126];
         pOut[172] = PmQf * pIn[172] + 3 * InvEtaABC2 * pIn[127];
         pOut[173] = PmQf * pIn[173] + 3 * InvEtaABC2 * pIn[128];
         pOut[174] = PmQf * pIn[174] + 3 * InvEtaABC2 * pIn[129];
         pOut[175] = PmQf * pIn[175] + InvEtaABC2 * pIn[130];
         pOut[176] = PmQf * pIn[176] + InvEtaABC2 * pIn[131];
         pOut[177] = PmQf * pIn[177] + InvEtaABC2 * pIn[132];
         pOut[178] = PmQf * pIn[178] + InvEtaABC2 * pIn[133];
         pOut[179] = PmQf * pIn[179] + InvEtaABC2 * pIn[134];
         pOut[180] = PmQf * pIn[180] + 8 * InvEtaABC2 * pIn[135];
         pOut[181] = PmQf * pIn[181] + 6 * InvEtaABC2 * pIn[136];
         pOut[182] = PmQf * pIn[182] + 6 * InvEtaABC2 * pIn[137];
         pOut[183] = PmQf * pIn[183] + 4 * InvEtaABC2 * pIn[138];
         pOut[184] = PmQf * pIn[184] + 4 * InvEtaABC2 * pIn[139];
         pOut[185] = PmQf * pIn[185] + 4 * InvEtaABC2 * pIn[140];
         pOut[186] = PmQf * pIn[186] + 2 * InvEtaABC2 * pIn[141];
         pOut[187] = PmQf * pIn[187] + 2 * InvEtaABC2 * pIn[142];
         pOut[188] = PmQf * pIn[188] + 2 * InvEtaABC2 * pIn[143];
         pOut[189] = PmQf * pIn[189] + 2 * InvEtaABC2 * pIn[144];
         pOut[190] = PmQf * pIn[190];
         pOut[191] = PmQf * pIn[191];
         pOut[192] = PmQf * pIn[192];
         pOut[193] = PmQf * pIn[193];
         pOut[194] = PmQf * pIn[194];
         pOut[195] = PmQf * pIn[195] + 8 * InvEtaABC2 * pIn[145];
         pOut[196] = PmQf * pIn[196] + 6 * InvEtaABC2 * pIn[146];
         pOut[197] = PmQf * pIn[197] + 6 * InvEtaABC2 * pIn[147];
         pOut[198] = PmQf * pIn[198] + 4 * InvEtaABC2 * pIn[148];
         pOut[199] = PmQf * pIn[199] + 4 * InvEtaABC2 * pIn[149];
         pOut[200] = PmQf * pIn[200] + 4 * InvEtaABC2 * pIn[150];
         pOut[201] = PmQf * pIn[201] + 2 * InvEtaABC2 * pIn[151];
         pOut[202] = PmQf * pIn[202] + 2 * InvEtaABC2 * pIn[152];
         pOut[203] = PmQf * pIn[203] + 2 * InvEtaABC2 * pIn[153];
         pOut[204] = PmQf * pIn[204] + 2 * InvEtaABC2 * pIn[154];
         pOut[205] = PmQf * pIn[205];
         pOut[206] = PmQf * pIn[206];
         pOut[207] = PmQf * pIn[207];
         pOut[208] = PmQf * pIn[208];
         pOut[209] = PmQf * pIn[209];
         pOut[210] = PmQf * pIn[210] + 7 * InvEtaABC2 * pIn[155];
         pOut[211] = PmQf * pIn[211] + 5 * InvEtaABC2 * pIn[156];
         pOut[212] = PmQf * pIn[212] + 5 * InvEtaABC2 * pIn[157];
         pOut[213] = PmQf * pIn[213] + 3 * InvEtaABC2 * pIn[158];
         pOut[214] = PmQf * pIn[214] + 3 * InvEtaABC2 * pIn[159];
         pOut[215] = PmQf * pIn[215] + 3 * InvEtaABC2 * pIn[160];
         pOut[216] = PmQf * pIn[216] + InvEtaABC2 * pIn[161];
         pOut[217] = PmQf * pIn[217] + InvEtaABC2 * pIn[162];
         pOut[218] = PmQf * pIn[218] + InvEtaABC2 * pIn[163];
         pOut[219] = PmQf * pIn[219] + InvEtaABC2 * pIn[164];
         if ( P.AngMomAB_Max == 9 ) break;

      #ifdef _DEBUG
         assert(0);
      #endif
         break;
      case (1*10 + 0):
         pOut[0] = PmQf * pIn[0];
         if ( P.AngMomAB_Max == 0 ) break;
      case (1*10 + 1):
         pOut[1] = PmQf * pIn[1];
         pOut[2] = PmQf * pIn[2] + InvEtaABC2 * pIn[0];
         pOut[3] = PmQf * pIn[3];
         if ( P.AngMomAB_Max == 1 ) break;
      case (1*10 + 2):
         pOut[4] = PmQf * pIn[4];
         pOut[5] = PmQf * pIn[5] + 2 * InvEtaABC2 * pIn[2];
         pOut[6] = PmQf * pIn[6];
         pOut[7] = PmQf * pIn[7] + InvEtaABC2 * pIn[1];
         pOut[8] = PmQf * pIn[8];
         pOut[9] = PmQf * pIn[9] + InvEtaABC2 * pIn[3];
         if ( P.AngMomAB_Max == 2 ) break;
      case (1*10 + 3):
         pOut[10] = PmQf * pIn[10];
         pOut[11] = PmQf * pIn[11] + 2 * InvEtaABC2 * pIn[7];
         pOut[12] = PmQf * pIn[12];
         pOut[13] = PmQf * pIn[13] + InvEtaABC2 * pIn[4];
         pOut[14] = PmQf * pIn[14] + 3 * InvEtaABC2 * pIn[5];
         pOut[15] = PmQf * pIn[15] + InvEtaABC2 * pIn[6];
         pOut[16] = PmQf * pIn[16];
         pOut[17] = PmQf * pIn[17] + 2 * InvEtaABC2 * pIn[9];
         pOut[18] = PmQf * pIn[18];
         pOut[19] = PmQf * pIn[19] + InvEtaABC2 * pIn[8];
         if ( P.AngMomAB_Max == 3 ) break;
      case (1*10 + 4):
         pOut[20] = PmQf * pIn[20];
         pOut[21] = PmQf * pIn[21] + 2 * InvEtaABC2 * pIn[13];
         pOut[22] = PmQf * pIn[22];
         pOut[23] = PmQf * pIn[23] + 4 * InvEtaABC2 * pIn[14];
         pOut[24] = PmQf * pIn[24] + 2 * InvEtaABC2 * pIn[15];
         pOut[25] = PmQf * pIn[25];
         pOut[26] = PmQf * pIn[26] + InvEtaABC2 * pIn[10];
         pOut[27] = PmQf * pIn[27] + 3 * InvEtaABC2 * pIn[11];
         pOut[28] = PmQf * pIn[28] + InvEtaABC2 * pIn[12];
         pOut[29] = PmQf * pIn[29];
         pOut[30] = PmQf * pIn[30] + 2 * InvEtaABC2 * pIn[19];
         pOut[31] = PmQf * pIn[31];
         pOut[32] = PmQf * pIn[32] + InvEtaABC2 * pIn[16];
         pOut[33] = PmQf * pIn[33] + 3 * InvEtaABC2 * pIn[17];
         pOut[34] = PmQf * pIn[34] + InvEtaABC2 * pIn[18];
         if ( P.AngMomAB_Max == 4 ) break;
      case (1*10 + 5):
         pOut[35] = PmQf * pIn[35];
         pOut[36] = PmQf * pIn[36] + 2 * InvEtaABC2 * pIn[26];
         pOut[37] = PmQf * pIn[37];
         pOut[38] = PmQf * pIn[38] + 4 * InvEtaABC2 * pIn[27];
         pOut[39] = PmQf * pIn[39] + 2 * InvEtaABC2 * pIn[28];
         pOut[40] = PmQf * pIn[40];
         pOut[41] = PmQf * pIn[41] + InvEtaABC2 * pIn[20];
         pOut[42] = PmQf * pIn[42] + 3 * InvEtaABC2 * pIn[21];
         pOut[43] = PmQf * pIn[43] + InvEtaABC2 * pIn[22];
         pOut[44] = PmQf * pIn[44] + 5 * InvEtaABC2 * pIn[23];
         pOut[45] = PmQf * pIn[45] + 3 * InvEtaABC2 * pIn[24];
         pOut[46] = PmQf * pIn[46] + InvEtaABC2 * pIn[25];
         pOut[47] = PmQf * pIn[47];
         pOut[48] = PmQf * pIn[48] + 2 * InvEtaABC2 * pIn[32];
         pOut[49] = PmQf * pIn[49];
         pOut[50] = PmQf * pIn[50] + 4 * InvEtaABC2 * pIn[33];
         pOut[51] = PmQf * pIn[51] + 2 * InvEtaABC2 * pIn[34];
         pOut[52] = PmQf * pIn[52];
         pOut[53] = PmQf * pIn[53] + InvEtaABC2 * pIn[29];
         pOut[54] = PmQf * pIn[54] + 3 * InvEtaABC2 * pIn[30];
         pOut[55] = PmQf * pIn[55] + InvEtaABC2 * pIn[31];
         if ( P.AngMomAB_Max == 5 ) break;
      case (1*10 + 6):
         pOut[56] = PmQf * pIn[56];
         pOut[57] = PmQf * pIn[57] + 2 * InvEtaABC2 * pIn[41];
         pOut[58] = PmQf * pIn[58];
         pOut[59] = PmQf * pIn[59] + 4 * InvEtaABC2 * pIn[42];
         pOut[60] = PmQf * pIn[60] + 2 * InvEtaABC2 * pIn[43];
         pOut[61] = PmQf * pIn[61];
         pOut[62] = PmQf * pIn[62] + 6 * InvEtaABC2 * pIn[44];
         pOut[63] = PmQf * pIn[63] + 4 * InvEtaABC2 * pIn[45];
         pOut[64] = PmQf * pIn[64] + 2 * InvEtaABC2 * pIn[46];
         pOut[65] = PmQf * pIn[65];
         pOut[66] = PmQf * pIn[66] + InvEtaABC2 * pIn[35];
         pOut[67] = PmQf * pIn[67] + 3 * InvEtaABC2 * pIn[36];
         pOut[68] = PmQf * pIn[68] + InvEtaABC2 * pIn[37];
         pOut[69] = PmQf * pIn[69] + 5 * InvEtaABC2 * pIn[38];
         pOut[70] = PmQf * pIn[70] + 3 * InvEtaABC2 * pIn[39];
         pOut[71] = PmQf * pIn[71] + InvEtaABC2 * pIn[40];
         pOut[72] = PmQf * pIn[72];
         pOut[73] = PmQf * pIn[73] + 2 * InvEtaABC2 * pIn[53];
         pOut[74] = PmQf * pIn[74];
         pOut[75] = PmQf * pIn[75] + 4 * InvEtaABC2 * pIn[54];
         pOut[76] = PmQf * pIn[76] + 2 * InvEtaABC2 * pIn[55];
         pOut[77] = PmQf * pIn[77];
         pOut[78] = PmQf * pIn[78] + InvEtaABC2 * pIn[47];
         pOut[79] = PmQf * pIn[79] + 3 * InvEtaABC2 * pIn[48];
         pOut[80] = PmQf * pIn[80] + InvEtaABC2 * pIn[49];
         pOut[81] = PmQf * pIn[81] + 5 * InvEtaABC2 * pIn[50];
         pOut[82] = PmQf * pIn[82] + 3 * InvEtaABC2 * pIn[51];
         pOut[83] = PmQf * pIn[83] + InvEtaABC2 * pIn[52];
         if ( P.AngMomAB_Max == 6 ) break;
      case (1*10 + 7):
         pOut[84] = PmQf * pIn[84];
         pOut[85] = PmQf * pIn[85] + 2 * InvEtaABC2 * pIn[66];
         pOut[86] = PmQf * pIn[86];
         pOut[87] = PmQf * pIn[87] + 4 * InvEtaABC2 * pIn[67];
         pOut[88] = PmQf * pIn[88] + 2 * InvEtaABC2 * pIn[68];
         pOut[89] = PmQf * pIn[89];
         pOut[90] = PmQf * pIn[90] + 6 * InvEtaABC2 * pIn[69];
         pOut[91] = PmQf * pIn[91] + 4 * InvEtaABC2 * pIn[70];
         pOut[92] = PmQf * pIn[92] + 2 * InvEtaABC2 * pIn[71];
         pOut[93] = PmQf * pIn[93];
         pOut[94] = PmQf * pIn[94] + InvEtaABC2 * pIn[56];
         pOut[95] = PmQf * pIn[95] + 3 * InvEtaABC2 * pIn[57];
         pOut[96] = PmQf * pIn[96] + InvEtaABC2 * pIn[58];
         pOut[97] = PmQf * pIn[97] + 5 * InvEtaABC2 * pIn[59];
         pOut[98] = PmQf * pIn[98] + 3 * InvEtaABC2 * pIn[60];
         pOut[99] = PmQf * pIn[99] + InvEtaABC2 * pIn[61];
         pOut[100] = PmQf * pIn[100] + 7 * InvEtaABC2 * pIn[62];
         pOut[101] = PmQf * pIn[101] + 5 * InvEtaABC2 * pIn[63];
         pOut[102] = PmQf * pIn[102] + 3 * InvEtaABC2 * pIn[64];
         pOut[103] = PmQf * pIn[103] + InvEtaABC2 * pIn[65];
         pOut[104] = PmQf * pIn[104];
         pOut[105] = PmQf * pIn[105] + 2 * InvEtaABC2 * pIn[78];
         pOut[106] = PmQf * pIn[106];
         pOut[107] = PmQf * pIn[107] + 4 * InvEtaABC2 * pIn[79];
         pOut[108] = PmQf * pIn[108] + 2 * InvEtaABC2 * pIn[80];
         pOut[109] = PmQf * pIn[109];
         pOut[110] = PmQf * pIn[110] + 6 * InvEtaABC2 * pIn[81];
         pOut[111] = PmQf * pIn[111] + 4 * InvEtaABC2 * pIn[82];
         pOut[112] = PmQf * pIn[112] + 2 * InvEtaABC2 * pIn[83];
         pOut[113] = PmQf * pIn[113];
         pOut[114] = PmQf * pIn[114] + InvEtaABC2 * pIn[72];
         pOut[115] = PmQf * pIn[115] + 3 * InvEtaABC2 * pIn[73];
         pOut[116] = PmQf * pIn[116] + InvEtaABC2 * pIn[74];
         pOut[117] = PmQf * pIn[117] + 5 * InvEtaABC2 * pIn[75];
         pOut[118] = PmQf * pIn[118] + 3 * InvEtaABC2 * pIn[76];
         pOut[119] = PmQf * pIn[119] + InvEtaABC2 * pIn[77];
         if ( P.AngMomAB_Max == 7 ) break;
      case (1*10 + 8):
         pOut[120] = PmQf * pIn[120];
         pOut[121] = PmQf * pIn[121] + 2 * InvEtaABC2 * pIn[94];
         pOut[122] = PmQf * pIn[122];
         pOut[123] = PmQf * pIn[123] + 4 * InvEtaABC2 * pIn[95];
         pOut[124] = PmQf * pIn[124] + 2 * InvEtaABC2 * pIn[96];
         pOut[125] = PmQf * pIn[125];
         pOut[126] = PmQf * pIn[126] + 6 * InvEtaABC2 * pIn[97];
         pOut[127] = PmQf * pIn[127] + 4 * InvEtaABC2 * pIn[98];
         pOut[128] = PmQf * pIn[128] + 2 * InvEtaABC2 * pIn[99];
         pOut[129] = PmQf * pIn[129];
         pOut[130] = PmQf * pIn[130] + 8 * InvEtaABC2 * pIn[100];
         pOut[131] = PmQf * pIn[131] + 6 * InvEtaABC2 * pIn[101];
         pOut[132] = PmQf * pIn[132] + 4 * InvEtaABC2 * pIn[102];
         pOut[133] = PmQf * pIn[133] + 2 * InvEtaABC2 * pIn[103];
         pOut[134] = PmQf * pIn[134];
         pOut[135] = PmQf * pIn[135] + InvEtaABC2 * pIn[84];
         pOut[136] = PmQf * pIn[136] + 3 * InvEtaABC2 * pIn[85];
         pOut[137] = PmQf * pIn[137] + InvEtaABC2 * pIn[86];
         pOut[138] = PmQf * pIn[138] + 5 * InvEtaABC2 * pIn[87];
         pOut[139] = PmQf * pIn[139] + 3 * InvEtaABC2 * pIn[88];
         pOut[140] = PmQf * pIn[140] + InvEtaABC2 * pIn[89];
         pOut[141] = PmQf * pIn[141] + 7 * InvEtaABC2 * pIn[90];
         pOut[142] = PmQf * pIn[142] + 5 * InvEtaABC2 * pIn[91];
         pOut[143] = PmQf * pIn[143] + 3 * InvEtaABC2 * pIn[92];
         pOut[144] = PmQf * pIn[144] + InvEtaABC2 * pIn[93];
         pOut[145] = PmQf * pIn[145];
         pOut[146] = PmQf * pIn[146] + 2 * InvEtaABC2 * pIn[114];
         pOut[147] = PmQf * pIn[147];
         pOut[148] = PmQf * pIn[148] + 4 * InvEtaABC2 * pIn[115];
         pOut[149] = PmQf * pIn[149] + 2 * InvEtaABC2 * pIn[116];
         pOut[150] = PmQf * pIn[150];
         pOut[151] = PmQf * pIn[151] + 6 * InvEtaABC2 * pIn[117];
         pOut[152] = PmQf * pIn[152] + 4 * InvEtaABC2 * pIn[118];
         pOut[153] = PmQf * pIn[153] + 2 * InvEtaABC2 * pIn[119];
         pOut[154] = PmQf * pIn[154];
         pOut[155] = PmQf * pIn[155] + InvEtaABC2 * pIn[104];
         pOut[156] = PmQf * pIn[156] + 3 * InvEtaABC2 * pIn[105];
         pOut[157] = PmQf * pIn[157] + InvEtaABC2 * pIn[106];
         pOut[158] = PmQf * pIn[158] + 5 * InvEtaABC2 * pIn[107];
         pOut[159] = PmQf * pIn[159] + 3 * InvEtaABC2 * pIn[108];
         pOut[160] = PmQf * pIn[160] + InvEtaABC2 * pIn[109];
         pOut[161] = PmQf * pIn[161] + 7 * InvEtaABC2 * pIn[110];
         pOut[162] = PmQf * pIn[162] + 5 * InvEtaABC2 * pIn[111];
         pOut[163] = PmQf * pIn[163] + 3 * InvEtaABC2 * pIn[112];
         pOut[164] = PmQf * pIn[164] + InvEtaABC2 * pIn[113];
         if ( P.AngMomAB_Max == 8 ) break;
      case (1*10 + 9):
         pOut[165] = PmQf * pIn[165];
         pOut[166] = PmQf * pIn[166] + 2 * InvEtaABC2 * pIn[135];
         pOut[167] = PmQf * pIn[167];
         pOut[168] = PmQf * pIn[168] + 4 * InvEtaABC2 * pIn[136];
         pOut[169] = PmQf * pIn[169] + 2 * InvEtaABC2 * pIn[137];
         pOut[170] = PmQf * pIn[170];
         pOut[171] = PmQf * pIn[171] + 6 * InvEtaABC2 * pIn[138];
         pOut[172] = PmQf * pIn[172] + 4 * InvEtaABC2 * pIn[139];
         pOut[173] = PmQf * pIn[173] + 2 * InvEtaABC2 * pIn[140];
         pOut[174] = PmQf * pIn[174];
         pOut[175] = PmQf * pIn[175] + 8 * InvEtaABC2 * pIn[141];
         pOut[176] = PmQf * pIn[176] + 6 * InvEtaABC2 * pIn[142];
         pOut[177] = PmQf * pIn[177] + 4 * InvEtaABC2 * pIn[143];
         pOut[178] = PmQf * pIn[178] + 2 * InvEtaABC2 * pIn[144];
         pOut[179] = PmQf * pIn[179];
         pOut[180] = PmQf * pIn[180] + InvEtaABC2 * pIn[120];
         pOut[181] = PmQf * pIn[181] + 3 * InvEtaABC2 * pIn[121];
         pOut[182] = PmQf * pIn[182] + InvEtaABC2 * pIn[122];
         pOut[183] = PmQf * pIn[183] + 5 * InvEtaABC2 * pIn[123];
         pOut[184] = PmQf * pIn[184] + 3 * InvEtaABC2 * pIn[124];
         pOut[185] = PmQf * pIn[185] + InvEtaABC2 * pIn[125];
         pOut[186] = PmQf * pIn[186] + 7 * InvEtaABC2 * pIn[126];
         pOut[187] = PmQf * pIn[187] + 5 * InvEtaABC2 * pIn[127];
         pOut[188] = PmQf * pIn[188] + 3 * InvEtaABC2 * pIn[128];
         pOut[189] = PmQf * pIn[189] + InvEtaABC2 * pIn[129];
         pOut[190] = PmQf * pIn[190] + 9 * InvEtaABC2 * pIn[130];
         pOut[191] = PmQf * pIn[191] + 7 * InvEtaABC2 * pIn[131];
         pOut[192] = PmQf * pIn[192] + 5 * InvEtaABC2 * pIn[132];
         pOut[193] = PmQf * pIn[193] + 3 * InvEtaABC2 * pIn[133];
         pOut[194] = PmQf * pIn[194] + InvEtaABC2 * pIn[134];
         pOut[195] = PmQf * pIn[195];
         pOut[196] = PmQf * pIn[196] + 2 * InvEtaABC2 * pIn[155];
         pOut[197] = PmQf * pIn[197];
         pOut[198] = PmQf * pIn[198] + 4 * InvEtaABC2 * pIn[156];
         pOut[199] = PmQf * pIn[199] + 2 * InvEtaABC2 * pIn[157];
         pOut[200] = PmQf * pIn[200];
         pOut[201] = PmQf * pIn[201] + 6 * InvEtaABC2 * pIn[158];
         pOut[202] = PmQf * pIn[202] + 4 * InvEtaABC2 * pIn[159];
         pOut[203] = PmQf * pIn[203] + 2 * InvEtaABC2 * pIn[160];
         pOut[204] = PmQf * pIn[204];
         pOut[205] = PmQf * pIn[205] + 8 * InvEtaABC2 * pIn[161];
         pOut[206] = PmQf * pIn[206] + 6 * InvEtaABC2 * pIn[162];
         pOut[207] = PmQf * pIn[207] + 4 * InvEtaABC2 * pIn[163];
         pOut[208] = PmQf * pIn[208] + 2 * InvEtaABC2 * pIn[164];
         pOut[209] = PmQf * pIn[209];
         pOut[210] = PmQf * pIn[210] + InvEtaABC2 * pIn[145];
         pOut[211] = PmQf * pIn[211] + 3 * InvEtaABC2 * pIn[146];
         pOut[212] = PmQf * pIn[212] + InvEtaABC2 * pIn[147];
         pOut[213] = PmQf * pIn[213] + 5 * InvEtaABC2 * pIn[148];
         pOut[214] = PmQf * pIn[214] + 3 * InvEtaABC2 * pIn[149];
         pOut[215] = PmQf * pIn[215] + InvEtaABC2 * pIn[150];
         pOut[216] = PmQf * pIn[216] + 7 * InvEtaABC2 * pIn[151];
         pOut[217] = PmQf * pIn[217] + 5 * InvEtaABC2 * pIn[152];
         pOut[218] = PmQf * pIn[218] + 3 * InvEtaABC2 * pIn[153];
         pOut[219] = PmQf * pIn[219] + InvEtaABC2 * pIn[154];
         if ( P.AngMomAB_Max == 9 ) break;

      #ifdef _DEBUG
         assert(0);
      #endif
         break;
      case (2*10 + 0):
         pOut[0] = PmQf * pIn[0];
         if ( P.AngMomAB_Max == 0 ) break;
      case (2*10 + 1):
         pOut[1] = PmQf * pIn[1];
         pOut[2] = PmQf * pIn[2];
         pOut[3] = PmQf * pIn[3] + InvEtaABC2 * pIn[0];
         if ( P.AngMomAB_Max == 1 ) break;
      case (2*10 + 2):
         pOut[4] = PmQf * pIn[4];
         pOut[5] = PmQf * pIn[5];
         pOut[6] = PmQf * pIn[6] + 2 * InvEtaABC2 * pIn[3];
         pOut[7] = PmQf * pIn[7];
         pOut[8] = PmQf * pIn[8] + InvEtaABC2 * pIn[1];
         pOut[9] = PmQf * pIn[9] + InvEtaABC2 * pIn[2];
         if ( P.AngMomAB_Max == 2 ) break;
      case (2*10 + 3):
         pOut[10] = PmQf * pIn[10];
         pOut[11] = PmQf * pIn[11];
         pOut[12] = PmQf * pIn[12] + 2 * InvEtaABC2 * pIn[8];
         pOut[13] = PmQf * pIn[13];
         pOut[14] = PmQf * pIn[14];
         pOut[15] = PmQf * pIn[15] + 2 * InvEtaABC2 * pIn[9];
         pOut[16] = PmQf * pIn[16] + InvEtaABC2 * pIn[4];
         pOut[17] = PmQf * pIn[17] + InvEtaABC2 * pIn[5];
         pOut[18] = PmQf * pIn[18] + 3 * InvEtaABC2 * pIn[6];
         pOut[19] = PmQf * pIn[19] + InvEtaABC2 * pIn[7];
         if ( P.AngMomAB_Max == 3 ) break;
      case (2*10 + 4):
         pOut[20] = PmQf * pIn[20];
         pOut[21] = PmQf * pIn[21];
         pOut[22] = PmQf * pIn[22] + 2 * InvEtaABC2 * pIn[16];
         pOut[23] = PmQf * pIn[23];
         pOut[24] = PmQf * pIn[24] + 2 * InvEtaABC2 * pIn[17];
         pOut[25] = PmQf * pIn[25] + 4 * InvEtaABC2 * pIn[18];
         pOut[26] = PmQf * pIn[26];
         pOut[27] = PmQf * pIn[27];
         pOut[28] = PmQf * pIn[28] + 2 * InvEtaABC2 * pIn[19];
         pOut[29] = PmQf * pIn[29] + InvEtaABC2 * pIn[10];
         pOut[30] = PmQf * pIn[30] + InvEtaABC2 * pIn[11];
         pOut[31] = PmQf * pIn[31] + 3 * InvEtaABC2 * pIn[12];
         pOut[32] = PmQf * pIn[32] + InvEtaABC2 * pIn[13];
         pOut[33] = PmQf * pIn[33] + InvEtaABC2 * pIn[14];
         pOut[34] = PmQf * pIn[34] + 3 * InvEtaABC2 * pIn[15];
         if ( P.AngMomAB_Max == 4 ) break;
      case (2*10 + 5):
         pOut[35] = PmQf * pIn[35];
         pOut[36] = PmQf * pIn[36];
         pOut[37] = PmQf * pIn[37] + 2 * InvEtaABC2 * pIn[29];
         pOut[38] = PmQf * pIn[38];
         pOut[39] = PmQf * pIn[39] + 2 * InvEtaABC2 * pIn[30];
         pOut[40] = PmQf * pIn[40] + 4 * InvEtaABC2 * pIn[31];
         pOut[41] = PmQf * pIn[41];
         pOut[42] = PmQf * pIn[42];
         pOut[43] = PmQf * pIn[43] + 2 * InvEtaABC2 * pIn[32];
         pOut[44] = PmQf * pIn[44];
         pOut[45] = PmQf * pIn[45] + 2 * InvEtaABC2 * pIn[33];
         pOut[46] = PmQf * pIn[46] + 4 * InvEtaABC2 * pIn[34];
         pOut[47] = PmQf * pIn[47] + InvEtaABC2 * pIn[20];
         pOut[48] = PmQf * pIn[48] + InvEtaABC2 * pIn[21];
         pOut[49] = PmQf * pIn[49] + 3 * InvEtaABC2 * pIn[22];
         pOut[50] = PmQf * pIn[50] + InvEtaABC2 * pIn[23];
         pOut[51] = PmQf * pIn[51] + 3 * InvEtaABC2 * pIn[24];
         pOut[52] = PmQf * pIn[52] + 5 * InvEtaABC2 * pIn[25];
         pOut[53] = PmQf * pIn[53] + InvEtaABC2 * pIn[26];
         pOut[54] = PmQf * pIn[54] + InvEtaABC2 * pIn[27];
         pOut[55] = PmQf * pIn[55] + 3 * InvEtaABC2 * pIn[28];
         if ( P.AngMomAB_Max == 5 ) break;
      case (2*10 + 6):
         pOut[56] = PmQf * pIn[56];
         pOut[57] = PmQf * pIn[57];
         pOut[58] = PmQf * pIn[58] + 2 * InvEtaABC2 * pIn[47];
         pOut[59] = PmQf * pIn[59];
         pOut[60] = PmQf * pIn[60] + 2 * InvEtaABC2 * pIn[48];
         pOut[61] = PmQf * pIn[61] + 4 * InvEtaABC2 * pIn[49];
         pOut[62] = PmQf * pIn[62];
         pOut[63] = PmQf * pIn[63] + 2 * InvEtaABC2 * pIn[50];
         pOut[64] = PmQf * pIn[64] + 4 * InvEtaABC2 * pIn[51];
         pOut[65] = PmQf * pIn[65] + 6 * InvEtaABC2 * pIn[52];
         pOut[66] = PmQf * pIn[66];
         pOut[67] = PmQf * pIn[67];
         pOut[68] = PmQf * pIn[68] + 2 * InvEtaABC2 * pIn[53];
         pOut[69] = PmQf * pIn[69];
         pOut[70] = PmQf * pIn[70] + 2 * InvEtaABC2 * pIn[54];
         pOut[71] = PmQf * pIn[71] + 4 * InvEtaABC2 * pIn[55];
         pOut[72] = PmQf * pIn[72] + InvEtaABC2 * pIn[35];
         pOut[73] = PmQf * pIn[73] + InvEtaABC2 * pIn[36];
         pOut[74] = PmQf * pIn[74] + 3 * InvEtaABC2 * pIn[37];
         pOut[75] = PmQf * pIn[75] + InvEtaABC2 * pIn[38];
         pOut[76] = PmQf * pIn[76] + 3 * InvEtaABC2 * pIn[39];
         pOut[77] = PmQf * pIn[77] + 5 * InvEtaABC2 * pIn[40];
         pOut[78] = PmQf * pIn[78] + InvEtaABC2 * pIn[41];
         pOut[79] = PmQf * pIn[79] + InvEtaABC2 * pIn[42];
         pOut[80] = PmQf * pIn[80] + 3 * InvEtaABC2 * pIn[43];
         pOut[81] = PmQf * pIn[81] + InvEtaABC2 * pIn[44];
         pOut[82] = PmQf * pIn[82] + 3 * InvEtaABC2 * pIn[45];
         pOut[83] = PmQf * pIn[83] + 5 * InvEtaABC2 * pIn[46];
         if ( P.AngMomAB_Max == 6 ) break;
      case (2*10 + 7):
         pOut[84] = PmQf * pIn[84];
         pOut[85] = PmQf * pIn[85];
         pOut[86] = PmQf * pIn[86] + 2 * InvEtaABC2 * pIn[72];
         pOut[87] = PmQf * pIn[87];
         pOut[88] = PmQf * pIn[88] + 2 * InvEtaABC2 * pIn[73];
         pOut[89] = PmQf * pIn[89] + 4 * InvEtaABC2 * pIn[74];
         pOut[90] = PmQf * pIn[90];
         pOut[91] = PmQf * pIn[91] + 2 * InvEtaABC2 * pIn[75];
         pOut[92] = PmQf * pIn[92] + 4 * InvEtaABC2 * pIn[76];
         pOut[93] = PmQf * pIn[93] + 6 * InvEtaABC2 * pIn[77];
         pOut[94] = PmQf * pIn[94];
         pOut[95] = PmQf * pIn[95];
         pOut[96] = PmQf * pIn[96] + 2 * InvEtaABC2 * pIn[78];
         pOut[97] = PmQf * pIn[97];
         pOut[98] = PmQf * pIn[98] + 2 * InvEtaABC2 * pIn[79];
         pOut[99] = PmQf * pIn[99] + 4 * InvEtaABC2 * pIn[80];
         pOut[100] = PmQf * pIn[100];
         pOut[101] = PmQf * pIn[101] + 2 * InvEtaABC2 * pIn[81];
         pOut[102] = PmQf * pIn[102] + 4 * InvEtaABC2 * pIn[82];
         pOut[103] = PmQf * pIn[103] + 6 * InvEtaABC2 * pIn[83];
         pOut[104] = PmQf * pIn[104] + InvEtaABC2 * pIn[56];
         pOut[105] = PmQf * pIn[105] + InvEtaABC2 * pIn[57];
         pOut[106] = PmQf * pIn[106] + 3 * InvEtaABC2 * pIn[58];
         pOut[107] = PmQf * pIn[107] + InvEtaABC2 * pIn[59];
         pOut[108] = PmQf * pIn[108] + 3 * InvEtaABC2 * pIn[60];
         pOut[109] = PmQf * pIn[109] + 5 * InvEtaABC2 * pIn[61];
         pOut[110] = PmQf * pIn[110] + InvEtaABC2 * pIn[62];
         pOut[111] = PmQf * pIn[111] + 3 * InvEtaABC2 * pIn[63];
         pOut[112] = PmQf * pIn[112] + 5 * InvEtaABC2 * pIn[64];
         pOut[113] = PmQf * pIn[113] + 7 * InvEtaABC2 * pIn[65];
         pOut[114] = PmQf * pIn[114] + InvEtaABC2 * pIn[66];
         pOut[115] = PmQf * pIn[115] + InvEtaABC2 * pIn[67];
         pOut[116] = PmQf * pIn[116] + 3 * InvEtaABC2 * pIn[68];
         pOut[117] = PmQf * pIn[117] + InvEtaABC2 * pIn[69];
         pOut[118] = PmQf * pIn[118] + 3 * InvEtaABC2 * pIn[70];
         pOut[119] = PmQf * pIn[119] + 5 * InvEtaABC2 * pIn[71];
         if ( P.AngMomAB_Max == 7 ) break;
      case (2*10 + 8):
         pOut[120] = PmQf * pIn[120];
         pOut[121] = PmQf * pIn[121];
         pOut[122] = PmQf * pIn[122] + 2 * InvEtaABC2 * pIn[104];
         pOut[123] = PmQf * pIn[123];
         pOut[124] = PmQf * pIn[124] + 2 * InvEtaABC2 * pIn[105];
         pOut[125] = PmQf * pIn[125] + 4 * InvEtaABC2 * pIn[106];
         pOut[126] = PmQf * pIn[126];
         pOut[127] = PmQf * pIn[127] + 2 * InvEtaABC2 * pIn[107];
         pOut[128] = PmQf * pIn[128] + 4 * InvEtaABC2 * pIn[108];
         pOut[129] = PmQf * pIn[129] + 6 * InvEtaABC2 * pIn[109];
         pOut[130] = PmQf * pIn[130];
         pOut[131] = PmQf * pIn[131] + 2 * InvEtaABC2 * pIn[110];
         pOut[132] = PmQf * pIn[132] + 4 * InvEtaABC2 * pIn[111];
         pOut[133] = PmQf * pIn[133] + 6 * InvEtaABC2 * pIn[112];
         pOut[134] = PmQf * pIn[134] + 8 * InvEtaABC2 * pIn[113];
         pOut[135] = PmQf * pIn[135];
         pOut[136] = PmQf * pIn[136];
         pOut[137] = PmQf * pIn[137] + 2 * InvEtaABC2 * pIn[114];
         pOut[138] = PmQf * pIn[138];
         pOut[139] = PmQf * pIn[139] + 2 * InvEtaABC2 * pIn[115];
         pOut[140] = PmQf * pIn[140] + 4 * InvEtaABC2 * pIn[116];
         pOut[141] = PmQf * pIn[141];
         pOut[142] = PmQf * pIn[142] + 2 * InvEtaABC2 * pIn[117];
         pOut[143] = PmQf * pIn[143] + 4 * InvEtaABC2 * pIn[118];
         pOut[144] = PmQf * pIn[144] + 6 * InvEtaABC2 * pIn[119];
         pOut[145] = PmQf * pIn[145] + InvEtaABC2 * pIn[84];
         pOut[146] = PmQf * pIn[146] + InvEtaABC2 * pIn[85];
         pOut[147] = PmQf * pIn[147] + 3 * InvEtaABC2 * pIn[86];
         pOut[148] = PmQf * pIn[148] + InvEtaABC2 * pIn[87];
         pOut[149] = PmQf * pIn[149] + 3 * InvEtaABC2 * pIn[88];
         pOut[150] = PmQf * pIn[150] + 5 * InvEtaABC2 * pIn[89];
         pOut[151] = PmQf * pIn[151] + InvEtaABC2 * pIn[90];
         pOut[152] = PmQf * pIn[152] + 3 * InvEtaABC2 * pIn[91];
         pOut[153] = PmQf * pIn[153] + 5 * InvEtaABC2 * pIn[92];
         pOut[154] = PmQf * pIn[154] + 7 * InvEtaABC2 * pIn[93];
         pOut[155] = PmQf * pIn[155] + InvEtaABC2 * pIn[94];
         pOut[156] = PmQf * pIn[156] + InvEtaABC2 * pIn[95];
         pOut[157] = PmQf * pIn[157] + 3 * InvEtaABC2 * pIn[96];
         pOut[158] = PmQf * pIn[158] + InvEtaABC2 * pIn[97];
         pOut[159] = PmQf * pIn[159] + 3 * InvEtaABC2 * pIn[98];
         pOut[160] = PmQf * pIn[160] + 5 * InvEtaABC2 * pIn[99];
         pOut[161] = PmQf * pIn[161] + InvEtaABC2 * pIn[100];
         pOut[162] = PmQf * pIn[162] + 3 * InvEtaABC2 * pIn[101];
         pOut[163] = PmQf * pIn[163] + 5 * InvEtaABC2 * pIn[102];
         pOut[164] = PmQf * pIn[164] + 7 * InvEtaABC2 * pIn[103];
         if ( P.AngMomAB_Max == 8 ) break;
      case (2*10 + 9):
         pOut[165] = PmQf * pIn[165];
         pOut[166] = PmQf * pIn[166];
         pOut[167] = PmQf * pIn[167] + 2 * InvEtaABC2 * pIn[145];
         pOut[168] = PmQf * pIn[168];
         pOut[169] = PmQf * pIn[169] + 2 * InvEtaABC2 * pIn[146];
         pOut[170] = PmQf * pIn[170] + 4 * InvEtaABC2 * pIn[147];
         pOut[171] = PmQf * pIn[171];
         pOut[172] = PmQf * pIn[172] + 2 * InvEtaABC2 * pIn[148];
         pOut[173] = PmQf * pIn[173] + 4 * InvEtaABC2 * pIn[149];
         pOut[174] = PmQf * pIn[174] + 6 * InvEtaABC2 * pIn[150];
         pOut[175] = PmQf * pIn[175];
         pOut[176] = PmQf * pIn[176] + 2 * InvEtaABC2 * pIn[151];
         pOut[177] = PmQf * pIn[177] + 4 * InvEtaABC2 * pIn[152];
         pOut[178] = PmQf * pIn[178] + 6 * InvEtaABC2 * pIn[153];
         pOut[179] = PmQf * pIn[179] + 8 * InvEtaABC2 * pIn[154];
         pOut[180] = PmQf * pIn[180];
         pOut[181] = PmQf * pIn[181];
         pOut[182] = PmQf * pIn[182] + 2 * InvEtaABC2 * pIn[155];
         pOut[183] = PmQf * pIn[183];
         pOut[184] = PmQf * pIn[184] + 2 * InvEtaABC2 * pIn[156];
         pOut[185] = PmQf * pIn[185] + 4 * InvEtaABC2 * pIn[157];
         pOut[186] = PmQf * pIn[186];
         pOut[187] = PmQf * pIn[187] + 2 * InvEtaABC2 * pIn[158];
         pOut[188] = PmQf * pIn[188] + 4 * InvEtaABC2 * pIn[159];
         pOut[189] = PmQf * pIn[189] + 6 * InvEtaABC2 * pIn[160];
         pOut[190] = PmQf * pIn[190];
         pOut[191] = PmQf * pIn[191] + 2 * InvEtaABC2 * pIn[161];
         pOut[192] = PmQf * pIn[192] + 4 * InvEtaABC2 * pIn[162];
         pOut[193] = PmQf * pIn[193] + 6 * InvEtaABC2 * pIn[163];
         pOut[194] = PmQf * pIn[194] + 8 * InvEtaABC2 * pIn[164];
         pOut[195] = PmQf * pIn[195] + InvEtaABC2 * pIn[120];
         pOut[196] = PmQf * pIn[196] + InvEtaABC2 * pIn[121];
         pOut[197] = PmQf * pIn[197] + 3 * InvEtaABC2 * pIn[122];
         pOut[198] = PmQf * pIn[198] + InvEtaABC2 * pIn[123];
         pOut[199] = PmQf * pIn[199] + 3 * InvEtaABC2 * pIn[124];
         pOut[200] = PmQf * pIn[200] + 5 * InvEtaABC2 * pIn[125];
         pOut[201] = PmQf * pIn[201] + InvEtaABC2 * pIn[126];
         pOut[202] = PmQf * pIn[202] + 3 * InvEtaABC2 * pIn[127];
         pOut[203] = PmQf * pIn[203] + 5 * InvEtaABC2 * pIn[128];
         pOut[204] = PmQf * pIn[204] + 7 * InvEtaABC2 * pIn[129];
         pOut[205] = PmQf * pIn[205] + InvEtaABC2 * pIn[130];
         pOut[206] = PmQf * pIn[206] + 3 * InvEtaABC2 * pIn[131];
         pOut[207] = PmQf * pIn[207] + 5 * InvEtaABC2 * pIn[132];
         pOut[208] = PmQf * pIn[208] + 7 * InvEtaABC2 * pIn[133];
         pOut[209] = PmQf * pIn[209] + 9 * InvEtaABC2 * pIn[134];
         pOut[210] = PmQf * pIn[210] + InvEtaABC2 * pIn[135];
         pOut[211] = PmQf * pIn[211] + InvEtaABC2 * pIn[136];
         pOut[212] = PmQf * pIn[212] + 3 * InvEtaABC2 * pIn[137];
         pOut[213] = PmQf * pIn[213] + InvEtaABC2 * pIn[138];
         pOut[214] = PmQf * pIn[214] + 3 * InvEtaABC2 * pIn[139];
         pOut[215] = PmQf * pIn[215] + 5 * InvEtaABC2 * pIn[140];
         pOut[216] = PmQf * pIn[216] + InvEtaABC2 * pIn[141];
         pOut[217] = PmQf * pIn[217] + 3 * InvEtaABC2 * pIn[142];
         pOut[218] = PmQf * pIn[218] + 5 * InvEtaABC2 * pIn[143];
         pOut[219] = PmQf * pIn[219] + 7 * InvEtaABC2 * pIn[144];
         if ( P.AngMomAB_Max == 9 ) break;

      #ifdef _DEBUG
         assert(0);
      #endif
         break;
   }

}

// In: Gm[0 .. 0] (inclusive), Out: [r]^0, ordered as AngularComps(0)
// (moment 0, 1 entries)
static void OsrrA0( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const PmA[3], FReal1 const PmQ[3], double rho, double InvEta ) AIC_NO_THROW
{
   pOut[0] = pIn[0];
   // 0 flops, 2 mops, 0.00kb stack
}


inline void ShTr0( FReal0 pOut[1], FReal0 r000 ) AIC_NO_THROW
{
   pOut[0] = r000;
   // 0 flops, 2 mops
}


static void OsrrC00( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   pOut[0*sa + 0*sb] += pIn[0];
   // 0 flops, 0 mops
}

// In: Gm[0 .. 1] (inclusive), Out: [r]^0, ordered as AngularComps(1)
// (moment 1, 3 entries)
static void OsrrA1( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const PmA[3], FReal1 const PmQ[3], double rho, double InvEta ) AIC_NO_THROW
{
   double
      riz = rho * InvEta,
      PmQf[3] = { riz*PmQ[0], riz*PmQ[1], riz*PmQ[2] };
   pOut[0] = pIn[0];
   // L = 1
   pOut[1] = PmA[0] * pIn[0] - PmQf[0] * pIn[1];
   pOut[2] = PmA[1] * pIn[0] - PmQf[1] * pIn[1];
   pOut[3] = PmA[2] * pIn[0] - PmQf[2] * pIn[1];
   // 9 flops, 8 mops, 0.00kb stack
}




inline void ShTr1( FReal0 pOut[3], FReal0 r100, FReal0 r010, FReal0 r001 ) AIC_NO_THROW
{
   pOut[0] = r100;
   pOut[1] = r010;
   pOut[2] = r001;
   // 0 flops, 6 mops
}


static void OsrrC10( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   pOut[0*sa + 0*sb] += pIn[1];
   pOut[1*sa + 0*sb] += pIn[2];
   pOut[2*sa + 0*sb] += pIn[3];
   // 0 flops, 0 mops
}

// In: Gm[0 .. 2] (inclusive), Out: [r]^0, ordered as AngularComps(2)
// (moment 2, 6 entries)
static void OsrrA2( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const PmA[3], FReal1 const PmQ[3], double rho, double InvEta ) AIC_NO_THROW
{
   double
      riz = rho * InvEta,
      iz2 = 0.5 * InvEta,
      PmQf[3] = { riz*PmQ[0], riz*PmQ[1], riz*PmQ[2] };
   pOut[0] = pIn[0];
   // L = 1
   pOut[1] = PmA[0] * pIn[0] - PmQf[0] * pIn[1];
   FReal2 r_100_1 = PmA[0] * pIn[1] - PmQf[0] * pIn[2];
   pOut[2] = PmA[1] * pIn[0] - PmQf[1] * pIn[1];
   FReal2 r_010_1 = PmA[1] * pIn[1] - PmQf[1] * pIn[2];
   pOut[3] = PmA[2] * pIn[0] - PmQf[2] * pIn[1];
   FReal2 r_001_1 = PmA[2] * pIn[1] - PmQf[2] * pIn[2];

   // L = 2
   pOut[4] = PmA[0] * pOut[1] - PmQf[0] * r_100_1 + iz2 * (pIn[0] - riz * pIn[1]);
   pOut[5] = PmA[1] * pOut[2] - PmQf[1] * r_010_1 + iz2 * (pIn[0] - riz * pIn[1]);
   pOut[6] = PmA[2] * pOut[3] - PmQf[2] * r_001_1 + iz2 * (pIn[0] - riz * pIn[1]);
   pOut[7] = PmA[0] * pOut[2] - PmQf[0] * r_010_1;
   pOut[8] = PmA[2] * pOut[1] - PmQf[2] * r_100_1;
   pOut[9] = PmA[1] * pOut[3] - PmQf[1] * r_001_1;
   // 45 flops, 32 mops, 0.02kb stack
}





static void OsrrC11( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   pOut[0*sa + 0*sb] += pIn[4] - BmA[0] * pIn[1];
   pOut[1*sa + 0*sb] += pIn[7] - BmA[0] * pIn[2];
   pOut[2*sa + 0*sb] += pIn[8] - BmA[0] * pIn[3];
   pOut[0*sa + 1*sb] += pIn[7] - BmA[1] * pIn[1];
   pOut[1*sa + 1*sb] += pIn[5] - BmA[1] * pIn[2];
   pOut[2*sa + 1*sb] += pIn[9] - BmA[1] * pIn[3];
   pOut[0*sa + 2*sb] += pIn[8] - BmA[2] * pIn[1];
   pOut[1*sa + 2*sb] += pIn[9] - BmA[2] * pIn[2];
   pOut[2*sa + 2*sb] += pIn[6] - BmA[2] * pIn[3];
   // 18 flops, 27 mops
}


inline void ShTr2( FReal0 pOut[5], FReal0 r200, FReal0 r020, FReal0 r002, FReal0 r110, FReal0 r101, FReal0 r011 ) AIC_NO_THROW
{
   const double c0 = 5.e-01;
   const double c1 = 1.7320508075688772;
   const double c2 = 8.660254037844386e-01;

   pOut[0] = -c0*r020 + r002 - c0*r200;
   pOut[1] = c1*r110;
   pOut[2] = c1*r101;
   pOut[3] = -c2*r020 + c2*r200;
   pOut[4] = c1*r011;
   // 11 flops, 21 mops
}


static void OsrrC20( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   const double c0 = 5.e-01;
   const double c1 = 1.7320508075688772;
   const double c2 = 8.660254037844386e-01;

   FReal2 a_0_000 = -c0*pIn[5] + pIn[6] - c0*pIn[4];
   FReal2 a_1_000 = c1*pIn[7];
   FReal2 a_2_000 = c1*pIn[8];
   FReal2 a_3_000 = -c2*pIn[5] + c2*pIn[4];
   FReal2 a_4_000 = c1*pIn[9];

   pOut[0*sa + 0*sb] += a_0_000;
   pOut[1*sa + 0*sb] += a_1_000;
   pOut[2*sa + 0*sb] += a_2_000;
   pOut[3*sa + 0*sb] += a_3_000;
   pOut[4*sa + 0*sb] += a_4_000;
   // 0 flops, 0 mops
}

// In: Gm[0 .. 3] (inclusive), Out: [r]^0, ordered as AngularComps(3)
// (moment 3, 10 entries)
static void OsrrA3( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const PmA[3], FReal1 const PmQ[3], double rho, double InvEta ) AIC_NO_THROW
{
   double
      riz = rho * InvEta,
      iz2 = 0.5 * InvEta,
      PmQf[3] = { riz*PmQ[0], riz*PmQ[1], riz*PmQ[2] };
   pOut[0] = pIn[0];
   // L = 1
   pOut[1] = PmA[0] * pIn[0] - PmQf[0] * pIn[1];
   FReal2 r_100_1 = PmA[0] * pIn[1] - PmQf[0] * pIn[2];
   FReal2 r_100_2 = PmA[0] * pIn[2] - PmQf[0] * pIn[3];
   pOut[2] = PmA[1] * pIn[0] - PmQf[1] * pIn[1];
   FReal2 r_010_1 = PmA[1] * pIn[1] - PmQf[1] * pIn[2];
   FReal2 r_010_2 = PmA[1] * pIn[2] - PmQf[1] * pIn[3];
   pOut[3] = PmA[2] * pIn[0] - PmQf[2] * pIn[1];
   FReal2 r_001_1 = PmA[2] * pIn[1] - PmQf[2] * pIn[2];
   FReal2 r_001_2 = PmA[2] * pIn[2] - PmQf[2] * pIn[3];

   // L = 2
   pOut[4] = PmA[0] * pOut[1] - PmQf[0] * r_100_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_200_1 = PmA[0] * r_100_1 - PmQf[0] * r_100_2 + iz2 * (pIn[1] - riz * pIn[2]);
   pOut[5] = PmA[1] * pOut[2] - PmQf[1] * r_010_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_020_1 = PmA[1] * r_010_1 - PmQf[1] * r_010_2 + iz2 * (pIn[1] - riz * pIn[2]);
   pOut[6] = PmA[2] * pOut[3] - PmQf[2] * r_001_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_002_1 = PmA[2] * r_001_1 - PmQf[2] * r_001_2 + iz2 * (pIn[1] - riz * pIn[2]);
   pOut[7] = PmA[0] * pOut[2] - PmQf[0] * r_010_1;
   pOut[8] = PmA[2] * pOut[1] - PmQf[2] * r_100_1;
   FReal2 r_101_1 = PmA[2] * r_100_1 - PmQf[2] * r_100_2;
   pOut[9] = PmA[1] * pOut[3] - PmQf[1] * r_001_1;

   // L = 3
   pOut[10] = PmA[0] * pOut[4] - PmQf[0] * r_200_1 + 2 * iz2 * (pOut[1] - riz * r_100_1);
   pOut[11] = PmA[0] * pOut[5] - PmQf[0] * r_020_1;
   pOut[12] = PmA[0] * pOut[6] - PmQf[0] * r_002_1;
   pOut[13] = PmA[1] * pOut[4] - PmQf[1] * r_200_1;
   pOut[14] = PmA[1] * pOut[5] - PmQf[1] * r_020_1 + 2 * iz2 * (pOut[2] - riz * r_010_1);
   pOut[15] = PmA[1] * pOut[6] - PmQf[1] * r_002_1;
   pOut[16] = PmA[0] * pOut[8] - PmQf[0] * r_101_1 + iz2 * (pOut[3] - riz * r_001_1);
   pOut[17] = PmA[2] * pOut[5] - PmQf[2] * r_020_1;
   pOut[18] = PmA[2] * pOut[6] - PmQf[2] * r_002_1 + 2 * iz2 * (pOut[3] - riz * r_001_1);
   pOut[19] = PmA[1] * pOut[8] - PmQf[1] * r_101_1;
   // 120 flops, 80 mops, 0.08kb stack
}







static void OsrrC21( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   const double c0 = 5.e-01;
   const double c1 = 1.7320508075688772;
   const double c2 = 8.660254037844386e-01;

   FReal2 a_0_000 = -c0*pIn[5] + pIn[6] - c0*pIn[4];
   FReal2 a_1_000 = c1*pIn[7];
   FReal2 a_2_000 = c1*pIn[8];
   FReal2 a_3_000 = -c2*pIn[5] + c2*pIn[4];
   FReal2 a_4_000 = c1*pIn[9];
   FReal2 a_0_100 = -c0*pIn[11] + pIn[12] - c0*pIn[10];
   FReal2 a_1_100 = c1*pIn[13];
   FReal2 a_2_100 = c1*pIn[16];
   FReal2 a_3_100 = -c2*pIn[11] + c2*pIn[10];
   FReal2 a_4_100 = c1*pIn[19];
   FReal2 a_0_010 = -c0*pIn[14] + pIn[15] - c0*pIn[13];
   FReal2 a_1_010 = c1*pIn[11];
   FReal2 a_2_010 = c1*pIn[19];
   FReal2 a_3_010 = -c2*pIn[14] + c2*pIn[13];
   FReal2 a_4_010 = c1*pIn[17];
   FReal2 a_0_001 = -c0*pIn[17] + pIn[18] - c0*pIn[16];
   FReal2 a_1_001 = c1*pIn[19];
   FReal2 a_2_001 = c1*pIn[12];
   FReal2 a_3_001 = -c2*pIn[17] + c2*pIn[16];
   FReal2 a_4_001 = c1*pIn[15];

   pOut[0*sa + 0*sb] += a_0_100 - BmA[0] * a_0_000;
   pOut[1*sa + 0*sb] += a_1_100 - BmA[0] * a_1_000;
   pOut[2*sa + 0*sb] += a_2_100 - BmA[0] * a_2_000;
   pOut[3*sa + 0*sb] += a_3_100 - BmA[0] * a_3_000;
   pOut[4*sa + 0*sb] += a_4_100 - BmA[0] * a_4_000;
   pOut[0*sa + 1*sb] += a_0_010 - BmA[1] * a_0_000;
   pOut[1*sa + 1*sb] += a_1_010 - BmA[1] * a_1_000;
   pOut[2*sa + 1*sb] += a_2_010 - BmA[1] * a_2_000;
   pOut[3*sa + 1*sb] += a_3_010 - BmA[1] * a_3_000;
   pOut[4*sa + 1*sb] += a_4_010 - BmA[1] * a_4_000;
   pOut[0*sa + 2*sb] += a_0_001 - BmA[2] * a_0_000;
   pOut[1*sa + 2*sb] += a_1_001 - BmA[2] * a_1_000;
   pOut[2*sa + 2*sb] += a_2_001 - BmA[2] * a_2_000;
   pOut[3*sa + 2*sb] += a_3_001 - BmA[2] * a_3_000;
   pOut[4*sa + 2*sb] += a_4_001 - BmA[2] * a_4_000;
   // 30 flops, 45 mops
}


static void ShTr3( FReal0 pOut[7], FReal0 r300, FReal0 r120, FReal0 r102, FReal0 r210, FReal0 r030, FReal0 r012, FReal0 r201, FReal0 r021, FReal0 r003, FReal0 r111 ) AIC_NO_THROW
{
   const double c0 = 6.1237243569579447e-01;
   const double c1 = 2.4494897427831779;
   const double c2 = 1.5;
   const double c3 = 2.3717082451262841;
   const double c4 = 7.9056941504209477e-01;
   const double c5 = 3.8729833462074166;
   const double c6 = 1.9364916731037083;

   pOut[0] = -c0*r300 + c1*r102 - c0*r120;
   pOut[1] = c1*r012 - c0*r030 - c0*r210;
   pOut[2] = -c2*r021 - c2*r201 + r003;
   pOut[3] = -c3*r120 + c4*r300;
   pOut[4] = c5*r111;
   pOut[5] = -c4*r030 + c3*r210;
   pOut[6] = -c6*r021 + c6*r201;
   // 25 flops, 39 mops
}


static void OsrrC30( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   FReal2 a_000[7]; ShTr3(a_000, pIn[10], pIn[11], pIn[12], pIn[13], pIn[14], pIn[15], pIn[16], pIn[17], pIn[18], pIn[19]);
   // ShTr: [25 flops, 49 mops]

   for ( unsigned ica = 0; ica < 7; ++ ica ) {
      pOut[ica*sa + 0*sb] += a_000[ica];
      // 7 * [0 flops, 0 mops]
   }
   // total: [25 flops, 49 mops]
}

// In: Gm[0 .. 4] (inclusive), Out: [r]^0, ordered as AngularComps(4)
// (moment 4, 15 entries)
static void OsrrA4( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const PmA[3], FReal1 const PmQ[3], double rho, double InvEta ) AIC_NO_THROW
{
   double
      riz = rho * InvEta,
      iz2 = 0.5 * InvEta,
      PmQf[3] = { riz*PmQ[0], riz*PmQ[1], riz*PmQ[2] };
   pOut[0] = pIn[0];
   // L = 1
   pOut[1] = PmA[0] * pIn[0] - PmQf[0] * pIn[1];
   FReal2 r_100_1 = PmA[0] * pIn[1] - PmQf[0] * pIn[2];
   FReal2 r_100_2 = PmA[0] * pIn[2] - PmQf[0] * pIn[3];
   FReal2 r_100_3 = PmA[0] * pIn[3] - PmQf[0] * pIn[4];
   pOut[2] = PmA[1] * pIn[0] - PmQf[1] * pIn[1];
   FReal2 r_010_1 = PmA[1] * pIn[1] - PmQf[1] * pIn[2];
   FReal2 r_010_2 = PmA[1] * pIn[2] - PmQf[1] * pIn[3];
   FReal2 r_010_3 = PmA[1] * pIn[3] - PmQf[1] * pIn[4];
   pOut[3] = PmA[2] * pIn[0] - PmQf[2] * pIn[1];
   FReal2 r_001_1 = PmA[2] * pIn[1] - PmQf[2] * pIn[2];
   FReal2 r_001_2 = PmA[2] * pIn[2] - PmQf[2] * pIn[3];
   FReal2 r_001_3 = PmA[2] * pIn[3] - PmQf[2] * pIn[4];

   // L = 2
   pOut[4] = PmA[0] * pOut[1] - PmQf[0] * r_100_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_200_1 = PmA[0] * r_100_1 - PmQf[0] * r_100_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_200_2 = PmA[0] * r_100_2 - PmQf[0] * r_100_3 + iz2 * (pIn[2] - riz * pIn[3]);
   pOut[5] = PmA[1] * pOut[2] - PmQf[1] * r_010_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_020_1 = PmA[1] * r_010_1 - PmQf[1] * r_010_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_020_2 = PmA[1] * r_010_2 - PmQf[1] * r_010_3 + iz2 * (pIn[2] - riz * pIn[3]);
   pOut[6] = PmA[2] * pOut[3] - PmQf[2] * r_001_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_002_1 = PmA[2] * r_001_1 - PmQf[2] * r_001_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_002_2 = PmA[2] * r_001_2 - PmQf[2] * r_001_3 + iz2 * (pIn[2] - riz * pIn[3]);
   pOut[7] = PmA[0] * pOut[2] - PmQf[0] * r_010_1;
   pOut[8] = PmA[2] * pOut[1] - PmQf[2] * r_100_1;
   FReal2 r_101_1 = PmA[2] * r_100_1 - PmQf[2] * r_100_2;
   pOut[9] = PmA[1] * pOut[3] - PmQf[1] * r_001_1;

   // L = 3
   pOut[10] = PmA[0] * pOut[4] - PmQf[0] * r_200_1 + 2 * iz2 * (pOut[1] - riz * r_100_1);
   FReal2 r_300_1 = PmA[0] * r_200_1 - PmQf[0] * r_200_2 + 2 * iz2 * (r_100_1 - riz * r_100_2);
   pOut[11] = PmA[0] * pOut[5] - PmQf[0] * r_020_1;
   FReal2 r_120_1 = PmA[0] * r_020_1 - PmQf[0] * r_020_2;
   pOut[12] = PmA[0] * pOut[6] - PmQf[0] * r_002_1;
   pOut[13] = PmA[1] * pOut[4] - PmQf[1] * r_200_1;
   pOut[14] = PmA[1] * pOut[5] - PmQf[1] * r_020_1 + 2 * iz2 * (pOut[2] - riz * r_010_1);
   FReal2 r_030_1 = PmA[1] * r_020_1 - PmQf[1] * r_020_2 + 2 * iz2 * (r_010_1 - riz * r_010_2);
   pOut[15] = PmA[1] * pOut[6] - PmQf[1] * r_002_1;
   FReal2 r_012_1 = PmA[1] * r_002_1 - PmQf[1] * r_002_2;
   pOut[16] = PmA[0] * pOut[8] - PmQf[0] * r_101_1 + iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_201_1 = PmA[2] * r_200_1 - PmQf[2] * r_200_2;
   pOut[17] = PmA[2] * pOut[5] - PmQf[2] * r_020_1;
   pOut[18] = PmA[2] * pOut[6] - PmQf[2] * r_002_1 + 2 * iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_003_1 = PmA[2] * r_002_1 - PmQf[2] * r_002_2 + 2 * iz2 * (r_001_1 - riz * r_001_2);
   pOut[19] = PmA[1] * pOut[8] - PmQf[1] * r_101_1;

   // L = 4
   pOut[20] = PmA[0] * pOut[10] - PmQf[0] * r_300_1 + 3 * iz2 * (pOut[4] - riz * r_200_1);
   pOut[21] = PmA[0] * pOut[11] - PmQf[0] * r_120_1 + iz2 * (pOut[5] - riz * r_020_1);
   pOut[22] = PmA[2] * pOut[16] - PmQf[2] * r_201_1 + iz2 * (pOut[4] - riz * r_200_1);
   pOut[23] = PmA[1] * pOut[14] - PmQf[1] * r_030_1 + 3 * iz2 * (pOut[5] - riz * r_020_1);
   pOut[24] = PmA[1] * pOut[15] - PmQf[1] * r_012_1 + iz2 * (pOut[6] - riz * r_002_1);
   pOut[25] = PmA[2] * pOut[18] - PmQf[2] * r_003_1 + 3 * iz2 * (pOut[6] - riz * r_002_1);
   pOut[26] = PmA[1] * pOut[10] - PmQf[1] * r_300_1;
   pOut[27] = PmA[0] * pOut[14] - PmQf[0] * r_030_1;
   pOut[28] = PmA[0] * pOut[15] - PmQf[0] * r_012_1;
   pOut[29] = PmA[0] * pOut[16] - PmQf[0] * r_201_1 + 2 * iz2 * (pOut[8] - riz * r_101_1);
   pOut[30] = PmA[2] * pOut[11] - PmQf[2] * r_120_1;
   pOut[31] = PmA[0] * pOut[18] - PmQf[0] * r_003_1;
   pOut[32] = PmA[1] * pOut[16] - PmQf[1] * r_201_1;
   pOut[33] = PmA[2] * pOut[14] - PmQf[2] * r_030_1;
   pOut[34] = PmA[1] * pOut[18] - PmQf[1] * r_003_1;
   // 247 flops, 160 mops, 0.17kb stack
}







static void OsrrC22( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   const double c0 = 5.e-01;
   const double c1 = 1.7320508075688772;
   const double c2 = 8.660254037844386e-01;

   FReal2 a_000[5]; ShTr2(a_000, pIn[4], pIn[5], pIn[6], pIn[7], pIn[8], pIn[9]);
   FReal2 a_100[5]; ShTr2(a_100, pIn[10], pIn[11], pIn[12], pIn[13], pIn[16], pIn[19]);
   FReal2 a_010[5]; ShTr2(a_010, pIn[13], pIn[14], pIn[15], pIn[11], pIn[19], pIn[17]);
   FReal2 a_001[5]; ShTr2(a_001, pIn[16], pIn[17], pIn[18], pIn[19], pIn[12], pIn[15]);
   FReal2 a_200[5]; ShTr2(a_200, pIn[20], pIn[21], pIn[22], pIn[26], pIn[29], pIn[32]);
   FReal2 a_020[5]; ShTr2(a_020, pIn[21], pIn[23], pIn[24], pIn[27], pIn[30], pIn[33]);
   FReal2 a_002[5]; ShTr2(a_002, pIn[22], pIn[24], pIn[25], pIn[28], pIn[31], pIn[34]);
   FReal2 a_110[5]; ShTr2(a_110, pIn[26], pIn[27], pIn[28], pIn[21], pIn[32], pIn[30]);
   FReal2 a_101[5]; ShTr2(a_101, pIn[29], pIn[30], pIn[31], pIn[32], pIn[22], pIn[28]);
   FReal2 a_011[5]; ShTr2(a_011, pIn[32], pIn[33], pIn[34], pIn[30], pIn[28], pIn[24]);
   // ShTr: [110 flops, 270 mops]

   for ( unsigned ica = 0; ica < 5; ++ ica ) {
      FReal2 r_000_100 = a_100[ica] - BmA[0] * a_000[ica];
      FReal2 r_100_100 = a_200[ica] - BmA[0] * a_100[ica];
      FReal2 r_001_100 = a_101[ica] - BmA[0] * a_001[ica];
      FReal2 r_000_010 = a_010[ica] - BmA[1] * a_000[ica];
      FReal2 r_100_010 = a_110[ica] - BmA[1] * a_100[ica];
      FReal2 r_010_010 = a_020[ica] - BmA[1] * a_010[ica];
      FReal2 r_000_001 = a_001[ica] - BmA[2] * a_000[ica];
      FReal2 r_010_001 = a_011[ica] - BmA[2] * a_010[ica];
      FReal2 r_001_001 = a_002[ica] - BmA[2] * a_001[ica];
      FReal2 b_200 = r_100_100 - BmA[0] * r_000_100;
      FReal2 b_020 = r_010_010 - BmA[1] * r_000_010;
      FReal2 b_002 = r_001_001 - BmA[2] * r_000_001;
      FReal2 b_110 = r_100_010 - BmA[0] * r_000_010;
      FReal2 b_101 = r_001_100 - BmA[2] * r_000_100;
      FReal2 b_011 = r_010_001 - BmA[1] * r_000_001;
      pOut[ica*sa + 0*sb] += -c0*b_020 + b_002 - c0*b_200;
      pOut[ica*sa + 1*sb] += c1*b_110;
      pOut[ica*sa + 2*sb] += c1*b_101;
      pOut[ica*sa + 3*sb] += -c2*b_020 + c2*b_200;
      pOut[ica*sa + 4*sb] += c1*b_011;
      // 5 * [30 flops, 45 mops]
   }
   // total: [260 flops, 495 mops]
}



static void OsrrC31( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   FReal2 a_000[7]; ShTr3(a_000, pIn[10], pIn[11], pIn[12], pIn[13], pIn[14], pIn[15], pIn[16], pIn[17], pIn[18], pIn[19]);
   FReal2 a_100[7]; ShTr3(a_100, pIn[20], pIn[21], pIn[22], pIn[26], pIn[27], pIn[28], pIn[29], pIn[30], pIn[31], pIn[32]);
   FReal2 a_010[7]; ShTr3(a_010, pIn[26], pIn[27], pIn[28], pIn[21], pIn[23], pIn[24], pIn[32], pIn[33], pIn[34], pIn[30]);
   FReal2 a_001[7]; ShTr3(a_001, pIn[29], pIn[30], pIn[31], pIn[32], pIn[33], pIn[34], pIn[22], pIn[24], pIn[25], pIn[28]);
   // ShTr: [100 flops, 196 mops]

   for ( unsigned ica = 0; ica < 7; ++ ica ) {
      FReal2 b_100 = a_100[ica] - BmA[0] * a_000[ica];
      FReal2 b_010 = a_010[ica] - BmA[1] * a_000[ica];
      FReal2 b_001 = a_001[ica] - BmA[2] * a_000[ica];
      pOut[ica*sa + 0*sb] += b_100;
      pOut[ica*sa + 1*sb] += b_010;
      pOut[ica*sa + 2*sb] += b_001;
      // 7 * [6 flops, 9 mops]
   }
   // total: [142 flops, 259 mops]
}


static void ShTr4( FReal0 pOut[9], FReal0 r400, FReal0 r220, FReal0 r202, FReal0 r040, FReal0 r022, FReal0 r004, FReal0 r310, FReal0 r130, FReal0 r112, FReal0 r301, FReal0 r121, FReal0 r103, FReal0 r211, FReal0 r031, FReal0 r013 ) AIC_NO_THROW
{
   const double c0 = 3.;
   const double c1 = 3.75e-01;
   const double c2 = 7.5e-01;
   const double c3 = 1.1180339887498947;
   const double c4 = 6.7082039324993685;
   const double c5 = 2.3717082451262841;
   const double c6 = 3.1622776601683791;
   const double c7 = 7.3950997288745191e-01;
   const double c8 = 4.4370598373247114;
   const double c9 = 3.3541019662496843;
   const double ca = 5.5901699437494734e-01;
   const double cb = 2.9580398915498076;
   const double cc = 6.2749501990055663;
   const double cd = 2.0916500663351885;

   pOut[0] = -c0*r022 + c1*r400 - c0*r202 + c2*r220 + r004 + c1*r040;
   pOut[1] = -c3*r310 - c3*r130 + c4*r112;
   pOut[2] = -c5*r301 - c5*r121 + c6*r103;
   pOut[3] = c7*r400 + c7*r040 - c8*r220;
   pOut[4] = -c5*r031 - c5*r211 + c6*r013;
   pOut[5] = c9*r202 - ca*r400 + ca*r040 - c9*r022;
   pOut[6] = cb*r310 - cb*r130;
   pOut[7] = -cc*r121 + cd*r301;
   pOut[8] = cc*r211 - cd*r031;
   // 47 flops, 65 mops
}


static void OsrrC40( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   FReal2 a_000[9]; ShTr4(a_000, pIn[20], pIn[21], pIn[22], pIn[23], pIn[24], pIn[25], pIn[26], pIn[27], pIn[28], pIn[29], pIn[30], pIn[31], pIn[32], pIn[33], pIn[34]);
   // ShTr: [47 flops, 80 mops]

   for ( unsigned ica = 0; ica < 9; ++ ica ) {
      pOut[ica*sa + 0*sb] += a_000[ica];
      // 9 * [0 flops, 0 mops]
   }
   // total: [47 flops, 80 mops]
}

// In: Gm[0 .. 5] (inclusive), Out: [r]^0, ordered as AngularComps(5)
// (moment 5, 21 entries)
static void OsrrA5( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const PmA[3], FReal1 const PmQ[3], double rho, double InvEta ) AIC_NO_THROW
{
   double
      riz = rho * InvEta,
      iz2 = 0.5 * InvEta,
      PmQf[3] = { riz*PmQ[0], riz*PmQ[1], riz*PmQ[2] };
   pOut[0] = pIn[0];
   // L = 1
   pOut[1] = PmA[0] * pIn[0] - PmQf[0] * pIn[1];
   FReal2 r_100_1 = PmA[0] * pIn[1] - PmQf[0] * pIn[2];
   FReal2 r_100_2 = PmA[0] * pIn[2] - PmQf[0] * pIn[3];
   FReal2 r_100_3 = PmA[0] * pIn[3] - PmQf[0] * pIn[4];
   FReal2 r_100_4 = PmA[0] * pIn[4] - PmQf[0] * pIn[5];
   pOut[2] = PmA[1] * pIn[0] - PmQf[1] * pIn[1];
   FReal2 r_010_1 = PmA[1] * pIn[1] - PmQf[1] * pIn[2];
   FReal2 r_010_2 = PmA[1] * pIn[2] - PmQf[1] * pIn[3];
   FReal2 r_010_3 = PmA[1] * pIn[3] - PmQf[1] * pIn[4];
   FReal2 r_010_4 = PmA[1] * pIn[4] - PmQf[1] * pIn[5];
   pOut[3] = PmA[2] * pIn[0] - PmQf[2] * pIn[1];
   FReal2 r_001_1 = PmA[2] * pIn[1] - PmQf[2] * pIn[2];
   FReal2 r_001_2 = PmA[2] * pIn[2] - PmQf[2] * pIn[3];
   FReal2 r_001_3 = PmA[2] * pIn[3] - PmQf[2] * pIn[4];
   FReal2 r_001_4 = PmA[2] * pIn[4] - PmQf[2] * pIn[5];

   // L = 2
   pOut[4] = PmA[0] * pOut[1] - PmQf[0] * r_100_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_200_1 = PmA[0] * r_100_1 - PmQf[0] * r_100_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_200_2 = PmA[0] * r_100_2 - PmQf[0] * r_100_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_200_3 = PmA[0] * r_100_3 - PmQf[0] * r_100_4 + iz2 * (pIn[3] - riz * pIn[4]);
   pOut[5] = PmA[1] * pOut[2] - PmQf[1] * r_010_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_020_1 = PmA[1] * r_010_1 - PmQf[1] * r_010_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_020_2 = PmA[1] * r_010_2 - PmQf[1] * r_010_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_020_3 = PmA[1] * r_010_3 - PmQf[1] * r_010_4 + iz2 * (pIn[3] - riz * pIn[4]);
   pOut[6] = PmA[2] * pOut[3] - PmQf[2] * r_001_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_002_1 = PmA[2] * r_001_1 - PmQf[2] * r_001_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_002_2 = PmA[2] * r_001_2 - PmQf[2] * r_001_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_002_3 = PmA[2] * r_001_3 - PmQf[2] * r_001_4 + iz2 * (pIn[3] - riz * pIn[4]);
   pOut[7] = PmA[0] * pOut[2] - PmQf[0] * r_010_1;
   pOut[8] = PmA[2] * pOut[1] - PmQf[2] * r_100_1;
   FReal2 r_101_1 = PmA[2] * r_100_1 - PmQf[2] * r_100_2;
   pOut[9] = PmA[1] * pOut[3] - PmQf[1] * r_001_1;

   // L = 3
   pOut[10] = PmA[0] * pOut[4] - PmQf[0] * r_200_1 + 2 * iz2 * (pOut[1] - riz * r_100_1);
   FReal2 r_300_1 = PmA[0] * r_200_1 - PmQf[0] * r_200_2 + 2 * iz2 * (r_100_1 - riz * r_100_2);
   FReal2 r_300_2 = PmA[0] * r_200_2 - PmQf[0] * r_200_3 + 2 * iz2 * (r_100_2 - riz * r_100_3);
   pOut[11] = PmA[0] * pOut[5] - PmQf[0] * r_020_1;
   FReal2 r_120_1 = PmA[0] * r_020_1 - PmQf[0] * r_020_2;
   FReal2 r_120_2 = PmA[0] * r_020_2 - PmQf[0] * r_020_3;
   pOut[12] = PmA[0] * pOut[6] - PmQf[0] * r_002_1;
   pOut[13] = PmA[1] * pOut[4] - PmQf[1] * r_200_1;
   pOut[14] = PmA[1] * pOut[5] - PmQf[1] * r_020_1 + 2 * iz2 * (pOut[2] - riz * r_010_1);
   FReal2 r_030_1 = PmA[1] * r_020_1 - PmQf[1] * r_020_2 + 2 * iz2 * (r_010_1 - riz * r_010_2);
   FReal2 r_030_2 = PmA[1] * r_020_2 - PmQf[1] * r_020_3 + 2 * iz2 * (r_010_2 - riz * r_010_3);
   pOut[15] = PmA[1] * pOut[6] - PmQf[1] * r_002_1;
   FReal2 r_012_1 = PmA[1] * r_002_1 - PmQf[1] * r_002_2;
   FReal2 r_012_2 = PmA[1] * r_002_2 - PmQf[1] * r_002_3;
   pOut[16] = PmA[0] * pOut[8] - PmQf[0] * r_101_1 + iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_201_1 = PmA[2] * r_200_1 - PmQf[2] * r_200_2;
   FReal2 r_201_2 = PmA[2] * r_200_2 - PmQf[2] * r_200_3;
   pOut[17] = PmA[2] * pOut[5] - PmQf[2] * r_020_1;
   pOut[18] = PmA[2] * pOut[6] - PmQf[2] * r_002_1 + 2 * iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_003_1 = PmA[2] * r_002_1 - PmQf[2] * r_002_2 + 2 * iz2 * (r_001_1 - riz * r_001_2);
   FReal2 r_003_2 = PmA[2] * r_002_2 - PmQf[2] * r_002_3 + 2 * iz2 * (r_001_2 - riz * r_001_3);
   pOut[19] = PmA[1] * pOut[8] - PmQf[1] * r_101_1;

   // L = 4
   pOut[20] = PmA[0] * pOut[10] - PmQf[0] * r_300_1 + 3 * iz2 * (pOut[4] - riz * r_200_1);
   FReal2 r_400_1 = PmA[0] * r_300_1 - PmQf[0] * r_300_2 + 3 * iz2 * (r_200_1 - riz * r_200_2);
   pOut[21] = PmA[0] * pOut[11] - PmQf[0] * r_120_1 + iz2 * (pOut[5] - riz * r_020_1);
   FReal2 r_220_1 = PmA[0] * r_120_1 - PmQf[0] * r_120_2 + iz2 * (r_020_1 - riz * r_020_2);
   pOut[22] = PmA[2] * pOut[16] - PmQf[2] * r_201_1 + iz2 * (pOut[4] - riz * r_200_1);
   FReal2 r_202_1 = PmA[2] * r_201_1 - PmQf[2] * r_201_2 + iz2 * (r_200_1 - riz * r_200_2);
   pOut[23] = PmA[1] * pOut[14] - PmQf[1] * r_030_1 + 3 * iz2 * (pOut[5] - riz * r_020_1);
   FReal2 r_040_1 = PmA[1] * r_030_1 - PmQf[1] * r_030_2 + 3 * iz2 * (r_020_1 - riz * r_020_2);
   pOut[24] = PmA[1] * pOut[15] - PmQf[1] * r_012_1 + iz2 * (pOut[6] - riz * r_002_1);
   FReal2 r_022_1 = PmA[1] * r_012_1 - PmQf[1] * r_012_2 + iz2 * (r_002_1 - riz * r_002_2);
   pOut[25] = PmA[2] * pOut[18] - PmQf[2] * r_003_1 + 3 * iz2 * (pOut[6] - riz * r_002_1);
   FReal2 r_004_1 = PmA[2] * r_003_1 - PmQf[2] * r_003_2 + 3 * iz2 * (r_002_1 - riz * r_002_2);
   pOut[26] = PmA[1] * pOut[10] - PmQf[1] * r_300_1;
   pOut[27] = PmA[0] * pOut[14] - PmQf[0] * r_030_1;
   FReal2 r_130_1 = PmA[0] * r_030_1 - PmQf[0] * r_030_2;
   pOut[28] = PmA[0] * pOut[15] - PmQf[0] * r_012_1;
   pOut[29] = PmA[0] * pOut[16] - PmQf[0] * r_201_1 + 2 * iz2 * (pOut[8] - riz * r_101_1);
   FReal2 r_301_1 = PmA[2] * r_300_1 - PmQf[2] * r_300_2;
   pOut[30] = PmA[2] * pOut[11] - PmQf[2] * r_120_1;
   pOut[31] = PmA[0] * pOut[18] - PmQf[0] * r_003_1;
   pOut[32] = PmA[1] * pOut[16] - PmQf[1] * r_201_1;
   pOut[33] = PmA[2] * pOut[14] - PmQf[2] * r_030_1;
   pOut[34] = PmA[1] * pOut[18] - PmQf[1] * r_003_1;
   FReal2 r_013_1 = PmA[1] * r_003_1 - PmQf[1] * r_003_2;

   // L = 5
   pOut[35] = PmA[0] * pOut[20] - PmQf[0] * r_400_1 + 4 * iz2 * (pOut[10] - riz * r_300_1);
   pOut[36] = PmA[0] * pOut[21] - PmQf[0] * r_220_1 + 2 * iz2 * (pOut[11] - riz * r_120_1);
   pOut[37] = PmA[2] * pOut[29] - PmQf[2] * r_301_1 + iz2 * (pOut[10] - riz * r_300_1);
   pOut[38] = PmA[0] * pOut[23] - PmQf[0] * r_040_1;
   pOut[39] = PmA[0] * pOut[24] - PmQf[0] * r_022_1;
   pOut[40] = PmA[0] * pOut[25] - PmQf[0] * r_004_1;
   pOut[41] = PmA[1] * pOut[20] - PmQf[1] * r_400_1;
   pOut[42] = PmA[0] * pOut[27] - PmQf[0] * r_130_1 + iz2 * (pOut[14] - riz * r_030_1);
   pOut[43] = PmA[1] * pOut[22] - PmQf[1] * r_202_1;
   pOut[44] = PmA[1] * pOut[23] - PmQf[1] * r_040_1 + 4 * iz2 * (pOut[14] - riz * r_030_1);
   pOut[45] = PmA[1] * pOut[24] - PmQf[1] * r_022_1 + 2 * iz2 * (pOut[15] - riz * r_012_1);
   pOut[46] = PmA[2] * pOut[34] - PmQf[2] * r_013_1 + 3 * iz2 * (pOut[15] - riz * r_012_1);
   pOut[47] = PmA[0] * pOut[29] - PmQf[0] * r_301_1 + 3 * iz2 * (pOut[16] - riz * r_201_1);
   pOut[48] = PmA[2] * pOut[21] - PmQf[2] * r_220_1;
   pOut[49] = PmA[2] * pOut[22] - PmQf[2] * r_202_1 + 2 * iz2 * (pOut[16] - riz * r_201_1);
   pOut[50] = PmA[2] * pOut[23] - PmQf[2] * r_040_1;
   pOut[51] = PmA[1] * pOut[34] - PmQf[1] * r_013_1 + iz2 * (pOut[18] - riz * r_003_1);
   pOut[52] = PmA[2] * pOut[25] - PmQf[2] * r_004_1 + 4 * iz2 * (pOut[18] - riz * r_003_1);
   pOut[53] = PmA[1] * pOut[29] - PmQf[1] * r_301_1;
   pOut[54] = PmA[2] * pOut[27] - PmQf[2] * r_130_1;
   pOut[55] = PmA[0] * pOut[34] - PmQf[0] * r_013_1;
   // 456 flops, 290 mops, 0.34kb stack
}









static void OsrrC32( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   const double c0 = 5.e-01;
   const double c1 = 1.7320508075688772;
   const double c2 = 8.660254037844386e-01;

   FReal2 a_000[7]; ShTr3(a_000, pIn[10], pIn[11], pIn[12], pIn[13], pIn[14], pIn[15], pIn[16], pIn[17], pIn[18], pIn[19]);
   FReal2 a_100[7]; ShTr3(a_100, pIn[20], pIn[21], pIn[22], pIn[26], pIn[27], pIn[28], pIn[29], pIn[30], pIn[31], pIn[32]);
   FReal2 a_010[7]; ShTr3(a_010, pIn[26], pIn[27], pIn[28], pIn[21], pIn[23], pIn[24], pIn[32], pIn[33], pIn[34], pIn[30]);
   FReal2 a_001[7]; ShTr3(a_001, pIn[29], pIn[30], pIn[31], pIn[32], pIn[33], pIn[34], pIn[22], pIn[24], pIn[25], pIn[28]);
   FReal2 a_200[7]; ShTr3(a_200, pIn[35], pIn[36], pIn[37], pIn[41], pIn[42], pIn[43], pIn[47], pIn[48], pIn[49], pIn[53]);
   FReal2 a_020[7]; ShTr3(a_020, pIn[36], pIn[38], pIn[39], pIn[42], pIn[44], pIn[45], pIn[48], pIn[50], pIn[51], pIn[54]);
   FReal2 a_002[7]; ShTr3(a_002, pIn[37], pIn[39], pIn[40], pIn[43], pIn[45], pIn[46], pIn[49], pIn[51], pIn[52], pIn[55]);
   FReal2 a_110[7]; ShTr3(a_110, pIn[41], pIn[42], pIn[43], pIn[36], pIn[38], pIn[39], pIn[53], pIn[54], pIn[55], pIn[48]);
   FReal2 a_101[7]; ShTr3(a_101, pIn[47], pIn[48], pIn[49], pIn[53], pIn[54], pIn[55], pIn[37], pIn[39], pIn[40], pIn[43]);
   FReal2 a_011[7]; ShTr3(a_011, pIn[53], pIn[54], pIn[55], pIn[48], pIn[50], pIn[51], pIn[43], pIn[45], pIn[46], pIn[39]);
   // ShTr: [250 flops, 490 mops]

   for ( unsigned ica = 0; ica < 7; ++ ica ) {
      FReal2 r_000_100 = a_100[ica] - BmA[0] * a_000[ica];
      FReal2 r_100_100 = a_200[ica] - BmA[0] * a_100[ica];
      FReal2 r_001_100 = a_101[ica] - BmA[0] * a_001[ica];
      FReal2 r_000_010 = a_010[ica] - BmA[1] * a_000[ica];
      FReal2 r_100_010 = a_110[ica] - BmA[1] * a_100[ica];
      FReal2 r_010_010 = a_020[ica] - BmA[1] * a_010[ica];
      FReal2 r_000_001 = a_001[ica] - BmA[2] * a_000[ica];
      FReal2 r_010_001 = a_011[ica] - BmA[2] * a_010[ica];
      FReal2 r_001_001 = a_002[ica] - BmA[2] * a_001[ica];
      FReal2 b_200 = r_100_100 - BmA[0] * r_000_100;
      FReal2 b_020 = r_010_010 - BmA[1] * r_000_010;
      FReal2 b_002 = r_001_001 - BmA[2] * r_000_001;
      FReal2 b_110 = r_100_010 - BmA[0] * r_000_010;
      FReal2 b_101 = r_001_100 - BmA[2] * r_000_100;
      FReal2 b_011 = r_010_001 - BmA[1] * r_000_001;
      pOut[ica*sa + 0*sb] += -c0*b_020 + b_002 - c0*b_200;
      pOut[ica*sa + 1*sb] += c1*b_110;
      pOut[ica*sa + 2*sb] += c1*b_101;
      pOut[ica*sa + 3*sb] += -c2*b_020 + c2*b_200;
      pOut[ica*sa + 4*sb] += c1*b_011;
      // 7 * [30 flops, 45 mops]
   }
   // total: [460 flops, 805 mops]
}



static void OsrrC41( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   FReal2 a_000[9]; ShTr4(a_000, pIn[20], pIn[21], pIn[22], pIn[23], pIn[24], pIn[25], pIn[26], pIn[27], pIn[28], pIn[29], pIn[30], pIn[31], pIn[32], pIn[33], pIn[34]);
   FReal2 a_100[9]; ShTr4(a_100, pIn[35], pIn[36], pIn[37], pIn[38], pIn[39], pIn[40], pIn[41], pIn[42], pIn[43], pIn[47], pIn[48], pIn[49], pIn[53], pIn[54], pIn[55]);
   FReal2 a_010[9]; ShTr4(a_010, pIn[41], pIn[42], pIn[43], pIn[44], pIn[45], pIn[46], pIn[36], pIn[38], pIn[39], pIn[53], pIn[54], pIn[55], pIn[48], pIn[50], pIn[51]);
   FReal2 a_001[9]; ShTr4(a_001, pIn[47], pIn[48], pIn[49], pIn[50], pIn[51], pIn[52], pIn[53], pIn[54], pIn[55], pIn[37], pIn[39], pIn[40], pIn[43], pIn[45], pIn[46]);
   // ShTr: [188 flops, 320 mops]

   for ( unsigned ica = 0; ica < 9; ++ ica ) {
      FReal2 b_100 = a_100[ica] - BmA[0] * a_000[ica];
      FReal2 b_010 = a_010[ica] - BmA[1] * a_000[ica];
      FReal2 b_001 = a_001[ica] - BmA[2] * a_000[ica];
      pOut[ica*sa + 0*sb] += b_100;
      pOut[ica*sa + 1*sb] += b_010;
      pOut[ica*sa + 2*sb] += b_001;
      // 9 * [6 flops, 9 mops]
   }
   // total: [242 flops, 401 mops]
}

// In: Gm[0 .. 6] (inclusive), Out: [r]^0, ordered as AngularComps(6)
// (moment 6, 28 entries)
static void OsrrA6( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const PmA[3], FReal1 const PmQ[3], double rho, double InvEta ) AIC_NO_THROW
{
   double
      riz = rho * InvEta,
      iz2 = 0.5 * InvEta,
      PmQf[3] = { riz*PmQ[0], riz*PmQ[1], riz*PmQ[2] };
   pOut[0] = pIn[0];
   // L = 1
   pOut[1] = PmA[0] * pIn[0] - PmQf[0] * pIn[1];
   FReal2 r_100_1 = PmA[0] * pIn[1] - PmQf[0] * pIn[2];
   FReal2 r_100_2 = PmA[0] * pIn[2] - PmQf[0] * pIn[3];
   FReal2 r_100_3 = PmA[0] * pIn[3] - PmQf[0] * pIn[4];
   FReal2 r_100_4 = PmA[0] * pIn[4] - PmQf[0] * pIn[5];
   FReal2 r_100_5 = PmA[0] * pIn[5] - PmQf[0] * pIn[6];
   pOut[2] = PmA[1] * pIn[0] - PmQf[1] * pIn[1];
   FReal2 r_010_1 = PmA[1] * pIn[1] - PmQf[1] * pIn[2];
   FReal2 r_010_2 = PmA[1] * pIn[2] - PmQf[1] * pIn[3];
   FReal2 r_010_3 = PmA[1] * pIn[3] - PmQf[1] * pIn[4];
   FReal2 r_010_4 = PmA[1] * pIn[4] - PmQf[1] * pIn[5];
   FReal2 r_010_5 = PmA[1] * pIn[5] - PmQf[1] * pIn[6];
   pOut[3] = PmA[2] * pIn[0] - PmQf[2] * pIn[1];
   FReal2 r_001_1 = PmA[2] * pIn[1] - PmQf[2] * pIn[2];
   FReal2 r_001_2 = PmA[2] * pIn[2] - PmQf[2] * pIn[3];
   FReal2 r_001_3 = PmA[2] * pIn[3] - PmQf[2] * pIn[4];
   FReal2 r_001_4 = PmA[2] * pIn[4] - PmQf[2] * pIn[5];
   FReal2 r_001_5 = PmA[2] * pIn[5] - PmQf[2] * pIn[6];

   // L = 2
   pOut[4] = PmA[0] * pOut[1] - PmQf[0] * r_100_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_200_1 = PmA[0] * r_100_1 - PmQf[0] * r_100_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_200_2 = PmA[0] * r_100_2 - PmQf[0] * r_100_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_200_3 = PmA[0] * r_100_3 - PmQf[0] * r_100_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_200_4 = PmA[0] * r_100_4 - PmQf[0] * r_100_5 + iz2 * (pIn[4] - riz * pIn[5]);
   pOut[5] = PmA[1] * pOut[2] - PmQf[1] * r_010_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_020_1 = PmA[1] * r_010_1 - PmQf[1] * r_010_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_020_2 = PmA[1] * r_010_2 - PmQf[1] * r_010_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_020_3 = PmA[1] * r_010_3 - PmQf[1] * r_010_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_020_4 = PmA[1] * r_010_4 - PmQf[1] * r_010_5 + iz2 * (pIn[4] - riz * pIn[5]);
   pOut[6] = PmA[2] * pOut[3] - PmQf[2] * r_001_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_002_1 = PmA[2] * r_001_1 - PmQf[2] * r_001_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_002_2 = PmA[2] * r_001_2 - PmQf[2] * r_001_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_002_3 = PmA[2] * r_001_3 - PmQf[2] * r_001_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_002_4 = PmA[2] * r_001_4 - PmQf[2] * r_001_5 + iz2 * (pIn[4] - riz * pIn[5]);
   pOut[7] = PmA[0] * pOut[2] - PmQf[0] * r_010_1;
   pOut[8] = PmA[2] * pOut[1] - PmQf[2] * r_100_1;
   FReal2 r_101_1 = PmA[2] * r_100_1 - PmQf[2] * r_100_2;
   pOut[9] = PmA[1] * pOut[3] - PmQf[1] * r_001_1;

   // L = 3
   pOut[10] = PmA[0] * pOut[4] - PmQf[0] * r_200_1 + 2 * iz2 * (pOut[1] - riz * r_100_1);
   FReal2 r_300_1 = PmA[0] * r_200_1 - PmQf[0] * r_200_2 + 2 * iz2 * (r_100_1 - riz * r_100_2);
   FReal2 r_300_2 = PmA[0] * r_200_2 - PmQf[0] * r_200_3 + 2 * iz2 * (r_100_2 - riz * r_100_3);
   FReal2 r_300_3 = PmA[0] * r_200_3 - PmQf[0] * r_200_4 + 2 * iz2 * (r_100_3 - riz * r_100_4);
   pOut[11] = PmA[0] * pOut[5] - PmQf[0] * r_020_1;
   FReal2 r_120_1 = PmA[0] * r_020_1 - PmQf[0] * r_020_2;
   FReal2 r_120_2 = PmA[0] * r_020_2 - PmQf[0] * r_020_3;
   pOut[12] = PmA[0] * pOut[6] - PmQf[0] * r_002_1;
   pOut[13] = PmA[1] * pOut[4] - PmQf[1] * r_200_1;
   pOut[14] = PmA[1] * pOut[5] - PmQf[1] * r_020_1 + 2 * iz2 * (pOut[2] - riz * r_010_1);
   FReal2 r_030_1 = PmA[1] * r_020_1 - PmQf[1] * r_020_2 + 2 * iz2 * (r_010_1 - riz * r_010_2);
   FReal2 r_030_2 = PmA[1] * r_020_2 - PmQf[1] * r_020_3 + 2 * iz2 * (r_010_2 - riz * r_010_3);
   FReal2 r_030_3 = PmA[1] * r_020_3 - PmQf[1] * r_020_4 + 2 * iz2 * (r_010_3 - riz * r_010_4);
   pOut[15] = PmA[1] * pOut[6] - PmQf[1] * r_002_1;
   FReal2 r_012_1 = PmA[1] * r_002_1 - PmQf[1] * r_002_2;
   FReal2 r_012_2 = PmA[1] * r_002_2 - PmQf[1] * r_002_3;
   pOut[16] = PmA[0] * pOut[8] - PmQf[0] * r_101_1 + iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_201_1 = PmA[2] * r_200_1 - PmQf[2] * r_200_2;
   FReal2 r_201_2 = PmA[2] * r_200_2 - PmQf[2] * r_200_3;
   FReal2 r_201_3 = PmA[2] * r_200_3 - PmQf[2] * r_200_4;
   pOut[17] = PmA[2] * pOut[5] - PmQf[2] * r_020_1;
   pOut[18] = PmA[2] * pOut[6] - PmQf[2] * r_002_1 + 2 * iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_003_1 = PmA[2] * r_002_1 - PmQf[2] * r_002_2 + 2 * iz2 * (r_001_1 - riz * r_001_2);
   FReal2 r_003_2 = PmA[2] * r_002_2 - PmQf[2] * r_002_3 + 2 * iz2 * (r_001_2 - riz * r_001_3);
   FReal2 r_003_3 = PmA[2] * r_002_3 - PmQf[2] * r_002_4 + 2 * iz2 * (r_001_3 - riz * r_001_4);
   pOut[19] = PmA[1] * pOut[8] - PmQf[1] * r_101_1;

   // L = 4
   pOut[20] = PmA[0] * pOut[10] - PmQf[0] * r_300_1 + 3 * iz2 * (pOut[4] - riz * r_200_1);
   FReal2 r_400_1 = PmA[0] * r_300_1 - PmQf[0] * r_300_2 + 3 * iz2 * (r_200_1 - riz * r_200_2);
   FReal2 r_400_2 = PmA[0] * r_300_2 - PmQf[0] * r_300_3 + 3 * iz2 * (r_200_2 - riz * r_200_3);
   pOut[21] = PmA[0] * pOut[11] - PmQf[0] * r_120_1 + iz2 * (pOut[5] - riz * r_020_1);
   FReal2 r_220_1 = PmA[0] * r_120_1 - PmQf[0] * r_120_2 + iz2 * (r_020_1 - riz * r_020_2);
   pOut[22] = PmA[2] * pOut[16] - PmQf[2] * r_201_1 + iz2 * (pOut[4] - riz * r_200_1);
   FReal2 r_202_1 = PmA[2] * r_201_1 - PmQf[2] * r_201_2 + iz2 * (r_200_1 - riz * r_200_2);
   FReal2 r_202_2 = PmA[2] * r_201_2 - PmQf[2] * r_201_3 + iz2 * (r_200_2 - riz * r_200_3);
   pOut[23] = PmA[1] * pOut[14] - PmQf[1] * r_030_1 + 3 * iz2 * (pOut[5] - riz * r_020_1);
   FReal2 r_040_1 = PmA[1] * r_030_1 - PmQf[1] * r_030_2 + 3 * iz2 * (r_020_1 - riz * r_020_2);
   FReal2 r_040_2 = PmA[1] * r_030_2 - PmQf[1] * r_030_3 + 3 * iz2 * (r_020_2 - riz * r_020_3);
   pOut[24] = PmA[1] * pOut[15] - PmQf[1] * r_012_1 + iz2 * (pOut[6] - riz * r_002_1);
   FReal2 r_022_1 = PmA[1] * r_012_1 - PmQf[1] * r_012_2 + iz2 * (r_002_1 - riz * r_002_2);
   pOut[25] = PmA[2] * pOut[18] - PmQf[2] * r_003_1 + 3 * iz2 * (pOut[6] - riz * r_002_1);
   FReal2 r_004_1 = PmA[2] * r_003_1 - PmQf[2] * r_003_2 + 3 * iz2 * (r_002_1 - riz * r_002_2);
   FReal2 r_004_2 = PmA[2] * r_003_2 - PmQf[2] * r_003_3 + 3 * iz2 * (r_002_2 - riz * r_002_3);
   pOut[26] = PmA[1] * pOut[10] - PmQf[1] * r_300_1;
   FReal2 r_310_1 = PmA[1] * r_300_1 - PmQf[1] * r_300_2;
   FReal2 r_310_2 = PmA[1] * r_300_2 - PmQf[1] * r_300_3;
   pOut[27] = PmA[0] * pOut[14] - PmQf[0] * r_030_1;
   FReal2 r_130_1 = PmA[0] * r_030_1 - PmQf[0] * r_030_2;
   FReal2 r_130_2 = PmA[0] * r_030_2 - PmQf[0] * r_030_3;
   pOut[28] = PmA[0] * pOut[15] - PmQf[0] * r_012_1;
   pOut[29] = PmA[0] * pOut[16] - PmQf[0] * r_201_1 + 2 * iz2 * (pOut[8] - riz * r_101_1);
   FReal2 r_301_1 = PmA[2] * r_300_1 - PmQf[2] * r_300_2;
   FReal2 r_301_2 = PmA[2] * r_300_2 - PmQf[2] * r_300_3;
   pOut[30] = PmA[2] * pOut[11] - PmQf[2] * r_120_1;
   pOut[31] = PmA[0] * pOut[18] - PmQf[0] * r_003_1;
   pOut[32] = PmA[1] * pOut[16] - PmQf[1] * r_201_1;
   pOut[33] = PmA[2] * pOut[14] - PmQf[2] * r_030_1;
   FReal2 r_031_1 = PmA[2] * r_030_1 - PmQf[2] * r_030_2;
   FReal2 r_031_2 = PmA[2] * r_030_2 - PmQf[2] * r_030_3;
   pOut[34] = PmA[1] * pOut[18] - PmQf[1] * r_003_1;
   FReal2 r_013_1 = PmA[1] * r_003_1 - PmQf[1] * r_003_2;
   FReal2 r_013_2 = PmA[1] * r_003_2 - PmQf[1] * r_003_3;

   // L = 5
   pOut[35] = PmA[0] * pOut[20] - PmQf[0] * r_400_1 + 4 * iz2 * (pOut[10] - riz * r_300_1);
   FReal2 r_500_1 = PmA[0] * r_400_1 - PmQf[0] * r_400_2 + 4 * iz2 * (r_300_1 - riz * r_300_2);
   pOut[36] = PmA[0] * pOut[21] - PmQf[0] * r_220_1 + 2 * iz2 * (pOut[11] - riz * r_120_1);
   FReal2 r_320_1 = PmA[1] * r_310_1 - PmQf[1] * r_310_2 + iz2 * (r_300_1 - riz * r_300_2);
   pOut[37] = PmA[2] * pOut[29] - PmQf[2] * r_301_1 + iz2 * (pOut[10] - riz * r_300_1);
   FReal2 r_302_1 = PmA[2] * r_301_1 - PmQf[2] * r_301_2 + iz2 * (r_300_1 - riz * r_300_2);
   pOut[38] = PmA[1] * pOut[27] - PmQf[1] * r_130_1 + 3 * iz2 * (pOut[11] - riz * r_120_1);
   FReal2 r_140_1 = PmA[1] * r_130_1 - PmQf[1] * r_130_2 + 3 * iz2 * (r_120_1 - riz * r_120_2);
   pOut[39] = PmA[0] * pOut[24] - PmQf[0] * r_022_1;
   pOut[40] = PmA[0] * pOut[25] - PmQf[0] * r_004_1;
   pOut[41] = PmA[1] * pOut[20] - PmQf[1] * r_400_1;
   pOut[42] = PmA[0] * pOut[27] - PmQf[0] * r_130_1 + iz2 * (pOut[14] - riz * r_030_1);
   FReal2 r_230_1 = PmA[0] * r_130_1 - PmQf[0] * r_130_2 + iz2 * (r_030_1 - riz * r_030_2);
   pOut[43] = PmA[1] * pOut[22] - PmQf[1] * r_202_1;
   FReal2 r_212_1 = PmA[1] * r_202_1 - PmQf[1] * r_202_2;
   pOut[44] = PmA[1] * pOut[23] - PmQf[1] * r_040_1 + 4 * iz2 * (pOut[14] - riz * r_030_1);
   FReal2 r_050_1 = PmA[1] * r_040_1 - PmQf[1] * r_040_2 + 4 * iz2 * (r_030_1 - riz * r_030_2);
   pOut[45] = PmA[1] * pOut[24] - PmQf[1] * r_022_1 + 2 * iz2 * (pOut[15] - riz * r_012_1);
   FReal2 r_032_1 = PmA[2] * r_031_1 - PmQf[2] * r_031_2 + iz2 * (r_030_1 - riz * r_030_2);
   pOut[46] = PmA[2] * pOut[34] - PmQf[2] * r_013_1 + 3 * iz2 * (pOut[15] - riz * r_012_1);
   FReal2 r_014_1 = PmA[2] * r_013_1 - PmQf[2] * r_013_2 + 3 * iz2 * (r_012_1 - riz * r_012_2);
   pOut[47] = PmA[0] * pOut[29] - PmQf[0] * r_301_1 + 3 * iz2 * (pOut[16] - riz * r_201_1);
   FReal2 r_401_1 = PmA[0] * r_301_1 - PmQf[0] * r_301_2 + 3 * iz2 * (r_201_1 - riz * r_201_2);
   pOut[48] = PmA[2] * pOut[21] - PmQf[2] * r_220_1;
   pOut[49] = PmA[2] * pOut[22] - PmQf[2] * r_202_1 + 2 * iz2 * (pOut[16] - riz * r_201_1);
   FReal2 r_203_1 = PmA[2] * r_202_1 - PmQf[2] * r_202_2 + 2 * iz2 * (r_201_1 - riz * r_201_2);
   pOut[50] = PmA[2] * pOut[23] - PmQf[2] * r_040_1;
   pOut[51] = PmA[1] * pOut[34] - PmQf[1] * r_013_1 + iz2 * (pOut[18] - riz * r_003_1);
   FReal2 r_023_1 = PmA[1] * r_013_1 - PmQf[1] * r_013_2 + iz2 * (r_003_1 - riz * r_003_2);
   pOut[52] = PmA[2] * pOut[25] - PmQf[2] * r_004_1 + 4 * iz2 * (pOut[18] - riz * r_003_1);
   FReal2 r_005_1 = PmA[2] * r_004_1 - PmQf[2] * r_004_2 + 4 * iz2 * (r_003_1 - riz * r_003_2);
   pOut[53] = PmA[1] * pOut[29] - PmQf[1] * r_301_1;
   pOut[54] = PmA[2] * pOut[27] - PmQf[2] * r_130_1;
   pOut[55] = PmA[0] * pOut[34] - PmQf[0] * r_013_1;

   // L = 6
   pOut[56] = PmA[0] * pOut[35] - PmQf[0] * r_500_1 + 5 * iz2 * (pOut[20] - riz * r_400_1);
   pOut[57] = PmA[0] * pOut[36] - PmQf[0] * r_320_1 + 3 * iz2 * (pOut[21] - riz * r_220_1);
   pOut[58] = PmA[0] * pOut[37] - PmQf[0] * r_302_1 + 3 * iz2 * (pOut[22] - riz * r_202_1);
   pOut[59] = PmA[0] * pOut[38] - PmQf[0] * r_140_1 + iz2 * (pOut[23] - riz * r_040_1);
   pOut[60] = PmA[1] * pOut[43] - PmQf[1] * r_212_1 + iz2 * (pOut[22] - riz * r_202_1);
   pOut[61] = PmA[2] * pOut[49] - PmQf[2] * r_203_1 + 3 * iz2 * (pOut[22] - riz * r_202_1);
   pOut[62] = PmA[1] * pOut[44] - PmQf[1] * r_050_1 + 5 * iz2 * (pOut[23] - riz * r_040_1);
   pOut[63] = PmA[1] * pOut[45] - PmQf[1] * r_032_1 + 3 * iz2 * (pOut[24] - riz * r_022_1);
   pOut[64] = PmA[2] * pOut[51] - PmQf[2] * r_023_1 + 3 * iz2 * (pOut[24] - riz * r_022_1);
   pOut[65] = PmA[2] * pOut[52] - PmQf[2] * r_005_1 + 5 * iz2 * (pOut[25] - riz * r_004_1);
   pOut[66] = PmA[1] * pOut[35] - PmQf[1] * r_500_1;
   pOut[67] = PmA[0] * pOut[42] - PmQf[0] * r_230_1 + 2 * iz2 * (pOut[27] - riz * r_130_1);
   pOut[68] = PmA[1] * pOut[37] - PmQf[1] * r_302_1;
   pOut[69] = PmA[0] * pOut[44] - PmQf[0] * r_050_1;
   pOut[70] = PmA[0] * pOut[45] - PmQf[0] * r_032_1;
   pOut[71] = PmA[0] * pOut[46] - PmQf[0] * r_014_1;
   pOut[72] = PmA[0] * pOut[47] - PmQf[0] * r_401_1 + 4 * iz2 * (pOut[29] - riz * r_301_1);
   pOut[73] = PmA[2] * pOut[36] - PmQf[2] * r_320_1;
   pOut[74] = PmA[2] * pOut[37] - PmQf[2] * r_302_1 + 2 * iz2 * (pOut[29] - riz * r_301_1);
   pOut[75] = PmA[2] * pOut[38] - PmQf[2] * r_140_1;
   pOut[76] = PmA[0] * pOut[51] - PmQf[0] * r_023_1;
   pOut[77] = PmA[0] * pOut[52] - PmQf[0] * r_005_1;
   pOut[78] = PmA[1] * pOut[47] - PmQf[1] * r_401_1;
   pOut[79] = PmA[2] * pOut[42] - PmQf[2] * r_230_1;
   pOut[80] = PmA[1] * pOut[49] - PmQf[1] * r_203_1;
   pOut[81] = PmA[2] * pOut[44] - PmQf[2] * r_050_1;
   pOut[82] = PmA[1] * pOut[51] - PmQf[1] * r_023_1 + 2 * iz2 * (pOut[34] - riz * r_013_1);
   pOut[83] = PmA[2] * pOut[46] - PmQf[2] * r_014_1 + 4 * iz2 * (pOut[34] - riz * r_013_1);
   // 783 flops, 490 mops, 0.60kb stack
}









static void OsrrC33( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   const double c0 = 6.1237243569579447e-01;
   const double c1 = 2.4494897427831779;
   const double c2 = 1.5;
   const double c3 = 2.3717082451262841;
   const double c4 = 7.9056941504209477e-01;
   const double c5 = 3.8729833462074166;
   const double c6 = 1.9364916731037083;

   FReal2 a_000[7]; ShTr3(a_000, pIn[10], pIn[11], pIn[12], pIn[13], pIn[14], pIn[15], pIn[16], pIn[17], pIn[18], pIn[19]);
   FReal2 a_100[7]; ShTr3(a_100, pIn[20], pIn[21], pIn[22], pIn[26], pIn[27], pIn[28], pIn[29], pIn[30], pIn[31], pIn[32]);
   FReal2 a_010[7]; ShTr3(a_010, pIn[26], pIn[27], pIn[28], pIn[21], pIn[23], pIn[24], pIn[32], pIn[33], pIn[34], pIn[30]);
   FReal2 a_001[7]; ShTr3(a_001, pIn[29], pIn[30], pIn[31], pIn[32], pIn[33], pIn[34], pIn[22], pIn[24], pIn[25], pIn[28]);
   FReal2 a_200[7]; ShTr3(a_200, pIn[35], pIn[36], pIn[37], pIn[41], pIn[42], pIn[43], pIn[47], pIn[48], pIn[49], pIn[53]);
   FReal2 a_020[7]; ShTr3(a_020, pIn[36], pIn[38], pIn[39], pIn[42], pIn[44], pIn[45], pIn[48], pIn[50], pIn[51], pIn[54]);
   FReal2 a_002[7]; ShTr3(a_002, pIn[37], pIn[39], pIn[40], pIn[43], pIn[45], pIn[46], pIn[49], pIn[51], pIn[52], pIn[55]);
   FReal2 a_110[7]; ShTr3(a_110, pIn[41], pIn[42], pIn[43], pIn[36], pIn[38], pIn[39], pIn[53], pIn[54], pIn[55], pIn[48]);
   FReal2 a_101[7]; ShTr3(a_101, pIn[47], pIn[48], pIn[49], pIn[53], pIn[54], pIn[55], pIn[37], pIn[39], pIn[40], pIn[43]);
   FReal2 a_011[7]; ShTr3(a_011, pIn[53], pIn[54], pIn[55], pIn[48], pIn[50], pIn[51], pIn[43], pIn[45], pIn[46], pIn[39]);
   FReal2 a_300[7]; ShTr3(a_300, pIn[56], pIn[57], pIn[58], pIn[66], pIn[67], pIn[68], pIn[72], pIn[73], pIn[74], pIn[78]);
   FReal2 a_120[7]; ShTr3(a_120, pIn[57], pIn[59], pIn[60], pIn[67], pIn[69], pIn[70], pIn[73], pIn[75], pIn[76], pIn[79]);
   FReal2 a_102[7]; ShTr3(a_102, pIn[58], pIn[60], pIn[61], pIn[68], pIn[70], pIn[71], pIn[74], pIn[76], pIn[77], pIn[80]);
   FReal2 a_210[7]; ShTr3(a_210, pIn[66], pIn[67], pIn[68], pIn[57], pIn[59], pIn[60], pIn[78], pIn[79], pIn[80], pIn[73]);
   FReal2 a_030[7]; ShTr3(a_030, pIn[67], pIn[69], pIn[70], pIn[59], pIn[62], pIn[63], pIn[79], pIn[81], pIn[82], pIn[75]);
   FReal2 a_012[7]; ShTr3(a_012, pIn[68], pIn[70], pIn[71], pIn[60], pIn[63], pIn[64], pIn[80], pIn[82], pIn[83], pIn[76]);
   FReal2 a_201[7]; ShTr3(a_201, pIn[72], pIn[73], pIn[74], pIn[78], pIn[79], pIn[80], pIn[58], pIn[60], pIn[61], pIn[68]);
   FReal2 a_021[7]; ShTr3(a_021, pIn[73], pIn[75], pIn[76], pIn[79], pIn[81], pIn[82], pIn[60], pIn[63], pIn[64], pIn[70]);
   FReal2 a_003[7]; ShTr3(a_003, pIn[74], pIn[76], pIn[77], pIn[80], pIn[82], pIn[83], pIn[61], pIn[64], pIn[65], pIn[71]);
   FReal2 a_111[7]; ShTr3(a_111, pIn[78], pIn[79], pIn[80], pIn[73], pIn[75], pIn[76], pIn[68], pIn[70], pIn[71], pIn[60]);
   // ShTr: [500 flops, 980 mops]

   for ( unsigned ica = 0; ica < 7; ++ ica ) {
      FReal2 r_000_100 = a_100[ica] - BmA[0] * a_000[ica];
      FReal2 r_100_100 = a_200[ica] - BmA[0] * a_100[ica];
      FReal2 r_010_100 = a_110[ica] - BmA[0] * a_010[ica];
      FReal2 r_001_100 = a_101[ica] - BmA[0] * a_001[ica];
      FReal2 r_200_100 = a_300[ica] - BmA[0] * a_200[ica];
      FReal2 r_110_100 = a_210[ica] - BmA[0] * a_110[ica];
      FReal2 r_101_100 = a_201[ica] - BmA[0] * a_101[ica];
      FReal2 r_011_100 = a_111[ica] - BmA[0] * a_011[ica];
      FReal2 r_000_010 = a_010[ica] - BmA[1] * a_000[ica];
      FReal2 r_100_010 = a_110[ica] - BmA[1] * a_100[ica];
      FReal2 r_010_010 = a_020[ica] - BmA[1] * a_010[ica];
      FReal2 r_001_010 = a_011[ica] - BmA[1] * a_001[ica];
      FReal2 r_020_010 = a_030[ica] - BmA[1] * a_020[ica];
      FReal2 r_110_010 = a_120[ica] - BmA[1] * a_110[ica];
      FReal2 r_011_010 = a_021[ica] - BmA[1] * a_011[ica];
      FReal2 r_000_001 = a_001[ica] - BmA[2] * a_000[ica];
      FReal2 r_100_001 = a_101[ica] - BmA[2] * a_100[ica];
      FReal2 r_010_001 = a_011[ica] - BmA[2] * a_010[ica];
      FReal2 r_001_001 = a_002[ica] - BmA[2] * a_001[ica];
      FReal2 r_002_001 = a_003[ica] - BmA[2] * a_002[ica];
      FReal2 r_101_001 = a_102[ica] - BmA[2] * a_101[ica];
      FReal2 r_011_001 = a_012[ica] - BmA[2] * a_011[ica];
      FReal2 r_000_200 = r_100_100 - BmA[0] * r_000_100;
      FReal2 r_100_200 = r_200_100 - BmA[0] * r_100_100;
      FReal2 r_010_200 = r_110_100 - BmA[0] * r_010_100;
      FReal2 r_001_200 = r_101_100 - BmA[0] * r_001_100;
      FReal2 r_000_020 = r_010_010 - BmA[1] * r_000_010;
      FReal2 r_100_020 = r_110_010 - BmA[1] * r_100_010;
      FReal2 r_010_020 = r_020_010 - BmA[1] * r_010_010;
      FReal2 r_001_020 = r_011_010 - BmA[1] * r_001_010;
      FReal2 r_000_002 = r_001_001 - BmA[2] * r_000_001;
      FReal2 r_100_002 = r_101_001 - BmA[2] * r_100_001;
      FReal2 r_010_002 = r_011_001 - BmA[2] * r_010_001;
      FReal2 r_001_002 = r_002_001 - BmA[2] * r_001_001;
      FReal2 r_000_101 = r_001_100 - BmA[2] * r_000_100;
      FReal2 r_010_101 = r_011_100 - BmA[2] * r_010_100;
      FReal2 b_300 = r_100_200 - BmA[0] * r_000_200;
      FReal2 b_120 = r_100_020 - BmA[0] * r_000_020;
      FReal2 b_102 = r_100_002 - BmA[0] * r_000_002;
      FReal2 b_210 = r_010_200 - BmA[1] * r_000_200;
      FReal2 b_030 = r_010_020 - BmA[1] * r_000_020;
      FReal2 b_012 = r_010_002 - BmA[1] * r_000_002;
      FReal2 b_201 = r_001_200 - BmA[2] * r_000_200;
      FReal2 b_021 = r_001_020 - BmA[2] * r_000_020;
      FReal2 b_003 = r_001_002 - BmA[2] * r_000_002;
      FReal2 b_111 = r_010_101 - BmA[1] * r_000_101;
      pOut[ica*sa + 0*sb] += -c0*b_300 + c1*b_102 - c0*b_120;
      pOut[ica*sa + 1*sb] += c1*b_012 - c0*b_030 - c0*b_210;
      pOut[ica*sa + 2*sb] += -c2*b_021 - c2*b_201 + b_003;
      pOut[ica*sa + 3*sb] += -c3*b_120 + c4*b_300;
      pOut[ica*sa + 4*sb] += c5*b_111;
      pOut[ica*sa + 5*sb] += -c4*b_030 + c3*b_210;
      pOut[ica*sa + 6*sb] += -c6*b_021 + c6*b_201;
      // 7 * [92 flops, 138 mops]
   }
   // total: [1144 flops, 1946 mops]
}



static void OsrrC42( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   const double c0 = 5.e-01;
   const double c1 = 1.7320508075688772;
   const double c2 = 8.660254037844386e-01;

   FReal2 a_000[9]; ShTr4(a_000, pIn[20], pIn[21], pIn[22], pIn[23], pIn[24], pIn[25], pIn[26], pIn[27], pIn[28], pIn[29], pIn[30], pIn[31], pIn[32], pIn[33], pIn[34]);
   FReal2 a_100[9]; ShTr4(a_100, pIn[35], pIn[36], pIn[37], pIn[38], pIn[39], pIn[40], pIn[41], pIn[42], pIn[43], pIn[47], pIn[48], pIn[49], pIn[53], pIn[54], pIn[55]);
   FReal2 a_010[9]; ShTr4(a_010, pIn[41], pIn[42], pIn[43], pIn[44], pIn[45], pIn[46], pIn[36], pIn[38], pIn[39], pIn[53], pIn[54], pIn[55], pIn[48], pIn[50], pIn[51]);
   FReal2 a_001[9]; ShTr4(a_001, pIn[47], pIn[48], pIn[49], pIn[50], pIn[51], pIn[52], pIn[53], pIn[54], pIn[55], pIn[37], pIn[39], pIn[40], pIn[43], pIn[45], pIn[46]);
   FReal2 a_200[9]; ShTr4(a_200, pIn[56], pIn[57], pIn[58], pIn[59], pIn[60], pIn[61], pIn[66], pIn[67], pIn[68], pIn[72], pIn[73], pIn[74], pIn[78], pIn[79], pIn[80]);
   FReal2 a_020[9]; ShTr4(a_020, pIn[57], pIn[59], pIn[60], pIn[62], pIn[63], pIn[64], pIn[67], pIn[69], pIn[70], pIn[73], pIn[75], pIn[76], pIn[79], pIn[81], pIn[82]);
   FReal2 a_002[9]; ShTr4(a_002, pIn[58], pIn[60], pIn[61], pIn[63], pIn[64], pIn[65], pIn[68], pIn[70], pIn[71], pIn[74], pIn[76], pIn[77], pIn[80], pIn[82], pIn[83]);
   FReal2 a_110[9]; ShTr4(a_110, pIn[66], pIn[67], pIn[68], pIn[69], pIn[70], pIn[71], pIn[57], pIn[59], pIn[60], pIn[78], pIn[79], pIn[80], pIn[73], pIn[75], pIn[76]);
   FReal2 a_101[9]; ShTr4(a_101, pIn[72], pIn[73], pIn[74], pIn[75], pIn[76], pIn[77], pIn[78], pIn[79], pIn[80], pIn[58], pIn[60], pIn[61], pIn[68], pIn[70], pIn[71]);
   FReal2 a_011[9]; ShTr4(a_011, pIn[78], pIn[79], pIn[80], pIn[81], pIn[82], pIn[83], pIn[73], pIn[75], pIn[76], pIn[68], pIn[70], pIn[71], pIn[60], pIn[63], pIn[64]);
   // ShTr: [470 flops, 800 mops]

   for ( unsigned ica = 0; ica < 9; ++ ica ) {
      FReal2 r_000_100 = a_100[ica] - BmA[0] * a_000[ica];
      FReal2 r_100_100 = a_200[ica] - BmA[0] * a_100[ica];
      FReal2 r_001_100 = a_101[ica] - BmA[0] * a_001[ica];
      FReal2 r_000_010 = a_010[ica] - BmA[1] * a_000[ica];
      FReal2 r_100_010 = a_110[ica] - BmA[1] * a_100[ica];
      FReal2 r_010_010 = a_020[ica] - BmA[1] * a_010[ica];
      FReal2 r_000_001 = a_001[ica] - BmA[2] * a_000[ica];
      FReal2 r_010_001 = a_011[ica] - BmA[2] * a_010[ica];
      FReal2 r_001_001 = a_002[ica] - BmA[2] * a_001[ica];
      FReal2 b_200 = r_100_100 - BmA[0] * r_000_100;
      FReal2 b_020 = r_010_010 - BmA[1] * r_000_010;
      FReal2 b_002 = r_001_001 - BmA[2] * r_000_001;
      FReal2 b_110 = r_100_010 - BmA[0] * r_000_010;
      FReal2 b_101 = r_001_100 - BmA[2] * r_000_100;
      FReal2 b_011 = r_010_001 - BmA[1] * r_000_001;
      pOut[ica*sa + 0*sb] += -c0*b_020 + b_002 - c0*b_200;
      pOut[ica*sa + 1*sb] += c1*b_110;
      pOut[ica*sa + 2*sb] += c1*b_101;
      pOut[ica*sa + 3*sb] += -c2*b_020 + c2*b_200;
      pOut[ica*sa + 4*sb] += c1*b_011;
      // 9 * [30 flops, 45 mops]
   }
   // total: [740 flops, 1205 mops]
}

// In: Gm[0 .. 7] (inclusive), Out: [r]^0, ordered as AngularComps(7)
// (moment 7, 36 entries)
static void OsrrA7( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const PmA[3], FReal1 const PmQ[3], double rho, double InvEta ) AIC_NO_THROW
{
   double
      riz = rho * InvEta,
      iz2 = 0.5 * InvEta,
      PmQf[3] = { riz*PmQ[0], riz*PmQ[1], riz*PmQ[2] };
   pOut[0] = pIn[0];
   // L = 1
   pOut[1] = PmA[0] * pIn[0] - PmQf[0] * pIn[1];
   FReal2 r_100_1 = PmA[0] * pIn[1] - PmQf[0] * pIn[2];
   FReal2 r_100_2 = PmA[0] * pIn[2] - PmQf[0] * pIn[3];
   FReal2 r_100_3 = PmA[0] * pIn[3] - PmQf[0] * pIn[4];
   FReal2 r_100_4 = PmA[0] * pIn[4] - PmQf[0] * pIn[5];
   FReal2 r_100_5 = PmA[0] * pIn[5] - PmQf[0] * pIn[6];
   FReal2 r_100_6 = PmA[0] * pIn[6] - PmQf[0] * pIn[7];
   pOut[2] = PmA[1] * pIn[0] - PmQf[1] * pIn[1];
   FReal2 r_010_1 = PmA[1] * pIn[1] - PmQf[1] * pIn[2];
   FReal2 r_010_2 = PmA[1] * pIn[2] - PmQf[1] * pIn[3];
   FReal2 r_010_3 = PmA[1] * pIn[3] - PmQf[1] * pIn[4];
   FReal2 r_010_4 = PmA[1] * pIn[4] - PmQf[1] * pIn[5];
   FReal2 r_010_5 = PmA[1] * pIn[5] - PmQf[1] * pIn[6];
   FReal2 r_010_6 = PmA[1] * pIn[6] - PmQf[1] * pIn[7];
   pOut[3] = PmA[2] * pIn[0] - PmQf[2] * pIn[1];
   FReal2 r_001_1 = PmA[2] * pIn[1] - PmQf[2] * pIn[2];
   FReal2 r_001_2 = PmA[2] * pIn[2] - PmQf[2] * pIn[3];
   FReal2 r_001_3 = PmA[2] * pIn[3] - PmQf[2] * pIn[4];
   FReal2 r_001_4 = PmA[2] * pIn[4] - PmQf[2] * pIn[5];
   FReal2 r_001_5 = PmA[2] * pIn[5] - PmQf[2] * pIn[6];
   FReal2 r_001_6 = PmA[2] * pIn[6] - PmQf[2] * pIn[7];

   // L = 2
   pOut[4] = PmA[0] * pOut[1] - PmQf[0] * r_100_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_200_1 = PmA[0] * r_100_1 - PmQf[0] * r_100_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_200_2 = PmA[0] * r_100_2 - PmQf[0] * r_100_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_200_3 = PmA[0] * r_100_3 - PmQf[0] * r_100_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_200_4 = PmA[0] * r_100_4 - PmQf[0] * r_100_5 + iz2 * (pIn[4] - riz * pIn[5]);
   FReal2 r_200_5 = PmA[0] * r_100_5 - PmQf[0] * r_100_6 + iz2 * (pIn[5] - riz * pIn[6]);
   pOut[5] = PmA[1] * pOut[2] - PmQf[1] * r_010_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_020_1 = PmA[1] * r_010_1 - PmQf[1] * r_010_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_020_2 = PmA[1] * r_010_2 - PmQf[1] * r_010_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_020_3 = PmA[1] * r_010_3 - PmQf[1] * r_010_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_020_4 = PmA[1] * r_010_4 - PmQf[1] * r_010_5 + iz2 * (pIn[4] - riz * pIn[5]);
   FReal2 r_020_5 = PmA[1] * r_010_5 - PmQf[1] * r_010_6 + iz2 * (pIn[5] - riz * pIn[6]);
   pOut[6] = PmA[2] * pOut[3] - PmQf[2] * r_001_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_002_1 = PmA[2] * r_001_1 - PmQf[2] * r_001_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_002_2 = PmA[2] * r_001_2 - PmQf[2] * r_001_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_002_3 = PmA[2] * r_001_3 - PmQf[2] * r_001_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_002_4 = PmA[2] * r_001_4 - PmQf[2] * r_001_5 + iz2 * (pIn[4] - riz * pIn[5]);
   FReal2 r_002_5 = PmA[2] * r_001_5 - PmQf[2] * r_001_6 + iz2 * (pIn[5] - riz * pIn[6]);
   pOut[7] = PmA[0] * pOut[2] - PmQf[0] * r_010_1;
   pOut[8] = PmA[2] * pOut[1] - PmQf[2] * r_100_1;
   FReal2 r_101_1 = PmA[2] * r_100_1 - PmQf[2] * r_100_2;
   pOut[9] = PmA[1] * pOut[3] - PmQf[1] * r_001_1;

   // L = 3
   pOut[10] = PmA[0] * pOut[4] - PmQf[0] * r_200_1 + 2 * iz2 * (pOut[1] - riz * r_100_1);
   FReal2 r_300_1 = PmA[0] * r_200_1 - PmQf[0] * r_200_2 + 2 * iz2 * (r_100_1 - riz * r_100_2);
   FReal2 r_300_2 = PmA[0] * r_200_2 - PmQf[0] * r_200_3 + 2 * iz2 * (r_100_2 - riz * r_100_3);
   FReal2 r_300_3 = PmA[0] * r_200_3 - PmQf[0] * r_200_4 + 2 * iz2 * (r_100_3 - riz * r_100_4);
   FReal2 r_300_4 = PmA[0] * r_200_4 - PmQf[0] * r_200_5 + 2 * iz2 * (r_100_4 - riz * r_100_5);
   pOut[11] = PmA[0] * pOut[5] - PmQf[0] * r_020_1;
   FReal2 r_120_1 = PmA[0] * r_020_1 - PmQf[0] * r_020_2;
   FReal2 r_120_2 = PmA[0] * r_020_2 - PmQf[0] * r_020_3;
   FReal2 r_120_3 = PmA[0] * r_020_3 - PmQf[0] * r_020_4;
   FReal2 r_120_4 = PmA[0] * r_020_4 - PmQf[0] * r_020_5;
   pOut[12] = PmA[0] * pOut[6] - PmQf[0] * r_002_1;
   pOut[13] = PmA[1] * pOut[4] - PmQf[1] * r_200_1;
   pOut[14] = PmA[1] * pOut[5] - PmQf[1] * r_020_1 + 2 * iz2 * (pOut[2] - riz * r_010_1);
   FReal2 r_030_1 = PmA[1] * r_020_1 - PmQf[1] * r_020_2 + 2 * iz2 * (r_010_1 - riz * r_010_2);
   FReal2 r_030_2 = PmA[1] * r_020_2 - PmQf[1] * r_020_3 + 2 * iz2 * (r_010_2 - riz * r_010_3);
   FReal2 r_030_3 = PmA[1] * r_020_3 - PmQf[1] * r_020_4 + 2 * iz2 * (r_010_3 - riz * r_010_4);
   FReal2 r_030_4 = PmA[1] * r_020_4 - PmQf[1] * r_020_5 + 2 * iz2 * (r_010_4 - riz * r_010_5);
   pOut[15] = PmA[1] * pOut[6] - PmQf[1] * r_002_1;
   FReal2 r_012_1 = PmA[1] * r_002_1 - PmQf[1] * r_002_2;
   FReal2 r_012_2 = PmA[1] * r_002_2 - PmQf[1] * r_002_3;
   pOut[16] = PmA[0] * pOut[8] - PmQf[0] * r_101_1 + iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_201_1 = PmA[2] * r_200_1 - PmQf[2] * r_200_2;
   FReal2 r_201_2 = PmA[2] * r_200_2 - PmQf[2] * r_200_3;
   FReal2 r_201_3 = PmA[2] * r_200_3 - PmQf[2] * r_200_4;
   pOut[17] = PmA[2] * pOut[5] - PmQf[2] * r_020_1;
   pOut[18] = PmA[2] * pOut[6] - PmQf[2] * r_002_1 + 2 * iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_003_1 = PmA[2] * r_002_1 - PmQf[2] * r_002_2 + 2 * iz2 * (r_001_1 - riz * r_001_2);
   FReal2 r_003_2 = PmA[2] * r_002_2 - PmQf[2] * r_002_3 + 2 * iz2 * (r_001_2 - riz * r_001_3);
   FReal2 r_003_3 = PmA[2] * r_002_3 - PmQf[2] * r_002_4 + 2 * iz2 * (r_001_3 - riz * r_001_4);
   FReal2 r_003_4 = PmA[2] * r_002_4 - PmQf[2] * r_002_5 + 2 * iz2 * (r_001_4 - riz * r_001_5);
   pOut[19] = PmA[1] * pOut[8] - PmQf[1] * r_101_1;

   // L = 4
   pOut[20] = PmA[0] * pOut[10] - PmQf[0] * r_300_1 + 3 * iz2 * (pOut[4] - riz * r_200_1);
   FReal2 r_400_1 = PmA[0] * r_300_1 - PmQf[0] * r_300_2 + 3 * iz2 * (r_200_1 - riz * r_200_2);
   FReal2 r_400_2 = PmA[0] * r_300_2 - PmQf[0] * r_300_3 + 3 * iz2 * (r_200_2 - riz * r_200_3);
   FReal2 r_400_3 = PmA[0] * r_300_3 - PmQf[0] * r_300_4 + 3 * iz2 * (r_200_3 - riz * r_200_4);
   pOut[21] = PmA[0] * pOut[11] - PmQf[0] * r_120_1 + iz2 * (pOut[5] - riz * r_020_1);
   FReal2 r_220_1 = PmA[0] * r_120_1 - PmQf[0] * r_120_2 + iz2 * (r_020_1 - riz * r_020_2);
   FReal2 r_220_2 = PmA[0] * r_120_2 - PmQf[0] * r_120_3 + iz2 * (r_020_2 - riz * r_020_3);
   FReal2 r_220_3 = PmA[0] * r_120_3 - PmQf[0] * r_120_4 + iz2 * (r_020_3 - riz * r_020_4);
   pOut[22] = PmA[2] * pOut[16] - PmQf[2] * r_201_1 + iz2 * (pOut[4] - riz * r_200_1);
   FReal2 r_202_1 = PmA[2] * r_201_1 - PmQf[2] * r_201_2 + iz2 * (r_200_1 - riz * r_200_2);
   FReal2 r_202_2 = PmA[2] * r_201_2 - PmQf[2] * r_201_3 + iz2 * (r_200_2 - riz * r_200_3);
   pOut[23] = PmA[1] * pOut[14] - PmQf[1] * r_030_1 + 3 * iz2 * (pOut[5] - riz * r_020_1);
   FReal2 r_040_1 = PmA[1] * r_030_1 - PmQf[1] * r_030_2 + 3 * iz2 * (r_020_1 - riz * r_020_2);
   FReal2 r_040_2 = PmA[1] * r_030_2 - PmQf[1] * r_030_3 + 3 * iz2 * (r_020_2 - riz * r_020_3);
   FReal2 r_040_3 = PmA[1] * r_030_3 - PmQf[1] * r_030_4 + 3 * iz2 * (r_020_3 - riz * r_020_4);
   pOut[24] = PmA[1] * pOut[15] - PmQf[1] * r_012_1 + iz2 * (pOut[6] - riz * r_002_1);
   FReal2 r_022_1 = PmA[1] * r_012_1 - PmQf[1] * r_012_2 + iz2 * (r_002_1 - riz * r_002_2);
   pOut[25] = PmA[2] * pOut[18] - PmQf[2] * r_003_1 + 3 * iz2 * (pOut[6] - riz * r_002_1);
   FReal2 r_004_1 = PmA[2] * r_003_1 - PmQf[2] * r_003_2 + 3 * iz2 * (r_002_1 - riz * r_002_2);
   FReal2 r_004_2 = PmA[2] * r_003_2 - PmQf[2] * r_003_3 + 3 * iz2 * (r_002_2 - riz * r_002_3);
   FReal2 r_004_3 = PmA[2] * r_003_3 - PmQf[2] * r_003_4 + 3 * iz2 * (r_002_3 - riz * r_002_4);
   pOut[26] = PmA[1] * pOut[10] - PmQf[1] * r_300_1;
   pOut[27] = PmA[0] * pOut[14] - PmQf[0] * r_030_1;
   FReal2 r_130_1 = PmA[0] * r_030_1 - PmQf[0] * r_030_2;
   FReal2 r_130_2 = PmA[0] * r_030_2 - PmQf[0] * r_030_3;
   FReal2 r_130_3 = PmA[0] * r_030_3 - PmQf[0] * r_030_4;
   pOut[28] = PmA[0] * pOut[15] - PmQf[0] * r_012_1;
   pOut[29] = PmA[0] * pOut[16] - PmQf[0] * r_201_1 + 2 * iz2 * (pOut[8] - riz * r_101_1);
   FReal2 r_301_1 = PmA[2] * r_300_1 - PmQf[2] * r_300_2;
   FReal2 r_301_2 = PmA[2] * r_300_2 - PmQf[2] * r_300_3;
   FReal2 r_301_3 = PmA[2] * r_300_3 - PmQf[2] * r_300_4;
   pOut[30] = PmA[2] * pOut[11] - PmQf[2] * r_120_1;
   pOut[31] = PmA[0] * pOut[18] - PmQf[0] * r_003_1;
   pOut[32] = PmA[1] * pOut[16] - PmQf[1] * r_201_1;
   pOut[33] = PmA[2] * pOut[14] - PmQf[2] * r_030_1;
   FReal2 r_031_1 = PmA[2] * r_030_1 - PmQf[2] * r_030_2;
   FReal2 r_031_2 = PmA[2] * r_030_2 - PmQf[2] * r_030_3;
   pOut[34] = PmA[1] * pOut[18] - PmQf[1] * r_003_1;
   FReal2 r_013_1 = PmA[1] * r_003_1 - PmQf[1] * r_003_2;
   FReal2 r_013_2 = PmA[1] * r_003_2 - PmQf[1] * r_003_3;
   FReal2 r_013_3 = PmA[1] * r_003_3 - PmQf[1] * r_003_4;

   // L = 5
   pOut[35] = PmA[0] * pOut[20] - PmQf[0] * r_400_1 + 4 * iz2 * (pOut[10] - riz * r_300_1);
   FReal2 r_500_1 = PmA[0] * r_400_1 - PmQf[0] * r_400_2 + 4 * iz2 * (r_300_1 - riz * r_300_2);
   FReal2 r_500_2 = PmA[0] * r_400_2 - PmQf[0] * r_400_3 + 4 * iz2 * (r_300_2 - riz * r_300_3);
   pOut[36] = PmA[0] * pOut[21] - PmQf[0] * r_220_1 + 2 * iz2 * (pOut[11] - riz * r_120_1);
   FReal2 r_320_1 = PmA[0] * r_220_1 - PmQf[0] * r_220_2 + 2 * iz2 * (r_120_1 - riz * r_120_2);
   FReal2 r_320_2 = PmA[0] * r_220_2 - PmQf[0] * r_220_3 + 2 * iz2 * (r_120_2 - riz * r_120_3);
   pOut[37] = PmA[2] * pOut[29] - PmQf[2] * r_301_1 + iz2 * (pOut[10] - riz * r_300_1);
   FReal2 r_302_1 = PmA[2] * r_301_1 - PmQf[2] * r_301_2 + iz2 * (r_300_1 - riz * r_300_2);
   FReal2 r_302_2 = PmA[2] * r_301_2 - PmQf[2] * r_301_3 + iz2 * (r_300_2 - riz * r_300_3);
   pOut[38] = PmA[1] * pOut[27] - PmQf[1] * r_130_1 + 3 * iz2 * (pOut[11] - riz * r_120_1);
   FReal2 r_140_1 = PmA[1] * r_130_1 - PmQf[1] * r_130_2 + 3 * iz2 * (r_120_1 - riz * r_120_2);
   pOut[39] = PmA[0] * pOut[24] - PmQf[0] * r_022_1;
   pOut[40] = PmA[0] * pOut[25] - PmQf[0] * r_004_1;
   FReal2 r_104_1 = PmA[0] * r_004_1 - PmQf[0] * r_004_2;
   FReal2 r_104_2 = PmA[0] * r_004_2 - PmQf[0] * r_004_3;
   pOut[41] = PmA[1] * pOut[20] - PmQf[1] * r_400_1;
   pOut[42] = PmA[0] * pOut[27] - PmQf[0] * r_130_1 + iz2 * (pOut[14] - riz * r_030_1);
   FReal2 r_230_1 = PmA[0] * r_130_1 - PmQf[0] * r_130_2 + iz2 * (r_030_1 - riz * r_030_2);
   FReal2 r_230_2 = PmA[0] * r_130_2 - PmQf[0] * r_130_3 + iz2 * (r_030_2 - riz * r_030_3);
   pOut[43] = PmA[1] * pOut[22] - PmQf[1] * r_202_1;
   FReal2 r_212_1 = PmA[1] * r_202_1 - PmQf[1] * r_202_2;
   pOut[44] = PmA[1] * pOut[23] - PmQf[1] * r_040_1 + 4 * iz2 * (pOut[14] - riz * r_030_1);
   FReal2 r_050_1 = PmA[1] * r_040_1 - PmQf[1] * r_040_2 + 4 * iz2 * (r_030_1 - riz * r_030_2);
   FReal2 r_050_2 = PmA[1] * r_040_2 - PmQf[1] * r_040_3 + 4 * iz2 * (r_030_2 - riz * r_030_3);
   pOut[45] = PmA[1] * pOut[24] - PmQf[1] * r_022_1 + 2 * iz2 * (pOut[15] - riz * r_012_1);
   FReal2 r_032_1 = PmA[2] * r_031_1 - PmQf[2] * r_031_2 + iz2 * (r_030_1 - riz * r_030_2);
   pOut[46] = PmA[2] * pOut[34] - PmQf[2] * r_013_1 + 3 * iz2 * (pOut[15] - riz * r_012_1);
   FReal2 r_014_1 = PmA[2] * r_013_1 - PmQf[2] * r_013_2 + 3 * iz2 * (r_012_1 - riz * r_012_2);
   FReal2 r_014_2 = PmA[1] * r_004_2 - PmQf[1] * r_004_3;
   pOut[47] = PmA[0] * pOut[29] - PmQf[0] * r_301_1 + 3 * iz2 * (pOut[16] - riz * r_201_1);
   FReal2 r_401_1 = PmA[0] * r_301_1 - PmQf[0] * r_301_2 + 3 * iz2 * (r_201_1 - riz * r_201_2);
   pOut[48] = PmA[2] * pOut[21] - PmQf[2] * r_220_1;
   FReal2 r_221_1 = PmA[2] * r_220_1 - PmQf[2] * r_220_2;
   FReal2 r_221_2 = PmA[2] * r_220_2 - PmQf[2] * r_220_3;
   pOut[49] = PmA[2] * pOut[22] - PmQf[2] * r_202_1 + 2 * iz2 * (pOut[16] - riz * r_201_1);
   FReal2 r_203_1 = PmA[2] * r_202_1 - PmQf[2] * r_202_2 + 2 * iz2 * (r_201_1 - riz * r_201_2);
   pOut[50] = PmA[2] * pOut[23] - PmQf[2] * r_040_1;
   FReal2 r_041_1 = PmA[2] * r_040_1 - PmQf[2] * r_040_2;
   FReal2 r_041_2 = PmA[2] * r_040_2 - PmQf[2] * r_040_3;
   pOut[51] = PmA[1] * pOut[34] - PmQf[1] * r_013_1 + iz2 * (pOut[18] - riz * r_003_1);
   FReal2 r_023_1 = PmA[1] * r_013_1 - PmQf[1] * r_013_2 + iz2 * (r_003_1 - riz * r_003_2);
   FReal2 r_023_2 = PmA[1] * r_013_2 - PmQf[1] * r_013_3 + iz2 * (r_003_2 - riz * r_003_3);
   pOut[52] = PmA[2] * pOut[25] - PmQf[2] * r_004_1 + 4 * iz2 * (pOut[18] - riz * r_003_1);
   FReal2 r_005_1 = PmA[2] * r_004_1 - PmQf[2] * r_004_2 + 4 * iz2 * (r_003_1 - riz * r_003_2);
   FReal2 r_005_2 = PmA[2] * r_004_2 - PmQf[2] * r_004_3 + 4 * iz2 * (r_003_2 - riz * r_003_3);
   pOut[53] = PmA[1] * pOut[29] - PmQf[1] * r_301_1;
   pOut[54] = PmA[2] * pOut[27] - PmQf[2] * r_130_1;
   pOut[55] = PmA[0] * pOut[34] - PmQf[0] * r_013_1;

   // L = 6
   pOut[56] = PmA[0] * pOut[35] - PmQf[0] * r_500_1 + 5 * iz2 * (pOut[20] - riz * r_400_1);
   FReal2 r_600_1 = PmA[0] * r_500_1 - PmQf[0] * r_500_2 + 5 * iz2 * (r_400_1 - riz * r_400_2);
   pOut[57] = PmA[0] * pOut[36] - PmQf[0] * r_320_1 + 3 * iz2 * (pOut[21] - riz * r_220_1);
   FReal2 r_420_1 = PmA[0] * r_320_1 - PmQf[0] * r_320_2 + 3 * iz2 * (r_220_1 - riz * r_220_2);
   pOut[58] = PmA[0] * pOut[37] - PmQf[0] * r_302_1 + 3 * iz2 * (pOut[22] - riz * r_202_1);
   pOut[59] = PmA[0] * pOut[38] - PmQf[0] * r_140_1 + iz2 * (pOut[23] - riz * r_040_1);
   FReal2 r_240_1 = PmA[1] * r_230_1 - PmQf[1] * r_230_2 + 3 * iz2 * (r_220_1 - riz * r_220_2);
   pOut[60] = PmA[1] * pOut[43] - PmQf[1] * r_212_1 + iz2 * (pOut[22] - riz * r_202_1);
   FReal2 r_222_1 = PmA[2] * r_221_1 - PmQf[2] * r_221_2 + iz2 * (r_220_1 - riz * r_220_2);
   pOut[61] = PmA[2] * pOut[49] - PmQf[2] * r_203_1 + 3 * iz2 * (pOut[22] - riz * r_202_1);
   FReal2 r_204_1 = PmA[0] * r_104_1 - PmQf[0] * r_104_2 + iz2 * (r_004_1 - riz * r_004_2);
   pOut[62] = PmA[1] * pOut[44] - PmQf[1] * r_050_1 + 5 * iz2 * (pOut[23] - riz * r_040_1);
   FReal2 r_060_1 = PmA[1] * r_050_1 - PmQf[1] * r_050_2 + 5 * iz2 * (r_040_1 - riz * r_040_2);
   pOut[63] = PmA[1] * pOut[45] - PmQf[1] * r_032_1 + 3 * iz2 * (pOut[24] - riz * r_022_1);
   FReal2 r_042_1 = PmA[2] * r_041_1 - PmQf[2] * r_041_2 + iz2 * (r_040_1 - riz * r_040_2);
   pOut[64] = PmA[2] * pOut[51] - PmQf[2] * r_023_1 + 3 * iz2 * (pOut[24] - riz * r_022_1);
   FReal2 r_024_1 = PmA[1] * r_014_1 - PmQf[1] * r_014_2 + iz2 * (r_004_1 - riz * r_004_2);
   pOut[65] = PmA[2] * pOut[52] - PmQf[2] * r_005_1 + 5 * iz2 * (pOut[25] - riz * r_004_1);
   FReal2 r_006_1 = PmA[2] * r_005_1 - PmQf[2] * r_005_2 + 5 * iz2 * (r_004_1 - riz * r_004_2);
   pOut[66] = PmA[1] * pOut[35] - PmQf[1] * r_500_1;
   pOut[67] = PmA[0] * pOut[42] - PmQf[0] * r_230_1 + 2 * iz2 * (pOut[27] - riz * r_130_1);
   FReal2 r_330_1 = PmA[0] * r_230_1 - PmQf[0] * r_230_2 + 2 * iz2 * (r_130_1 - riz * r_130_2);
   pOut[68] = PmA[1] * pOut[37] - PmQf[1] * r_302_1;
   FReal2 r_312_1 = PmA[1] * r_302_1 - PmQf[1] * r_302_2;
   pOut[69] = PmA[1] * pOut[38] - PmQf[1] * r_140_1 + 4 * iz2 * (pOut[27] - riz * r_130_1);
   FReal2 r_150_1 = PmA[0] * r_050_1 - PmQf[0] * r_050_2;
   pOut[70] = PmA[0] * pOut[45] - PmQf[0] * r_032_1;
   pOut[71] = PmA[0] * pOut[46] - PmQf[0] * r_014_1;
   pOut[72] = PmA[0] * pOut[47] - PmQf[0] * r_401_1 + 4 * iz2 * (pOut[29] - riz * r_301_1);
   FReal2 r_501_1 = PmA[2] * r_500_1 - PmQf[2] * r_500_2;
   pOut[73] = PmA[2] * pOut[36] - PmQf[2] * r_320_1;
   pOut[74] = PmA[2] * pOut[37] - PmQf[2] * r_302_1 + 2 * iz2 * (pOut[29] - riz * r_301_1);
   FReal2 r_303_1 = PmA[2] * r_302_1 - PmQf[2] * r_302_2 + 2 * iz2 * (r_301_1 - riz * r_301_2);
   pOut[75] = PmA[2] * pOut[38] - PmQf[2] * r_140_1;
   pOut[76] = PmA[0] * pOut[51] - PmQf[0] * r_023_1;
   FReal2 r_123_1 = PmA[0] * r_023_1 - PmQf[0] * r_023_2;
   pOut[77] = PmA[0] * pOut[52] - PmQf[0] * r_005_1;
   pOut[78] = PmA[1] * pOut[47] - PmQf[1] * r_401_1;
   pOut[79] = PmA[2] * pOut[42] - PmQf[2] * r_230_1;
   pOut[80] = PmA[1] * pOut[49] - PmQf[1] * r_203_1;
   pOut[81] = PmA[2] * pOut[44] - PmQf[2] * r_050_1;
   pOut[82] = PmA[1] * pOut[51] - PmQf[1] * r_023_1 + 2 * iz2 * (pOut[34] - riz * r_013_1);
   FReal2 r_033_1 = PmA[1] * r_023_1 - PmQf[1] * r_023_2 + 2 * iz2 * (r_013_1 - riz * r_013_2);
   pOut[83] = PmA[2] * pOut[46] - PmQf[2] * r_014_1 + 4 * iz2 * (pOut[34] - riz * r_013_1);
   FReal2 r_015_1 = PmA[2] * r_014_1 - PmQf[2] * r_014_2 + 4 * iz2 * (r_013_1 - riz * r_013_2);

   // L = 7
   pOut[84] = PmA[0] * pOut[56] - PmQf[0] * r_600_1 + 6 * iz2 * (pOut[35] - riz * r_500_1);
   pOut[85] = PmA[0] * pOut[57] - PmQf[0] * r_420_1 + 4 * iz2 * (pOut[36] - riz * r_320_1);
   pOut[86] = PmA[2] * pOut[72] - PmQf[2] * r_501_1 + iz2 * (pOut[35] - riz * r_500_1);
   pOut[87] = PmA[0] * pOut[59] - PmQf[0] * r_240_1 + 2 * iz2 * (pOut[38] - riz * r_140_1);
   pOut[88] = PmA[1] * pOut[68] - PmQf[1] * r_312_1 + iz2 * (pOut[37] - riz * r_302_1);
   pOut[89] = PmA[0] * pOut[61] - PmQf[0] * r_204_1 + 2 * iz2 * (pOut[40] - riz * r_104_1);
   pOut[90] = PmA[0] * pOut[62] - PmQf[0] * r_060_1;
   pOut[91] = PmA[0] * pOut[63] - PmQf[0] * r_042_1;
   pOut[92] = PmA[0] * pOut[64] - PmQf[0] * r_024_1;
   pOut[93] = PmA[0] * pOut[65] - PmQf[0] * r_006_1;
   pOut[94] = PmA[1] * pOut[56] - PmQf[1] * r_600_1;
   pOut[95] = PmA[0] * pOut[67] - PmQf[0] * r_330_1 + 3 * iz2 * (pOut[42] - riz * r_230_1);
   pOut[96] = PmA[0] * pOut[68] - PmQf[0] * r_312_1 + 3 * iz2 * (pOut[43] - riz * r_212_1);
   pOut[97] = PmA[0] * pOut[69] - PmQf[0] * r_150_1 + iz2 * (pOut[44] - riz * r_050_1);
   pOut[98] = PmA[1] * pOut[60] - PmQf[1] * r_222_1 + 2 * iz2 * (pOut[43] - riz * r_212_1);
   pOut[99] = PmA[1] * pOut[61] - PmQf[1] * r_204_1;
   pOut[100] = PmA[1] * pOut[62] - PmQf[1] * r_060_1 + 6 * iz2 * (pOut[44] - riz * r_050_1);
   pOut[101] = PmA[1] * pOut[63] - PmQf[1] * r_042_1 + 4 * iz2 * (pOut[45] - riz * r_032_1);
   pOut[102] = PmA[1] * pOut[64] - PmQf[1] * r_024_1 + 2 * iz2 * (pOut[46] - riz * r_014_1);
   pOut[103] = PmA[2] * pOut[83] - PmQf[2] * r_015_1 + 5 * iz2 * (pOut[46] - riz * r_014_1);
   pOut[104] = PmA[0] * pOut[72] - PmQf[0] * r_501_1 + 5 * iz2 * (pOut[47] - riz * r_401_1);
   pOut[105] = PmA[2] * pOut[57] - PmQf[2] * r_420_1;
   pOut[106] = PmA[0] * pOut[74] - PmQf[0] * r_303_1 + 3 * iz2 * (pOut[49] - riz * r_203_1);
   pOut[107] = PmA[2] * pOut[59] - PmQf[2] * r_240_1;
   pOut[108] = PmA[0] * pOut[76] - PmQf[0] * r_123_1 + iz2 * (pOut[51] - riz * r_023_1);
   pOut[109] = PmA[2] * pOut[61] - PmQf[2] * r_204_1 + 4 * iz2 * (pOut[49] - riz * r_203_1);
   pOut[110] = PmA[2] * pOut[62] - PmQf[2] * r_060_1;
   pOut[111] = PmA[2] * pOut[63] - PmQf[2] * r_042_1 + 2 * iz2 * (pOut[50] - riz * r_041_1);
   pOut[112] = PmA[2] * pOut[64] - PmQf[2] * r_024_1 + 4 * iz2 * (pOut[51] - riz * r_023_1);
   pOut[113] = PmA[2] * pOut[65] - PmQf[2] * r_006_1 + 6 * iz2 * (pOut[52] - riz * r_005_1);
   pOut[114] = PmA[1] * pOut[72] - PmQf[1] * r_501_1;
   pOut[115] = PmA[2] * pOut[67] - PmQf[2] * r_330_1;
   pOut[116] = PmA[1] * pOut[74] - PmQf[1] * r_303_1;
   pOut[117] = PmA[2] * pOut[69] - PmQf[2] * r_150_1;
   pOut[118] = PmA[0] * pOut[82] - PmQf[0] * r_033_1;
   pOut[119] = PmA[0] * pOut[83] - PmQf[0] * r_015_1;
   // 1232 flops, 764 mops, 0.98kb stack
}











static void OsrrC43( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   const double c0 = 6.1237243569579447e-01;
   const double c1 = 2.4494897427831779;
   const double c2 = 1.5;
   const double c3 = 2.3717082451262841;
   const double c4 = 7.9056941504209477e-01;
   const double c5 = 3.8729833462074166;
   const double c6 = 1.9364916731037083;

   FReal2 a_000[9]; ShTr4(a_000, pIn[20], pIn[21], pIn[22], pIn[23], pIn[24], pIn[25], pIn[26], pIn[27], pIn[28], pIn[29], pIn[30], pIn[31], pIn[32], pIn[33], pIn[34]);
   FReal2 a_100[9]; ShTr4(a_100, pIn[35], pIn[36], pIn[37], pIn[38], pIn[39], pIn[40], pIn[41], pIn[42], pIn[43], pIn[47], pIn[48], pIn[49], pIn[53], pIn[54], pIn[55]);
   FReal2 a_010[9]; ShTr4(a_010, pIn[41], pIn[42], pIn[43], pIn[44], pIn[45], pIn[46], pIn[36], pIn[38], pIn[39], pIn[53], pIn[54], pIn[55], pIn[48], pIn[50], pIn[51]);
   FReal2 a_001[9]; ShTr4(a_001, pIn[47], pIn[48], pIn[49], pIn[50], pIn[51], pIn[52], pIn[53], pIn[54], pIn[55], pIn[37], pIn[39], pIn[40], pIn[43], pIn[45], pIn[46]);
   FReal2 a_200[9]; ShTr4(a_200, pIn[56], pIn[57], pIn[58], pIn[59], pIn[60], pIn[61], pIn[66], pIn[67], pIn[68], pIn[72], pIn[73], pIn[74], pIn[78], pIn[79], pIn[80]);
   FReal2 a_020[9]; ShTr4(a_020, pIn[57], pIn[59], pIn[60], pIn[62], pIn[63], pIn[64], pIn[67], pIn[69], pIn[70], pIn[73], pIn[75], pIn[76], pIn[79], pIn[81], pIn[82]);
   FReal2 a_002[9]; ShTr4(a_002, pIn[58], pIn[60], pIn[61], pIn[63], pIn[64], pIn[65], pIn[68], pIn[70], pIn[71], pIn[74], pIn[76], pIn[77], pIn[80], pIn[82], pIn[83]);
   FReal2 a_110[9]; ShTr4(a_110, pIn[66], pIn[67], pIn[68], pIn[69], pIn[70], pIn[71], pIn[57], pIn[59], pIn[60], pIn[78], pIn[79], pIn[80], pIn[73], pIn[75], pIn[76]);
   FReal2 a_101[9]; ShTr4(a_101, pIn[72], pIn[73], pIn[74], pIn[75], pIn[76], pIn[77], pIn[78], pIn[79], pIn[80], pIn[58], pIn[60], pIn[61], pIn[68], pIn[70], pIn[71]);
   FReal2 a_011[9]; ShTr4(a_011, pIn[78], pIn[79], pIn[80], pIn[81], pIn[82], pIn[83], pIn[73], pIn[75], pIn[76], pIn[68], pIn[70], pIn[71], pIn[60], pIn[63], pIn[64]);
   FReal2 a_300[9]; ShTr4(a_300, pIn[84], pIn[85], pIn[86], pIn[87], pIn[88], pIn[89], pIn[94], pIn[95], pIn[96], pIn[104], pIn[105], pIn[106], pIn[114], pIn[115], pIn[116]);
   FReal2 a_120[9]; ShTr4(a_120, pIn[85], pIn[87], pIn[88], pIn[90], pIn[91], pIn[92], pIn[95], pIn[97], pIn[98], pIn[105], pIn[107], pIn[108], pIn[115], pIn[117], pIn[118]);
   FReal2 a_102[9]; ShTr4(a_102, pIn[86], pIn[88], pIn[89], pIn[91], pIn[92], pIn[93], pIn[96], pIn[98], pIn[99], pIn[106], pIn[108], pIn[109], pIn[116], pIn[118], pIn[119]);
   FReal2 a_210[9]; ShTr4(a_210, pIn[94], pIn[95], pIn[96], pIn[97], pIn[98], pIn[99], pIn[85], pIn[87], pIn[88], pIn[114], pIn[115], pIn[116], pIn[105], pIn[107], pIn[108]);
   FReal2 a_030[9]; ShTr4(a_030, pIn[95], pIn[97], pIn[98], pIn[100], pIn[101], pIn[102], pIn[87], pIn[90], pIn[91], pIn[115], pIn[117], pIn[118], pIn[107], pIn[110], pIn[111]);
   FReal2 a_012[9]; ShTr4(a_012, pIn[96], pIn[98], pIn[99], pIn[101], pIn[102], pIn[103], pIn[88], pIn[91], pIn[92], pIn[116], pIn[118], pIn[119], pIn[108], pIn[111], pIn[112]);
   FReal2 a_201[9]; ShTr4(a_201, pIn[104], pIn[105], pIn[106], pIn[107], pIn[108], pIn[109], pIn[114], pIn[115], pIn[116], pIn[86], pIn[88], pIn[89], pIn[96], pIn[98], pIn[99]);
   FReal2 a_021[9]; ShTr4(a_021, pIn[105], pIn[107], pIn[108], pIn[110], pIn[111], pIn[112], pIn[115], pIn[117], pIn[118], pIn[88], pIn[91], pIn[92], pIn[98], pIn[101], pIn[102]);
   FReal2 a_003[9]; ShTr4(a_003, pIn[106], pIn[108], pIn[109], pIn[111], pIn[112], pIn[113], pIn[116], pIn[118], pIn[119], pIn[89], pIn[92], pIn[93], pIn[99], pIn[102], pIn[103]);
   FReal2 a_111[9]; ShTr4(a_111, pIn[114], pIn[115], pIn[116], pIn[117], pIn[118], pIn[119], pIn[105], pIn[107], pIn[108], pIn[96], pIn[98], pIn[99], pIn[88], pIn[91], pIn[92]);
   // ShTr: [940 flops, 1600 mops]

   for ( unsigned ica = 0; ica < 9; ++ ica ) {
      FReal2 r_000_100 = a_100[ica] - BmA[0] * a_000[ica];
      FReal2 r_100_100 = a_200[ica] - BmA[0] * a_100[ica];
      FReal2 r_010_100 = a_110[ica] - BmA[0] * a_010[ica];
      FReal2 r_001_100 = a_101[ica] - BmA[0] * a_001[ica];
      FReal2 r_200_100 = a_300[ica] - BmA[0] * a_200[ica];
      FReal2 r_110_100 = a_210[ica] - BmA[0] * a_110[ica];
      FReal2 r_101_100 = a_201[ica] - BmA[0] * a_101[ica];
      FReal2 r_011_100 = a_111[ica] - BmA[0] * a_011[ica];
      FReal2 r_000_010 = a_010[ica] - BmA[1] * a_000[ica];
      FReal2 r_100_010 = a_110[ica] - BmA[1] * a_100[ica];
      FReal2 r_010_010 = a_020[ica] - BmA[1] * a_010[ica];
      FReal2 r_001_010 = a_011[ica] - BmA[1] * a_001[ica];
      FReal2 r_020_010 = a_030[ica] - BmA[1] * a_020[ica];
      FReal2 r_110_010 = a_120[ica] - BmA[1] * a_110[ica];
      FReal2 r_011_010 = a_021[ica] - BmA[1] * a_011[ica];
      FReal2 r_000_001 = a_001[ica] - BmA[2] * a_000[ica];
      FReal2 r_100_001 = a_101[ica] - BmA[2] * a_100[ica];
      FReal2 r_010_001 = a_011[ica] - BmA[2] * a_010[ica];
      FReal2 r_001_001 = a_002[ica] - BmA[2] * a_001[ica];
      FReal2 r_002_001 = a_003[ica] - BmA[2] * a_002[ica];
      FReal2 r_101_001 = a_102[ica] - BmA[2] * a_101[ica];
      FReal2 r_011_001 = a_012[ica] - BmA[2] * a_011[ica];
      FReal2 r_000_200 = r_100_100 - BmA[0] * r_000_100;
      FReal2 r_100_200 = r_200_100 - BmA[0] * r_100_100;
      FReal2 r_010_200 = r_110_100 - BmA[0] * r_010_100;
      FReal2 r_001_200 = r_101_100 - BmA[0] * r_001_100;
      FReal2 r_000_020 = r_010_010 - BmA[1] * r_000_010;
      FReal2 r_100_020 = r_110_010 - BmA[1] * r_100_010;
      FReal2 r_010_020 = r_020_010 - BmA[1] * r_010_010;
      FReal2 r_001_020 = r_011_010 - BmA[1] * r_001_010;
      FReal2 r_000_002 = r_001_001 - BmA[2] * r_000_001;
      FReal2 r_100_002 = r_101_001 - BmA[2] * r_100_001;
      FReal2 r_010_002 = r_011_001 - BmA[2] * r_010_001;
      FReal2 r_001_002 = r_002_001 - BmA[2] * r_001_001;
      FReal2 r_000_101 = r_001_100 - BmA[2] * r_000_100;
      FReal2 r_010_101 = r_011_100 - BmA[2] * r_010_100;
      FReal2 b_300 = r_100_200 - BmA[0] * r_000_200;
      FReal2 b_120 = r_100_020 - BmA[0] * r_000_020;
      FReal2 b_102 = r_100_002 - BmA[0] * r_000_002;
      FReal2 b_210 = r_010_200 - BmA[1] * r_000_200;
      FReal2 b_030 = r_010_020 - BmA[1] * r_000_020;
      FReal2 b_012 = r_010_002 - BmA[1] * r_000_002;
      FReal2 b_201 = r_001_200 - BmA[2] * r_000_200;
      FReal2 b_021 = r_001_020 - BmA[2] * r_000_020;
      FReal2 b_003 = r_001_002 - BmA[2] * r_000_002;
      FReal2 b_111 = r_010_101 - BmA[1] * r_000_101;
      pOut[ica*sa + 0*sb] += -c0*b_300 + c1*b_102 - c0*b_120;
      pOut[ica*sa + 1*sb] += c1*b_012 - c0*b_030 - c0*b_210;
      pOut[ica*sa + 2*sb] += -c2*b_021 - c2*b_201 + b_003;
      pOut[ica*sa + 3*sb] += -c3*b_120 + c4*b_300;
      pOut[ica*sa + 4*sb] += c5*b_111;
      pOut[ica*sa + 5*sb] += -c4*b_030 + c3*b_210;
      pOut[ica*sa + 6*sb] += -c6*b_021 + c6*b_201;
      // 9 * [92 flops, 138 mops]
   }
   // total: [1768 flops, 2842 mops]
}

// In: Gm[0 .. 8] (inclusive), Out: [r]^0, ordered as AngularComps(8)
// (moment 8, 45 entries)
static void OsrrA8( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const PmA[3], FReal1 const PmQ[3], double rho, double InvEta ) AIC_NO_THROW
{
   double
      riz = rho * InvEta,
      iz2 = 0.5 * InvEta,
      PmQf[3] = { riz*PmQ[0], riz*PmQ[1], riz*PmQ[2] };
   pOut[0] = pIn[0];
   // L = 1
   pOut[1] = PmA[0] * pIn[0] - PmQf[0] * pIn[1];
   FReal2 r_100_1 = PmA[0] * pIn[1] - PmQf[0] * pIn[2];
   FReal2 r_100_2 = PmA[0] * pIn[2] - PmQf[0] * pIn[3];
   FReal2 r_100_3 = PmA[0] * pIn[3] - PmQf[0] * pIn[4];
   FReal2 r_100_4 = PmA[0] * pIn[4] - PmQf[0] * pIn[5];
   FReal2 r_100_5 = PmA[0] * pIn[5] - PmQf[0] * pIn[6];
   FReal2 r_100_6 = PmA[0] * pIn[6] - PmQf[0] * pIn[7];
   FReal2 r_100_7 = PmA[0] * pIn[7] - PmQf[0] * pIn[8];
   pOut[2] = PmA[1] * pIn[0] - PmQf[1] * pIn[1];
   FReal2 r_010_1 = PmA[1] * pIn[1] - PmQf[1] * pIn[2];
   FReal2 r_010_2 = PmA[1] * pIn[2] - PmQf[1] * pIn[3];
   FReal2 r_010_3 = PmA[1] * pIn[3] - PmQf[1] * pIn[4];
   FReal2 r_010_4 = PmA[1] * pIn[4] - PmQf[1] * pIn[5];
   FReal2 r_010_5 = PmA[1] * pIn[5] - PmQf[1] * pIn[6];
   FReal2 r_010_6 = PmA[1] * pIn[6] - PmQf[1] * pIn[7];
   FReal2 r_010_7 = PmA[1] * pIn[7] - PmQf[1] * pIn[8];
   pOut[3] = PmA[2] * pIn[0] - PmQf[2] * pIn[1];
   FReal2 r_001_1 = PmA[2] * pIn[1] - PmQf[2] * pIn[2];
   FReal2 r_001_2 = PmA[2] * pIn[2] - PmQf[2] * pIn[3];
   FReal2 r_001_3 = PmA[2] * pIn[3] - PmQf[2] * pIn[4];
   FReal2 r_001_4 = PmA[2] * pIn[4] - PmQf[2] * pIn[5];
   FReal2 r_001_5 = PmA[2] * pIn[5] - PmQf[2] * pIn[6];
   FReal2 r_001_6 = PmA[2] * pIn[6] - PmQf[2] * pIn[7];
   FReal2 r_001_7 = PmA[2] * pIn[7] - PmQf[2] * pIn[8];

   // L = 2
   pOut[4] = PmA[0] * pOut[1] - PmQf[0] * r_100_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_200_1 = PmA[0] * r_100_1 - PmQf[0] * r_100_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_200_2 = PmA[0] * r_100_2 - PmQf[0] * r_100_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_200_3 = PmA[0] * r_100_3 - PmQf[0] * r_100_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_200_4 = PmA[0] * r_100_4 - PmQf[0] * r_100_5 + iz2 * (pIn[4] - riz * pIn[5]);
   FReal2 r_200_5 = PmA[0] * r_100_5 - PmQf[0] * r_100_6 + iz2 * (pIn[5] - riz * pIn[6]);
   FReal2 r_200_6 = PmA[0] * r_100_6 - PmQf[0] * r_100_7 + iz2 * (pIn[6] - riz * pIn[7]);
   pOut[5] = PmA[1] * pOut[2] - PmQf[1] * r_010_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_020_1 = PmA[1] * r_010_1 - PmQf[1] * r_010_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_020_2 = PmA[1] * r_010_2 - PmQf[1] * r_010_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_020_3 = PmA[1] * r_010_3 - PmQf[1] * r_010_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_020_4 = PmA[1] * r_010_4 - PmQf[1] * r_010_5 + iz2 * (pIn[4] - riz * pIn[5]);
   FReal2 r_020_5 = PmA[1] * r_010_5 - PmQf[1] * r_010_6 + iz2 * (pIn[5] - riz * pIn[6]);
   FReal2 r_020_6 = PmA[1] * r_010_6 - PmQf[1] * r_010_7 + iz2 * (pIn[6] - riz * pIn[7]);
   pOut[6] = PmA[2] * pOut[3] - PmQf[2] * r_001_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_002_1 = PmA[2] * r_001_1 - PmQf[2] * r_001_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_002_2 = PmA[2] * r_001_2 - PmQf[2] * r_001_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_002_3 = PmA[2] * r_001_3 - PmQf[2] * r_001_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_002_4 = PmA[2] * r_001_4 - PmQf[2] * r_001_5 + iz2 * (pIn[4] - riz * pIn[5]);
   FReal2 r_002_5 = PmA[2] * r_001_5 - PmQf[2] * r_001_6 + iz2 * (pIn[5] - riz * pIn[6]);
   FReal2 r_002_6 = PmA[2] * r_001_6 - PmQf[2] * r_001_7 + iz2 * (pIn[6] - riz * pIn[7]);
   pOut[7] = PmA[0] * pOut[2] - PmQf[0] * r_010_1;
   pOut[8] = PmA[2] * pOut[1] - PmQf[2] * r_100_1;
   FReal2 r_101_1 = PmA[2] * r_100_1 - PmQf[2] * r_100_2;
   pOut[9] = PmA[1] * pOut[3] - PmQf[1] * r_001_1;

   // L = 3
   pOut[10] = PmA[0] * pOut[4] - PmQf[0] * r_200_1 + 2 * iz2 * (pOut[1] - riz * r_100_1);
   FReal2 r_300_1 = PmA[0] * r_200_1 - PmQf[0] * r_200_2 + 2 * iz2 * (r_100_1 - riz * r_100_2);
   FReal2 r_300_2 = PmA[0] * r_200_2 - PmQf[0] * r_200_3 + 2 * iz2 * (r_100_2 - riz * r_100_3);
   FReal2 r_300_3 = PmA[0] * r_200_3 - PmQf[0] * r_200_4 + 2 * iz2 * (r_100_3 - riz * r_100_4);
   FReal2 r_300_4 = PmA[0] * r_200_4 - PmQf[0] * r_200_5 + 2 * iz2 * (r_100_4 - riz * r_100_5);
   FReal2 r_300_5 = PmA[0] * r_200_5 - PmQf[0] * r_200_6 + 2 * iz2 * (r_100_5 - riz * r_100_6);
   pOut[11] = PmA[0] * pOut[5] - PmQf[0] * r_020_1;
   FReal2 r_120_1 = PmA[0] * r_020_1 - PmQf[0] * r_020_2;
   FReal2 r_120_2 = PmA[0] * r_020_2 - PmQf[0] * r_020_3;
   FReal2 r_120_3 = PmA[0] * r_020_3 - PmQf[0] * r_020_4;
   FReal2 r_120_4 = PmA[0] * r_020_4 - PmQf[0] * r_020_5;
   pOut[12] = PmA[0] * pOut[6] - PmQf[0] * r_002_1;
   pOut[13] = PmA[1] * pOut[4] - PmQf[1] * r_200_1;
   pOut[14] = PmA[1] * pOut[5] - PmQf[1] * r_020_1 + 2 * iz2 * (pOut[2] - riz * r_010_1);
   FReal2 r_030_1 = PmA[1] * r_020_1 - PmQf[1] * r_020_2 + 2 * iz2 * (r_010_1 - riz * r_010_2);
   FReal2 r_030_2 = PmA[1] * r_020_2 - PmQf[1] * r_020_3 + 2 * iz2 * (r_010_2 - riz * r_010_3);
   FReal2 r_030_3 = PmA[1] * r_020_3 - PmQf[1] * r_020_4 + 2 * iz2 * (r_010_3 - riz * r_010_4);
   FReal2 r_030_4 = PmA[1] * r_020_4 - PmQf[1] * r_020_5 + 2 * iz2 * (r_010_4 - riz * r_010_5);
   FReal2 r_030_5 = PmA[1] * r_020_5 - PmQf[1] * r_020_6 + 2 * iz2 * (r_010_5 - riz * r_010_6);
   pOut[15] = PmA[1] * pOut[6] - PmQf[1] * r_002_1;
   FReal2 r_012_1 = PmA[1] * r_002_1 - PmQf[1] * r_002_2;
   FReal2 r_012_2 = PmA[1] * r_002_2 - PmQf[1] * r_002_3;
   FReal2 r_012_3 = PmA[1] * r_002_3 - PmQf[1] * r_002_4;
   pOut[16] = PmA[0] * pOut[8] - PmQf[0] * r_101_1 + iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_201_1 = PmA[2] * r_200_1 - PmQf[2] * r_200_2;
   FReal2 r_201_2 = PmA[2] * r_200_2 - PmQf[2] * r_200_3;
   FReal2 r_201_3 = PmA[2] * r_200_3 - PmQf[2] * r_200_4;
   pOut[17] = PmA[2] * pOut[5] - PmQf[2] * r_020_1;
   pOut[18] = PmA[2] * pOut[6] - PmQf[2] * r_002_1 + 2 * iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_003_1 = PmA[2] * r_002_1 - PmQf[2] * r_002_2 + 2 * iz2 * (r_001_1 - riz * r_001_2);
   FReal2 r_003_2 = PmA[2] * r_002_2 - PmQf[2] * r_002_3 + 2 * iz2 * (r_001_2 - riz * r_001_3);
   FReal2 r_003_3 = PmA[2] * r_002_3 - PmQf[2] * r_002_4 + 2 * iz2 * (r_001_3 - riz * r_001_4);
   FReal2 r_003_4 = PmA[2] * r_002_4 - PmQf[2] * r_002_5 + 2 * iz2 * (r_001_4 - riz * r_001_5);
   FReal2 r_003_5 = PmA[2] * r_002_5 - PmQf[2] * r_002_6 + 2 * iz2 * (r_001_5 - riz * r_001_6);
   pOut[19] = PmA[1] * pOut[8] - PmQf[1] * r_101_1;

   // L = 4
   pOut[20] = PmA[0] * pOut[10] - PmQf[0] * r_300_1 + 3 * iz2 * (pOut[4] - riz * r_200_1);
   FReal2 r_400_1 = PmA[0] * r_300_1 - PmQf[0] * r_300_2 + 3 * iz2 * (r_200_1 - riz * r_200_2);
   FReal2 r_400_2 = PmA[0] * r_300_2 - PmQf[0] * r_300_3 + 3 * iz2 * (r_200_2 - riz * r_200_3);
   FReal2 r_400_3 = PmA[0] * r_300_3 - PmQf[0] * r_300_4 + 3 * iz2 * (r_200_3 - riz * r_200_4);
   FReal2 r_400_4 = PmA[0] * r_300_4 - PmQf[0] * r_300_5 + 3 * iz2 * (r_200_4 - riz * r_200_5);
   pOut[21] = PmA[0] * pOut[11] - PmQf[0] * r_120_1 + iz2 * (pOut[5] - riz * r_020_1);
   FReal2 r_220_1 = PmA[0] * r_120_1 - PmQf[0] * r_120_2 + iz2 * (r_020_1 - riz * r_020_2);
   FReal2 r_220_2 = PmA[0] * r_120_2 - PmQf[0] * r_120_3 + iz2 * (r_020_2 - riz * r_020_3);
   FReal2 r_220_3 = PmA[0] * r_120_3 - PmQf[0] * r_120_4 + iz2 * (r_020_3 - riz * r_020_4);
   pOut[22] = PmA[2] * pOut[16] - PmQf[2] * r_201_1 + iz2 * (pOut[4] - riz * r_200_1);
   FReal2 r_202_1 = PmA[2] * r_201_1 - PmQf[2] * r_201_2 + iz2 * (r_200_1 - riz * r_200_2);
   FReal2 r_202_2 = PmA[2] * r_201_2 - PmQf[2] * r_201_3 + iz2 * (r_200_2 - riz * r_200_3);
   pOut[23] = PmA[1] * pOut[14] - PmQf[1] * r_030_1 + 3 * iz2 * (pOut[5] - riz * r_020_1);
   FReal2 r_040_1 = PmA[1] * r_030_1 - PmQf[1] * r_030_2 + 3 * iz2 * (r_020_1 - riz * r_020_2);
   FReal2 r_040_2 = PmA[1] * r_030_2 - PmQf[1] * r_030_3 + 3 * iz2 * (r_020_2 - riz * r_020_3);
   FReal2 r_040_3 = PmA[1] * r_030_3 - PmQf[1] * r_030_4 + 3 * iz2 * (r_020_3 - riz * r_020_4);
   FReal2 r_040_4 = PmA[1] * r_030_4 - PmQf[1] * r_030_5 + 3 * iz2 * (r_020_4 - riz * r_020_5);
   pOut[24] = PmA[1] * pOut[15] - PmQf[1] * r_012_1 + iz2 * (pOut[6] - riz * r_002_1);
   FReal2 r_022_1 = PmA[1] * r_012_1 - PmQf[1] * r_012_2 + iz2 * (r_002_1 - riz * r_002_2);
   FReal2 r_022_2 = PmA[1] * r_012_2 - PmQf[1] * r_012_3 + iz2 * (r_002_2 - riz * r_002_3);
   pOut[25] = PmA[2] * pOut[18] - PmQf[2] * r_003_1 + 3 * iz2 * (pOut[6] - riz * r_002_1);
   FReal2 r_004_1 = PmA[2] * r_003_1 - PmQf[2] * r_003_2 + 3 * iz2 * (r_002_1 - riz * r_002_2);
   FReal2 r_004_2 = PmA[2] * r_003_2 - PmQf[2] * r_003_3 + 3 * iz2 * (r_002_2 - riz * r_002_3);
   FReal2 r_004_3 = PmA[2] * r_003_3 - PmQf[2] * r_003_4 + 3 * iz2 * (r_002_3 - riz * r_002_4);
   FReal2 r_004_4 = PmA[2] * r_003_4 - PmQf[2] * r_003_5 + 3 * iz2 * (r_002_4 - riz * r_002_5);
   pOut[26] = PmA[1] * pOut[10] - PmQf[1] * r_300_1;
   pOut[27] = PmA[0] * pOut[14] - PmQf[0] * r_030_1;
   FReal2 r_130_1 = PmA[0] * r_030_1 - PmQf[0] * r_030_2;
   FReal2 r_130_2 = PmA[0] * r_030_2 - PmQf[0] * r_030_3;
   FReal2 r_130_3 = PmA[0] * r_030_3 - PmQf[0] * r_030_4;
   FReal2 r_130_4 = PmA[0] * r_030_4 - PmQf[0] * r_030_5;
   pOut[28] = PmA[0] * pOut[15] - PmQf[0] * r_012_1;
   pOut[29] = PmA[0] * pOut[16] - PmQf[0] * r_201_1 + 2 * iz2 * (pOut[8] - riz * r_101_1);
   FReal2 r_301_1 = PmA[2] * r_300_1 - PmQf[2] * r_300_2;
   FReal2 r_301_2 = PmA[2] * r_300_2 - PmQf[2] * r_300_3;
   FReal2 r_301_3 = PmA[2] * r_300_3 - PmQf[2] * r_300_4;
   FReal2 r_301_4 = PmA[2] * r_300_4 - PmQf[2] * r_300_5;
   pOut[30] = PmA[2] * pOut[11] - PmQf[2] * r_120_1;
   pOut[31] = PmA[0] * pOut[18] - PmQf[0] * r_003_1;
   FReal2 r_103_2 = PmA[0] * r_003_2 - PmQf[0] * r_003_3;
   FReal2 r_103_3 = PmA[0] * r_003_3 - PmQf[0] * r_003_4;
   pOut[32] = PmA[1] * pOut[16] - PmQf[1] * r_201_1;
   pOut[33] = PmA[2] * pOut[14] - PmQf[2] * r_030_1;
   FReal2 r_031_2 = PmA[2] * r_030_2 - PmQf[2] * r_030_3;
   FReal2 r_031_3 = PmA[2] * r_030_3 - PmQf[2] * r_030_4;
   pOut[34] = PmA[1] * pOut[18] - PmQf[1] * r_003_1;
   FReal2 r_013_1 = PmA[1] * r_003_1 - PmQf[1] * r_003_2;
   FReal2 r_013_2 = PmA[1] * r_003_2 - PmQf[1] * r_003_3;
   FReal2 r_013_3 = PmA[1] * r_003_3 - PmQf[1] * r_003_4;
   FReal2 r_013_4 = PmA[1] * r_003_4 - PmQf[1] * r_003_5;

   // L = 5
   pOut[35] = PmA[0] * pOut[20] - PmQf[0] * r_400_1 + 4 * iz2 * (pOut[10] - riz * r_300_1);
   FReal2 r_500_1 = PmA[0] * r_400_1 - PmQf[0] * r_400_2 + 4 * iz2 * (r_300_1 - riz * r_300_2);
   FReal2 r_500_2 = PmA[0] * r_400_2 - PmQf[0] * r_400_3 + 4 * iz2 * (r_300_2 - riz * r_300_3);
   FReal2 r_500_3 = PmA[0] * r_400_3 - PmQf[0] * r_400_4 + 4 * iz2 * (r_300_3 - riz * r_300_4);
   pOut[36] = PmA[0] * pOut[21] - PmQf[0] * r_220_1 + 2 * iz2 * (pOut[11] - riz * r_120_1);
   FReal2 r_320_1 = PmA[0] * r_220_1 - PmQf[0] * r_220_2 + 2 * iz2 * (r_120_1 - riz * r_120_2);
   FReal2 r_320_2 = PmA[0] * r_220_2 - PmQf[0] * r_220_3 + 2 * iz2 * (r_120_2 - riz * r_120_3);
   pOut[37] = PmA[2] * pOut[29] - PmQf[2] * r_301_1 + iz2 * (pOut[10] - riz * r_300_1);
   FReal2 r_302_1 = PmA[2] * r_301_1 - PmQf[2] * r_301_2 + iz2 * (r_300_1 - riz * r_300_2);
   FReal2 r_302_2 = PmA[2] * r_301_2 - PmQf[2] * r_301_3 + iz2 * (r_300_2 - riz * r_300_3);
   FReal2 r_302_3 = PmA[2] * r_301_3 - PmQf[2] * r_301_4 + iz2 * (r_300_3 - riz * r_300_4);
   pOut[38] = PmA[1] * pOut[27] - PmQf[1] * r_130_1 + 3 * iz2 * (pOut[11] - riz * r_120_1);
   FReal2 r_140_1 = PmA[1] * r_130_1 - PmQf[1] * r_130_2 + 3 * iz2 * (r_120_1 - riz * r_120_2);
   pOut[39] = PmA[0] * pOut[24] - PmQf[0] * r_022_1;
   pOut[40] = PmA[0] * pOut[25] - PmQf[0] * r_004_1;
   pOut[41] = PmA[1] * pOut[20] - PmQf[1] * r_400_1;
   pOut[42] = PmA[0] * pOut[27] - PmQf[0] * r_130_1 + iz2 * (pOut[14] - riz * r_030_1);
   FReal2 r_230_1 = PmA[0] * r_130_1 - PmQf[0] * r_130_2 + iz2 * (r_030_1 - riz * r_030_2);
   FReal2 r_230_2 = PmA[0] * r_130_2 - PmQf[0] * r_130_3 + iz2 * (r_030_2 - riz * r_030_3);
   FReal2 r_230_3 = PmA[0] * r_130_3 - PmQf[0] * r_130_4 + iz2 * (r_030_3 - riz * r_030_4);
   pOut[43] = PmA[1] * pOut[22] - PmQf[1] * r_202_1;
   FReal2 r_212_1 = PmA[1] * r_202_1 - PmQf[1] * r_202_2;
   pOut[44] = PmA[1] * pOut[23] - PmQf[1] * r_040_1 + 4 * iz2 * (pOut[14] - riz * r_030_1);
   FReal2 r_050_1 = PmA[1] * r_040_1 - PmQf[1] * r_040_2 + 4 * iz2 * (r_030_1 - riz * r_030_2);
   FReal2 r_050_2 = PmA[1] * r_040_2 - PmQf[1] * r_040_3 + 4 * iz2 * (r_030_2 - riz * r_030_3);
   FReal2 r_050_3 = PmA[1] * r_040_3 - PmQf[1] * r_040_4 + 4 * iz2 * (r_030_3 - riz * r_030_4);
   pOut[45] = PmA[1] * pOut[24] - PmQf[1] * r_022_1 + 2 * iz2 * (pOut[15] - riz * r_012_1);
   FReal2 r_032_1 = PmA[1] * r_022_1 - PmQf[1] * r_022_2 + 2 * iz2 * (r_012_1 - riz * r_012_2);
   FReal2 r_032_2 = PmA[2] * r_031_2 - PmQf[2] * r_031_3 + iz2 * (r_030_2 - riz * r_030_3);
   pOut[46] = PmA[2] * pOut[34] - PmQf[2] * r_013_1 + 3 * iz2 * (pOut[15] - riz * r_012_1);
   FReal2 r_014_1 = PmA[2] * r_013_1 - PmQf[2] * r_013_2 + 3 * iz2 * (r_012_1 - riz * r_012_2);
   pOut[47] = PmA[0] * pOut[29] - PmQf[0] * r_301_1 + 3 * iz2 * (pOut[16] - riz * r_201_1);
   FReal2 r_401_1 = PmA[0] * r_301_1 - PmQf[0] * r_301_2 + 3 * iz2 * (r_201_1 - riz * r_201_2);
   pOut[48] = PmA[2] * pOut[21] - PmQf[2] * r_220_1;
   FReal2 r_221_1 = PmA[2] * r_220_1 - PmQf[2] * r_220_2;
   FReal2 r_221_2 = PmA[2] * r_220_2 - PmQf[2] * r_220_3;
   pOut[49] = PmA[2] * pOut[22] - PmQf[2] * r_202_1 + 2 * iz2 * (pOut[16] - riz * r_201_1);
   FReal2 r_203_1 = PmA[2] * r_202_1 - PmQf[2] * r_202_2 + 2 * iz2 * (r_201_1 - riz * r_201_2);
   FReal2 r_203_2 = PmA[0] * r_103_2 - PmQf[0] * r_103_3 + iz2 * (r_003_2 - riz * r_003_3);
   pOut[50] = PmA[2] * pOut[23] - PmQf[2] * r_040_1;
   pOut[51] = PmA[1] * pOut[34] - PmQf[1] * r_013_1 + iz2 * (pOut[18] - riz * r_003_1);
   FReal2 r_023_1 = PmA[1] * r_013_1 - PmQf[1] * r_013_2 + iz2 * (r_003_1 - riz * r_003_2);
   FReal2 r_023_2 = PmA[1] * r_013_2 - PmQf[1] * r_013_3 + iz2 * (r_003_2 - riz * r_003_3);
   FReal2 r_023_3 = PmA[1] * r_013_3 - PmQf[1] * r_013_4 + iz2 * (r_003_3 - riz * r_003_4);
   pOut[52] = PmA[2] * pOut[25] - PmQf[2] * r_004_1 + 4 * iz2 * (pOut[18] - riz * r_003_1);
   FReal2 r_005_1 = PmA[2] * r_004_1 - PmQf[2] * r_004_2 + 4 * iz2 * (r_003_1 - riz * r_003_2);
   FReal2 r_005_2 = PmA[2] * r_004_2 - PmQf[2] * r_004_3 + 4 * iz2 * (r_003_2 - riz * r_003_3);
   FReal2 r_005_3 = PmA[2] * r_004_3 - PmQf[2] * r_004_4 + 4 * iz2 * (r_003_3 - riz * r_003_4);
   pOut[53] = PmA[1] * pOut[29] - PmQf[1] * r_301_1;
   pOut[54] = PmA[2] * pOut[27] - PmQf[2] * r_130_1;
   pOut[55] = PmA[0] * pOut[34] - PmQf[0] * r_013_1;

   // L = 6
   pOut[56] = PmA[0] * pOut[35] - PmQf[0] * r_500_1 + 5 * iz2 * (pOut[20] - riz * r_400_1);
   FReal2 r_600_1 = PmA[0] * r_500_1 - PmQf[0] * r_500_2 + 5 * iz2 * (r_400_1 - riz * r_400_2);
   FReal2 r_600_2 = PmA[0] * r_500_2 - PmQf[0] * r_500_3 + 5 * iz2 * (r_400_2 - riz * r_400_3);
   pOut[57] = PmA[0] * pOut[36] - PmQf[0] * r_320_1 + 3 * iz2 * (pOut[21] - riz * r_220_1);
   pOut[58] = PmA[0] * pOut[37] - PmQf[0] * r_302_1 + 3 * iz2 * (pOut[22] - riz * r_202_1);
   FReal2 r_402_1 = PmA[0] * r_302_1 - PmQf[0] * r_302_2 + 3 * iz2 * (r_202_1 - riz * r_202_2);
   pOut[59] = PmA[0] * pOut[38] - PmQf[0] * r_140_1 + iz2 * (pOut[23] - riz * r_040_1);
   FReal2 r_240_1 = PmA[1] * r_230_1 - PmQf[1] * r_230_2 + 3 * iz2 * (r_220_1 - riz * r_220_2);
   pOut[60] = PmA[1] * pOut[43] - PmQf[1] * r_212_1 + iz2 * (pOut[22] - riz * r_202_1);
   FReal2 r_222_1 = PmA[2] * r_221_1 - PmQf[2] * r_221_2 + iz2 * (r_220_1 - riz * r_220_2);
   pOut[61] = PmA[2] * pOut[49] - PmQf[2] * r_203_1 + 3 * iz2 * (pOut[22] - riz * r_202_1);
   FReal2 r_204_1 = PmA[2] * r_203_1 - PmQf[2] * r_203_2 + 3 * iz2 * (r_202_1 - riz * r_202_2);
   pOut[62] = PmA[1] * pOut[44] - PmQf[1] * r_050_1 + 5 * iz2 * (pOut[23] - riz * r_040_1);
   FReal2 r_060_1 = PmA[1] * r_050_1 - PmQf[1] * r_050_2 + 5 * iz2 * (r_040_1 - riz * r_040_2);
   FReal2 r_060_2 = PmA[1] * r_050_2 - PmQf[1] * r_050_3 + 5 * iz2 * (r_040_2 - riz * r_040_3);
   pOut[63] = PmA[1] * pOut[45] - PmQf[1] * r_032_1 + 3 * iz2 * (pOut[24] - riz * r_022_1);
   FReal2 r_042_1 = PmA[1] * r_032_1 - PmQf[1] * r_032_2 + 3 * iz2 * (r_022_1 - riz * r_022_2);
   pOut[64] = PmA[2] * pOut[51] - PmQf[2] * r_023_1 + 3 * iz2 * (pOut[24] - riz * r_022_1);
   FReal2 r_024_1 = PmA[2] * r_023_1 - PmQf[2] * r_023_2 + 3 * iz2 * (r_022_1 - riz * r_022_2);
   pOut[65] = PmA[2] * pOut[52] - PmQf[2] * r_005_1 + 5 * iz2 * (pOut[25] - riz * r_004_1);
   FReal2 r_006_1 = PmA[2] * r_005_1 - PmQf[2] * r_005_2 + 5 * iz2 * (r_004_1 - riz * r_004_2);
   FReal2 r_006_2 = PmA[2] * r_005_2 - PmQf[2] * r_005_3 + 5 * iz2 * (r_004_2 - riz * r_004_3);
   pOut[66] = PmA[1] * pOut[35] - PmQf[1] * r_500_1;
   FReal2 r_510_1 = PmA[1] * r_500_1 - PmQf[1] * r_500_2;
   FReal2 r_510_2 = PmA[1] * r_500_2 - PmQf[1] * r_500_3;
   pOut[67] = PmA[0] * pOut[42] - PmQf[0] * r_230_1 + 2 * iz2 * (pOut[27] - riz * r_130_1);
   FReal2 r_330_1 = PmA[0] * r_230_1 - PmQf[0] * r_230_2 + 2 * iz2 * (r_130_1 - riz * r_130_2);
   FReal2 r_330_2 = PmA[0] * r_230_2 - PmQf[0] * r_230_3 + 2 * iz2 * (r_130_2 - riz * r_130_3);
   pOut[68] = PmA[1] * pOut[37] - PmQf[1] * r_302_1;
   FReal2 r_312_1 = PmA[1] * r_302_1 - PmQf[1] * r_302_2;
   FReal2 r_312_2 = PmA[1] * r_302_2 - PmQf[1] * r_302_3;
   pOut[69] = PmA[1] * pOut[38] - PmQf[1] * r_140_1 + 4 * iz2 * (pOut[27] - riz * r_130_1);
   FReal2 r_150_1 = PmA[0] * r_050_1 - PmQf[0] * r_050_2;
   FReal2 r_150_2 = PmA[0] * r_050_2 - PmQf[0] * r_050_3;
   pOut[70] = PmA[0] * pOut[45] - PmQf[0] * r_032_1;
   FReal2 r_132_1 = PmA[0] * r_032_1 - PmQf[0] * r_032_2;
   pOut[71] = PmA[0] * pOut[46] - PmQf[0] * r_014_1;
   pOut[72] = PmA[0] * pOut[47] - PmQf[0] * r_401_1 + 4 * iz2 * (pOut[29] - riz * r_301_1);
   FReal2 r_501_1 = PmA[2] * r_500_1 - PmQf[2] * r_500_2;
   FReal2 r_501_2 = PmA[2] * r_500_2 - PmQf[2] * r_500_3;
   pOut[73] = PmA[2] * pOut[36] - PmQf[2] * r_320_1;
   FReal2 r_321_1 = PmA[2] * r_320_1 - PmQf[2] * r_320_2;
   pOut[74] = PmA[2] * pOut[37] - PmQf[2] * r_302_1 + 2 * iz2 * (pOut[29] - riz * r_301_1);
   FReal2 r_303_1 = PmA[2] * r_302_1 - PmQf[2] * r_302_2 + 2 * iz2 * (r_301_1 - riz * r_301_2);
   FReal2 r_303_2 = PmA[2] * r_302_2 - PmQf[2] * r_302_3 + 2 * iz2 * (r_301_2 - riz * r_301_3);
   pOut[75] = PmA[2] * pOut[38] - PmQf[2] * r_140_1;
   pOut[76] = PmA[0] * pOut[51] - PmQf[0] * r_023_1;
   FReal2 r_123_1 = PmA[0] * r_023_1 - PmQf[0] * r_023_2;
   FReal2 r_123_2 = PmA[0] * r_023_2 - PmQf[0] * r_023_3;
   pOut[77] = PmA[0] * pOut[52] - PmQf[0] * r_005_1;
   FReal2 r_105_1 = PmA[0] * r_005_1 - PmQf[0] * r_005_2;
   FReal2 r_105_2 = PmA[0] * r_005_2 - PmQf[0] * r_005_3;
   pOut[78] = PmA[1] * pOut[47] - PmQf[1] * r_401_1;
   pOut[79] = PmA[2] * pOut[42] - PmQf[2] * r_230_1;
   FReal2 r_231_1 = PmA[2] * r_230_1 - PmQf[2] * r_230_2;
   FReal2 r_231_2 = PmA[2] * r_230_2 - PmQf[2] * r_230_3;
   pOut[80] = PmA[1] * pOut[49] - PmQf[1] * r_203_1;
   pOut[81] = PmA[2] * pOut[44] - PmQf[2] * r_050_1;
   FReal2 r_051_1 = PmA[2] * r_050_1 - PmQf[2] * r_050_2;
   FReal2 r_051_2 = PmA[2] * r_050_2 - PmQf[2] * r_050_3;
   pOut[82] = PmA[1] * pOut[51] - PmQf[1] * r_023_1 + 2 * iz2 * (pOut[34] - riz * r_013_1);
   FReal2 r_033_1 = PmA[1] * r_023_1 - PmQf[1] * r_023_2 + 2 * iz2 * (r_013_1 - riz * r_013_2);
   FReal2 r_033_2 = PmA[1] * r_023_2 - PmQf[1] * r_023_3 + 2 * iz2 * (r_013_2 - riz * r_013_3);
   pOut[83] = PmA[2] * pOut[46] - PmQf[2] * r_014_1 + 4 * iz2 * (pOut[34] - riz * r_013_1);
   FReal2 r_015_1 = PmA[1] * r_005_1 - PmQf[1] * r_005_2;
   FReal2 r_015_2 = PmA[1] * r_005_2 - PmQf[1] * r_005_3;

   // L = 7
   pOut[84] = PmA[0] * pOut[56] - PmQf[0] * r_600_1 + 6 * iz2 * (pOut[35] - riz * r_500_1);
   FReal2 r_700_1 = PmA[0] * r_600_1 - PmQf[0] * r_600_2 + 6 * iz2 * (r_500_1 - riz * r_500_2);
   pOut[85] = PmA[1] * pOut[66] - PmQf[1] * r_510_1 + iz2 * (pOut[35] - riz * r_500_1);
   FReal2 r_520_1 = PmA[1] * r_510_1 - PmQf[1] * r_510_2 + iz2 * (r_500_1 - riz * r_500_2);
   pOut[86] = PmA[2] * pOut[72] - PmQf[2] * r_501_1 + iz2 * (pOut[35] - riz * r_500_1);
   FReal2 r_502_1 = PmA[2] * r_501_1 - PmQf[2] * r_501_2 + iz2 * (r_500_1 - riz * r_500_2);
   pOut[87] = PmA[1] * pOut[67] - PmQf[1] * r_330_1 + 3 * iz2 * (pOut[36] - riz * r_320_1);
   FReal2 r_340_1 = PmA[1] * r_330_1 - PmQf[1] * r_330_2 + 3 * iz2 * (r_320_1 - riz * r_320_2);
   pOut[88] = PmA[1] * pOut[68] - PmQf[1] * r_312_1 + iz2 * (pOut[37] - riz * r_302_1);
   FReal2 r_322_1 = PmA[1] * r_312_1 - PmQf[1] * r_312_2 + iz2 * (r_302_1 - riz * r_302_2);
   pOut[89] = PmA[2] * pOut[74] - PmQf[2] * r_303_1 + 3 * iz2 * (pOut[37] - riz * r_302_1);
   FReal2 r_304_1 = PmA[2] * r_303_1 - PmQf[2] * r_303_2 + 3 * iz2 * (r_302_1 - riz * r_302_2);
   pOut[90] = PmA[1] * pOut[69] - PmQf[1] * r_150_1 + 5 * iz2 * (pOut[38] - riz * r_140_1);
   FReal2 r_160_1 = PmA[0] * r_060_1 - PmQf[0] * r_060_2;
   pOut[91] = PmA[0] * pOut[63] - PmQf[0] * r_042_1;
   pOut[92] = PmA[0] * pOut[64] - PmQf[0] * r_024_1;
   pOut[93] = PmA[0] * pOut[65] - PmQf[0] * r_006_1;
   pOut[94] = PmA[1] * pOut[56] - PmQf[1] * r_600_1;
   FReal2 r_610_1 = PmA[1] * r_600_1 - PmQf[1] * r_600_2;
   pOut[95] = PmA[0] * pOut[67] - PmQf[0] * r_330_1 + 3 * iz2 * (pOut[42] - riz * r_230_1);
   FReal2 r_430_1 = PmA[0] * r_330_1 - PmQf[0] * r_330_2 + 3 * iz2 * (r_230_1 - riz * r_230_2);
   pOut[96] = PmA[0] * pOut[68] - PmQf[0] * r_312_1 + 3 * iz2 * (pOut[43] - riz * r_212_1);
   pOut[97] = PmA[0] * pOut[69] - PmQf[0] * r_150_1 + iz2 * (pOut[44] - riz * r_050_1);
   FReal2 r_250_1 = PmA[0] * r_150_1 - PmQf[0] * r_150_2 + iz2 * (r_050_1 - riz * r_050_2);
   pOut[98] = PmA[1] * pOut[60] - PmQf[1] * r_222_1 + 2 * iz2 * (pOut[43] - riz * r_212_1);
   FReal2 r_232_1 = PmA[2] * r_231_1 - PmQf[2] * r_231_2 + iz2 * (r_230_1 - riz * r_230_2);
   pOut[99] = PmA[1] * pOut[61] - PmQf[1] * r_204_1;
   pOut[100] = PmA[1] * pOut[62] - PmQf[1] * r_060_1 + 6 * iz2 * (pOut[44] - riz * r_050_1);
   FReal2 r_070_1 = PmA[1] * r_060_1 - PmQf[1] * r_060_2 + 6 * iz2 * (r_050_1 - riz * r_050_2);
   pOut[101] = PmA[1] * pOut[63] - PmQf[1] * r_042_1 + 4 * iz2 * (pOut[45] - riz * r_032_1);
   FReal2 r_052_1 = PmA[2] * r_051_1 - PmQf[2] * r_051_2 + iz2 * (r_050_1 - riz * r_050_2);
   pOut[102] = PmA[1] * pOut[64] - PmQf[1] * r_024_1 + 2 * iz2 * (pOut[46] - riz * r_014_1);
   FReal2 r_034_1 = PmA[2] * r_033_1 - PmQf[2] * r_033_2 + 3 * iz2 * (r_032_1 - riz * r_032_2);
   pOut[103] = PmA[2] * pOut[83] - PmQf[2] * r_015_1 + 5 * iz2 * (pOut[46] - riz * r_014_1);
   FReal2 r_016_1 = PmA[1] * r_006_1 - PmQf[1] * r_006_2;
   pOut[104] = PmA[0] * pOut[72] - PmQf[0] * r_501_1 + 5 * iz2 * (pOut[47] - riz * r_401_1);
   FReal2 r_601_1 = PmA[2] * r_600_1 - PmQf[2] * r_600_2;
   pOut[105] = PmA[0] * pOut[73] - PmQf[0] * r_321_1 + 3 * iz2 * (pOut[48] - riz * r_221_1);
   pOut[106] = PmA[0] * pOut[74] - PmQf[0] * r_303_1 + 3 * iz2 * (pOut[49] - riz * r_203_1);
   FReal2 r_403_1 = PmA[0] * r_303_1 - PmQf[0] * r_303_2 + 3 * iz2 * (r_203_1 - riz * r_203_2);
   pOut[107] = PmA[1] * pOut[79] - PmQf[1] * r_231_1 + 3 * iz2 * (pOut[48] - riz * r_221_1);
   pOut[108] = PmA[0] * pOut[76] - PmQf[0] * r_123_1 + iz2 * (pOut[51] - riz * r_023_1);
   FReal2 r_223_1 = PmA[0] * r_123_1 - PmQf[0] * r_123_2 + iz2 * (r_023_1 - riz * r_023_2);
   pOut[109] = PmA[2] * pOut[61] - PmQf[2] * r_204_1 + 4 * iz2 * (pOut[49] - riz * r_203_1);
   FReal2 r_205_1 = PmA[0] * r_105_1 - PmQf[0] * r_105_2 + iz2 * (r_005_1 - riz * r_005_2);
   pOut[110] = PmA[2] * pOut[62] - PmQf[2] * r_060_1;
   pOut[111] = PmA[1] * pOut[82] - PmQf[1] * r_033_1 + 3 * iz2 * (pOut[51] - riz * r_023_1);
   pOut[112] = PmA[1] * pOut[83] - PmQf[1] * r_015_1 + iz2 * (pOut[52] - riz * r_005_1);
   FReal2 r_025_1 = PmA[1] * r_015_1 - PmQf[1] * r_015_2 + iz2 * (r_005_1 - riz * r_005_2);
   pOut[113] = PmA[2] * pOut[65] - PmQf[2] * r_006_1 + 6 * iz2 * (pOut[52] - riz * r_005_1);
   FReal2 r_007_1 = PmA[2] * r_006_1 - PmQf[2] * r_006_2 + 6 * iz2 * (r_005_1 - riz * r_005_2);
   pOut[114] = PmA[1] * pOut[72] - PmQf[1] * r_501_1;
   pOut[115] = PmA[2] * pOut[67] - PmQf[2] * r_330_1;
   pOut[116] = PmA[1] * pOut[74] - PmQf[1] * r_303_1;
   pOut[117] = PmA[2] * pOut[69] - PmQf[2] * r_150_1;
   pOut[118] = PmA[0] * pOut[82] - PmQf[0] * r_033_1;
   FReal2 r_133_1 = PmA[0] * r_033_1 - PmQf[0] * r_033_2;
   pOut[119] = PmA[0] * pOut[83] - PmQf[0] * r_015_1;

   // L = 8
   pOut[120] = PmA[0] * pOut[84] - PmQf[0] * r_700_1 + 7 * iz2 * (pOut[56] - riz * r_600_1);
   pOut[121] = PmA[1] * pOut[94] - PmQf[1] * r_610_1 + iz2 * (pOut[56] - riz * r_600_1);
   pOut[122] = PmA[2] * pOut[104] - PmQf[2] * r_601_1 + iz2 * (pOut[56] - riz * r_600_1);
   pOut[123] = PmA[0] * pOut[87] - PmQf[0] * r_340_1 + 3 * iz2 * (pOut[59] - riz * r_240_1);
   pOut[124] = PmA[0] * pOut[88] - PmQf[0] * r_322_1 + 3 * iz2 * (pOut[60] - riz * r_222_1);
   pOut[125] = PmA[2] * pOut[106] - PmQf[2] * r_403_1 + 3 * iz2 * (pOut[58] - riz * r_402_1);
   pOut[126] = PmA[0] * pOut[90] - PmQf[0] * r_160_1 + iz2 * (pOut[62] - riz * r_060_1);
   pOut[127] = PmA[1] * pOut[98] - PmQf[1] * r_232_1 + 3 * iz2 * (pOut[60] - riz * r_222_1);
   pOut[128] = PmA[2] * pOut[108] - PmQf[2] * r_223_1 + 3 * iz2 * (pOut[60] - riz * r_222_1);
   pOut[129] = PmA[2] * pOut[109] - PmQf[2] * r_205_1 + 5 * iz2 * (pOut[61] - riz * r_204_1);
   pOut[130] = PmA[1] * pOut[100] - PmQf[1] * r_070_1 + 7 * iz2 * (pOut[62] - riz * r_060_1);
   pOut[131] = PmA[1] * pOut[101] - PmQf[1] * r_052_1 + 5 * iz2 * (pOut[63] - riz * r_042_1);
   pOut[132] = PmA[1] * pOut[102] - PmQf[1] * r_034_1 + 3 * iz2 * (pOut[64] - riz * r_024_1);
   pOut[133] = PmA[1] * pOut[103] - PmQf[1] * r_016_1 + iz2 * (pOut[65] - riz * r_006_1);
   pOut[134] = PmA[2] * pOut[113] - PmQf[2] * r_007_1 + 7 * iz2 * (pOut[65] - riz * r_006_1);
   pOut[135] = PmA[0] * pOut[94] - PmQf[0] * r_610_1 + 6 * iz2 * (pOut[66] - riz * r_510_1);
   pOut[136] = PmA[1] * pOut[85] - PmQf[1] * r_520_1 + 2 * iz2 * (pOut[66] - riz * r_510_1);
   pOut[137] = PmA[1] * pOut[86] - PmQf[1] * r_502_1;
   pOut[138] = PmA[0] * pOut[97] - PmQf[0] * r_250_1 + 2 * iz2 * (pOut[69] - riz * r_150_1);
   pOut[139] = PmA[0] * pOut[98] - PmQf[0] * r_232_1 + 2 * iz2 * (pOut[70] - riz * r_132_1);
   pOut[140] = PmA[1] * pOut[89] - PmQf[1] * r_304_1;
   pOut[141] = PmA[0] * pOut[100] - PmQf[0] * r_070_1;
   pOut[142] = PmA[0] * pOut[101] - PmQf[0] * r_052_1;
   pOut[143] = PmA[2] * pOut[118] - PmQf[2] * r_133_1 + 3 * iz2 * (pOut[70] - riz * r_132_1);
   pOut[144] = PmA[0] * pOut[103] - PmQf[0] * r_016_1;
   pOut[145] = PmA[0] * pOut[104] - PmQf[0] * r_601_1 + 6 * iz2 * (pOut[72] - riz * r_501_1);
   pOut[146] = PmA[2] * pOut[85] - PmQf[2] * r_520_1;
   pOut[147] = PmA[2] * pOut[86] - PmQf[2] * r_502_1 + 2 * iz2 * (pOut[72] - riz * r_501_1);
   pOut[148] = PmA[2] * pOut[87] - PmQf[2] * r_340_1;
   pOut[149] = PmA[2] * pOut[88] - PmQf[2] * r_322_1 + 2 * iz2 * (pOut[73] - riz * r_321_1);
   pOut[150] = PmA[0] * pOut[109] - PmQf[0] * r_205_1 + 2 * iz2 * (pOut[77] - riz * r_105_1);
   pOut[151] = PmA[2] * pOut[90] - PmQf[2] * r_160_1;
   pOut[152] = PmA[1] * pOut[118] - PmQf[1] * r_133_1 + 3 * iz2 * (pOut[76] - riz * r_123_1);
   pOut[153] = PmA[0] * pOut[112] - PmQf[0] * r_025_1;
   pOut[154] = PmA[0] * pOut[113] - PmQf[0] * r_007_1;
   pOut[155] = PmA[1] * pOut[104] - PmQf[1] * r_601_1;
   pOut[156] = PmA[2] * pOut[95] - PmQf[2] * r_430_1;
   pOut[157] = PmA[1] * pOut[106] - PmQf[1] * r_403_1;
   pOut[158] = PmA[2] * pOut[97] - PmQf[2] * r_250_1;
   pOut[159] = PmA[0] * pOut[118] - PmQf[0] * r_133_1 + iz2 * (pOut[82] - riz * r_033_1);
   pOut[160] = PmA[1] * pOut[109] - PmQf[1] * r_205_1;
   pOut[161] = PmA[2] * pOut[100] - PmQf[2] * r_070_1;
   pOut[162] = PmA[2] * pOut[101] - PmQf[2] * r_052_1 + 2 * iz2 * (pOut[81] - riz * r_051_1);
   pOut[163] = PmA[1] * pOut[112] - PmQf[1] * r_025_1 + 2 * iz2 * (pOut[83] - riz * r_015_1);
   pOut[164] = PmA[2] * pOut[103] - PmQf[2] * r_016_1 + 6 * iz2 * (pOut[83] - riz * r_015_1);
   // 1829 flops, 1128 mops, 1.48kb stack
}











static void OsrrC44( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const BmA[3], unsigned sa, unsigned sb ) AIC_NO_THROW
{
   const double c0 = 3.;
   const double c1 = 3.75e-01;
   const double c2 = 7.5e-01;
   const double c3 = 1.1180339887498947;
   const double c4 = 6.7082039324993685;
   const double c5 = 2.3717082451262841;
   const double c6 = 3.1622776601683791;
   const double c7 = 7.3950997288745191e-01;
   const double c8 = 4.4370598373247114;
   const double c9 = 3.3541019662496843;
   const double ca = 5.5901699437494734e-01;
   const double cb = 2.9580398915498076;
   const double cc = 6.2749501990055663;
   const double cd = 2.0916500663351885;

   FReal2 a_000[9]; ShTr4(a_000, pIn[20], pIn[21], pIn[22], pIn[23], pIn[24], pIn[25], pIn[26], pIn[27], pIn[28], pIn[29], pIn[30], pIn[31], pIn[32], pIn[33], pIn[34]);
   FReal2 a_100[9]; ShTr4(a_100, pIn[35], pIn[36], pIn[37], pIn[38], pIn[39], pIn[40], pIn[41], pIn[42], pIn[43], pIn[47], pIn[48], pIn[49], pIn[53], pIn[54], pIn[55]);
   FReal2 a_010[9]; ShTr4(a_010, pIn[41], pIn[42], pIn[43], pIn[44], pIn[45], pIn[46], pIn[36], pIn[38], pIn[39], pIn[53], pIn[54], pIn[55], pIn[48], pIn[50], pIn[51]);
   FReal2 a_001[9]; ShTr4(a_001, pIn[47], pIn[48], pIn[49], pIn[50], pIn[51], pIn[52], pIn[53], pIn[54], pIn[55], pIn[37], pIn[39], pIn[40], pIn[43], pIn[45], pIn[46]);
   FReal2 a_200[9]; ShTr4(a_200, pIn[56], pIn[57], pIn[58], pIn[59], pIn[60], pIn[61], pIn[66], pIn[67], pIn[68], pIn[72], pIn[73], pIn[74], pIn[78], pIn[79], pIn[80]);
   FReal2 a_020[9]; ShTr4(a_020, pIn[57], pIn[59], pIn[60], pIn[62], pIn[63], pIn[64], pIn[67], pIn[69], pIn[70], pIn[73], pIn[75], pIn[76], pIn[79], pIn[81], pIn[82]);
   FReal2 a_002[9]; ShTr4(a_002, pIn[58], pIn[60], pIn[61], pIn[63], pIn[64], pIn[65], pIn[68], pIn[70], pIn[71], pIn[74], pIn[76], pIn[77], pIn[80], pIn[82], pIn[83]);
   FReal2 a_110[9]; ShTr4(a_110, pIn[66], pIn[67], pIn[68], pIn[69], pIn[70], pIn[71], pIn[57], pIn[59], pIn[60], pIn[78], pIn[79], pIn[80], pIn[73], pIn[75], pIn[76]);
   FReal2 a_101[9]; ShTr4(a_101, pIn[72], pIn[73], pIn[74], pIn[75], pIn[76], pIn[77], pIn[78], pIn[79], pIn[80], pIn[58], pIn[60], pIn[61], pIn[68], pIn[70], pIn[71]);
   FReal2 a_011[9]; ShTr4(a_011, pIn[78], pIn[79], pIn[80], pIn[81], pIn[82], pIn[83], pIn[73], pIn[75], pIn[76], pIn[68], pIn[70], pIn[71], pIn[60], pIn[63], pIn[64]);
   FReal2 a_300[9]; ShTr4(a_300, pIn[84], pIn[85], pIn[86], pIn[87], pIn[88], pIn[89], pIn[94], pIn[95], pIn[96], pIn[104], pIn[105], pIn[106], pIn[114], pIn[115], pIn[116]);
   FReal2 a_120[9]; ShTr4(a_120, pIn[85], pIn[87], pIn[88], pIn[90], pIn[91], pIn[92], pIn[95], pIn[97], pIn[98], pIn[105], pIn[107], pIn[108], pIn[115], pIn[117], pIn[118]);
   FReal2 a_102[9]; ShTr4(a_102, pIn[86], pIn[88], pIn[89], pIn[91], pIn[92], pIn[93], pIn[96], pIn[98], pIn[99], pIn[106], pIn[108], pIn[109], pIn[116], pIn[118], pIn[119]);
   FReal2 a_210[9]; ShTr4(a_210, pIn[94], pIn[95], pIn[96], pIn[97], pIn[98], pIn[99], pIn[85], pIn[87], pIn[88], pIn[114], pIn[115], pIn[116], pIn[105], pIn[107], pIn[108]);
   FReal2 a_030[9]; ShTr4(a_030, pIn[95], pIn[97], pIn[98], pIn[100], pIn[101], pIn[102], pIn[87], pIn[90], pIn[91], pIn[115], pIn[117], pIn[118], pIn[107], pIn[110], pIn[111]);
   FReal2 a_012[9]; ShTr4(a_012, pIn[96], pIn[98], pIn[99], pIn[101], pIn[102], pIn[103], pIn[88], pIn[91], pIn[92], pIn[116], pIn[118], pIn[119], pIn[108], pIn[111], pIn[112]);
   FReal2 a_201[9]; ShTr4(a_201, pIn[104], pIn[105], pIn[106], pIn[107], pIn[108], pIn[109], pIn[114], pIn[115], pIn[116], pIn[86], pIn[88], pIn[89], pIn[96], pIn[98], pIn[99]);
   FReal2 a_021[9]; ShTr4(a_021, pIn[105], pIn[107], pIn[108], pIn[110], pIn[111], pIn[112], pIn[115], pIn[117], pIn[118], pIn[88], pIn[91], pIn[92], pIn[98], pIn[101], pIn[102]);
   FReal2 a_003[9]; ShTr4(a_003, pIn[106], pIn[108], pIn[109], pIn[111], pIn[112], pIn[113], pIn[116], pIn[118], pIn[119], pIn[89], pIn[92], pIn[93], pIn[99], pIn[102], pIn[103]);
   FReal2 a_111[9]; ShTr4(a_111, pIn[114], pIn[115], pIn[116], pIn[117], pIn[118], pIn[119], pIn[105], pIn[107], pIn[108], pIn[96], pIn[98], pIn[99], pIn[88], pIn[91], pIn[92]);
   FReal2 a_400[9]; ShTr4(a_400, pIn[120], pIn[121], pIn[122], pIn[123], pIn[124], pIn[125], pIn[135], pIn[136], pIn[137], pIn[145], pIn[146], pIn[147], pIn[155], pIn[156], pIn[157]);
   FReal2 a_220[9]; ShTr4(a_220, pIn[121], pIn[123], pIn[124], pIn[126], pIn[127], pIn[128], pIn[136], pIn[138], pIn[139], pIn[146], pIn[148], pIn[149], pIn[156], pIn[158], pIn[159]);
   FReal2 a_202[9]; ShTr4(a_202, pIn[122], pIn[124], pIn[125], pIn[127], pIn[128], pIn[129], pIn[137], pIn[139], pIn[140], pIn[147], pIn[149], pIn[150], pIn[157], pIn[159], pIn[160]);
   FReal2 a_040[9]; ShTr4(a_040, pIn[123], pIn[126], pIn[127], pIn[130], pIn[131], pIn[132], pIn[138], pIn[141], pIn[142], pIn[148], pIn[151], pIn[152], pIn[158], pIn[161], pIn[162]);
   FReal2 a_022[9]; ShTr4(a_022, pIn[124], pIn[127], pIn[128], pIn[131], pIn[132], pIn[133], pIn[139], pIn[142], pIn[143], pIn[149], pIn[152], pIn[153], pIn[159], pIn[162], pIn[163]);
   FReal2 a_004[9]; ShTr4(a_004, pIn[125], pIn[128], pIn[129], pIn[132], pIn[133], pIn[134], pIn[140], pIn[143], pIn[144], pIn[150], pIn[153], pIn[154], pIn[160], pIn[163], pIn[164]);
   FReal2 a_310[9]; ShTr4(a_310, pIn[135], pIn[136], pIn[137], pIn[138], pIn[139], pIn[140], pIn[121], pIn[123], pIn[124], pIn[155], pIn[156], pIn[157], pIn[146], pIn[148], pIn[149]);
   FReal2 a_130[9]; ShTr4(a_130, pIn[136], pIn[138], pIn[139], pIn[141], pIn[142], pIn[143], pIn[123], pIn[126], pIn[127], pIn[156], pIn[158], pIn[159], pIn[148], pIn[151], pIn[152]);
   FReal2 a_112[9]; ShTr4(a_112, pIn[137], pIn[139], pIn[140], pIn[142], pIn[143], pIn[144], pIn[124], pIn[127], pIn[128], pIn[157], pIn[159], pIn[160], pIn[149], pIn[152], pIn[153]);
   FReal2 a_301[9]; ShTr4(a_301, pIn[145], pIn[146], pIn[147], pIn[148], pIn[149], pIn[150], pIn[155], pIn[156], pIn[157], pIn[122], pIn[124], pIn[125], pIn[137], pIn[139], pIn[140]);
   FReal2 a_121[9]; ShTr4(a_121, pIn[146], pIn[148], pIn[149], pIn[151], pIn[152], pIn[153], pIn[156], pIn[158], pIn[159], pIn[124], pIn[127], pIn[128], pIn[139], pIn[142], pIn[143]);
   FReal2 a_103[9]; ShTr4(a_103, pIn[147], pIn[149], pIn[150], pIn[152], pIn[153], pIn[154], pIn[157], pIn[159], pIn[160], pIn[125], pIn[128], pIn[129], pIn[140], pIn[143], pIn[144]);
   FReal2 a_211[9]; ShTr4(a_211, pIn[155], pIn[156], pIn[157], pIn[158], pIn[159], pIn[160], pIn[146], pIn[148], pIn[149], pIn[137], pIn[139], pIn[140], pIn[124], pIn[127], pIn[128]);
   FReal2 a_031[9]; ShTr4(a_031, pIn[156], pIn[158], pIn[159], pIn[161], pIn[162], pIn[163], pIn[148], pIn[151], pIn[152], pIn[139], pIn[142], pIn[143], pIn[127], pIn[131], pIn[132]);
   FReal2 a_013[9]; ShTr4(a_013, pIn[157], pIn[159], pIn[160], pIn[162], pIn[163], pIn[164], pIn[149], pIn[152], pIn[153], pIn[140], pIn[143], pIn[144], pIn[128], pIn[132], pIn[133]);
   // ShTr: [1645 flops, 2800 mops]

   for ( unsigned ica = 0; ica < 9; ++ ica ) {
      FReal2 r_000_100 = a_100[ica] - BmA[0] * a_000[ica];
      FReal2 r_100_100 = a_200[ica] - BmA[0] * a_100[ica];
      FReal2 r_010_100 = a_110[ica] - BmA[0] * a_010[ica];
      FReal2 r_001_100 = a_101[ica] - BmA[0] * a_001[ica];
      FReal2 r_200_100 = a_300[ica] - BmA[0] * a_200[ica];
      FReal2 r_002_100 = a_102[ica] - BmA[0] * a_002[ica];
      FReal2 r_110_100 = a_210[ica] - BmA[0] * a_110[ica];
      FReal2 r_101_100 = a_201[ica] - BmA[0] * a_101[ica];
      FReal2 r_011_100 = a_111[ica] - BmA[0] * a_011[ica];
      FReal2 r_300_100 = a_400[ica] - BmA[0] * a_300[ica];
      FReal2 r_102_100 = a_202[ica] - BmA[0] * a_102[ica];
      FReal2 r_210_100 = a_310[ica] - BmA[0] * a_210[ica];
      FReal2 r_201_100 = a_301[ica] - BmA[0] * a_201[ica];
      FReal2 r_111_100 = a_211[ica] - BmA[0] * a_111[ica];
      FReal2 r_000_010 = a_010[ica] - BmA[1] * a_000[ica];
      FReal2 r_100_010 = a_110[ica] - BmA[1] * a_100[ica];
      FReal2 r_010_010 = a_020[ica] - BmA[1] * a_010[ica];
      FReal2 r_001_010 = a_011[ica] - BmA[1] * a_001[ica];
      FReal2 r_200_010 = a_210[ica] - BmA[1] * a_200[ica];
      FReal2 r_020_010 = a_030[ica] - BmA[1] * a_020[ica];
      FReal2 r_110_010 = a_120[ica] - BmA[1] * a_110[ica];
      FReal2 r_101_010 = a_111[ica] - BmA[1] * a_101[ica];
      FReal2 r_011_010 = a_021[ica] - BmA[1] * a_011[ica];
      FReal2 r_120_010 = a_130[ica] - BmA[1] * a_120[ica];
      FReal2 r_210_010 = a_220[ica] - BmA[1] * a_210[ica];
      FReal2 r_030_010 = a_040[ica] - BmA[1] * a_030[ica];
      FReal2 r_021_010 = a_031[ica] - BmA[1] * a_021[ica];
      FReal2 r_111_010 = a_121[ica] - BmA[1] * a_111[ica];
      FReal2 r_000_001 = a_001[ica] - BmA[2] * a_000[ica];
      FReal2 r_100_001 = a_101[ica] - BmA[2] * a_100[ica];
      FReal2 r_010_001 = a_011[ica] - BmA[2] * a_010[ica];
      FReal2 r_001_001 = a_002[ica] - BmA[2] * a_001[ica];
      FReal2 r_020_001 = a_021[ica] - BmA[2] * a_020[ica];
      FReal2 r_002_001 = a_003[ica] - BmA[2] * a_002[ica];
      FReal2 r_110_001 = a_111[ica] - BmA[2] * a_110[ica];
      FReal2 r_101_001 = a_102[ica] - BmA[2] * a_101[ica];
      FReal2 r_011_001 = a_012[ica] - BmA[2] * a_011[ica];
      FReal2 r_102_001 = a_103[ica] - BmA[2] * a_102[ica];
      FReal2 r_012_001 = a_013[ica] - BmA[2] * a_012[ica];
      FReal2 r_021_001 = a_022[ica] - BmA[2] * a_021[ica];
      FReal2 r_003_001 = a_004[ica] - BmA[2] * a_003[ica];
      FReal2 r_111_001 = a_112[ica] - BmA[2] * a_111[ica];
      FReal2 r_000_200 = r_100_100 - BmA[0] * r_000_100;
      FReal2 r_100_200 = r_200_100 - BmA[0] * r_100_100;
      FReal2 r_010_200 = r_110_100 - BmA[0] * r_010_100;
      FReal2 r_001_200 = r_101_100 - BmA[0] * r_001_100;
      FReal2 r_200_200 = r_300_100 - BmA[0] * r_200_100;
      FReal2 r_002_200 = r_102_100 - BmA[0] * r_002_100;
      FReal2 r_110_200 = r_210_100 - BmA[0] * r_110_100;
      FReal2 r_101_200 = r_201_100 - BmA[0] * r_101_100;
      FReal2 r_011_200 = r_111_100 - BmA[0] * r_011_100;
      FReal2 r_000_020 = r_010_010 - BmA[1] * r_000_010;
      FReal2 r_100_020 = r_110_010 - BmA[1] * r_100_010;
      FReal2 r_010_020 = r_020_010 - BmA[1] * r_010_010;
      FReal2 r_001_020 = r_011_010 - BmA[1] * r_001_010;
      FReal2 r_200_020 = r_210_010 - BmA[1] * r_200_010;
      FReal2 r_020_020 = r_030_010 - BmA[1] * r_020_010;
      FReal2 r_110_020 = r_120_010 - BmA[1] * r_110_010;
      FReal2 r_101_020 = r_111_010 - BmA[1] * r_101_010;
      FReal2 r_011_020 = r_021_010 - BmA[1] * r_011_010;
      FReal2 r_000_002 = r_001_001 - BmA[2] * r_000_001;
      FReal2 r_100_002 = r_101_001 - BmA[2] * r_100_001;
      FReal2 r_010_002 = r_011_001 - BmA[2] * r_010_001;
      FReal2 r_001_002 = r_002_001 - BmA[2] * r_001_001;
      FReal2 r_020_002 = r_021_001 - BmA[2] * r_020_001;
      FReal2 r_002_002 = r_003_001 - BmA[2] * r_002_001;
      FReal2 r_110_002 = r_111_001 - BmA[2] * r_110_001;
      FReal2 r_101_002 = r_102_001 - BmA[2] * r_101_001;
      FReal2 r_011_002 = r_012_001 - BmA[2] * r_011_001;
      FReal2 r_000_300 = r_100_200 - BmA[0] * r_000_200;
      FReal2 r_100_300 = r_200_200 - BmA[0] * r_100_200;
      FReal2 r_010_300 = r_110_200 - BmA[0] * r_010_200;
      FReal2 r_001_300 = r_101_200 - BmA[0] * r_001_200;
      FReal2 r_000_120 = r_100_020 - BmA[0] * r_000_020;
      FReal2 r_100_120 = r_200_020 - BmA[0] * r_100_020;
      FReal2 r_001_120 = r_101_020 - BmA[0] * r_001_020;
      FReal2 r_000_030 = r_010_020 - BmA[1] * r_000_020;
      FReal2 r_100_030 = r_110_020 - BmA[1] * r_100_020;
      FReal2 r_010_030 = r_020_020 - BmA[1] * r_010_020;
      FReal2 r_001_030 = r_011_020 - BmA[1] * r_001_020;
      FReal2 r_000_012 = r_010_002 - BmA[1] * r_000_002;
      FReal2 r_100_012 = r_110_002 - BmA[1] * r_100_002;
      FReal2 r_010_012 = r_020_002 - BmA[1] * r_010_002;
      FReal2 r_000_201 = r_001_200 - BmA[2] * r_000_200;
      FReal2 r_010_201 = r_011_200 - BmA[2] * r_010_200;
      FReal2 r_001_201 = r_002_200 - BmA[2] * r_001_200;
      FReal2 r_000_003 = r_001_002 - BmA[2] * r_000_002;
      FReal2 r_100_003 = r_101_002 - BmA[2] * r_100_002;
      FReal2 r_010_003 = r_011_002 - BmA[2] * r_010_002;
      FReal2 r_001_003 = r_002_002 - BmA[2] * r_001_002;
      FReal2 b_400 = r_100_300 - BmA[0] * r_000_300;
      FReal2 b_220 = r_100_120 - BmA[0] * r_000_120;
      FReal2 b_202 = r_001_201 - BmA[2] * r_000_201;
      FReal2 b_040 = r_010_030 - BmA[1] * r_000_030;
      FReal2 b_022 = r_010_012 - BmA[1] * r_000_012;
      FReal2 b_004 = r_001_003 - BmA[2] * r_000_003;
      FReal2 b_310 = r_010_300 - BmA[1] * r_000_300;
      FReal2 b_130 = r_100_030 - BmA[0] * r_000_030;
      FReal2 b_112 = r_100_012 - BmA[0] * r_000_012;
      FReal2 b_301 = r_001_300 - BmA[2] * r_000_300;
      FReal2 b_121 = r_001_120 - BmA[2] * r_000_120;
      FReal2 b_103 = r_100_003 - BmA[0] * r_000_003;
      FReal2 b_211 = r_010_201 - BmA[1] * r_000_201;
      FReal2 b_031 = r_001_030 - BmA[2] * r_000_030;
      FReal2 b_013 = r_010_003 - BmA[1] * r_000_003;
      pOut[ica*sa + 0*sb] += -c0*b_022 + c1*b_400 - c0*b_202 + c2*b_220 + b_004 + c1*b_040;
      pOut[ica*sa + 1*sb] += -c3*b_310 - c3*b_130 + c4*b_112;
      pOut[ica*sa + 2*sb] += -c5*b_301 - c5*b_121 + c6*b_103;
      pOut[ica*sa + 3*sb] += c7*b_400 + c7*b_040 - c8*b_220;
      pOut[ica*sa + 4*sb] += -c5*b_031 - c5*b_211 + c6*b_013;
      pOut[ica*sa + 5*sb] += c9*b_202 - ca*b_400 + ca*b_040 - c9*b_022;
      pOut[ica*sa + 6*sb] += cb*b_310 - cb*b_130;
      pOut[ica*sa + 7*sb] += -cc*b_121 + cd*b_301;
      pOut[ica*sa + 8*sb] += cc*b_211 - cd*b_031;
      // 9 * [210 flops, 315 mops]
   }
   // total: [3535 flops, 5635 mops]
}

// In: Gm[0 .. 9] (inclusive), Out: [r]^0, ordered as AngularComps(9)
// (moment 9, 55 entries)
static void OsrrA9( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, FReal1 const PmA[3], FReal1 const PmQ[3], double rho, double InvEta ) AIC_NO_THROW
{
   double
      riz = rho * InvEta,
      iz2 = 0.5 * InvEta,
      PmQf[3] = { riz*PmQ[0], riz*PmQ[1], riz*PmQ[2] };
   pOut[0] = pIn[0];
   // L = 1
   pOut[1] = PmA[0] * pIn[0] - PmQf[0] * pIn[1];
   FReal2 r_100_1 = PmA[0] * pIn[1] - PmQf[0] * pIn[2];
   FReal2 r_100_2 = PmA[0] * pIn[2] - PmQf[0] * pIn[3];
   FReal2 r_100_3 = PmA[0] * pIn[3] - PmQf[0] * pIn[4];
   FReal2 r_100_4 = PmA[0] * pIn[4] - PmQf[0] * pIn[5];
   FReal2 r_100_5 = PmA[0] * pIn[5] - PmQf[0] * pIn[6];
   FReal2 r_100_6 = PmA[0] * pIn[6] - PmQf[0] * pIn[7];
   FReal2 r_100_7 = PmA[0] * pIn[7] - PmQf[0] * pIn[8];
   FReal2 r_100_8 = PmA[0] * pIn[8] - PmQf[0] * pIn[9];
   pOut[2] = PmA[1] * pIn[0] - PmQf[1] * pIn[1];
   FReal2 r_010_1 = PmA[1] * pIn[1] - PmQf[1] * pIn[2];
   FReal2 r_010_2 = PmA[1] * pIn[2] - PmQf[1] * pIn[3];
   FReal2 r_010_3 = PmA[1] * pIn[3] - PmQf[1] * pIn[4];
   FReal2 r_010_4 = PmA[1] * pIn[4] - PmQf[1] * pIn[5];
   FReal2 r_010_5 = PmA[1] * pIn[5] - PmQf[1] * pIn[6];
   FReal2 r_010_6 = PmA[1] * pIn[6] - PmQf[1] * pIn[7];
   FReal2 r_010_7 = PmA[1] * pIn[7] - PmQf[1] * pIn[8];
   FReal2 r_010_8 = PmA[1] * pIn[8] - PmQf[1] * pIn[9];
   pOut[3] = PmA[2] * pIn[0] - PmQf[2] * pIn[1];
   FReal2 r_001_1 = PmA[2] * pIn[1] - PmQf[2] * pIn[2];
   FReal2 r_001_2 = PmA[2] * pIn[2] - PmQf[2] * pIn[3];
   FReal2 r_001_3 = PmA[2] * pIn[3] - PmQf[2] * pIn[4];
   FReal2 r_001_4 = PmA[2] * pIn[4] - PmQf[2] * pIn[5];
   FReal2 r_001_5 = PmA[2] * pIn[5] - PmQf[2] * pIn[6];
   FReal2 r_001_6 = PmA[2] * pIn[6] - PmQf[2] * pIn[7];
   FReal2 r_001_7 = PmA[2] * pIn[7] - PmQf[2] * pIn[8];
   FReal2 r_001_8 = PmA[2] * pIn[8] - PmQf[2] * pIn[9];

   // L = 2
   pOut[4] = PmA[0] * pOut[1] - PmQf[0] * r_100_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_200_1 = PmA[0] * r_100_1 - PmQf[0] * r_100_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_200_2 = PmA[0] * r_100_2 - PmQf[0] * r_100_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_200_3 = PmA[0] * r_100_3 - PmQf[0] * r_100_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_200_4 = PmA[0] * r_100_4 - PmQf[0] * r_100_5 + iz2 * (pIn[4] - riz * pIn[5]);
   FReal2 r_200_5 = PmA[0] * r_100_5 - PmQf[0] * r_100_6 + iz2 * (pIn[5] - riz * pIn[6]);
   FReal2 r_200_6 = PmA[0] * r_100_6 - PmQf[0] * r_100_7 + iz2 * (pIn[6] - riz * pIn[7]);
   FReal2 r_200_7 = PmA[0] * r_100_7 - PmQf[0] * r_100_8 + iz2 * (pIn[7] - riz * pIn[8]);
   pOut[5] = PmA[1] * pOut[2] - PmQf[1] * r_010_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_020_1 = PmA[1] * r_010_1 - PmQf[1] * r_010_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_020_2 = PmA[1] * r_010_2 - PmQf[1] * r_010_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_020_3 = PmA[1] * r_010_3 - PmQf[1] * r_010_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_020_4 = PmA[1] * r_010_4 - PmQf[1] * r_010_5 + iz2 * (pIn[4] - riz * pIn[5]);
   FReal2 r_020_5 = PmA[1] * r_010_5 - PmQf[1] * r_010_6 + iz2 * (pIn[5] - riz * pIn[6]);
   FReal2 r_020_6 = PmA[1] * r_010_6 - PmQf[1] * r_010_7 + iz2 * (pIn[6] - riz * pIn[7]);
   FReal2 r_020_7 = PmA[1] * r_010_7 - PmQf[1] * r_010_8 + iz2 * (pIn[7] - riz * pIn[8]);
   pOut[6] = PmA[2] * pOut[3] - PmQf[2] * r_001_1 + iz2 * (pIn[0] - riz * pIn[1]);
   FReal2 r_002_1 = PmA[2] * r_001_1 - PmQf[2] * r_001_2 + iz2 * (pIn[1] - riz * pIn[2]);
   FReal2 r_002_2 = PmA[2] * r_001_2 - PmQf[2] * r_001_3 + iz2 * (pIn[2] - riz * pIn[3]);
   FReal2 r_002_3 = PmA[2] * r_001_3 - PmQf[2] * r_001_4 + iz2 * (pIn[3] - riz * pIn[4]);
   FReal2 r_002_4 = PmA[2] * r_001_4 - PmQf[2] * r_001_5 + iz2 * (pIn[4] - riz * pIn[5]);
   FReal2 r_002_5 = PmA[2] * r_001_5 - PmQf[2] * r_001_6 + iz2 * (pIn[5] - riz * pIn[6]);
   FReal2 r_002_6 = PmA[2] * r_001_6 - PmQf[2] * r_001_7 + iz2 * (pIn[6] - riz * pIn[7]);
   FReal2 r_002_7 = PmA[2] * r_001_7 - PmQf[2] * r_001_8 + iz2 * (pIn[7] - riz * pIn[8]);
   pOut[7] = PmA[0] * pOut[2] - PmQf[0] * r_010_1;
   pOut[8] = PmA[2] * pOut[1] - PmQf[2] * r_100_1;
   FReal2 r_101_1 = PmA[2] * r_100_1 - PmQf[2] * r_100_2;
   pOut[9] = PmA[1] * pOut[3] - PmQf[1] * r_001_1;

   // L = 3
   pOut[10] = PmA[0] * pOut[4] - PmQf[0] * r_200_1 + 2 * iz2 * (pOut[1] - riz * r_100_1);
   FReal2 r_300_1 = PmA[0] * r_200_1 - PmQf[0] * r_200_2 + 2 * iz2 * (r_100_1 - riz * r_100_2);
   FReal2 r_300_2 = PmA[0] * r_200_2 - PmQf[0] * r_200_3 + 2 * iz2 * (r_100_2 - riz * r_100_3);
   FReal2 r_300_3 = PmA[0] * r_200_3 - PmQf[0] * r_200_4 + 2 * iz2 * (r_100_3 - riz * r_100_4);
   FReal2 r_300_4 = PmA[0] * r_200_4 - PmQf[0] * r_200_5 + 2 * iz2 * (r_100_4 - riz * r_100_5);
   FReal2 r_300_5 = PmA[0] * r_200_5 - PmQf[0] * r_200_6 + 2 * iz2 * (r_100_5 - riz * r_100_6);
   FReal2 r_300_6 = PmA[0] * r_200_6 - PmQf[0] * r_200_7 + 2 * iz2 * (r_100_6 - riz * r_100_7);
   pOut[11] = PmA[0] * pOut[5] - PmQf[0] * r_020_1;
   FReal2 r_120_1 = PmA[0] * r_020_1 - PmQf[0] * r_020_2;
   FReal2 r_120_2 = PmA[0] * r_020_2 - PmQf[0] * r_020_3;
   FReal2 r_120_3 = PmA[0] * r_020_3 - PmQf[0] * r_020_4;
   FReal2 r_120_4 = PmA[0] * r_020_4 - PmQf[0] * r_020_5;
   pOut[12] = PmA[0] * pOut[6] - PmQf[0] * r_002_1;
   pOut[13] = PmA[1] * pOut[4] - PmQf[1] * r_200_1;
   pOut[14] = PmA[1] * pOut[5] - PmQf[1] * r_020_1 + 2 * iz2 * (pOut[2] - riz * r_010_1);
   FReal2 r_030_1 = PmA[1] * r_020_1 - PmQf[1] * r_020_2 + 2 * iz2 * (r_010_1 - riz * r_010_2);
   FReal2 r_030_2 = PmA[1] * r_020_2 - PmQf[1] * r_020_3 + 2 * iz2 * (r_010_2 - riz * r_010_3);
   FReal2 r_030_3 = PmA[1] * r_020_3 - PmQf[1] * r_020_4 + 2 * iz2 * (r_010_3 - riz * r_010_4);
   FReal2 r_030_4 = PmA[1] * r_020_4 - PmQf[1] * r_020_5 + 2 * iz2 * (r_010_4 - riz * r_010_5);
   FReal2 r_030_5 = PmA[1] * r_020_5 - PmQf[1] * r_020_6 + 2 * iz2 * (r_010_5 - riz * r_010_6);
   FReal2 r_030_6 = PmA[1] * r_020_6 - PmQf[1] * r_020_7 + 2 * iz2 * (r_010_6 - riz * r_010_7);
   pOut[15] = PmA[1] * pOut[6] - PmQf[1] * r_002_1;
   FReal2 r_012_1 = PmA[1] * r_002_1 - PmQf[1] * r_002_2;
   FReal2 r_012_2 = PmA[1] * r_002_2 - PmQf[1] * r_002_3;
   FReal2 r_012_3 = PmA[1] * r_002_3 - PmQf[1] * r_002_4;
   pOut[16] = PmA[0] * pOut[8] - PmQf[0] * r_101_1 + iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_201_1 = PmA[2] * r_200_1 - PmQf[2] * r_200_2;
   FReal2 r_201_2 = PmA[2] * r_200_2 - PmQf[2] * r_200_3;
   FReal2 r_201_3 = PmA[2] * r_200_3 - PmQf[2] * r_200_4;
   pOut[17] = PmA[2] * pOut[5] - PmQf[2] * r_020_1;
   pOut[18] = PmA[2] * pOut[6] - PmQf[2] * r_002_1 + 2 * iz2 * (pOut[3] - riz * r_001_1);
   FReal2 r_003_1 = PmA[2] * r_002_1 - PmQf[2] * r_002_2 + 2 * iz2 * (r_001_1 - riz * r_001_2);
   FReal2 r_003_2 = PmA[2] * r_002_2 - PmQf[2] * r_002_3 + 2 * iz2 * (r_001_2 - riz * r_001_3);
   FReal2 r_003_3 = PmA[2] * r_002_3 - PmQf[2] * r_002_4 + 2 * iz2 * (r_001_3 - riz * r_001_4);
   FReal2 r_003_4 = PmA[2] * r_002_4 - PmQf[2] * r_002_5 + 2 * iz2 * (r_001_4 - riz * r_001_5);
   FReal2 r_003_5 = PmA[2] * r_002_5 - PmQf[2] * r_002_6 + 2 * iz2 * (r_001_5 - riz * r_001_6);
   FReal2 r_003_6 = PmA[2] * r_002_6 - PmQf[2] * r_002_7 + 2 * iz2 * (r_001_6 - riz * r_001_7);
   pOut[19] = PmA[1] * pOut[8] - PmQf[1] * r_101_1;

   // L = 4
   pOut[20] = PmA[0] * pOut[10] - PmQf[0] * r_300_1 + 3 * iz2 * (pOut[4] - riz * r_200_1);
   FReal2 r_400_1 = PmA[0] * r_300_1 - PmQf[0] * r_300_2 + 3 * iz2 * (r_200_1 - riz * r_200_2);
   FReal2 r_400_2 = PmA[0] * r_300_2 - PmQf[0] * r_300_3 + 3 * iz2 * (r_200_2 - riz * r_200_3);
   FReal2 r_400_3 = PmA[0] * r_300_3 - PmQf[0] * r_300_4 + 3 * iz2 * (r_200_3 - riz * r_200_4);
   FReal2 r_400_4 = PmA[0] * r_300_4 - PmQf[0] * r_300_5 + 3 * iz2 * (r_200_4 - riz * r_200_5);
   FReal2 r_400_5 = PmA[0] * r_300_5 - PmQf[0] * r_300_6 + 3 * iz2 * (r_200_5 - riz * r_200_6);
   pOut[21] = PmA[0] * pOut[11] - PmQf[0] * r_120_1 + iz2 * (pOut[5] - riz * r_020_1);
   FReal2 r_220_1 = PmA[0] * r_120_1 - PmQf[0] * r_120_2 + iz2 * (r_020_1 - riz * r_020_2);
   FReal2 r_220_2 = PmA[0] * r_120_2 - PmQf[0] * r_120_3 + iz2 * (r_020_2 - riz * r_020_3);
   FReal2 r_220_3 = PmA[0] * r_120_3 - PmQf[0] * r_120_4 + iz2 * (r_020_3 - riz * r_020_4);
   pOut[22] = PmA[2] * pOut[16] - PmQf[2] * r_201_1 + iz2 * (pOut[4] - riz * r_200_1);
   FReal2 r_202_1 = PmA[2] * r_201_1 - PmQf[2] * r_201_2 + iz2 * (r_200_1 - riz * r_200_2);
   FReal2 r_202_2 = PmA[2] * r_201_2 - PmQf[2] * r_201_3 + iz2 * (r_200_2 - riz * r_200_3);
   pOut[23] = PmA[1] * pOut[14] - PmQf[1] * r_030_1 + 3 * iz2 * (pOut[5] - riz * r_020_1);
   FReal2 r_040_1 = PmA[1] * r_030_1 - PmQf[1] * r_030_2 + 3 * iz2 * (r_020_1 - riz * r_020_2);
   FReal2 r_040_2 = PmA[1] * r_030_2 - PmQf[1] * r_030_3 + 3 * iz2 * (r_020_2 - riz * r_020_3);
   FReal2 r_040_3 = PmA[1] * r_030_3 - PmQf[1] * r_030_4 + 3 * iz2 * (r_020_3 - riz * r_020_4);
   FReal2 r_040_4 = PmA[1] * r_030_4 - PmQf[1] * r_030_5 + 3 * iz2 * (r_020_4 - riz * r_020_5);
   FReal2 r_040_5 = PmA[1] * r_030_5 - PmQf[1] * r_030_6 + 3 * iz2 * (r_020_5 - riz * r_020_6);
   pOut[24] = PmA[1] * pOut[15] - PmQf[1] * r_012_1 + iz2 * (pOut[6] - riz * r_002_1);
   FReal2 r_022_1 = PmA[1] * r_012_1 - PmQf[1] * r_012_2 + iz2 * (r_002_1 - riz * r_002_2);
   FReal2 r_022_2 = PmA[1] * r_012_2 - PmQf[1] * r_012_3 + iz2 * (r_002_2 - riz * r_002_3);
   pOut[25] = PmA[2] * pOut[18] - PmQf[2] * r_003_1 + 3 * iz2 * (pOut[6] - riz * r_002_1);
   FReal2 r_004_1 = PmA[2] * r_003_1 - PmQf[2] * r_003_2 + 3 * iz2 * (r_002_1 - riz * r_002_2);
   FReal2 r_004_2 = PmA[2] * r_003_2 - PmQf[2] * r_003_3 + 3 * iz2 * (r_002_2 - riz * r_002_3);
   FReal2 r_004_3 = PmA[2] * r_003_3 - PmQf[2] * r_003_4 + 3 * iz2 * (r_002_3 - riz * r_002_4);
   FReal2 r_004_4 = PmA[2] * r_003_4 - PmQf[2] * r_003_5 + 3 * iz2 * (r_002_4 - riz * r_002_5);
   FReal2 r_004_5 = PmA[2] * r_003_5 - PmQf[2] * r_003_6 + 3 * iz2 * (r_002_5 - riz * r_002_6);
   pOut[26] = PmA[1] * pOut[10] - PmQf[1] * r_300_1;
   pOut[27] = PmA[0] * pOut[14] - PmQf[0] * r_030_1;
   FReal2 r_130_1 = PmA[0] * r_030_1 - PmQf[0] * r_030_2;
   FReal2 r_130_2 = PmA[0] * r_030_2 - PmQf[0] * r_030_3;
   FReal2 r_130_3 = PmA[0] * r_030_3 - PmQf[0] * r_030_4;
   FReal2 r_130_4 = PmA[0] * r_030_4 - PmQf[0] * r_030_5;
   FReal2 r_130_5 = PmA[0] * r_030_5 - PmQf[0] * r_030_6;
   pOut[28] = PmA[0] * pOut[15] - PmQf[0] * r_012_1;
   pOut[29] = PmA[0] * pOut[16] - PmQf[0] * r_201_1 + 2 * iz2 * (pOut[8] - riz * r_101_1);
   FReal2 r_301_1 = PmA[2] * r_300_1 - PmQf[2] * r_300_2;
   FReal2 r_301_2 = PmA[2] * r_300_2 - PmQf[2] * r_300_3;
   FReal2 r_301_3 = PmA[2] * r_300_3 - PmQf[2] * r_300_4;
   FReal2 r_301_4 = PmA[2] * r_300_4 - PmQf[2] * r_300_5;
   FReal2 r_301_5 = PmA[2] * r_300_5 - PmQf[2] * r_300_6;
   pOut[30] = PmA[2] * pOut[11] - PmQf[2] * r_120_1;
   pOut[31] = PmA[0] * pOut[18] - PmQf[0] * r_003_1;
   pOut[32] = PmA[1] * pOut[16] - PmQf[1] * r_201_1;
   pOut[33] = PmA[2] * pOut[14] - PmQf[2] * r_030_1;
   FReal2 r_031_2 = PmA[2] * r_030_2 - PmQf[2] * r_030_3;
   FReal2 r_031_3 = PmA[2] * r_030_3 - PmQf[2] * r_030_4;
   FReal2 r_031_4 = PmA[2] * r_030_4 - PmQf[2] * r_030_5;
   pOut[34] = PmA[1] * pOut[18] - PmQf[1] * r_003_1;
   FReal2 r_013_1 = PmA[1] * r_003_1 - PmQf[1] * r_003_2;
   FReal2 r_013_2 = PmA[1] * r_003_2 - PmQf[1] * r_003_3;
   FReal2 r_013_3 = PmA[1] * r_003_3 - PmQf[1] * r_003_4;
   FReal2 r_013_4 = PmA[1] * r_003_4 - PmQf[1] * r_003_5;
   FReal2 r_013_5 = PmA[1] * r_003_5 - PmQf[1] * r_003_6;

   // L = 5
   pOut[35] = PmA[0] * pOut[20] - PmQf[0] * r_400_1 + 4 * iz2 * (pOut[10] - riz * r_300_1);
   FReal2 r_500_1 = PmA[0] * r_400_1 - PmQf[0] * r_400_2 + 4 * iz2 * (r_300_1 - riz * r_300_2);
   FReal2 r_500_2 = PmA[0] * r_400_2 - PmQf[0] * r_400_3 + 4 * iz2 * (r_300_2 - riz * r_300_3);
   FReal2 r_500_3 = PmA[0] * r_400_3 - PmQf[0] * r_400_4 + 4 * iz2 * (r_300_3 - riz * r_300_4);
   FReal2 r_500_4 = PmA[0] * r_400_4 - PmQf[0] * r_400_5 + 4 * iz2 * (r_300_4 - riz * r_300_5);
   pOut[36] = PmA[0] * pOut[21] - PmQf[0] * r_220_1 + 2 * iz2 * (pOut[11] - riz * r_120_1);
   FReal2 r_320_1 = PmA[0] * r_220_1 - PmQf[0] * r_220_2 + 2 * iz2 * (r_120_1 - riz * r_120_2);
   pOut[37] = PmA[2] * pOut[29] - PmQf[2] * r_301_1 + iz2 * (pOut[10] - riz * r_300_1);
   FReal2 r_302_1 = PmA[2] * r_301_1 - PmQf[2] * r_301_2 + iz2 * (r_300_1 - riz * r_300_2);
   FReal2 r_302_2 = PmA[2] * r_301_2 - PmQf[2] * r_301_3 + iz2 * (r_300_2 - riz * r_300_3);
   FReal2 r_302_3 = PmA[2] * r_301_3 - PmQf[2] * r_301_4 + iz2 * (r_300_3 - riz * r_300_4);
   FReal2 r_302_4 = PmA[2] * r_301_4 - PmQf[2] * r_301_5 + iz2 * (r_300_4 - riz * r_300_5);
   pOut[38] = PmA[1] * pOut[27] - PmQf[1] * r_130_1 + 3 * iz2 * (pOut[11] - riz * r_120_1);
   FReal2 r_140_1 = PmA[1] * r_130_1 - PmQf[1] * r_130_2 + 3 * iz2 * (r_120_1 - riz * r_120_2);
   FReal2 r_140_2 = PmA[1] * r_130_2 - PmQf[1] * r_130_3 + 3 * iz2 * (r_120_2 - riz * r_120_3);
   FReal2 r_140_3 = PmA[1] * r_130_3 - PmQf[1] * r_130_4 + 3 * iz2 * (r_120_3 - riz * r_120_4);
   FReal2 r_140_4 = PmA[0] * r_040_4 - PmQf[0] * r_040_5;
   pOut[39] = PmA[0] * pOut[24] - PmQf[0] * r_022_1;
   pOut[40] = PmA[0] * pOut[25] - PmQf[0] * r_004_1;
   FReal2 r_104_1 = PmA[0] * r_004_1 - PmQf[0] * r_004_2;
   FReal2 r_104_2 = PmA[0] * r_004_2 - PmQf[0] * r_004_3;
   pOut[41] = PmA[1] * pOut[20] - PmQf[1] * r_400_1;
   FReal2 r_410_1 = PmA[1] * r_400_1 - PmQf[1] * r_400_2;
   FReal2 r_410_2 = PmA[1] * r_400_2 - PmQf[1] * r_400_3;
   FReal2 r_410_3 = PmA[1] * r_400_3 - PmQf[1] * r_400_4;
   pOut[42] = PmA[0] * pOut[27] - PmQf[0] * r_130_1 + iz2 * (pOut[14] - riz * r_030_1);
   FReal2 r_230_1 = PmA[0] * r_130_1 - PmQf[0] * r_130_2 + iz2 * (r_030_1 - riz * r_030_2);
   FReal2 r_230_2 = PmA[0] * r_130_2 - PmQf[0] * r_130_3 + iz2 * (r_030_2 - riz * r_030_3);
   FReal2 r_230_3 = PmA[0] * r_130_3 - PmQf[0] * r_130_4 + iz2 * (r_030_3 - riz * r_030_4);
   FReal2 r_230_4 = PmA[0] * r_130_4 - PmQf[0] * r_130_5 + iz2 * (r_030_4 - riz * r_030_5);
   pOut[43] = PmA[1] * pOut[22] - PmQf[1] * r_202_1;
   FReal2 r_212_1 = PmA[1] * r_202_1 - PmQf[1] * r_202_2;
   pOut[44] = PmA[1] * pOut[23] - PmQf[1] * r_040_1 + 4 * iz2 * (pOut[14] - riz * r_030_1);
   FReal2 r_050_1 = PmA[1] * r_040_1 - PmQf[1] * r_040_2 + 4 * iz2 * (r_030_1 - riz * r_030_2);
   FReal2 r_050_2 = PmA[1] * r_040_2 - PmQf[1] * r_040_3 + 4 * iz2 * (r_030_2 - riz * r_030_3);
   FReal2 r_050_3 = PmA[1] * r_040_3 - PmQf[1] * r_040_4 + 4 * iz2 * (r_030_3 - riz * r_030_4);
   FReal2 r_050_4 = PmA[1] * r_040_4 - PmQf[1] * r_040_5 + 4 * iz2 * (r_030_4 - riz * r_030_5);
   pOut[45] = PmA[1] * pOut[24] - PmQf[1] * r_022_1 + 2 * iz2 * (pOut[15] - riz * r_012_1);
   FReal2 r_032_1 = PmA[1] * r_022_1 - PmQf[1] * r_022_2 + 2 * iz2 * (r_012_1 - riz * r_012_2);
   FReal2 r_032_2 = PmA[2] * r_031_2 - PmQf[2] * r_031_3 + iz2 * (r_030_2 - riz * r_030_3);
   FReal2 r_032_3 = PmA[2] * r_031_3 - PmQf[2] * r_031_4 + iz2 * (r_030_3 - riz * r_030_4);
   pOut[46] = PmA[2] * pOut[34] - PmQf[2] * r_013_1 + 3 * iz2 * (pOut[15] - riz * r_012_1);
   FReal2 r_014_1 = PmA[2] * r_013_1 - PmQf[2] * r_013_2 + 3 * iz2 * (r_012_1 - riz * r_012_2);
   FReal2 r_014_2 = PmA[2] * r_013_2 - PmQf[2] * r_013_3 + 3 * iz2 * (r_012_2 - riz * r_012_3);
   FReal2 r_014_3 = PmA[1] * r_004_3 - PmQf[1] * r_004_4;
   FReal2 r_014_4 = PmA[1] * r_004_4 - PmQf[1] * r_004_5;
   pOut[47] = PmA[0] * pOut[29] - PmQf[0] * r_301_1 + 3 * iz2 * (pOut[16] - riz * r_201_1);
   FReal2 r_401_1 = PmA[0] * r_301_1 - PmQf[0] * r_301_2 + 3 * iz2 * (r_201_1 - riz * r_201_2);
   FReal2 r_401_2 = PmA[0] * r_301_2 - PmQf[0] * r_301_3 + 3 * iz2 * (r_201_2 - riz * r_201_3);
   FReal2 r_401_3 = PmA[2] * r_400_3 - PmQf[2] * r_400_4;
   FReal2 r_401_4 = PmA[2] * r_400_4 - PmQf[2] * r_400_5;
   pOut[48] = PmA[2] * pOut[21] - PmQf[2] * r_220_1;
   FReal2 r_221_1 = PmA[2] * r_220_1 - PmQf[2] * r_220_2;
   FReal2 r_221_2 = PmA[2] * r_220_2 - PmQf[2] * r_220_3;
   pOut[49] = PmA[2] * pOut[22] - PmQf[2] * r_202_1 + 2 * iz2 * (pOut[16] - riz * r_201_1);
   FReal2 r_203_1 = PmA[2] * r_202_1 - PmQf[2] * r_202_2 + 2 * iz2 * (r_201_1 - riz * r_201_2);
   pOut[50] = PmA[2] * pOut[23] - PmQf[2] * r_040_1;
   FReal2 r_041_1 = PmA[2] * r_040_1 - PmQf[2] * r_040_2;
   FReal2 r_041_2 = PmA[2] * r_040_2 - PmQf[2] * r_040_3;
   FReal2 r_041_3 = PmA[2] * r_040_3 - PmQf[2] * r_040_4;
   FReal2 r_041_4 = PmA[2] * r_040_4 - PmQf[2] * r_040_5;
   pOut[51] = PmA[1] * pOut[34] - PmQf[1] * r_013_1 + iz2 * (pOut[18] - riz * r_003_1);
   FReal2 r_023_1 = PmA[1] * r_013_1 - PmQf[1] * r_013_2 + iz2 * (r_003_1 - riz * r_003_2);
   FReal2 r_023_2 = PmA[1] * r_013_2 - PmQf[1] * r_013_3 + iz2 * (r_003_2 - riz * r_003_3);
   FReal2 r_023_3 = PmA[1] * r_013_3 - PmQf[1] * r_013_4 + iz2 * (r_003_3 - riz * r_003_4);
   FReal2 r_023_4 = PmA[1] * r_013_4 - PmQf[1] * r_013_5 + iz2 * (r_003_4 - riz * r_003_5);
   pOut[52] = PmA[2] * pOut[25] - PmQf[2] * r_004_1 + 4 * iz2 * (pOut[18] - riz * r_003_1);
   FReal2 r_005_1 = PmA[2] * r_004_1 - PmQf[2] * r_004_2 + 4 * iz2 * (r_003_1 - riz * r_003_2);
   FReal2 r_005_2 = PmA[2] * r_004_2 - PmQf[2] * r_004_3 + 4 * iz2 * (r_003_2 - riz * r_003_3);
   FReal2 r_005_3 = PmA[2] * r_004_3 - PmQf[2] * r_004_4 + 4 * iz2 * (r_003_3 - riz * r_003_4);
   FReal2 r_005_4 = PmA[2] * r_004_4 - PmQf[2] * r_004_5 + 4 * iz2 * (r_003_4 - riz * r_003_5);
   pOut[53] = PmA[1] * pOut[29] - PmQf[1] * r_301_1;
   pOut[54] = PmA[2] * pOut[27] - PmQf[2] * r_130_1;
   pOut[55] = PmA[0] * pOut[34] - PmQf[0] * r_013_1;

   // L = 6
   pOut[56] = PmA[0] * pOut[35] - PmQf[0] * r_500_1 + 5 * iz2 * (pOut[20] - riz * r_400_1);
   FReal2 r_600_1 = PmA[0] * r_500_1 - PmQf[0] * r_500_2 + 5 * iz2 * (r_400_1 - riz * r_400_2);
   FReal2 r_600_2 = PmA[0] * r_500_2 - PmQf[0] * r_500_3 + 5 * iz2 * (r_400_2 - riz * r_400_3);
   FReal2 r_600_3 = PmA[0] * r_500_3 - PmQf[0] * r_500_4 + 5 * iz2 * (r_400_3 - riz * r_400_4);
   pOut[57] = PmA[0] * pOut[36] - PmQf[0] * r_320_1 + 3 * iz2 * (pOut[21] - riz * r_220_1);
   FReal2 r_420_1 = PmA[1] * r_410_1 - PmQf[1] * r_410_2 + iz2 * (r_400_1 - riz * r_400_2);
   FReal2 r_420_2 = PmA[1] * r_410_2 - PmQf[1] * r_410_3 + iz2 * (r_400_2 - riz * r_400_3);
   pOut[58] = PmA[0] * pOut[37] - PmQf[0] * r_302_1 + 3 * iz2 * (pOut[22] - riz * r_202_1);
   FReal2 r_402_1 = PmA[0] * r_302_1 - PmQf[0] * r_302_2 + 3 * iz2 * (r_202_1 - riz * r_202_2);
   FReal2 r_402_2 = PmA[2] * r_401_2 - PmQf[2] * r_401_3 + iz2 * (r_400_2 - riz * r_400_3);
   FReal2 r_402_3 = PmA[2] * r_401_3 - PmQf[2] * r_401_4 + iz2 * (r_400_3 - riz * r_400_4);
   pOut[59] = PmA[1] * pOut[42] - PmQf[1] * r_230_1 + 3 * iz2 * (pOut[21] - riz * r_220_1);
   FReal2 r_240_1 = PmA[1] * r_230_1 - PmQf[1] * r_230_2 + 3 * iz2 * (r_220_1 - riz * r_220_2);
   FReal2 r_240_2 = PmA[1] * r_230_2 - PmQf[1] * r_230_3 + 3 * iz2 * (r_220_2 - riz * r_220_3);
   FReal2 r_240_3 = PmA[0] * r_140_3 - PmQf[0] * r_140_4 + iz2 * (r_040_3 - riz * r_040_4);
   pOut[60] = PmA[1] * pOut[43] - PmQf[1] * r_212_1 + iz2 * (pOut[22] - riz * r_202_1);
   FReal2 r_222_1 = PmA[2] * r_221_1 - PmQf[2] * r_221_2 + iz2 * (r_220_1 - riz * r_220_2);
   pOut[61] = PmA[2] * pOut[49] - PmQf[2] * r_203_1 + 3 * iz2 * (pOut[22] - riz * r_202_1);
   FReal2 r_204_1 = PmA[0] * r_104_1 - PmQf[0] * r_104_2 + iz2 * (r_004_1 - riz * r_004_2);
   pOut[62] = PmA[1] * pOut[44] - PmQf[1] * r_050_1 + 5 * iz2 * (pOut[23] - riz * r_040_1);
   FReal2 r_060_1 = PmA[1] * r_050_1 - PmQf[1] * r_050_2 + 5 * iz2 * (r_040_1 - riz * r_040_2);
   FReal2 r_060_2 = PmA[1] * r_050_2 - PmQf[1] * r_050_3 + 5 * iz2 * (r_040_2 - riz * r_040_3);
   FReal2 r_060_3 = PmA[1] * r_050_3 - PmQf[1] * r_050_4 + 5 * iz2 * (r_040_3 - riz * r_040_4);
   pOut[63] = PmA[1] * pOut[45] - PmQf[1] * r_032_1 + 3 * iz2 * (pOut[24] - riz * r_022_1);
   FReal2 r_042_1 = PmA[1] * r_032_1 - PmQf[1] * r_032_2 + 3 * iz2 * (r_022_1 - riz * r_022_2);
   FReal2 r_042_2 = PmA[2] * r_041_2 - PmQf[2] * r_041_3 + iz2 * (r_040_2 - riz * r_040_3);
   FReal2 r_042_3 = PmA[2] * r_041_3 - PmQf[2] * r_041_4 + iz2 * (r_040_3 - riz * r_040_4);
   pOut[64] = PmA[2] * pOut[51] - PmQf[2] * r_023_1 + 3 * iz2 * (pOut[24] - riz * r_022_1);
   FReal2 r_024_1 = PmA[2] * r_023_1 - PmQf[2] * r_023_2 + 3 * iz2 * (r_022_1 - riz * r_022_2);
   FReal2 r_024_2 = PmA[1] * r_014_2 - PmQf[1] * r_014_3 + iz2 * (r_004_2 - riz * r_004_3);
   FReal2 r_024_3 = PmA[1] * r_014_3 - PmQf[1] * r_014_4 + iz2 * (r_004_3 - riz * r_004_4);
   pOut[65] = PmA[2] * pOut[52] - PmQf[2] * r_005_1 + 5 * iz2 * (pOut[25] - riz * r_004_1);
   FReal2 r_006_1 = PmA[2] * r_005_1 - PmQf[2] * r_005_2 + 5 * iz2 * (r_004_1 - riz * r_004_2);
   FReal2 r_006_2 = PmA[2] * r_005_2 - PmQf[2] * r_005_3 + 5 * iz2 * (r_004_2 - riz * r_004_3);
   FReal2 r_006_3 = PmA[2] * r_005_3 - PmQf[2] * r_005_4 + 5 * iz2 * (r_004_3 - riz * r_004_4);
   pOut[66] = PmA[1] * pOut[35] - PmQf[1] * r_500_1;
   FReal2 r_510_1 = PmA[1] * r_500_1 - PmQf[1] * r_500_2;
   FReal2 r_510_2 = PmA[1] * r_500_2 - PmQf[1] * r_500_3;
   FReal2 r_510_3 = PmA[1] * r_500_3 - PmQf[1] * r_500_4;
   pOut[67] = PmA[0] * pOut[42] - PmQf[0] * r_230_1 + 2 * iz2 * (pOut[27] - riz * r_130_1);
   FReal2 r_330_1 = PmA[0] * r_230_1 - PmQf[0] * r_230_2 + 2 * iz2 * (r_130_1 - riz * r_130_2);
   pOut[68] = PmA[1] * pOut[37] - PmQf[1] * r_302_1;
   FReal2 r_312_1 = PmA[1] * r_302_1 - PmQf[1] * r_302_2;
   FReal2 r_312_2 = PmA[1] * r_302_2 - PmQf[1] * r_302_3;
   pOut[69] = PmA[1] * pOut[38] - PmQf[1] * r_140_1 + 4 * iz2 * (pOut[27] - riz * r_130_1);
   FReal2 r_150_1 = PmA[1] * r_140_1 - PmQf[1] * r_140_2 + 4 * iz2 * (r_130_1 - riz * r_130_2);
   FReal2 r_150_2 = PmA[1] * r_140_2 - PmQf[1] * r_140_3 + 4 * iz2 * (r_130_2 - riz * r_130_3);
   pOut[70] = PmA[0] * pOut[45] - PmQf[0] * r_032_1;
   FReal2 r_132_1 = PmA[0] * r_032_1 - PmQf[0] * r_032_2;
   FReal2 r_132_2 = PmA[0] * r_032_2 - PmQf[0] * r_032_3;
   pOut[71] = PmA[0] * pOut[46] - PmQf[0] * r_014_1;
   pOut[72] = PmA[0] * pOut[47] - PmQf[0] * r_401_1 + 4 * iz2 * (pOut[29] - riz * r_301_1);
   FReal2 r_501_1 = PmA[0] * r_401_1 - PmQf[0] * r_401_2 + 4 * iz2 * (r_301_1 - riz * r_301_2);
   pOut[73] = PmA[2] * pOut[36] - PmQf[2] * r_320_1;
   pOut[74] = PmA[2] * pOut[37] - PmQf[2] * r_302_1 + 2 * iz2 * (pOut[29] - riz * r_301_1);
   FReal2 r_303_1 = PmA[2] * r_302_1 - PmQf[2] * r_302_2 + 2 * iz2 * (r_301_1 - riz * r_301_2);
   FReal2 r_303_2 = PmA[2] * r_302_2 - PmQf[2] * r_302_3 + 2 * iz2 * (r_301_2 - riz * r_301_3);
   FReal2 r_303_3 = PmA[2] * r_302_3 - PmQf[2] * r_302_4 + 2 * iz2 * (r_301_3 - riz * r_301_4);
   pOut[75] = PmA[2] * pOut[38] - PmQf[2] * r_140_1;
   pOut[76] = PmA[0] * pOut[51] - PmQf[0] * r_023_1;
   FReal2 r_123_1 = PmA[0] * r_023_1 - PmQf[0] * r_023_2;
   FReal2 r_123_2 = PmA[0] * r_023_2 - PmQf[0] * r_023_3;
   pOut[77] = PmA[0] * pOut[52] - PmQf[0] * r_005_1;
   FReal2 r_105_1 = PmA[0] * r_005_1 - PmQf[0] * r_005_2;
   FReal2 r_105_2 = PmA[0] * r_005_2 - PmQf[0] * r_005_3;
   pOut[78] = PmA[1] * pOut[47] - PmQf[1] * r_401_1;
   pOut[79] = PmA[2] * pOut[42] - PmQf[2] * r_230_1;
   FReal2 r_231_1 = PmA[2] * r_230_1 - PmQf[2] * r_230_2;
   FReal2 r_231_2 = PmA[2] * r_230_2 - PmQf[2] * r_230_3;
   FReal2 r_231_3 = PmA[2] * r_230_3 - PmQf[2] * r_230_4;
   pOut[80] = PmA[1] * pOut[49] - PmQf[1] * r_203_1;
   pOut[81] = PmA[2] * pOut[44] - PmQf[2] * r_050_1;
   pOut[82] = PmA[1] * pOut[51] - PmQf[1] * r_023_1 + 2 * iz2 * (pOut[34] - riz * r_013_1);
   FReal2 r_033_1 = PmA[1] * r_023_1 - PmQf[1] * r_023_2 + 2 * iz2 * (r_013_1 - riz * r_013_2);
   FReal2 r_033_2 = PmA[1] * r_023_2 - PmQf[1] * r_023_3 + 2 * iz2 * (r_013_2 - riz * r_013_3);
   FReal2 r_033_3 = PmA[1] * r_023_3 - PmQf[1] * r_023_4 + 2 * iz2 * (r_013_3 - riz * r_013_4);
   pOut[83] = PmA[2] * pOut[46] - PmQf[2] * r_014_1 + 4 * iz2 * (pOut[34] - riz * r_013_1);
   FReal2 r_015_1 = PmA[2] * r_014_1 - PmQf[2] * r_014_2 + 4 * iz2 * (r_013_1 - riz * r_013_2);

   // L = 7
   pOut[84] = PmA[0] * pOut[56] - PmQf[0] * r_600_1 + 6 * iz2 * (pOut[35] - riz * r_500_1);
   FReal2 r_700_1 = PmA[0] * r_600_1 - PmQf[0] * r_600_2 + 6 * iz2 * (r_500_1 - riz * r_500_2);
   FReal2 r_700_2 = PmA[0] * r_600_2 - PmQf[0] * r_600_3 + 6 * iz2 * (r_500_2 - riz * r_500_3);
   pOut[85] = PmA[0] * pOut[57] - PmQf[0] * r_420_1 + 4 * iz2 * (pOut[36] - riz * r_320_1);
   FReal2 r_520_1 = PmA[1] * r_510_1 - PmQf[1] * r_510_2 + iz2 * (r_500_1 - riz * r_500_2);
   FReal2 r_520_2 = PmA[1] * r_510_2 - PmQf[1] * r_510_3 + iz2 * (r_500_2 - riz * r_500_3);
   pOut[86] = PmA[2] * pOut[72] - PmQf[2] * r_501_1 + iz2 * (pOut[35] - riz * r_500_1);
   pOut[87] = PmA[1] * pOut[67] - PmQf[1] * r_330_1 + 3 * iz2 * (pOut[36] - riz * r_320_1);
   FReal2 r_340_1 = PmA[0] * r_240_1 - PmQf[0] * r_240_2 + 2 * iz2 * (r_140_1 - riz * r_140_2);
   FReal2 r_340_2 = PmA[0] * r_240_2 - PmQf[0] * r_240_3 + 2 * iz2 * (r_140_2 - riz * r_140_3);
   pOut[88] = PmA[1] * pOut[68] - PmQf[1] * r_312_1 + iz2 * (pOut[37] - riz * r_302_1);
   FReal2 r_322_1 = PmA[1] * r_312_1 - PmQf[1] * r_312_2 + iz2 * (r_302_1 - riz * r_302_2);
   pOut[89] = PmA[2] * pOut[74] - PmQf[2] * r_303_1 + 3 * iz2 * (pOut[37] - riz * r_302_1);
   FReal2 r_304_1 = PmA[2] * r_303_1 - PmQf[2] * r_303_2 + 3 * iz2 * (r_302_1 - riz * r_302_2);
   FReal2 r_304_2 = PmA[2] * r_303_2 - PmQf[2] * r_303_3 + 3 * iz2 * (r_302_2 - riz * r_302_3);
   pOut[90] = PmA[1] * pOut[69] - PmQf[1] * r_150_1 + 5 * iz2 * (pOut[38] - riz * r_140_1);
   FReal2 r_160_1 = PmA[1] * r_150_1 - PmQf[1] * r_150_2 + 5 * iz2 * (r_140_1 - riz * r_140_2);
   pOut[91] = PmA[0] * pOut[63] - PmQf[0] * r_042_1;
   pOut[92] = PmA[0] * pOut[64] - PmQf[0] * r_024_1;
   FReal2 r_124_1 = PmA[0] * r_024_1 - PmQf[0] * r_024_2;
   FReal2 r_124_2 = PmA[0] * r_024_2 - PmQf[0] * r_024_3;
   pOut[93] = PmA[0] * pOut[65] - PmQf[0] * r_006_1;
   FReal2 r_106_1 = PmA[0] * r_006_1 - PmQf[0] * r_006_2;
   FReal2 r_106_2 = PmA[0] * r_006_2 - PmQf[0] * r_006_3;
   pOut[94] = PmA[1] * pOut[56] - PmQf[1] * r_600_1;
   pOut[95] = PmA[1] * pOut[57] - PmQf[1] * r_420_1 + 2 * iz2 * (pOut[41] - riz * r_410_1);
   FReal2 r_430_1 = PmA[1] * r_420_1 - PmQf[1] * r_420_2 + 2 * iz2 * (r_410_1 - riz * r_410_2);
   pOut[96] = PmA[0] * pOut[68] - PmQf[0] * r_312_1 + 3 * iz2 * (pOut[43] - riz * r_212_1);
   FReal2 r_412_1 = PmA[1] * r_402_1 - PmQf[1] * r_402_2;
   FReal2 r_412_2 = PmA[1] * r_402_2 - PmQf[1] * r_402_3;
   pOut[97] = PmA[0] * pOut[69] - PmQf[0] * r_150_1 + iz2 * (pOut[44] - riz * r_050_1);
   FReal2 r_250_1 = PmA[1] * r_240_1 - PmQf[1] * r_240_2 + 4 * iz2 * (r_230_1 - riz * r_230_2);
   FReal2 r_250_2 = PmA[1] * r_240_2 - PmQf[1] * r_240_3 + 4 * iz2 * (r_230_2 - riz * r_230_3);
   pOut[98] = PmA[1] * pOut[60] - PmQf[1] * r_222_1 + 2 * iz2 * (pOut[43] - riz * r_212_1);
   FReal2 r_232_1 = PmA[2] * r_231_1 - PmQf[2] * r_231_2 + iz2 * (r_230_1 - riz * r_230_2);
   FReal2 r_232_2 = PmA[2] * r_231_2 - PmQf[2] * r_231_3 + iz2 * (r_230_2 - riz * r_230_3);
   pOut[99] = PmA[1] * pOut[61] - PmQf[1] * r_204_1;
   pOut[100] = PmA[1] * pOut[62] - PmQf[1] * r_060_1 + 6 * iz2 * (pOut[44] - riz * r_050_1);
   FReal2 r_070_1 = PmA[1] * r_060_1 - PmQf[1] * r_060_2 + 6 * iz2 * (r_050_1 - riz * r_050_2);
   FReal2 r_070_2 = PmA[1] * r_060_2 - PmQf[1] * r_060_3 + 6 * iz2 * (r_050_2 - riz * r_050_3);
   pOut[101] = PmA[1] * pOut[63] - PmQf[1] * r_042_1 + 4 * iz2 * (pOut[45] - riz * r_032_1);
   FReal2 r_052_1 = PmA[1] * r_042_1 - PmQf[1] * r_042_2 + 4 * iz2 * (r_032_1 - riz * r_032_2);
   FReal2 r_052_2 = PmA[1] * r_042_2 - PmQf[1] * r_042_3 + 4 * iz2 * (r_032_2 - riz * r_032_3);
   pOut[102] = PmA[1] * pOut[64] - PmQf[1] * r_024_1 + 2 * iz2 * (pOut[46] - riz * r_014_1);
   FReal2 r_034_1 = PmA[2] * r_033_1 - PmQf[2] * r_033_2 + 3 * iz2 * (r_032_1 - riz * r_032_2);
   FReal2 r_034_2 = PmA[2] * r_033_2 - PmQf[2] * r_033_3 + 3 * iz2 * (r_032_2 - riz * r_032_3);
   pOut[103] = PmA[2] * pOut[83] - PmQf[2] * r_015_1 + 5 * iz2 * (pOut[46] - riz * r_014_1);
   FReal2 r_016_1 = PmA[1] * r_006_1 - PmQf[1] * r_006_2;
   FReal2 r_016_2 = PmA[1] * r_006_2 - PmQf[1] * r_006_3;
   pOut[104] = PmA[0] * pOut[72] - PmQf[0] * r_501_1 + 5 * iz2 * (pOut[47] - riz * r_401_1);
   FReal2 r_601_1 = PmA[2] * r_600_1 - PmQf[2] * r_600_2;
   FReal2 r_601_2 = PmA[2] * r_600_2 - PmQf[2] * r_600_3;
   pOut[105] = PmA[2] * pOut[57] - PmQf[2] * r_420_1;
   pOut[106] = PmA[0] * pOut[74] - PmQf[0] * r_303_1 + 3 * iz2 * (pOut[49] - riz * r_203_1);
   FReal2 r_403_1 = PmA[2] * r_402_1 - PmQf[2] * r_402_2 + 2 * iz2 * (r_401_1 - riz * r_401_2);
   FReal2 r_403_2 = PmA[2] * r_402_2 - PmQf[2] * r_402_3 + 2 * iz2 * (r_401_2 - riz * r_401_3);
   pOut[107] = PmA[2] * pOut[59] - PmQf[2] * r_240_1;
   FReal2 r_241_1 = PmA[2] * r_240_1 - PmQf[2] * r_240_2;
   FReal2 r_241_2 = PmA[2] * r_240_2 - PmQf[2] * r_240_3;
   pOut[108] = PmA[0] * pOut[76] - PmQf[0] * r_123_1 + iz2 * (pOut[51] - riz * r_023_1);
   FReal2 r_223_1 = PmA[0] * r_123_1 - PmQf[0] * r_123_2 + iz2 * (r_023_1 - riz * r_023_2);
   pOut[109] = PmA[2] * pOut[61] - PmQf[2] * r_204_1 + 4 * iz2 * (pOut[49] - riz * r_203_1);
   FReal2 r_205_1 = PmA[0] * r_105_1 - PmQf[0] * r_105_2 + iz2 * (r_005_1 - riz * r_005_2);
   pOut[110] = PmA[2] * pOut[62] - PmQf[2] * r_060_1;
   pOut[111] = PmA[1] * pOut[82] - PmQf[1] * r_033_1 + 3 * iz2 * (pOut[51] - riz * r_023_1);
   FReal2 r_043_1 = PmA[2] * r_042_1 - PmQf[2] * r_042_2 + 2 * iz2 * (r_041_1 - riz * r_041_2);
   FReal2 r_043_2 = PmA[1] * r_033_2 - PmQf[1] * r_033_3 + 3 * iz2 * (r_023_2 - riz * r_023_3);
   pOut[112] = PmA[1] * pOut[83] - PmQf[1] * r_015_1 + iz2 * (pOut[52] - riz * r_005_1);
   pOut[113] = PmA[2] * pOut[65] - PmQf[2] * r_006_1 + 6 * iz2 * (pOut[52] - riz * r_005_1);
   FReal2 r_007_1 = PmA[2] * r_006_1 - PmQf[2] * r_006_2 + 6 * iz2 * (r_005_1 - riz * r_005_2);
   FReal2 r_007_2 = PmA[2] * r_006_2 - PmQf[2] * r_006_3 + 6 * iz2 * (r_005_2 - riz * r_005_3);
   pOut[114] = PmA[1] * pOut[72] - PmQf[1] * r_501_1;
   pOut[115] = PmA[2] * pOut[67] - PmQf[2] * r_330_1;
   pOut[116] = PmA[1] * pOut[74] - PmQf[1] * r_303_1;
   FReal2 r_313_1 = PmA[1] * r_303_1 - PmQf[1] * r_303_2;
   FReal2 r_313_2 = PmA[1] * r_303_2 - PmQf[1] * r_303_3;
   pOut[117] = PmA[2] * pOut[69] - PmQf[2] * r_150_1;
   pOut[118] = PmA[0] * pOut[82] - PmQf[0] * r_033_1;
   FReal2 r_133_1 = PmA[0] * r_033_1 - PmQf[0] * r_033_2;
   FReal2 r_133_2 = PmA[0] * r_033_2 - PmQf[0] * r_033_3;
   pOut[119] = PmA[0] * pOut[83] - PmQf[0] * r_015_1;

   // L = 8
   pOut[120] = PmA[0] * pOut[84] - PmQf[0] * r_700_1 + 7 * iz2 * (pOut[56] - riz * r_600_1);
   FReal2 r_800_1 = PmA[0] * r_700_1 - PmQf[0] * r_700_2 + 7 * iz2 * (r_600_1 - riz * r_600_2);
   pOut[121] = PmA[0] * pOut[85] - PmQf[0] * r_520_1 + 5 * iz2 * (pOut[57] - riz * r_420_1);
   FReal2 r_620_1 = PmA[0] * r_520_1 - PmQf[0] * r_520_2 + 5 * iz2 * (r_420_1 - riz * r_420_2);
   pOut[122] = PmA[2] * pOut[104] - PmQf[2] * r_601_1 + iz2 * (pOut[56] - riz * r_600_1);
   FReal2 r_602_1 = PmA[2] * r_601_1 - PmQf[2] * r_601_2 + iz2 * (r_600_1 - riz * r_600_2);
   pOut[123] = PmA[0] * pOut[87] - PmQf[0] * r_340_1 + 3 * iz2 * (pOut[59] - riz * r_240_1);
   FReal2 r_440_1 = PmA[0] * r_340_1 - PmQf[0] * r_340_2 + 3 * iz2 * (r_240_1 - riz * r_240_2);
   pOut[124] = PmA[0] * pOut[88] - PmQf[0] * r_322_1 + 3 * iz2 * (pOut[60] - riz * r_222_1);
   FReal2 r_422_1 = PmA[1] * r_412_1 - PmQf[1] * r_412_2 + iz2 * (r_402_1 - riz * r_402_2);
   pOut[125] = PmA[2] * pOut[106] - PmQf[2] * r_403_1 + 3 * iz2 * (pOut[58] - riz * r_402_1);
   FReal2 r_404_1 = PmA[2] * r_403_1 - PmQf[2] * r_403_2 + 3 * iz2 * (r_402_1 - riz * r_402_2);
   pOut[126] = PmA[1] * pOut[97] - PmQf[1] * r_250_1 + 5 * iz2 * (pOut[59] - riz * r_240_1);
   FReal2 r_260_1 = PmA[1] * r_250_1 - PmQf[1] * r_250_2 + 5 * iz2 * (r_240_1 - riz * r_240_2);
   pOut[127] = PmA[1] * pOut[98] - PmQf[1] * r_232_1 + 3 * iz2 * (pOut[60] - riz * r_222_1);
   FReal2 r_242_1 = PmA[2] * r_241_1 - PmQf[2] * r_241_2 + iz2 * (r_240_1 - riz * r_240_2);
   pOut[128] = PmA[2] * pOut[108] - PmQf[2] * r_223_1 + 3 * iz2 * (pOut[60] - riz * r_222_1);
   FReal2 r_224_1 = PmA[0] * r_124_1 - PmQf[0] * r_124_2 + iz2 * (r_024_1 - riz * r_024_2);
   pOut[129] = PmA[2] * pOut[109] - PmQf[2] * r_205_1 + 5 * iz2 * (pOut[61] - riz * r_204_1);
   FReal2 r_206_1 = PmA[0] * r_106_1 - PmQf[0] * r_106_2 + iz2 * (r_006_1 - riz * r_006_2);
   pOut[130] = PmA[1] * pOut[100] - PmQf[1] * r_070_1 + 7 * iz2 * (pOut[62] - riz * r_060_1);
   FReal2 r_080_1 = PmA[1] * r_070_1 - PmQf[1] * r_070_2 + 7 * iz2 * (r_060_1 - riz * r_060_2);
   pOut[131] = PmA[1] * pOut[101] - PmQf[1] * r_052_1 + 5 * iz2 * (pOut[63] - riz * r_042_1);
   FReal2 r_062_1 = PmA[1] * r_052_1 - PmQf[1] * r_052_2 + 5 * iz2 * (r_042_1 - riz * r_042_2);
   pOut[132] = PmA[1] * pOut[102] - PmQf[1] * r_034_1 + 3 * iz2 * (pOut[64] - riz * r_024_1);
   FReal2 r_044_1 = PmA[1] * r_034_1 - PmQf[1] * r_034_2 + 3 * iz2 * (r_024_1 - riz * r_024_2);
   pOut[133] = PmA[1] * pOut[103] - PmQf[1] * r_016_1 + iz2 * (pOut[65] - riz * r_006_1);
   FReal2 r_026_1 = PmA[1] * r_016_1 - PmQf[1] * r_016_2 + iz2 * (r_006_1 - riz * r_006_2);
   pOut[134] = PmA[2] * pOut[113] - PmQf[2] * r_007_1 + 7 * iz2 * (pOut[65] - riz * r_006_1);
   FReal2 r_008_1 = PmA[2] * r_007_1 - PmQf[2] * r_007_2 + 7 * iz2 * (r_006_1 - riz * r_006_2);
   pOut[135] = PmA[1] * pOut[84] - PmQf[1] * r_700_1;
   pOut[136] = PmA[1] * pOut[85] - PmQf[1] * r_520_1 + 2 * iz2 * (pOut[66] - riz * r_510_1);
   FReal2 r_530_1 = PmA[1] * r_520_1 - PmQf[1] * r_520_2 + 2 * iz2 * (r_510_1 - riz * r_510_2);
   pOut[137] = PmA[0] * pOut[96] - PmQf[0] * r_412_1 + 4 * iz2 * (pOut[68] - riz * r_312_1);
   pOut[138] = PmA[1] * pOut[87] - PmQf[1] * r_340_1 + 4 * iz2 * (pOut[67] - riz * r_330_1);
   FReal2 r_350_1 = PmA[0] * r_250_1 - PmQf[0] * r_250_2 + 2 * iz2 * (r_150_1 - riz * r_150_2);
   pOut[139] = PmA[0] * pOut[98] - PmQf[0] * r_232_1 + 2 * iz2 * (pOut[70] - riz * r_132_1);
   FReal2 r_332_1 = PmA[0] * r_232_1 - PmQf[0] * r_232_2 + 2 * iz2 * (r_132_1 - riz * r_132_2);
   pOut[140] = PmA[2] * pOut[116] - PmQf[2] * r_313_1 + 3 * iz2 * (pOut[68] - riz * r_312_1);
   FReal2 r_314_1 = PmA[1] * r_304_1 - PmQf[1] * r_304_2;
   pOut[141] = PmA[1] * pOut[90] - PmQf[1] * r_160_1 + 6 * iz2 * (pOut[69] - riz * r_150_1);
   FReal2 r_170_1 = PmA[0] * r_070_1 - PmQf[0] * r_070_2;
   pOut[142] = PmA[0] * pOut[101] - PmQf[0] * r_052_1;
   pOut[143] = PmA[2] * pOut[118] - PmQf[2] * r_133_1 + 3 * iz2 * (pOut[70] - riz * r_132_1);
   pOut[144] = PmA[0] * pOut[103] - PmQf[0] * r_016_1;
   pOut[145] = PmA[0] * pOut[104] - PmQf[0] * r_601_1 + 6 * iz2 * (pOut[72] - riz * r_501_1);
   FReal2 r_701_1 = PmA[2] * r_700_1 - PmQf[2] * r_700_2;
   pOut[146] = PmA[2] * pOut[85] - PmQf[2] * r_520_1;
   pOut[147] = PmA[0] * pOut[106] - PmQf[0] * r_403_1 + 4 * iz2 * (pOut[74] - riz * r_303_1);
   FReal2 r_503_1 = PmA[0] * r_403_1 - PmQf[0] * r_403_2 + 4 * iz2 * (r_303_1 - riz * r_303_2);
   pOut[148] = PmA[2] * pOut[87] - PmQf[2] * r_340_1;
   FReal2 r_341_1 = PmA[2] * r_340_1 - PmQf[2] * r_340_2;
   pOut[149] = PmA[1] * pOut[116] - PmQf[1] * r_313_1 + iz2 * (pOut[74] - riz * r_303_1);
   FReal2 r_323_1 = PmA[1] * r_313_1 - PmQf[1] * r_313_2 + iz2 * (r_303_1 - riz * r_303_2);
   pOut[150] = PmA[0] * pOut[109] - PmQf[0] * r_205_1 + 2 * iz2 * (pOut[77] - riz * r_105_1);
   pOut[151] = PmA[2] * pOut[90] - PmQf[2] * r_160_1;
   pOut[152] = PmA[1] * pOut[118] - PmQf[1] * r_133_1 + 3 * iz2 * (pOut[76] - riz * r_123_1);
   pOut[153] = PmA[2] * pOut[92] - PmQf[2] * r_124_1 + 4 * iz2 * (pOut[76] - riz * r_123_1);
   pOut[154] = PmA[0] * pOut[113] - PmQf[0] * r_007_1;
   pOut[155] = PmA[1] * pOut[104] - PmQf[1] * r_601_1;
   pOut[156] = PmA[2] * pOut[95] - PmQf[2] * r_430_1;
   pOut[157] = PmA[1] * pOut[106] - PmQf[1] * r_403_1;
   pOut[158] = PmA[1] * pOut[107] - PmQf[1] * r_241_1 + 4 * iz2 * (pOut[79] - riz * r_231_1);
   pOut[159] = PmA[0] * pOut[118] - PmQf[0] * r_133_1 + iz2 * (pOut[82] - riz * r_033_1);
   FReal2 r_233_1 = PmA[0] * r_133_1 - PmQf[0] * r_133_2 + iz2 * (r_033_1 - riz * r_033_2);
   pOut[160] = PmA[1] * pOut[109] - PmQf[1] * r_205_1;
   pOut[161] = PmA[2] * pOut[100] - PmQf[2] * r_070_1;
   pOut[162] = PmA[1] * pOut[111] - PmQf[1] * r_043_1 + 4 * iz2 * (pOut[82] - riz * r_033_1);
   FReal2 r_053_1 = PmA[1] * r_043_1 - PmQf[1] * r_043_2 + 4 * iz2 * (r_033_1 - riz * r_033_2);
   pOut[163] = PmA[2] * pOut[102] - PmQf[2] * r_034_1 + 4 * iz2 * (pOut[82] - riz * r_033_1);
   FReal2 r_035_1 = PmA[2] * r_034_1 - PmQf[2] * r_034_2 + 4 * iz2 * (r_033_1 - riz * r_033_2);
   pOut[164] = PmA[2] * pOut[103] - PmQf[2] * r_016_1 + 6 * iz2 * (pOut[83] - riz * r_015_1);
   FReal2 r_017_1 = PmA[1] * r_007_1 - PmQf[1] * r_007_2;

   // L = 9
   pOut[165] = PmA[0] * pOut[120] - PmQf[0] * r_800_1 + 8 * iz2 * (pOut[84] - riz * r_700_1);
   pOut[166] = PmA[0] * pOut[121] - PmQf[0] * r_620_1 + 6 * iz2 * (pOut[85] - riz * r_520_1);
   pOut[167] = PmA[2] * pOut[145] - PmQf[2] * r_701_1 + iz2 * (pOut[84] - riz * r_700_1);
   pOut[168] = PmA[1] * pOut[136] - PmQf[1] * r_530_1 + 3 * iz2 * (pOut[85] - riz * r_520_1);
   pOut[169] = PmA[0] * pOut[124] - PmQf[0] * r_422_1 + 4 * iz2 * (pOut[88] - riz * r_322_1);
   pOut[170] = PmA[0] * pOut[125] - PmQf[0] * r_404_1 + 4 * iz2 * (pOut[89] - riz * r_304_1);
   pOut[171] = PmA[0] * pOut[126] - PmQf[0] * r_260_1 + 2 * iz2 * (pOut[90] - riz * r_160_1);
   pOut[172] = PmA[2] * pOut[148] - PmQf[2] * r_341_1 + iz2 * (pOut[87] - riz * r_340_1);
   pOut[173] = PmA[1] * pOut[140] - PmQf[1] * r_314_1 + iz2 * (pOut[89] - riz * r_304_1);
   pOut[174] = PmA[0] * pOut[129] - PmQf[0] * r_206_1 + 2 * iz2 * (pOut[93] - riz * r_106_1);
   pOut[175] = PmA[0] * pOut[130] - PmQf[0] * r_080_1;
   pOut[176] = PmA[0] * pOut[131] - PmQf[0] * r_062_1;
   pOut[177] = PmA[0] * pOut[132] - PmQf[0] * r_044_1;
   pOut[178] = PmA[0] * pOut[133] - PmQf[0] * r_026_1;
   pOut[179] = PmA[0] * pOut[134] - PmQf[0] * r_008_1;
   pOut[180] = PmA[1] * pOut[120] - PmQf[1] * r_800_1;
   pOut[181] = PmA[0] * pOut[136] - PmQf[0] * r_530_1 + 5 * iz2 * (pOut[95] - riz * r_430_1);
   pOut[182] = PmA[1] * pOut[122] - PmQf[1] * r_602_1;
   pOut[183] = PmA[0] * pOut[138] - PmQf[0] * r_350_1 + 3 * iz2 * (pOut[97] - riz * r_250_1);
   pOut[184] = PmA[0] * pOut[139] - PmQf[0] * r_332_1 + 3 * iz2 * (pOut[98] - riz * r_232_1);
   pOut[185] = PmA[1] * pOut[125] - PmQf[1] * r_404_1;
   pOut[186] = PmA[0] * pOut[141] - PmQf[0] * r_170_1 + iz2 * (pOut[100] - riz * r_070_1);
   pOut[187] = PmA[1] * pOut[127] - PmQf[1] * r_242_1 + 4 * iz2 * (pOut[98] - riz * r_232_1);
   pOut[188] = PmA[2] * pOut[159] - PmQf[2] * r_233_1 + 3 * iz2 * (pOut[98] - riz * r_232_1);
   pOut[189] = PmA[1] * pOut[129] - PmQf[1] * r_206_1;
   pOut[190] = PmA[1] * pOut[130] - PmQf[1] * r_080_1 + 8 * iz2 * (pOut[100] - riz * r_070_1);
   pOut[191] = PmA[1] * pOut[131] - PmQf[1] * r_062_1 + 6 * iz2 * (pOut[101] - riz * r_052_1);
   pOut[192] = PmA[2] * pOut[162] - PmQf[2] * r_053_1 + 3 * iz2 * (pOut[101] - riz * r_052_1);
   pOut[193] = PmA[1] * pOut[133] - PmQf[1] * r_026_1 + 2 * iz2 * (pOut[103] - riz * r_016_1);
   pOut[194] = PmA[2] * pOut[164] - PmQf[2] * r_017_1 + 7 * iz2 * (pOut[103] - riz * r_016_1);
   pOut[195] = PmA[0] * pOut[145] - PmQf[0] * r_701_1 + 7 * iz2 * (pOut[104] - riz * r_601_1);
   pOut[196] = PmA[2] * pOut[121] - PmQf[2] * r_620_1;
   pOut[197] = PmA[2] * pOut[122] - PmQf[2] * r_602_1 + 2 * iz2 * (pOut[104] - riz * r_601_1);
   pOut[198] = PmA[2] * pOut[123] - PmQf[2] * r_440_1;
   pOut[199] = PmA[0] * pOut[149] - PmQf[0] * r_323_1 + 3 * iz2 * (pOut[108] - riz * r_223_1);
   pOut[200] = PmA[2] * pOut[125] - PmQf[2] * r_404_1 + 4 * iz2 * (pOut[106] - riz * r_403_1);
   pOut[201] = PmA[2] * pOut[126] - PmQf[2] * r_260_1;
   pOut[202] = PmA[1] * pOut[159] - PmQf[1] * r_233_1 + 3 * iz2 * (pOut[108] - riz * r_223_1);
   pOut[203] = PmA[2] * pOut[128] - PmQf[2] * r_224_1 + 4 * iz2 * (pOut[108] - riz * r_223_1);
   pOut[204] = PmA[2] * pOut[129] - PmQf[2] * r_206_1 + 6 * iz2 * (pOut[109] - riz * r_205_1);
   pOut[205] = PmA[2] * pOut[130] - PmQf[2] * r_080_1;
   pOut[206] = PmA[1] * pOut[162] - PmQf[1] * r_053_1 + 5 * iz2 * (pOut[111] - riz * r_043_1);
   pOut[207] = PmA[2] * pOut[132] - PmQf[2] * r_044_1 + 4 * iz2 * (pOut[111] - riz * r_043_1);
   pOut[208] = PmA[1] * pOut[164] - PmQf[1] * r_017_1 + iz2 * (pOut[113] - riz * r_007_1);
   pOut[209] = PmA[2] * pOut[134] - PmQf[2] * r_008_1 + 8 * iz2 * (pOut[113] - riz * r_007_1);
   pOut[210] = PmA[1] * pOut[145] - PmQf[1] * r_701_1;
   pOut[211] = PmA[2] * pOut[136] - PmQf[2] * r_530_1;
   pOut[212] = PmA[1] * pOut[147] - PmQf[1] * r_503_1;
   pOut[213] = PmA[2] * pOut[138] - PmQf[2] * r_350_1;
   pOut[214] = PmA[1] * pOut[149] - PmQf[1] * r_323_1 + 2 * iz2 * (pOut[116] - riz * r_313_1);
   pOut[215] = PmA[2] * pOut[140] - PmQf[2] * r_314_1 + 4 * iz2 * (pOut[116] - riz * r_313_1);
   pOut[216] = PmA[2] * pOut[141] - PmQf[2] * r_170_1;
   pOut[217] = PmA[0] * pOut[162] - PmQf[0] * r_053_1;
   pOut[218] = PmA[0] * pOut[163] - PmQf[0] * r_035_1;
   pOut[219] = PmA[0] * pOut[164] - PmQf[0] * r_017_1;
   // 2698 flops, 1652 mops, 2.27kb stack
}

// Expand solid harmonic S[l=0,m](r-B) centered at B into cartesians of degree 0..0 centered at A.
static void OsrrRx0( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, unsigned si, FReal1 const AmB[3] ) AIC_NO_THROW
{
   pOut[0] = pIn[0*si];
}

// Expand solid harmonic S[l=1,m](r-B) centered at B into cartesians of degree 0..1 centered at A.
static void OsrrRx1( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, unsigned si, FReal1 const AmB[3] ) AIC_NO_THROW
{
   pOut[1] = pIn[0*si];
   pOut[2] = pIn[1*si];
   pOut[3] = pIn[2*si];
   pOut[0] =   AmB[0]*pOut[1] +   AmB[1]*pOut[2] +   AmB[2]*pOut[3];
}

// Expand solid harmonic S[l=2,m](r-B) centered at B into cartesians of degree 0..2 centered at A.
static void OsrrRx2( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, unsigned si, FReal1 const AmB[3] ) AIC_NO_THROW
{
   CaTrN1(&pOut[4], &pIn[0], si, 2);

   if (AmB[0] == 0 && AmB[1] == 0 && AmB[2] == 0) {
      for ( unsigned i = 0; i < 4; ++ i )
         pOut[i] = 0;
      return;
   } else {
      pOut[1] = 2*AmB[0]*pOut[4] +   AmB[1]*pOut[7] +   AmB[2]*pOut[8];
      pOut[2] =   AmB[0]*pOut[7] + 2*AmB[1]*pOut[5] +   AmB[2]*pOut[9];
      pOut[3] =   AmB[0]*pOut[8] +   AmB[1]*pOut[9] + 2*AmB[2]*pOut[6];
      double const pf2 = 1.0/2.0;
      pOut[0] = pf2 * (  AmB[0]*pOut[1] +   AmB[1]*pOut[2] +   AmB[2]*pOut[3]);
   }
}

// Expand solid harmonic S[l=3,m](r-B) centered at B into cartesians of degree 0..3 centered at A.
static void OsrrRx3( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, unsigned si, FReal1 const AmB[3] ) AIC_NO_THROW
{
   CaTrN1(&pOut[10], &pIn[0], si, 3);

   if (AmB[0] == 0 && AmB[1] == 0 && AmB[2] == 0) {
      for ( unsigned i = 0; i < 10; ++ i )
         pOut[i] = 0;
      return;
   } else {
      pOut[4] = 3*AmB[0]*pOut[10] +   AmB[1]*pOut[13] +   AmB[2]*pOut[16];
      pOut[5] =   AmB[0]*pOut[11] + 3*AmB[1]*pOut[14] +   AmB[2]*pOut[17];
      pOut[6] =   AmB[0]*pOut[12] +   AmB[1]*pOut[15] + 3*AmB[2]*pOut[18];
      pOut[7] = 2*AmB[0]*pOut[13] + 2*AmB[1]*pOut[11] +   AmB[2]*pOut[19];
      pOut[8] = 2*AmB[0]*pOut[16] +   AmB[1]*pOut[19] + 2*AmB[2]*pOut[12];
      pOut[9] =   AmB[0]*pOut[19] + 2*AmB[1]*pOut[17] + 2*AmB[2]*pOut[15];
      double const pf2 = 1.0/2.0;
      pOut[1] = pf2 * (2*AmB[0]*pOut[4] +   AmB[1]*pOut[7] +   AmB[2]*pOut[8]);
      pOut[2] = pf2 * (  AmB[0]*pOut[7] + 2*AmB[1]*pOut[5] +   AmB[2]*pOut[9]);
      pOut[3] = pf2 * (  AmB[0]*pOut[8] +   AmB[1]*pOut[9] + 2*AmB[2]*pOut[6]);
      double const pf3 = 1.0/3.0;
      pOut[0] = pf3 * (  AmB[0]*pOut[1] +   AmB[1]*pOut[2] +   AmB[2]*pOut[3]);
   }
}

// Expand solid harmonic S[l=4,m](r-B) centered at B into cartesians of degree 0..4 centered at A.
static void OsrrRx4( FReal0 *AIC_RP pOut, FReal0 const *AIC_RP pIn, unsigned si, FReal1 const AmB[3] ) AIC_NO_THROW
{
   CaTrN1(&pOut[20], &pIn[0], si, 4);

   if (AmB[0] == 0 && AmB[1] == 0 && AmB[2] == 0) {
      for ( unsigned i = 0; i < 20; ++ i )
         pOut[i] = 0;
      return;
   } else {
      pOut[10] = 4*AmB[0]*pOut[20] +   AmB[1]*pOut[26] +   AmB[2]*pOut[29];
      pOut[11] = 2*AmB[0]*pOut[21] + 3*AmB[1]*pOut[27] +   AmB[2]*pOut[30];
      pOut[12] = 2*AmB[0]*pOut[22] +   AmB[1]*pOut[28] + 3*AmB[2]*pOut[31];
      pOut[13] = 3*AmB[0]*pOut[26] + 2*AmB[1]*pOut[21] +   AmB[2]*pOut[32];
      pOut[14] =   AmB[0]*pOut[27] + 4*AmB[1]*pOut[23] +   AmB[2]*pOut[33];
      pOut[15] =   AmB[0]*pOut[28] + 2*AmB[1]*pOut[24] + 3*AmB[2]*pOut[34];
      pOut[16] = 3*AmB[0]*pOut[29] +   AmB[1]*pOut[32] + 2*AmB[2]*pOut[22];
      pOut[17] =   AmB[0]*pOut[30] + 3*AmB[1]*pOut[33] + 2*AmB[2]*pOut[24];
      pOut[18] =   AmB[0]*pOut[31] +   AmB[1]*pOut[34] + 4*AmB[2]*pOut[25];
      pOut[19] = 2*AmB[0]*pOut[32] + 2*AmB[1]*pOut[30] + 2*AmB[2]*pOut[28];
      double const pf2 = 1.0/2.0;
      pOut[4] = pf2 * (3*AmB[0]*pOut[10] +   AmB[1]*pOut[13] +   AmB[2]*pOut[16]);
      pOut[5] = pf2 * (  AmB[0]*pOut[11] + 3*AmB[1]*pOut[14] +   AmB[2]*pOut[17]);
      pOut[6] = pf2 * (  AmB[0]*pOut[12] +   AmB[1]*pOut[15] + 3*AmB[2]*pOut[18]);
      pOut[7] = pf2 * (2*AmB[0]*pOut[13] + 2*AmB[1]*pOut[11] +   AmB[2]*pOut[19]);
      pOut[8] = pf2 * (2*AmB[0]*pOut[16] +   AmB[1]*pOut[19] + 2*AmB[2]*pOut[12]);
      pOut[9] = pf2 * (  AmB[0]*pOut[19] + 2*AmB[1]*pOut[17] + 2*AmB[2]*pOut[15]);
      double const pf3 = 1.0/3.0;
      pOut[1] = pf3 * (2*AmB[0]*pOut[4] +   AmB[1]*pOut[7] +   AmB[2]*pOut[8]);
      pOut[2] = pf3 * (  AmB[0]*pOut[7] + 2*AmB[1]*pOut[5] +   AmB[2]*pOut[9]);
      pOut[3] = pf3 * (  AmB[0]*pOut[8] +   AmB[1]*pOut[9] + 2*AmB[2]*pOut[6]);
      double const pf4 = 1.0/4.0;
      pOut[0] = pf4 * (  AmB[0]*pOut[1] +   AmB[1]*pOut[2] +   AmB[2]*pOut[3]);
   }
}


FOsrrFnA const
   OsrrA[MaxLab+1] = { &OsrrA0, &OsrrA1, &OsrrA2, &OsrrA3, &OsrrA4, &OsrrA5, &OsrrA6, &OsrrA7, &OsrrA8, &OsrrA9 };

FOsrrFnC const
   OsrrC[MaxL1a+1][MaxL1b+1] = {
      {&OsrrC00, 0, 0, 0, 0},
      {&OsrrC10, &OsrrC11, 0, 0, 0},
      {&OsrrC20, &OsrrC21, &OsrrC22, 0, 0},
      {&OsrrC30, &OsrrC31, &OsrrC32, &OsrrC33, 0},
      {&OsrrC40, &OsrrC41, &OsrrC42, &OsrrC43, &OsrrC44}
};

// nCartY(0) x nCartX(4) = 1 x 35 array of cartesian indices into nCartX(4) (35 entries)
unsigned short iCartProd_la0[35] = {
     0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
    32,  33,  34
};


// nCartY(1) x nCartX(4) = 3 x 35 array of cartesian indices into nCartX(5) (56 entries)
unsigned short iCartProd_la1[105] = {
     1,   2,   3,   4,   7,   8,   7,   5,   9,   8,   9,   6,  10,  13,  16,  11,
    14,  17,  12,  15,  18,  13,  11,  19,  16,  19,  12,  19,  17,  15,  20,  26,
    29,  21,  27,  30,  22,  28,  31,  26,  21,  32,  27,  23,  33,  28,  24,  34,
    29,  32,  22,  30,  33,  24,  31,  34,  25,  32,  30,  28,  35,  41,  47,  36,
    42,  48,  37,  43,  49,  38,  44,  50,  39,  45,  51,  40,  46,  52,  41,  36,
    53,  42,  38,  54,  43,  39,  55,  47,  53,  37,  48,  54,  39,  49,  55,  40,
    53,  48,  43,  54,  50,  45,  55,  51,  46
};


// nCartY(2) x nCartX(4) = 6 x 35 array of cartesian indices into nCartX(6) (84 entries)
unsigned short iCartProd_la2[210] = {
     4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  16,  19,  13,  14,  15,  11,
    19,  17,  16,  17,  18,  19,  12,  15,  20,  21,  22,  26,  29,  32,  21,  23,
    24,  27,  30,  33,  22,  24,  25,  28,  31,  34,  26,  27,  28,  21,  32,  30,
    29,  30,  31,  32,  22,  28,  32,  33,  34,  30,  28,  24,  35,  36,  37,  41,
    47,  53,  36,  38,  39,  42,  48,  54,  37,  39,  40,  43,  49,  55,  41,  42,
    43,  36,  53,  48,  42,  44,  45,  38,  54,  50,  43,  45,  46,  39,  55,  51,
    47,  48,  49,  53,  37,  43,  48,  50,  51,  54,  39,  45,  49,  51,  52,  55,
    40,  46,  53,  54,  55,  48,  43,  39,  56,  57,  58,  66,  72,  78,  57,  59,
    60,  67,  73,  79,  58,  60,  61,  68,  74,  80,  59,  62,  63,  69,  75,  81,
    60,  63,  64,  70,  76,  82,  61,  64,  65,  71,  77,  83,  66,  67,  68,  57,
    78,  73,  67,  69,  70,  59,  79,  75,  68,  70,  71,  60,  80,  76,  72,  73,
    74,  78,  58,  68,  73,  75,  76,  79,  60,  70,  74,  76,  77,  80,  61,  71,
    78,  79,  80,  73,  68,  60,  79,  81,  82,  75,  70,  63,  80,  82,  83,  76,
    71,  64
};


// nCartY(3) x nCartX(4) = 10 x 35 array of cartesian indices into nCartX(7) (120 entries)
unsigned short iCartProd_la3[350] = {
    10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  26,  27,  28,
    29,  30,  31,  32,  26,  27,  28,  21,  23,  24,  32,  33,  34,  30,  29,  30,
    31,  32,  33,  34,  22,  24,  25,  28,  35,  36,  37,  41,  42,  43,  47,  48,
    49,  53,  36,  38,  39,  42,  44,  45,  48,  50,  51,  54,  37,  39,  40,  43,
    45,  46,  49,  51,  52,  55,  41,  42,  43,  36,  38,  39,  53,  54,  55,  48,
    47,  48,  49,  53,  54,  55,  37,  39,  40,  43,  53,  54,  55,  48,  50,  51,
    43,  45,  46,  39,  56,  57,  58,  66,  67,  68,  72,  73,  74,  78,  57,  59,
    60,  67,  69,  70,  73,  75,  76,  79,  58,  60,  61,  68,  70,  71,  74,  76,
    77,  80,  66,  67,  68,  57,  59,  60,  78,  79,  80,  73,  67,  69,  70,  59,
    62,  63,  79,  81,  82,  75,  68,  70,  71,  60,  63,  64,  80,  82,  83,  76,
    72,  73,  74,  78,  79,  80,  58,  60,  61,  68,  73,  75,  76,  79,  81,  82,
    60,  63,  64,  70,  74,  76,  77,  80,  82,  83,  61,  64,  65,  71,  78,  79,
    80,  73,  75,  76,  68,  70,  71,  60,  84,  85,  86,  94,  95,  96, 104, 105,
   106, 114,  85,  87,  88,  95,  97,  98, 105, 107, 108, 115,  86,  88,  89,  96,
    98,  99, 106, 108, 109, 116,  87,  90,  91,  97, 100, 101, 107, 110, 111, 117,
    88,  91,  92,  98, 101, 102, 108, 111, 112, 118,  89,  92,  93,  99, 102, 103,
   109, 112, 113, 119,  94,  95,  96,  85,  87,  88, 114, 115, 116, 105,  95,  97,
    98,  87,  90,  91, 115, 117, 118, 107,  96,  98,  99,  88,  91,  92, 116, 118,
   119, 108, 104, 105, 106, 114, 115, 116,  86,  88,  89,  96, 105, 107, 108, 115,
   117, 118,  88,  91,  92,  98, 106, 108, 109, 116, 118, 119,  89,  92,  93,  99,
   114, 115, 116, 105, 107, 108,  96,  98,  99,  88, 115, 117, 118, 107, 110, 111,
    98, 101, 102,  91, 116, 118, 119, 108, 111, 112,  99, 102, 103,  92
};


// nCartY(4) x nCartX(4) = 15 x 35 array of cartesian indices into nCartX(8) (165 entries)
unsigned short iCartProd_la4[525] = {
    20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,
    36,  37,  38,  39,  40,  41,  42,  43,  47,  48,  49,  53,  54,  55,  41,  42,
    43,  44,  45,  46,  36,  38,  39,  53,  54,  55,  48,  50,  51,  47,  48,  49,
    50,  51,  52,  53,  54,  55,  37,  39,  40,  43,  45,  46,  56,  57,  58,  59,
    60,  61,  66,  67,  68,  72,  73,  74,  78,  79,  80,  57,  59,  60,  62,  63,
    64,  67,  69,  70,  73,  75,  76,  79,  81,  82,  58,  60,  61,  63,  64,  65,
    68,  70,  71,  74,  76,  77,  80,  82,  83,  66,  67,  68,  69,  70,  71,  57,
    59,  60,  78,  79,  80,  73,  75,  76,  72,  73,  74,  75,  76,  77,  78,  79,
    80,  58,  60,  61,  68,  70,  71,  78,  79,  80,  81,  82,  83,  73,  75,  76,
    68,  70,  71,  60,  63,  64,  84,  85,  86,  87,  88,  89,  94,  95,  96, 104,
   105, 106, 114, 115, 116,  85,  87,  88,  90,  91,  92,  95,  97,  98, 105, 107,
   108, 115, 117, 118,  86,  88,  89,  91,  92,  93,  96,  98,  99, 106, 108, 109,
   116, 118, 119,  94,  95,  96,  97,  98,  99,  85,  87,  88, 114, 115, 116, 105,
   107, 108,  95,  97,  98, 100, 101, 102,  87,  90,  91, 115, 117, 118, 107, 110,
   111,  96,  98,  99, 101, 102, 103,  88,  91,  92, 116, 118, 119, 108, 111, 112,
   104, 105, 106, 107, 108, 109, 114, 115, 116,  86,  88,  89,  96,  98,  99, 105,
   107, 108, 110, 111, 112, 115, 117, 118,  88,  91,  92,  98, 101, 102, 106, 108,
   109, 111, 112, 113, 116, 118, 119,  89,  92,  93,  99, 102, 103, 114, 115, 116,
   117, 118, 119, 105, 107, 108,  96,  98,  99,  88,  91,  92, 120, 121, 122, 123,
   124, 125, 135, 136, 137, 145, 146, 147, 155, 156, 157, 121, 123, 124, 126, 127,
   128, 136, 138, 139, 146, 148, 149, 156, 158, 159, 122, 124, 125, 127, 128, 129,
   137, 139, 140, 147, 149, 150, 157, 159, 160, 123, 126, 127, 130, 131, 132, 138,
   141, 142, 148, 151, 152, 158, 161, 162, 124, 127, 128, 131, 132, 133, 139, 142,
   143, 149, 152, 153, 159, 162, 163, 125, 128, 129, 132, 133, 134, 140, 143, 144,
   150, 153, 154, 160, 163, 164, 135, 136, 137, 138, 139, 140, 121, 123, 124, 155,
   156, 157, 146, 148, 149, 136, 138, 139, 141, 142, 143, 123, 126, 127, 156, 158,
   159, 148, 151, 152, 137, 139, 140, 142, 143, 144, 124, 127, 128, 157, 159, 160,
   149, 152, 153, 145, 146, 147, 148, 149, 150, 155, 156, 157, 122, 124, 125, 137,
   139, 140, 146, 148, 149, 151, 152, 153, 156, 158, 159, 124, 127, 128, 139, 142,
   143, 147, 149, 150, 152, 153, 154, 157, 159, 160, 125, 128, 129, 140, 143, 144,
   155, 156, 157, 158, 159, 160, 146, 148, 149, 137, 139, 140, 124, 127, 128, 156,
   158, 159, 161, 162, 163, 148, 151, 152, 139, 142, 143, 127, 131, 132, 157, 159,
   160, 162, 163, 164, 149, 152, 153, 140, 143, 144, 128, 132, 133
};


FOsrrFnR OsrrRx[5] = {
   OsrrRx0, OsrrRx1, OsrrRx2, OsrrRx3, OsrrRx4
};

unsigned const short * piCartProd[5] = {
   &iCartProd_la0[0], &iCartProd_la1[0], &iCartProd_la2[0], &iCartProd_la3[0], &iCartProd_la4[0]

};

FOsrrFnC_dA const
   OsrrC_dA[MaxL1a_Der1+1][MaxL1b_Der1+1] = {

};

FOsrrFnC_dB const
   OsrrC_dB[MaxL1a_Der1+1][MaxL1b_Der1+1] = {

};


} // namespace Osrr

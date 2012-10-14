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

/* AicOsrr.h v20121011 EST [charge, Gerald Knizia] */
#ifndef _AIC_OSRR_H
#define _AIC_OSRR_H

#include "AicDefs.h"

namespace aic{

typedef double
   FReal0;      // data input & output
typedef double
   FReal1;      // parameter input
typedef double
   FReal2;      // internal recursion intermediate

const unsigned
   MaxL1a = 4, // largest L at (a.|.)
   MaxL1b = 4, // largest L at (.b|.)
   MaxL1c = 6, // largest L at (..|c)
   MinL1 = 4,  // max(L1a,L1b,L1c) (i.e., largest L supported everywhere)
   MaxL1 = 6,  // max(L1a,L1b,L1c) (i.e., largest L supported anywhere)
   MaxLab = 9, // largest total L of (ab| for OsrrA and OsrrB
   MaxL1a_Der1 = -1, // largest L at (a.|.) for which first derivatives are available
   MaxL1b_Der1 = -1, // largest L at (.b|.) for which first derivatives are available
   MaxL1c_Der1 = MaxL1c - 1; // largest L at (..|c) for which first derivatives are available

struct FOsrrParamsB {
   FReal1
      InvEtaABC2,
      PmQf[3];
   int
      AngMomAB_Min, AngMomAB_Max;
};
void OsrrKrnB_3cen( FReal2 *AIC_RP pOut, FReal2 const *AIC_RP pIn, unsigned iDir, FOsrrParamsB const &P );

typedef void
   (*FOsrrFnA)( FReal0 *AIC_RP, FReal0 const *AIC_RP, FReal1 const [3], FReal1 const [3], double, double ),
   (*FOsrrFnC)( FReal0 *AIC_RP, FReal0 const *AIC_RP, FReal1 const [3], unsigned sa, unsigned sb ),
   (*FOsrrFnR)( FReal0 *AIC_RP, FReal0 const *AIC_RP, unsigned si, FReal1 const [3] ),
   (*FOsrrFnC_dA)( FReal0 *, FReal0 const *AIC_RP, FReal0 const *AIC_RP, FReal1 const [3], unsigned, unsigned, unsigned ),
   (*FOsrrFnC_dB)( FReal0 *, FReal0 const *AIC_RP, FReal0 const *AIC_RP, FReal1 const [3], unsigned, unsigned, unsigned );
extern FOsrrFnA const
   OsrrA[MaxLab+1];
extern FOsrrFnC const
   OsrrC[MaxL1a+1][MaxL1b+1];
extern FOsrrFnC_dA const
   OsrrC_dA[MaxL1a_Der1+1][MaxL1b_Der1+1];
extern FOsrrFnC_dB const
   OsrrC_dB[MaxL1a_Der1+1][MaxL1b_Der1+1];
extern FOsrrFnR
   OsrrRx[MaxL1b+1];
extern unsigned const short
   *piCartProd[MaxL1b+1];


} // namespace aic

#endif // _AIC_OSRR_H


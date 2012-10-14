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

/* AicSolidHarmonics.h v20121011 EST [charge, Gerald Knizia] */
#ifndef _AIC_SOLID_HARMONICS_2_H
#define _AIC_SOLID_HARMONICS_2_H

#include "AicDefs.h"
#ifndef assert
   #include <assert.h>
#endif

// Spherical component order (Molpro/Molcas):
// L=0: s0
// L=1: p1+ p1- p0
// L=2: d0  d2- d1+ d2+ d1-
// L=3: f1+ f1- f0  f3+ f2- f3- f2+
// L=4: g0  g2- g1+ g4+ g1- g2+ g4- g3+ g3-
// L=5: h1+ h1- h2+ h3+ h4- h3- h4+ h5- h0  h5+ h2-
// L=6: i6+ i2- i5+ i4+ i5- i2+ i6- i3+ i4- i0  i3- i1- i1+

namespace aic {
   // SlmSymSig[l*l + si(m)] is the bit pattern denoting the symmetry signature
   // of S(l,m). Axes: bit0: X, bit1: Y, bit2: Z.  Bit #i: If even, function is
   // symmetric wrt mirroring at axis #i, if odd, it is antisymmetric.
   extern unsigned char SlmSymSig[49];

   // transform one cartesian coordinate to solid harmonics.
   void ShTrN( double *AIC_RP pOut, unsigned StrideOut, double const *AIC_RP pIn, unsigned StrideIn, unsigned L, unsigned Count );

   // transform one solid harmonic coordinate to cartesians.
   void CaTrN1( double *AIC_RP pOut, double const *AIC_RP pIn, unsigned StrideIn, unsigned L );

   // evaluate all S^l_m(R) for l = 0 .. L. Results stored in Out[iSlmA(l,m)].
   void EvalSlmX( double *AIC_RP pOut, double const *AIC_RP R, int L );

   // evaluate all S^l_m(R) and d/dx_i S^l_m(R) for l = 0 .. L.
   // Results stored in Out[4*iSlmA(l,m) + iComp]. iComp = 0: value, =1,2,3: d/dx_i.
   void EvalSlmX_Grd1( double *AIC_RP Out, double const *AIC_RP R, int L );

   // evaluate all S^l_m(R), d/dx_i S^l_m(R), and [d^2/(dx_i dx_j)] S^l_m(R) for l = 0 .. L.
   // Results stored in Out[10*iSlmA(l,m) + iComp]. iComp = 0: value, =1,2,3: d/dx_i; 4..10: D(xx),D(yy),D(zz),D(xy),D(xz),D(zz)
   void EvalSlmX_Grd2( double *AIC_RP Out, double const *AIC_RP R, int L );
}

#endif // _AIC_SOLID_HARMONICS_2_H

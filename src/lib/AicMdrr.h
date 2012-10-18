#ifndef _AIC_MDRR_H
#define _AIC_MDRR_H

#include "AicDefs.h"

namespace aic{

typedef double
   FReal0;      // data input & output
typedef double
   FReal1;      // parameter input
typedef double
   FReal2;      // internal recursion intermediate

// McMurchie-Davidson-recurrence relations for given Ls:
//    [r]^m = R_i [r - 1_i]^(m+1) - (r_i - 1)[r - 2_i]^(m+1)
// where r are 3-vectors, 1_i and 2_i are unit/scaled vectors in the
// i-direction and [r - 2_i] is considered zero if it does not exist.
//
// Each function gets [0]^m for m in [0..L] (inclusive) and R
// as input, and outputs [r]^0 where r are the angular components of L
// (i.e., the cartesian monomials x^nx y^ny z^nz with nx+ny+nz=L
// in the order defined by AngularComps(L)).
//
// Effectively the MDRR (and therefore these functions) calculate
//   D^r f(t)
// from f^[m](t) = (D/Dt)^m f(t) and R, where t = R^2, D^r means
// \prod_i (D/D{R_i})^{r_x}, and f is some arbitrary scalar function.
// (of which you supply the m'th derivatives with respect to t as [0]^m).
//
// Note: In AIC, MDRR is used for two-center integrals, where the shell-
// MDRR implemented here (i.e., only [R] with sum(r)==L occur) allows for
// calculating two-electron two-center integrals directly from *contracted*
// Gm(rho,T) kernel functions, additionally to beeing exceedingly
// numerically stable.

const unsigned
   MaxL1_MDRR = 6,
   MaxL_MDRR = 14;

typedef void
   (*FMdrrFn)(FReal0 *AIC_RP, FReal0 const *AIC_RP, FReal1 const [3]);
extern FMdrrFn const
   ShellMdrr[MaxL_MDRR+1];

// [r1][r2]: pointers to tables nCartY(r1) x nCartY(r2) indexing into nCartY(r1+r2)
extern unsigned short *iCartPxx[15][7];

// Scatter 2-shell result: write a (2la + 1) x (2lc + 1) block of data at pIn
// into memory at pOut[ia*sa + ic*sc], overwriting the previous result.
void Scatter2e2c_O(FReal0 *AIC_RP pOut, unsigned sa, unsigned sc, FReal2 const *AIC_RP pIn, unsigned la, unsigned lc);
// Scatter 2-shell result: write a (2la + 1) x (2lc + 1) block of data at pIn
// into memory at pOut[ia*sa + ic*sc], adding to the previous result.
void Scatter2e2c_A(FReal0 *AIC_RP pOut, unsigned sa, unsigned sc, FReal2 const *AIC_RP pIn, unsigned la, unsigned lc);

typedef void
   (*FScatter2e2cFn)(FReal0 *, unsigned, unsigned, FReal2 const *, unsigned, unsigned);
extern FScatter2e2cFn const
   // [0] overwrites, [1] adds
   Scatter2e2c[2];


} // namespace aic

#endif // _AIC_MDRR_H


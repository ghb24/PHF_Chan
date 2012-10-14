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

#include <cmath>
#ifndef AIC_LOCAL_INCLUDES
   #define MACHINES_H_DEFINES_ONLY
   #include "util/machines.h" // for FORTINT and FORT_Extern
#else
   #include "CxFortranInt.h"
#endif
#include <boost/intrusive_ptr.hpp>
#include <stdint.h> // for uint32_t vs uint64_t (used for component presence flags)

#include "AicCommon.h"
#include "AicDefs.h"
#include "AicKernels.h"
#include "AicBoysFn.h"

namespace aic {

const uint
   MaxJ = 40;

uint FIntegralKernel::GetNumComponents() const
{
   return 1;
}


void FCoulombKernel::EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const
{
   BoysFn( pOut, MaxM, T, (2*M_PI)*Prefactor/rho );
}

static void EvalErfCoulombGm(double *pOut, double rho, double T, uint MaxM, double Prefactor, double Omega)
{
   if ( MaxM != 0 )
      assert_rt("!FIXME: implement erf screened coulomb with higher angular momenta.");
   // dx.doi.org/10.1039/b605188j eq.52, 2nd term.
   double f = Omega/std::sqrt(Omega*Omega + rho);
   BoysFn( pOut, MaxM, f*f*T, f*Prefactor );
};

void FErfCoulombKernel::EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const
{
   EvalErfCoulombGm(pOut, rho, T, MaxM, (2*M_PI)*Prefactor/rho, m_Omega);
};

void FErfcCoulombKernel::EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const
{
   // form regular coulomb kernel and subtract the long-range kernel.
   double
      Fm[MaxJ],
      PrefactorFm = (2 * M_PI) * Prefactor/rho;
   BoysFn( &Fm[0], MaxM, T, PrefactorFm );
   EvalErfCoulombGm(pOut, rho, T, MaxM, PrefactorFm, m_Omega);
   for ( uint i = 0; i <= MaxM; ++ i )
      pOut[i] -= Fm[i];
};



void FOverlapKernel::EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const
{
//    double
//       f = M_PI/rho;
//    printf("T=%f  rho=%f",T/rho,rho);
//    pOut[0] = Prefactor * f * std::sqrt(f) * std::exp(-T);
//    for ( uint i = 1; i <= MaxM; ++i )
//       pOut[i] = rho * pOut[i-1];

   pOut[0] = Prefactor * std::exp(-T);
   for ( uint i = 1; i <= MaxM; ++i )
      pOut[i] = pOut[i-1]; // (-d/dT)^n G(0). Note that T = rho |P-Q|^2 includes the rho!
}


// extern "C" {
// void gmgausskernel_(double const *Omega, double const *Coeff, long const &nExp, double const &rho, long const &MaxM, double const &T, double *Out, long const &iStride, double const &Prefactor);
// }
//
// // F_{12}(r_{12}) = \sum_i c_i exp(-zeta_i r_{12}^2)
// void FGaussKernel::EvalGm( double *pOut, double Rho, double T, uint MaxM, double Prefactor ) const
// {
//    gmgausskernel_(&m_pGaussExp->Exponents[0], &m_pGaussExp->Coeffs[0],
//       m_pGaussExp->Coeffs.size(), Rho, MaxM, T, pOut, 1, Prefactor);
//
//
//
// /*   assert( m_pGaussExp->Exponents.size() == m_pGaussExp->Coeffs.size() );
//    // based on:
//    //   [1]: PCCP 10 3390 (2008).
//    double const
//       *AIC_RP Omega = &m_pGaussExp->Exponents[0],
//       *AIC_RP Coeff = &m_pGaussExp->Coeffs[0];
//
//    for ( uint m = 0; m <= MaxM; ++ m )
//       pOut[m] = 0;
//
//    for ( uint iExp = 0; iExp < m_pGaussExp->Exponents.size(); ++ iExp ) {
//       double
//          f = 1.0/(Rho + Omega[iExp]),
//          RhoTilde = f * Omega[iExp],
//          ExpArg = RhoTilde * T;
// //       if ( ExpArg > 40 )
// //          continue; // exp(-40) = ~5e-18
//
//       // [1], eq. (24):
//       //    G_m(rho,T) = \sum_i c[i] RhoTilde[i]^m (pi/(rho+omega_i))^(3/2) exp(-RhoTilde[i] * T).
//       double
//          u = M_PI * f,
//          Pref = Prefactor * Coeff[iExp] * u * std::sqrt(u) * std::exp(-ExpArg),
//          RhoTildeM = 1.0;
//       for ( uint m = 0; m <= MaxM; ++ m ) {
//          pOut[m] += Pref * RhoTildeM;
//          RhoTildeM *= RhoTilde;
//       }
//    }*/
// };

// // F_{12}(r_{12}) = \sum_i c_i exp(-zeta_i r_{12}^2)
// void FGaussKernel::EvalGm( double *pOut, double Rho, double T, uint MaxM, double Prefactor ) const
// {
//    assert( m_pGaussExp->Exponents.size() == m_pGaussExp->Coeffs.size() );
//    // based on:
//    //   [1]: PCCP 10 3390 (2008).
//    double const
//       *AIC_RP Omega = &m_pGaussExp->Exponents[0],
//       *AIC_RP Coeff = &m_pGaussExp->Coeffs[0];
//
//    const uint
//       N = 25;
//    double
//       Pref[N],
//       RhoTilde[N],
//       RhoTildeM[N];
//    uint
//       nExp = m_pGaussExp->Exponents.size();
//
//    for ( uint iExp = 0; iExp < nExp; ++ iExp ) {
//       double
//          f = 1.0/(Rho + Omega[iExp]);
//       double
//          u = M_PI * f;
//       RhoTilde[iExp] = f * Omega[iExp];
//       Pref[iExp] = Prefactor * Coeff[iExp] * u * std::sqrt(u) * std::exp(-RhoTilde[iExp] * T);
//    }
//    for ( uint iExp = 0; iExp < nExp; ++ iExp )
//       RhoTildeM[iExp] = 1.0;
//    for ( uint m = 0; m <= MaxM; ++ m ) {
//       double o = 0;
//       for ( uint iExp = 0; iExp < nExp; ++ iExp )
//          o += Pref[iExp] * RhoTildeM[iExp];
//       for ( uint iExp = 0; iExp < nExp; ++ iExp )
//          RhoTildeM[iExp] *= RhoTilde[iExp];
//       pOut[m] = o;
//    }
// };


void FGaussKernel::EvalGm( double *pOut, double Rho, double T, uint MaxM, double Prefactor ) const
{
   assert( m_pGaussExp->Exponents.size() * m_pGaussExp->nCo == m_pGaussExp->Coeffs.size() );
   // based on:
   //   [1]: PCCP 10 3390 (2008).
   double const
      *AIC_RP Omega = &m_pGaussExp->Exponents[0],
      *AIC_RP Coeff = &m_pGaussExp->Coeffs[0];
   uint
      nExp = m_pGaussExp->Exponents.size(),
      nCo = m_pGaussExp->nCo, // number of contractions
      nSt = MaxM + 1; // output component stride

   for ( uint m = 0; m < nSt * nCo; ++ m )
      pOut[m] = 0;

   for ( uint iExp = 0; iExp < nExp; ++ iExp ) {
      double
         f = 1.0/(Rho + Omega[iExp]),
         RhoTilde = f * Omega[iExp];

      // [1], eq. (24):
      //    G_m(rho,T) = \sum_i c[i] RhoTilde[i]^m (pi/(rho+omega_i))^(3/2) exp(-RhoTilde[i] * T).
      double
         u = M_PI * f,
         PrefBase = Prefactor * u * std::sqrt(u) * std::exp(-RhoTilde * T);
      for ( uint iCo = 0; iCo < nCo; ++ iCo ) {
         double
            Pref = Coeff[iCo + nCo*iExp] * PrefBase,
            RhoTildeM = 1.0;
         for ( uint m = 0; m <= MaxM; ++ m ) {
            pOut[m + iCo*nSt] += Pref * RhoTildeM;
            RhoTildeM *= RhoTilde;
         }
      }
   }
}


#if 0
// we might want to tabulate this.
static std::size_t NoverK(std::size_t n_, std::size_t k_)
{
   // FIXME: FUNCTION IS HIGHLY DOUBIOUS!! INVESTIGATE THIS FIRST
   //        WHEN TF INTEGRALS DONT WORK!
   // should probably use boost::binomial_coefficient. But current
   // boost version installed here doesn't have it yet.
   // Also note the high chance of overflow..
   // FIXME: add checking code for high angular momenta here!!
   //
   if ( k_ == 0 || k_ == n_ )
      return 1;
   if ( k_ == 1 || k_ == n_-1 )
      return n_;
   std::size_t
      n = n_ + 1,
      k = k_ + 1;
//    if ( k > n )
//       return 0;
   std::size_t
      num = 1,
      denom = 1;
   for ( std::size_t i = 1; i < k; ++ i )
      denom *= i;
   for ( std::size_t i = n-k+1; i < n; ++ i )
      num *= i;
   return num/denom;
}

#else

// address with comb(N,k) = Tab[N*(N-1)/2 + k]
const uint MAX_NoverK = 20; // 6*3+2, should be sufficient for up to (ii|i)^2nd derivative shells.
static int s_NoverK_Tab[191] = {
    1, 1, 2, 1, 3, 3, 1, 4, 6, 4,
    1, 5, 10, 10, 5, 1, 6, 15, 20, 15,
    6, 1, 7, 21, 35, 35, 21, 7, 1, 8,
    28, 56, 70, 56, 28, 8, 1, 9, 36, 84,
    126, 126, 84, 36, 9, 1, 10, 45, 120, 210,
    252, 210, 120, 45, 10, 1, 11, 55, 165, 330,
    462, 462, 330, 165, 55, 11, 1, 12, 66, 220,
    495, 792, 924, 792, 495, 220, 66, 12, 1, 13,
    78, 286, 715, 1287, 1716, 1716, 1287, 715, 286, 78,
    13, 1, 14, 91, 364, 1001, 2002, 3003, 3432, 3003,
    2002, 1001, 364, 91, 14, 1, 15, 105, 455, 1365,
    3003, 5005, 6435, 6435, 5005, 3003, 1365, 455, 105, 15,
    1, 16, 120, 560, 1820, 4368, 8008, 11440, 12870, 11440,
    8008, 4368, 1820, 560, 120, 16, 1, 17, 136, 680,
    2380, 6188, 12376, 19448, 24310, 24310, 19448, 12376, 6188, 2380,
    680, 136, 17, 1, 18, 153, 816, 3060, 8568, 18564,
    31824, 43758, 48620, 43758, 31824, 18564, 8568, 3060, 816, 153,
    18, 1, 19, 171, 969, 3876, 11628, 27132, 50388, 75582,
    92378, 92378, 75582, 50388, 27132, 11628, 3876, 969, 171, 19, 1
}; // 0.76 kb.

inline int NoverK(uint n_, uint k_){
   assert_rt( n_ <= MAX_NoverK );
   assert( k_ <= n_ );
   return s_NoverK_Tab[n_*(n_-1)/2 + k_];
}
#endif


// F_{12}/r
void FGaussCoulombKernel::EvalGm( double *pOut, double Rho, double T, uint MaxM, double Prefactor ) const
{
   assert( m_pGaussExp->Exponents.size() * m_pGaussExp->nCo == m_pGaussExp->Coeffs.size() );
   // based on:
   //   [1]: PCCP 10 3390 (2008).
   double const
      *Omega = &m_pGaussExp->Exponents[0],
      *Coeff = &m_pGaussExp->Coeffs[0];
   uint
      nExp = m_pGaussExp->Exponents.size(),
      nCo = m_pGaussExp->nCo, // number of contractions
      nSt = MaxM + 1; // output component stride
   double
      Fm[MaxJ];
   assert_rt( MaxM < MaxJ );

   for ( uint m = 0; m < nSt * nCo; ++ m )
      pOut[m] = 0;

   for ( uint iExp = 0; iExp < nExp; ++ iExp ) {
      double
         f = 1.0/(Rho + Omega[iExp]),
         RhoTilde = f * Omega[iExp],
         RhoHat = f * Rho,
         InvRhoHat = 1.0/RhoHat;

      BoysFn( &Fm[0], MaxM, RhoHat * T, 1.0 );

      double
         Pref = Prefactor * (2.0 * M_PI) * f * std::exp(-RhoTilde * T),
         RhoHatM = 1.0; // RhoHat^m
      for ( uint m = 0; m <= MaxM; ++ m ) {
         double
            r = 0,
            t = RhoHatM,
            u = 1.0;
         // this is supposed to form
         //   r = \sum_k NoverK(m,k) * RhoTilde**k * RhoHat**(m-k) * Fm[m-k];
         //   (second line of [1] eq. (30))
         for ( uint k = 0; k <= m; ++ k ) {
            r += NoverK(m,k) * u * t * Fm[m-k];
            u *= RhoTilde;
            t *= InvRhoHat;
         }

         RhoHatM *= RhoHat;
         for ( uint iCo = 0; iCo < nCo; ++ iCo )
            pOut[m + nSt*iCo] += Coeff[iCo + nCo*iExp] * Pref * r;
      }
   }
}



// [grad F_{12}]^2
void FGaussKineticKernel::EvalGm( double *pOut, double Rho, double T, uint MaxM, double Prefactor ) const
{
   assert( m_pGaussExp->Exponents.size() * m_pGaussExp->nCo == m_pGaussExp->Coeffs.size() );
   // based on:
   //   [1]: PCCP 10 3390 (2008).
   double const
      *Omega = &m_pGaussExp->Exponents[0],
      *Coeff = &m_pGaussExp->Coeffs[0];
   uint
      nExp = m_pGaussExp->Exponents.size(),
      nCo = m_pGaussExp->nCo, // number of contractions
      nSt = MaxM + 1; // output component stride

   for ( uint m = 0; m < nSt * nCo; ++ m )
      pOut[m] = 0;

   // [1] eq. (32).
   for ( uint j = 0; j < nExp; ++ j ) {
      for ( uint i = j; i < nExp; ++ i ) {
         double
            wi = Omega[i],
            wj = Omega[j],
            wij = wi + wj,
            f = 1.0/(Rho + wij),
            RhoTilde = f * wij,
            RhoHat = f * Rho,
            Pref = Prefactor * 4.0 * wi * wj *
                   f*f*sqrt((M_PI*M_PI*M_PI)*f) *
                   exp(-RhoTilde * T);
         if ( i != j )
            Pref *= 2.0;
         double
            t = RhoTilde * (1.5 + RhoHat * T),
            PrefRhoTildeM1_ = Pref/RhoTilde;
         for ( uint iCo = 0; iCo < nCo; ++ iCo ) {
            double
               PrefRhoTildeM1 = Coeff[i*nCo + iCo] * Coeff[j*nCo + iCo] * PrefRhoTildeM1_;
            for ( uint m = 0; m <= MaxM; ++ m ) {
               pOut[m + nSt*iCo] += PrefRhoTildeM1 * (t - m * RhoHat);
               PrefRhoTildeM1 *= RhoTilde;
            }
         }
      }
   }
}


FIntegralKernel::~FIntegralKernel() {}
FCoulombKernel::~FCoulombKernel() {}
FOverlapKernel::~FOverlapKernel() {}
FGaussKernel::~FGaussKernel() {}
FGaussCoulombKernel::~FGaussCoulombKernel() {}
FGaussKineticKernel::~FGaussKineticKernel() {}
FErfCoulombKernel::~FErfCoulombKernel() {}
FErfcCoulombKernel::~FErfcCoulombKernel() {}


} // namespace aic

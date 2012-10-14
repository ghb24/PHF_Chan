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

#ifndef _AIC_INTKERNELS_H
#define _AIC_INTKERNELS_H

#include <boost/intrusive_ptr.hpp>
#include <vector>

#include "AicCommon.h"

namespace aic {

void InitBoysFnTable(); // AicBoysFn.cpp. Only actually does something when called for the first time.

struct FIntegralKernel : public FIntrusivePtrDest
{
   virtual char const *GetName() const = 0;

   // obtain number of related components of the integral (default: 1).
   // Used, for example, when one set of primitive Gauss geminals is
   // contracted into multiple different total functions.
   virtual uint GetNumComponents() const;
   // Calculates
   //   pOut[m + (1+MaxM)*iComp] = G(m) * Prefactor
   // for m = 0 .. MaxM (inclusive MaxM) and iComp = 0..nComps-1 (kernel components)
   virtual void EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const = 0;
   virtual ~FIntegralKernel();
};

typedef boost::intrusive_ptr<FIntegralKernel>
   FKernelPtr;

struct FCoulombKernel : public FIntegralKernel
{
   char const *GetName() const { return "r^-1"; }
   void EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const; // override
   FCoulombKernel() { InitBoysFnTable(); }
   ~FCoulombKernel(); // override
};

struct FOverlapKernel : public FIntegralKernel
{
   char const *GetName() const { return "id"; } // override
   void EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const; // override
   ~FOverlapKernel(); // override
};


// describes one or more general geminals expanded in terms of a
// common set primitive Gauss geminals.
struct FGaussGeminalExp
{
   typedef std::vector<double>
      FScalarArray;
   uint
      // number of contracted geminals
      nCo;
   FScalarArray
      Exponents,
      // nCo x nExp contraction matrix.
      Coeffs;

   inline double const &Coeff( uint iExp, uint iCo ){
      assert(iCo < nCo && Coeffs.size() == nCo * Exponents.size());
      return Coeffs[iCo + nCo * iExp];
   }

   void Add(double Exponent, double Coeff) {
      Exponents.push_back(Exponent);
      Coeffs.push_back(Coeff);
   }
   void Clear(uint nReserve = 0, uint nGem = 0) {
      Exponents.clear();
      Exponents.reserve(nReserve);
      Coeffs.clear();
      Coeffs.reserve(nReserve * nGem);
      nCo = 1;
   };
   void SetLength(uint nLength) {
      Exponents.resize(nLength);
      Coeffs.resize(nLength);
      nCo = 1;
   }
   void Set(double *pExponents, double *pCoeffs, uint nLength){
      Clear(nLength);
      for ( uint i = 0; i < nLength; ++ i )
         Add( pExponents[i], pCoeffs[i] );
      nCo = 1;
   }
   void Set(double *pExponents, double *pCoeffs, uint nExp, uint nGem){
      Clear(nExp, nGem);
      nCo = nGem;
      for ( uint iExp = 0; iExp < nExp; ++ iExp )
         Exponents.push_back(pExponents[iExp]);
      for ( uint iExp = 0; iExp < nExp; ++ iExp )
         for ( uint iCo = 0; iCo < nCo; ++ iCo )
            Coeffs.push_back(pCoeffs[iExp + nExp * iCo]);
   }
};

// F_{12}(r_{12}) = \sum_i c_i exp(-zeta_i r_{12}^2)
struct FGaussKernel : public FIntegralKernel
{
   explicit FGaussKernel( FGaussGeminalExp const *pGaussExp ) : m_pGaussExp(pGaussExp) {};
   char const *GetName() const { return "F12"; } // override
   virtual uint GetNumComponents() const { return m_pGaussExp->nCo; } // override
   void EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const; // override
   ~FGaussKernel();
private:
   FGaussGeminalExp const
      *m_pGaussExp;
};

// erf(omega r)/r  (long-range coulomb)
struct FErfCoulombKernel : public FIntegralKernel
{
   explicit FErfCoulombKernel( double Omega_ ) : m_Omega(Omega_) { InitBoysFnTable(); };
   char const *GetName() const { return "erf(w r)/r"; } // override
   virtual uint GetNumComponents() const { return 1; } // override
   void EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const; // override
   ~FErfCoulombKernel();
private:
   double
      m_Omega;
};

// erfc(omega r)/r  (short-range coulomb)
struct FErfcCoulombKernel : public FIntegralKernel
{
   explicit FErfcCoulombKernel( double Omega_ ) : m_Omega(Omega_) { InitBoysFnTable(); };
   char const *GetName() const { return "erfc(w r)/r"; } // override
   virtual uint GetNumComponents() const { return 1; } // override
   void EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const; // override
   ~FErfcCoulombKernel();
private:
   double
      m_Omega;
};


// F_{12}/r
struct FGaussCoulombKernel : public FIntegralKernel
{
   explicit FGaussCoulombKernel( FGaussGeminalExp const *pGaussExp ) : m_pGaussExp(pGaussExp) { InitBoysFnTable(); };
   char const *GetName() const { return "F12/r"; } // override
   virtual uint GetNumComponents() const { return m_pGaussExp->nCo; } // override
   void EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const; // override
   ~FGaussCoulombKernel();
private:
   FGaussGeminalExp const
      *m_pGaussExp;
};

// [grad F_{12}]^2
struct FGaussKineticKernel : public FIntegralKernel
{
   explicit FGaussKineticKernel( FGaussGeminalExp const *pGaussExp ) : m_pGaussExp(pGaussExp) {};
   char const *GetName() const { return "[[Del,F],F]"; } // override
   virtual uint GetNumComponents() const { return m_pGaussExp->nCo; } // override
   void EvalGm( double *pOut, double rho, double T, uint MaxM, double Prefactor ) const; // override
   ~FGaussKineticKernel();
private:
   FGaussGeminalExp const
      *m_pGaussExp;
};

} // namespace aic

#endif // _AIC_INT_KERNELS_H

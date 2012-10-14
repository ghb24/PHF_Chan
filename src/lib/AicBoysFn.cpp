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

// Code for evaluating the Boys function F_m(T) required in Coulomb integrals
// over Gaussian functions. Original version belongs to the AIC integral core.
// (Gerald Knizia)
#include <math.h>
#include <vector>

#include "AicDefs.h"
#include "AicBoysFn.h"

using namespace std;

typedef unsigned int
   uint;

namespace aic {

// in the near region (T<'TableEnd') of Fj(T), tabulate 'TalorOrder' order
// Taylor expansions every 'Step' Ts.
static const double
   Step = 0.1,
   InvStep = 1.0/Step,
   TableEnd = 16.0;
//    TableEnd = 12.0;
static const unsigned
   TaylorOrder = 8;

struct FBoysFnTable
{
   typedef std::vector<double>
      FTable;
   FTable
      Table, // will store nPoints * J scalars, representing
             // the F(t,j) for a finite grid of expansion points.
      Denominators; // hardly useful, but i'll give it a try.
   uint
      MaxJ,
      NumPoints,
      JsPerPoint; // number of J values for T

   inline double& Entry( uint IndexT, uint j ){
      assert( IndexT * JsPerPoint + j < Table.size() );
      return Table[IndexT * JsPerPoint + j];
   }

   static double FjtIterative( double T, int J );

   // cgk: inline because used only in one single place.
   inline void Eval( double *Result, int J, double T, double Prefactor );

   FBoysFnTable( int MaxJ );
   ~FBoysFnTable();
};


// calculate Fj(T) pseudoexactly by using a (slow but stable) iterative form.
double FBoysFnTable::FjtIterative( double T, int J )
{
   double
      Term = 1.0/(2*J + 1);
   volatile double
      c = 0, // Kahan summation.
      Sum = Term;
   uint
      j;
   for ( j = 1; j < 400000; ++ j )
   {  // This is the taylor expansion of Fj(T) at T=0.
      // note that all contributions to 'Sum' have the same sign,
      // therefore the order in which we add doesn't matter much.
      // However, for large T it requires a huge amount of terms to
      // converge, so it isn't totally stable there, too.
      Term *= (2*T) / (2*(J+j) + 1);
      volatile double
         y = Term - c,
         t = Sum + y;
      c = (t - Sum) - y;
      Sum = t;
      if ( Term/Sum < 1e-16 )
         break;
   }
//    if ( T == 16.0 )
//       cout << "iterative Fj(T): Number of sum terms for T == " << T <<  "  j =  " << j << std::endl;
   return Sum * exp( -T );
}


// calculate Fj(T) pseudoexactly by using a (slow but stable) iterative form.
// double FBoysFnTable::FjtIterative( double T, int J )
// {
//    double
//       Term = 1.0/(2*J + 1),
//       Sum = Term,
//       ExpT = exp(-T);
//    for ( uint j = 1; j < 400000; ++ j )
//    {
//       Term *= (2*T) / (2*(J+j) + 1);
//       Sum += Term;
//       if ( (Term/Sum) < 1e-15 )
//          break;
//    }
//    return Sum * ExpT;
// }


FBoysFnTable::FBoysFnTable( int MaxJ_ )
{
   this->MaxJ = static_cast<uint>(MaxJ_);
   JsPerPoint = MaxJ + TaylorOrder + 1;
   // ^- we need these additional values as derivatives of
   // the Boys function for our taylor expansions.
   NumPoints = ( 1 + static_cast<int>( ( TableEnd + 0.00001 ) / Step ) );
   Table.resize( JsPerPoint * NumPoints, 0.0 );


   Denominators.resize( 1 + JsPerPoint, 0.0 );
   Denominators[0] = 0.0;
   for ( uint i = 1; i <= JsPerPoint; ++i )
      Denominators[i] = 1.0/(2*i - 1);

   for ( uint i = 0; i < NumPoints; ++ i ){
      double
         T = Step * i;
      double
         ExpT = exp(-T);

      // calculate highest derivative with slow but stable iteration.
      Entry(i, JsPerPoint - 1) = FjtIterative(T, JsPerPoint - 1);

      // calculate remaining derivatives by downward recursion.
      for ( uint J = JsPerPoint - 1; J >= 1; -- J )
         Entry(i,J-1) = (Entry(i,J) * 2.0*T + ExpT)*Denominators[J];
   }

}

FBoysFnTable::~FBoysFnTable()
{}


inline void FBoysFnTable::Eval( double *Result_, int J, double T, double Prefactor )
{
   double
     *AIC_RP Result = Result_;
   double const
      SqrtPiHalf = 0.886226925452758; // sqrt(pi)/2

   assert( J <= (int)MaxJ );

//    if ( T < 30.0 )
//    {
//      Result[J] = Prefactor * FjtIterative( T, J );
//      ExpT = exp( -T ) * Prefactor;
//      for ( uint i = J; i > 0; --i ){
//         // Denom = Denom - 2.0;
//         Result[i-1] = ( 2*T * Result[i] + ExpT ) * Denominators[i];
//      }
//      return;
//    };

   uint
      iTable = static_cast<uint>( InvStep * T + 0.5 );
   // ^- table index.

   // outside of pretabulated range?
   if ( iTable >= NumPoints ){
      double
         InvT = 1.0/T,
         ExpT = 0,
         ErfFactor = 1;
      if ( T < 36 + 2*J ) {
         // large T asymptotic expansion:
         //    F_{0}(T) = sqrt(pi)/(2*sqrt(T)) * Erf(sqrt(T))
         //    F_{j+1}(T) = ((2j+1) * F_{j+1}(T) - exp(-T))/(2*T)
         // here assuming that exp(T) is 0 already (exp(-36) is about 1e-16)
         // and Erf(Sqrt(T)) is 1 -- but ExpT might still lead to non-zero
         // terms in the upward recursion (that's what the 2*J are for above)
         ExpT = exp(-T);
      }
      if ( T < 36 ) {
//       ErfFactor = erf(T*SqrtInvT);
//          ^- really slow.
         // F0[T_] = Sqrt[Pi]/2 * [  Erf[Sqrt[T]]/Sqrt[T] ]
         //        = Sqrt[Pi]/2 * [ (1 - Erfc[Sqrt[T])/Sqrt[T]] ]
         //        = Sqrt[Pi]/2 * (1/Sqrt[T] - Erfc[Sqrt[T]/Sqrt[T]])
         //        = Sqrt[Pi]/2 *(1/Sqrt[T]) * (1 - RationalA[T]*Exp[-T])

         // ErfcX is a minimax fit for Erfc[Sqrt[T]]/Exp[-T] for the region 10--20:
         // mma = GeneralMiniMaxApproximation[{T, Erfc[Sqrt[T]]/Exp[-T],1},
         //               {T, {10, 20}, 4, 4}, T, MaxItertations -> 500,
         //               Brake -> {100, 100}, WorkingPrecision -> 25]
         // For T > 10 it is good for a relative accuracy of <1e-18 for
         //     Erf[Sqrt[T]] \approx (1 - ErfcX * ExpT)
         double
            ErfcX =   (0.6718093055073827864155565 +
                        T*(0.4137432635707518305147398 +
                           T*(0.04197356405743738516935509 +
                              (0.0008197659391705324705371961 +
                                 1.560629841333219438084952e-6*T)*T)))/
                        (1. + T*(1.355035861352927307555204 +
                           T*(0.2919074786046937026719713 +
                              (0.01300581332979264635911687 +
                                 0.00009509476618141353796874409*T)*T)));
         ErfFactor = (1 - ErfcX * ExpT);
      }
      ExpT *= Prefactor;
      Result[0] = ErfFactor * SqrtPiHalf * sqrt(InvT) * Prefactor;

      // Compute the rest from Result[0] by upward recursion
      // (instable for small T).
      double
         InvTh = 0.5 * InvT;
      // upward recursion.
      for ( int i = 1; i <= J; ++ i )
         Result[i] = ((2*i-1) * Result[i-1] - ExpT) * InvTh;
   } else {
      // difference 'actual value - next point in table'
      double
         DeltaT = T - Step * iTable;

      // get value at reference point and evaluate Taylor expansion
      // for the rest
      double
         *AIC_RP Fjs = &Entry( iTable, J );

      // (-1)^n/(n!)
      double const
         Coeff2 = 1./2,
         Coeff3 = -1./(2*3),
         Coeff4 = 1./(2*3*4),
         Coeff5 = -1./(2*3*4*5),
         Coeff6 = 1./(2*3*4*5*6),
         Coeff7 = -1./(2*3*4*5*6*7),
         Coeff8 = 1./(2*3*4*5*6*7*8);

      Result[J] =
         Fjs[0] + DeltaT*( - Fjs[1]
                     + DeltaT * ( Coeff2 * Fjs[2]
                     + DeltaT * ( Coeff3 * Fjs[3]
                     + DeltaT * ( Coeff4 * Fjs[4]
                     + DeltaT * ( Coeff5 * Fjs[5]
                     + DeltaT * ( Coeff6 * Fjs[6]
                     + DeltaT * ( Coeff7 * Fjs[7]
                     + DeltaT * ( Coeff8 * Fjs[8] ))
               ))))));

      Result[J] *= Prefactor;

      if ( J == 0 )
         return;

      // Compute lower derivatives at target point by downward recursion
      double
         T2 = 2.0 * T,
         ExpT = exp(-T) * Prefactor;
      uint i = J;
      for ( ; i > 3; i -= 4 ){
         Result[i-1] = ( T2 * Result[i] + ExpT ) * Denominators[i];
         Result[i-2] = ( T2 * Result[i-1] + ExpT ) * Denominators[i-1];
         Result[i-3] = ( T2 * Result[i-2] + ExpT ) * Denominators[i-2];
         Result[i-4] = ( T2 * Result[i-3] + ExpT ) * Denominators[i-3];
      }
      for ( ; i != 0; --i )
         Result[i-1] = ( T2 * Result[i] + ExpT ) * Denominators[i];
   }
}

// extern "C" {
// 	void repulsion_f_(long const &nmax, double const &T, double *out, long const &istride, double const &fac);
// }


// static FBoysFnTable
//    BoysFnTable(30);
//
//
// // will write values to Target[0 .. J] (i.e. J+1 values).
// void BoysFn( double *Target, int J, double T, double Prefactor )
// {
//     BoysFnTable.Eval( Target, J, T, Prefactor );
// };


static FBoysFnTable
   *s_pBoysFnTable = 0;
void InitBoysFnTable()
{
    if ( s_pBoysFnTable == 0 )
        s_pBoysFnTable = new FBoysFnTable(30);
}

// will write values to Target[0 .. J] (i.e. J+1 values).
void BoysFn( double *Target, int J, double T, double Prefactor )
{
    assert(s_pBoysFnTable != 0); // this construction is a courtesy of the sunf90/gcc compiler
                                 // combination which does not initialize static objects.
    s_pBoysFnTable->Eval( Target, J, T, Prefactor );
// 	repulsion_f_(J, T, Target, 1, Prefactor);
}


} // namespace aic


// kate: space-indent on; tab-indent on; indent-width 3; mixedindent off; indent-mode normal;

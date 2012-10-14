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

#ifndef _AIC_PTR_H
#define _AIC_PTR_H

#include <string.h> // for memset
#include <vector>
#include <cmath>
#include <boost/intrusive_ptr.hpp>
#include "AicDefs.h"

#include "CxMemoryStack.h"
#include "CxTypes.h"

using std::size_t;

typedef unsigned int
   uint;

namespace aic {

typedef ct::FMemoryStack
   FMemoryStack;
typedef ct::FMemoryStack2
   FMemoryStack2;
typedef ct::FIntrusivePtrDest
   FIntrusivePtrDest;



template<class T>
inline T sqr(T const &a){
   return a*a;
}

template<class T>
struct TVector3
{
   typedef T value_type;

   TVector3() {}
   TVector3( T const &x_, T const &y_, T const &z_ ) {
      m[0] = x_; m[1] = y_; m[2] = z_;
   }
   explicit TVector3( T const (&xyz_)[3] ) {
      m[0] = xyz_[0]; m[1] = xyz_[1]; m[2] = xyz_[2];
   }

   T LengthSq() const {
      return m[0]*m[0] + m[1]*m[1] + m[2]*m[2];
   }

   T const &operator[] ( unsigned i ) const { return m[i]; }
   T &operator[] ( unsigned i ) { return m[i]; }

   inline bool operator == ( TVector3 const &other ) const;
   inline bool operator != ( TVector3 const &other ) const;
   inline void operator += ( TVector3 const &other );
   inline void operator -= ( TVector3 const &other );
   inline void operator *= ( T const &f );

   void Normalize() {
      (*this) *= 1./std::sqrt(this->LengthSq());
   }

   T const &x () const { return m[0]; }
   T const &y () const { return m[1]; }
   T const &z () const { return m[2]; }
   T &x () { return m[0]; }
   T &y () { return m[1]; }
   T &z () { return m[2]; }

   T m[3];
};

template<class T>
inline TVector3<T> operator + ( TVector3<T> const &a, TVector3<T> const &b ) {
   return TVector3<T>( a.m[0] + b.m[0], a.m[1] + b.m[1], a.m[2] + b.m[2] );
}
template<class T>
inline TVector3<T> operator - ( TVector3<T> const &a, TVector3<T> const &b ) {
   return TVector3<T>( a.m[0] - b.m[0], a.m[1] - b.m[1], a.m[2] - b.m[2] );
}
template<class T>
inline TVector3<T> operator * ( T const &f, TVector3<T> const &b ) {
   return TVector3<T>( f * b.m[0], f * b.m[1], f * b.m[2] );
}
template<class T>
inline T Dot( TVector3<T> const &a, TVector3<T> const &b ){
   return a.m[0] * b.m[0] + a.m[1] * b.m[1] + a.m[2] * b.m[2];
}
template<class T>
inline T DistanceSq( TVector3<T> const &a, TVector3<T> const &b ){
   return ::aic::sqr(a.m[0] - b.m[0]) + ::aic::sqr(a.m[1] - b.m[1]) + ::aic::sqr(a.m[2] - b.m[2]);
   // ^- HP compiler does not like sqr() without the namespace specification...
}
template<class T>
inline T Distance( TVector3<T> const &a, TVector3<T> const &b ){
   return std::sqrt(DistanceSq(a,b));
}
template<class T>
inline bool TVector3<T>::operator == ( TVector3 const &other ) const {
   double r = ::aic::sqr(m[0]-other.m[0]) + ::aic::sqr(m[1]-other.m[1]) + ::aic::sqr(m[2]-other.m[2]);
   return r == 0.0;
}
template<class T>
inline bool TVector3<T>::operator != ( TVector3 const &other ) const {
   return !this->operator == (other);
}
template<class T>
inline void TVector3<T>::operator += ( TVector3 const &other ) {
   m[0] += other.m[0]; m[1] += other.m[1]; m[2] += other.m[2];
}
template<class T>
inline void TVector3<T>::operator -= ( TVector3 const &other ) {
   m[0] -= other.m[0]; m[1] -= other.m[1]; m[2] -= other.m[2];
}
template<class T>
inline void TVector3<T>::operator *= ( T const &f ) {
   m[0] *= f; m[1] *= f; m[2] *= f;
}

template<class T>
static void Cross(TVector3<T> &Out, TVector3<T> const &a, TVector3<T> const &b)
{
   Out[0] = a[1]*b[2] - b[1] * a[2];
   Out[1] = a[2]*b[0] - b[2] * a[0];
   Out[2] = a[0]*b[1] - b[0] * a[1];
}

template<class T>
static TVector3<T> Cross(TVector3<T> const &a, TVector3<T> const &b)
{
   TVector3<T> Out;
   Cross(Out, a, b);
   return Out;
}

typedef TVector3<double>
   FVector3;


// r += f * x
template<class FScalar>
inline void Add( FScalar *AIC_RP r, FScalar const *AIC_RP x, FScalar f, std::size_t n )
{
   for ( std::size_t i = 0; i < n; ++ i )
      r[i] += f * x[i];
}


} // namespace aic


#endif // _AIC_PTR_H

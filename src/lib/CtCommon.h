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

#ifndef _CT_COMMON
#define _CT_COMMON

#include "AicDefs.h"
#include <stdint.h>
#define RESTRICT AIC_RP

#include "AicCommon.h"
using namespace aic;

#define _for_each(it,con) for ( (it) = (con).begin(); (it) != (con).end(); ++(it) )

namespace ct {
   typedef double
      FScalar;
   extern std::ostream
      &xerr, &xout;
   extern int
      Verbosity;
}


#ifndef assert_fast
    #define assert_fast assert
#endif
#ifndef assert_rt
    #define assert_rt assert
#endif

// note: these are macros and not real functions in order to allow the
// shadowing of xout and xout_ind by local variables. This design is by
// intention.
#define _xout_impl(level,x) \
    if ( Verbosity >= (level) ) \
    { xout << x << std::endl; }\
    else \
    {}
#define _xout(level,x) _xout_impl(level,x)
#define _xout0(x) _xout_impl(0,x)
#define _xout1(x) _xout_impl(1,x)
#define _xout2(x) _xout_impl(2,x)
#define _xout3(x) _xout_impl(3,x)
#define _xout4(x) _xout_impl(4,x)

// for "println debugging" in the absence of a good debugger...
#define _xdbgemit(x) "   " << #x << " = " << (x)
#define _xdbg(x) xout << _xdbgemit(x) << std::endl;
#define _xdbg2(x0,x1) xout << _xdbgemit(x0) << _xdbgemit(x1) << std::endl;
#define _xdbg3(x0,x1,x2) xout << _xdbgemit(x0) << _xdbgemit(x1) << _xdbgemit(x2) << std::endl;
#define _xdbg4(x0,x1,x2,x3) xout << _xdbgemit(x0) << _xdbgemit(x1) << _xdbgemit(x2) << _xdbgemit(x3) << std::endl;
#define _xdbg5(x0,x1,x2,x3,x4) xout << _xdbgemit(x0) << _xdbgemit(x1) << _xdbgemit(x2) << _xdbgemit(x3) << _xdbgemit(x4) << std::endl;
#define _xdbg6(x0,x1,x2,x3,x4,x5) xout << _xdbgemit(x0) << _xdbgemit(x1) << _xdbgemit(x2) << _xdbgemit(x3) << _xdbgemit(x4) << _xdbgemit(x5) << std::endl;
#define _xdbg7(x0,x1,x2,x3,x4,x5,x6) xout << _xdbgemit(x0) << _xdbgemit(x1) << _xdbgemit(x2) << _xdbgemit(x3) << _xdbgemit(x4) << _xdbgemit(x5) << _xdbgemit(x6) << std::endl;
#define _xdbg8(x0,x1,x2,x3,x4,x5,x6,x7) xout << _xdbgemit(x0) << _xdbgemit(x1) << _xdbgemit(x2) << _xdbgemit(x3) << _xdbgemit(x4) << _xdbgemit(x5) << _xdbgemit(x6) << _xdbgemit(x7) << std::endl;

// converts a stream output operation to a string.
#define _xstr(x) static_cast<const std::stringstream&>(std::stringstream() << x).str()


#ifdef __GNUC__
    // tell the compiler that this function will never return
    #define DECL_NORETURN __attribute__((noreturn))
    // 'insert software breakpoint here'.
    #define DEBUG_BREAK __asm__( "int $0x03" );
#elif _MSC_VER // microsoft c++. intel c++ may also understand this syntax.
    #define DECL_NORETURN __declspec(noreturn)
    #define DEBUG_BREAK __asm{ int 3 };
#else
    #define DECL_NORETURN
    #define DEBUG_BREAK
#endif


namespace ct {
   void FatalError( std::string const &Message,
       char const *pFromWhere, int nLine = -1 );// DECL_NORETURN;
}
#define ErrorExit(Why) ct::FatalError(Why, __FILE__, __LINE__)


namespace ct {
   // return 'true' if Key is mapped to something by Map.
   template< class TMap >
   inline bool exists( typename TMap::key_type const &Key, TMap const &Map ){
      return Map.end() != Map.find( Key );
   }

   // assert Key exists in map and return a const reference to its value.
   template< class TMap >
   inline typename TMap::mapped_type const&
   value( typename TMap::key_type const &Key, TMap const &Map ){
      typename TMap::const_iterator it = Map.find( Key );
      assert( Map.end() != it );
      return it->second;
   }

   // assert Key exists in map and return a reference to its value.
   template< class TMap >
   inline typename TMap::mapped_type &
   value( typename TMap::key_type const &Key, TMap &Map ){
      typename TMap::iterator it = Map.find( Key );
      assert( Map.end() != it );
      return it->second;
   }
}

#endif // _CT_COMMON

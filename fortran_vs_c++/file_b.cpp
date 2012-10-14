#include <iostream>
#include <boost/format.hpp>
using boost::format;
#include "CxFortranInt.h"
#include "CxPodArray.h"
using ct::TArray;

// a note: padding in structs might be different in C and Fortran.
// This will not be a problem when using 8byte integers and doubles,
// but it might be a problem in 32bit systems.
struct FIntsAndFloats{
   double a, b;
   FORTINT x;
   double c;
   FORTINT y;
};

extern "C" {
   void dummy_func_1_(FIntsAndFloats &obj1);
   void dummy_func_2_(TArray<FORTINT> &obj2);
   void dummy_func_3_(TArray<double> &obj3);
}

int main()
{
   FIntsAndFloats obj1 = {1.2, 2.3, 4, 5.6, 7};

   std::cout << "\nPASSING SIMPLE STRUCTURES\n";
   // works with 64bit integers.
   dummy_func_1_(obj1);
   std::cout << format("(C++ side) a = %f  b = %f  x = %i  c = %f  y = %i")
                  % obj1.a % obj1.b % obj1.x % obj1.c % obj1.y
             << std::endl;

   std::cout << "\nPASSING DYNAMIC ARRAYS: INTEGERS\n";
   // now for the interesting part...
   TArray<FORTINT> obj2;
   for ( uint i = 0; i < 8; ++ i )
      obj2.push_back(i*123+1012.);

//    std::cout << format("(C++ side) objs.pData = %i") % (uint64_t)&obj2[0] << std::endl;
//    std::cout << format("(C++ side) objs.pData[0] = %f") % (&obj2[0])[0] << std::endl;
   dummy_func_2_(obj2);
   std::cout << format("(C++ side) objs[] =  ");
   for ( uint i = 0; i < 8; ++ i )
      std::cout << format(" %10.4f") % obj2[i];
   std::cout << std::endl;

   std::cout << "\nPASSING DYNAMIC ARRAYS: FLOATS\n";
   TArray<double> obj3;
   for ( uint i = 0; i < 8; ++ i )
      obj3.push_back(i*123.23+1012.);

//    std::cout << format("(C++ side) objs.pData = %i") % (uint64_t)&obj3[0] << std::endl;
//    std::cout << format("(C++ side) objs.pData[0] = %f") % (&obj3[0])[0] << std::endl;
   dummy_func_3_(obj3);
   std::cout << format("(C++ side) objs[] =  ");
   for ( uint i = 0; i < 8; ++ i )
      std::cout << format(" %10.4f") % obj3[i];
   std::cout << std::endl;
}



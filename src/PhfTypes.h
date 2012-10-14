#ifndef PHF_TYPES_H
#define PHF_TYPES_H

#include "lib/CxTypes.h"

#include <complex>
typedef std::complex<double>
   complex_double;

#include <cstddef>
using std::size_t;
using std::ptrdiff_t;

#include "lib/CxFortranInt.h"
#include "lib/CxPodArray.h"
using ct::TArray;
using ct::TArrayRef;

#include "lib/AicCommon.h"
using aic::FVector3;
using aic::TVector3;
typedef TVector3<FORTINT>
   FVector3i;

#endif // PHF_TYPES_H

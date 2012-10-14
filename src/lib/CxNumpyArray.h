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

#ifndef CX_NUMPY_ARRAY_H
#define CX_NUMPY_ARRAY_H

#include <stdio.h> // for FILE*.
#include <stdexcept>
#include "CxPodArray.h"

// support for reading and writing array data in .npy format.
// (that's what numpy.save() and numpy.load() use)

namespace ct {

struct FIoExceptionNpy : public std::runtime_error {
   explicit FIoExceptionNpy(std::string const &Msg, std::string const &FileName = "");
};

typedef TArray<std::size_t>
   FShapeNpy;

// write the continuous array pData[i,j,k,..] to the given file in .npy
// format. pShape[i] gives the number of indices in dimension #i, nDim
// gives the total dimension of the array.
void WriteNpy(FILE *File, double const *pData, std::size_t const pShape[], std::size_t nDim);
void WriteNpy(std::string const &FileName, double const *pData, std::size_t const pShape[], std::size_t nDim);
void WriteNpy(FILE *File, double const *pData, FShapeNpy const &Shape);
void WriteNpy(std::string const &FileName, double const *pData, FShapeNpy const &Shape);

// convenience functions for making array shape objects.
FShapeNpy MakeShape(std::size_t i0);
FShapeNpy MakeShape(std::size_t i0, std::size_t i1);
FShapeNpy MakeShape(std::size_t i0, std::size_t i1, std::size_t i2);
FShapeNpy MakeShape(std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3);

struct FArrayNpy {
   TArray<double>
      Data;
   FShapeNpy
      Shape,
      Strides; // <- may be inverted when reading in stuff in C order.
   std::size_t Rank() const { return Shape.size(); }
};

// read array from File in .npy format.
void ReadNpy(FArrayNpy &Out, FILE *File);
void ReadNpy(FArrayNpy &Out, std::string const &FileName);

} // namespace ct

#endif // CX_NUMPY_ARRAY_H

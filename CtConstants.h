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

#ifndef CT_CONSTANTS_H
#define CT_CONSTANTS_H

namespace ct {
   static double const
      ToAng2006 = 0.529177249, // 2006 value.
      ToAng = 0.5291772108, // molpro default.
      ToAngCp2k = 0.529177211,
      ToEv2006 = 27.2113961,
      ToKcal = 627.5096;
} // namespace ct



#endif // CT_CONSTANTS_H

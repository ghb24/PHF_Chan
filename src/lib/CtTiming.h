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

// This file is released under the GNU General Public License ("GPL", version 2)
// as part of the CT8K program. Copying, modification, creating derivative works
// and redistribution of the file is allowed, but _only_ subject to the terms
// of that GPL. You should have received a version of this license along
// with this source code. The program comes "as is", without any kind of
// warranty.
//
// Authors/Copyright holders:
//  - Gerald Knizia, 2006 (tag: cgk, contact: cgk.d@gmx.net)


#ifndef _TIMING_H
#define _TIMING_H

#include <string>

namespace ct {

double Second();

struct FTimer
{
   FTimer() : Start(Second()) {};
   operator double () { return Second() - Start; };
private:
   double Start;
};

// Get time stamp counter (TSC) from current CPU core.
// WARNING: This is fragile and unportable! Use only for developer-side
// performance analysis and debugging!
// Notes:
//   - This counts CPU cycles on the current core
//   - Function should be fast (my guess: ~10-30 cycles), so it can be
//     used for timing small parts of code.
//   - TSCs on different CPU cores may differ. If the calling thread gets
//     re-distributed across cores between two successive calls, the TSC
//     may jump.
//   - In the presence of power management, espect varying cpu rates and
//     large jumps in the counter (suspended mode etc)
unsigned long long dbgGetTimeStampCounter();

} // namespace ct

#endif // _TIMING_H

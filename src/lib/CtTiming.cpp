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


#include <ctime>
#include <sys/time.h>

#include "CtTiming.h"

namespace ct {

double Second()
{
	timeval
		tv;
	gettimeofday( &tv, 0 );
	return tv.tv_sec + tv.tv_usec/1000000.0;
}

// get time stamp counter from current CPU. Fragile!
// Use only for debugging, if at all!
unsigned long long dbgGetTimeStampCounter()
{
	unsigned a, d;
	asm volatile("rdtsc" : "=a" (a), "=d" (d));

	return (((unsigned long long)a) | (((unsigned long long)d) << 32));
}

} // namespace ct

// kate: space-indent off; tab-indent on;

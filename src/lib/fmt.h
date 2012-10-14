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

#ifndef _CT8K_FMT_H
#define _CT8K_FMT_H

// provides some string-formatting adapters intended to make dealing
// with iostreams directly (without boost::format) less painfull.
//
// usage example:
//    void DoStuff{
//       using namespace fmt;
//       std::cout << fc("Some important energy:",20) << ff(76.02341234,14,6) << std::endl;
//    }
// This will print the string constant left-aligned with a width of 20 chars,
// and the float with with 14 and 6 decimal places.
//
// Supported adapters:
//    Floats:
//        ff(val,widh,prec)     fixed float, e.g. '12.54345'
//        fe(val,widh,prec)     exponential float, e.g. '1.254e+1'
//        fd(val,widh,prec)     "default float": fixed if reasonable value, exp. otherwise
//    Integers:
//        fi(val,width)
//    Strings:
//        fc(val,width)         accepts 'char const *' constants directly (0-terminated)
//        fs(val,width)         accepts std::strings.
//    Containers:
//        fl(cont,ewidth,prec,sep)  accepts any container, entry width, entry precision and entry separator
//        fl(first,last,ewidth,prec,sep)
//
// width and/or precission need not be specified.
// Contrary to FORTRAN formats, in C++'s IO streams 'width' specifies the
// _minimum_ number of characters printed. If the given value in the current
// setting does not fit into this number of characters, more will be printed
// instead of '******'.

#include <ostream>
#include <sstream>
#include <string>
// #include "assert.h"

namespace fmt{
    typedef std::ios_base io;

    template<class T>
    struct FFmtBase{
        typedef FFmtBase<T>
            FBase;
        // -1 for width or precision: nothing. Precision only works for floats.
        FFmtBase( io::fmtflags flags, T const &value_, int width_ = -1, int prec_ = -1 )
            : value(value_), width(width_), prec(prec_), ios_flags(flags)
        {}
        FFmtBase rebind( T const &new_value ){
            return FFmtBase(ios_flags,new_value,width,prec);
        }
        T const &value;
        int width;
        int prec;
        io::fmtflags ios_flags;
	private:
		void operator = ( FFmtBase const & ); // not implemented.
    };

    // integer
    struct fi : public FFmtBase<int>{ fi( int const &value, int width = -1 )
        : FBase(io::right, value, width) {}; };
    // fixed float
    struct ff : public FFmtBase<double>{ ff( double const &value, int width = -1, int prec = -1 )
        : FBase(io::right | io::fixed, value, width, prec) {}; };
    // exponential float
    struct fe : public FFmtBase<double>{ fe( double const &value, int width = -1, int prec = -1 )
        : FBase(io::right | io::scientific, value, width, prec) {}; };
    // default float (neither)
    struct fd : public FFmtBase<double>{ fd( double const &value, int width = -1, int prec = -1 )
        : FBase(io::right, value, width, prec) {}; };
    // string (should work for char const * and std::string)
    template<class FString>
    struct FFmtStr : public FFmtBase<FString>{ FFmtStr( FString const &value, int width = -1, int align = -1 )
        : FFmtBase<FString>((align==-1)?(io::left):(io::right),value,width) {}; };
    typedef FFmtStr<char const *> fc;
    typedef FFmtStr<std::string> fs;

    template<class T>
    std::ostream &operator << ( std::ostream &out, FFmtBase<T> const &f )
    {
        io::fmtflags
            old_flags = out.flags(f.ios_flags);
        int // <- gcc doesn't like ios_base::streamsize
            old_prec = static_cast<int>( out.precision() );
        if ( f.prec != -1 ) out.precision(f.prec);
        if ( f.width != -1 ) out.width(f.width);
        out << f.value;
        out.precision(old_prec);
        out.flags(old_flags);
        return out;
    };

    // not in class style because this would necessiate specifying the
    // concrete type of the container as template argument.
    template<class FIterator>
    std::string fl(FIterator first, FIterator last, int EntryWidth = -1, int EntryPrec = -1, char const *Separator = " ")
    {
        std::stringstream
            str;
        if ( EntryPrec != -1 )
            str.flags(std::ios_base::fixed);
        for ( FIterator it = first; it != last; ++it ){
            if ( it != first )
                str << Separator;
            if ( EntryWidth != -1 )
                str.width( EntryWidth );
            if ( EntryPrec != -1 )
                str.precision( EntryPrec );
            str << *it;
        }
        return str.str();
    };

    template<class FContainer>
    std::string fl(FContainer const &List, int EntryWidth = -1, int EntryPrec = -1, char const *Separator = " ")
    {
        return fl(List.begin(), List.end(), EntryWidth, EntryPrec, Separator);
    };



    // output float in full non-redundant format, but not longer than necessary.
    struct ffull{
        double v;
        ffull(double v_) : v(v_) {};
    };
    std::ostream &operator << ( std::ostream &out, ffull const &f );
}

#endif // _CT8K_FMT_H

// kate: space-indent on; tab-indent on; backspace-indent on; tab-width 4; indent-width 4; mixedindent off; indent-mode normal;

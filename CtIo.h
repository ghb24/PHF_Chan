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

#ifndef CT_IO_H
#define CT_IO_H

namespace ct {

extern char const
   *pResultFmt,
   *pResultFmtAnnoted,
   *pResultFmtI,
   *pResultFmtIAnnoted,
   *pTimingFmt;

} // namespace ct


#include "CtCommon.h"
#include "CxPodArray.h"


namespace ct {

// makes a string lower case. Will break with multibyte characters.
std::string tolower( std::string const &in );
// returns string 'in' with all spaces (but not tabs) removed.
std::string stripwhitespace( std::string const &in );

// makes a string to be output into the log as the first line of a major
// program component.
void MajorProgramIntro( std::ostream &out, const std::string &Name, const std::string &Version="" );

// will load an entire file into the stringstream str. Returns false if failed.
// if 0 != FileLength, *FileLength will receive the number of loaded chars.
bool LoadFileIntoMemory( TArray<char> &pFileContent,
        std::string const &FileName, uint *pFileLength = 0 );


enum FEnergyUnit{
    ENERGY_AtomicUnits, // = Hartree (1 Eh = 2 Ry = 27.211 3845(23) eV)
    ENERGY_Hartree = ENERGY_AtomicUnits,
    ENERGY_ElectronVolts,
    ENERGY_KjPerMol,
    ENERGY_KCalPerMol
};

extern FEnergyUnit
    g_StandardEnergyOutputUnit;

enum FDistanceUnit{
    DISTANCE_AtomicUnits,
    DISTANCE_Abohr = DISTANCE_AtomicUnits, // 0.5291772108 angstrom
    DISTANCE_Angstrom,
    DISTANCE_Picometer
};

extern FDistanceUnit
    g_StandardDistanceOutputUnit;

// controls the behavior of string outputting functions, namely Energy() and
// Distance(). (not supposed to do anything with internal values. These
// are still stored in a.u.)
void SetOutputUnits( FEnergyUnit EnergyUnit, FDistanceUnit DistanceUnit );

// this struct thing is used for the sole purpose of formatting energies
// when they are output into a stream by using the Energy() function;
// this allows us to transform internal energies into user specified formats.
struct FQuantityOutputAdapter{
    enum FQuantityType{
        QUANTITY_Energy,
        QUANTITY_Distance
    };

    FQuantityType
        QuantityType;
    double
        QuantityInAu; // energy/distance this thing represents in atomic units.
    bool
        ShowUnit;
    FQuantityOutputAdapter( FQuantityType Type, double QuantiyInAu_, bool ShowUnit_ = true )
        : QuantityType(Type), QuantityInAu( QuantiyInAu_ ), ShowUnit( ShowUnit_ )
    {}
};

std::ostream& operator << ( std::ostream& str, FQuantityOutputAdapter Quantiy );

// intended use:  cout << Energy( 20.00 ) << std::endl; Will convert energy
// to format set by g_StandardEnergyOutputUnit.
inline FQuantityOutputAdapter Energy( double EnergyInAu, bool ShowUnit = true )
{   return FQuantityOutputAdapter( FQuantityOutputAdapter::QUANTITY_Energy, EnergyInAu, ShowUnit ); }
std::string EnergyUnit(); // returns a string for energy unit which is currently set.
FScalar ConvertAuToOutputEnergy( FScalar EnergyInAu );

// intended use:  cout << Distance( 20.00 ) << std::endl; Will convert energy
// to format set by g_StandardDistanceOutputUnit.
inline FQuantityOutputAdapter Distance( double DistanceInAu, bool ShowUnit = true )
{   return FQuantityOutputAdapter( FQuantityOutputAdapter::QUANTITY_Distance, DistanceInAu, ShowUnit ); }
std::string DistanceUnit(); // returns a string for energy unit which is currently set.
FScalar ConvertAuToOutputDistance( FScalar DistanceInAu );

} // namespace ct.


#endif // CT_IO_H

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

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>

#include "CtCommon.h"
#include "CtConstants.h"
#include "CtIo.h"

#define INDENTSTREAM_IMPL
#include "CxIndentStream.h"
namespace ct {

char const
   *pResultFmt = "%-32s%18.12f\n",
   *pResultFmtAnnoted = "%-32s%18.12f  (%s)\n",
   *pResultFmtI = "%-32s%5i\n",
   *pResultFmtIAnnoted = "%-32s%5i  (%s)\n",
   *pTimingFmt = "Time for %s:%40t%10.2f sec\n";

int
   Verbosity = 0;

static fmt::FIndentStream1
    s_OutputStreamAdp(std::cout,true,1);
std::ostream
    &xout = s_OutputStreamAdp.stream,
    &xerr = std::cerr;

void FatalError( std::string const &Message,
    char const *pFromWhere, int nLine )
{
   std::stringstream str;
   str << Message;
   if ( pFromWhere ) {
      str << " at " << pFromWhere << ":" << nLine;
   };
   throw std::runtime_error(str.str());
};



FEnergyUnit
    g_StandardEnergyOutputUnit = ENERGY_AtomicUnits;
//  g_StandardEnergyOutputUnit = ENERGY_ElectronVolts;
//  g_StandardEnergyOutputUnit = ENERGY_KjPerMol;
FDistanceUnit
    g_StandardDistanceOutputUnit = DISTANCE_Angstrom;


void SetOutputUnits( FEnergyUnit EnergyUnit, FDistanceUnit DistanceUnit )
{
    g_StandardEnergyOutputUnit = EnergyUnit;
    g_StandardDistanceOutputUnit = DistanceUnit;
}


// makes a string to be output into the log as the first line of a major
// program component.
void MajorProgramIntro( std::ostream &out, const std::string &Name, const std::string &Version )
{
//  std::stringstream s;
    out << fmt::unind();
    out << "\n*** " << Name;
    if ( Version != "" )
        out << " [Ver. " << Version << "]";
    out << " ***\n";
    out << fmt::eind();
};

// makes a string lower case. Will break with localized multibyte characters,
// but we don't such.
std::string tolower( std::string const &in )
{
    std::string
        Result(in);
    std::string::iterator
        it;
    for( it = Result.begin(); it != Result.end(); ++it )
        *it = std::tolower(*it);
    return Result;
};

// returns string 'in' with all spaces (but not tabs) removed.
std::string stripwhitespace( std::string const &in )
{
    std::string
        Result;
    Result.reserve( in.size() );
    for ( size_t i = 0; i < in.size(); ++i )
        if ( ' ' != in[i] )
            Result.push_back( in[i] );
    return Result;
};

// returns a string for energy unit which is currently set.
std::string EnergyUnit()
{
    switch( g_StandardEnergyOutputUnit )
    {
        case ENERGY_AtomicUnits: return "Eh";
        case ENERGY_ElectronVolts: return "eV";
        case ENERGY_KjPerMol: return "kJ/mol";

        case ENERGY_KCalPerMol: return "kcal/mol";
        default:
            assert(0);
    };
    return "unk";
};


// converts the numerical value of the energy given to it in a.u. into the
// numerical value of the same energy in the current energy output format.
FScalar ConvertAuToOutputEnergy( FScalar EnergyInAu )
{
    switch( g_StandardEnergyOutputUnit )
    {
        case ENERGY_AtomicUnits:
            return EnergyInAu;
        case ENERGY_ElectronVolts: // all constants are taken from google.
            // 1 eV = 3.674 932 45(31)e-2 Hartree
            return EnergyInAu / 0.0367493245;
        case ENERGY_KjPerMol:
            // 1eV = 1.60217646e-19J   N_A = 6.02214199e23
            // E/[kJ/mol] = E/[eV]*N_A*1eV
            return EnergyInAu / 0.0367493245 * 6.02214199e23 * 1.60217646e-19 / 1000.0;
        case ENERGY_KCalPerMol:
            //  1 calorie = 4.18400 joules
            return EnergyInAu / 0.0367493245 * 6.02214199e23 * 1.60217646e-19 / 1000.0 / 4.18400;
        default:
            assert_rt(0);
            return 0.0;
    };
}

// returns a string for distance unit which is currently set.
std::string DistanceUnit()
{
    switch( g_StandardDistanceOutputUnit )
    {
        case DISTANCE_Abohr: return "abohr";
        case DISTANCE_Angstrom: return "A";
        case DISTANCE_Picometer: return "pm";
        default:
            assert(0);
    };
    return "unk";
};


// converts the numerical value of the distance given to it in a.u. (=abohr)
// into the numerical value of the same distance in the current energy output
// format.
FScalar ConvertAuToOutputDistance( FScalar DistanceInAu )
{
    switch( g_StandardDistanceOutputUnit )
    {
        case DISTANCE_Abohr:
            return DistanceInAu;
        case DISTANCE_Angstrom:
            return DistanceInAu * ToAng;
        case DISTANCE_Picometer:
            return DistanceInAu * 100.0 * ToAng;
        default:
            assert_rt(0);
            return 0.0;
    };
}


std::ostream& operator << ( std::ostream& str, FQuantityOutputAdapter Adapter )
{
    switch( Adapter.QuantityType )
    {
        case FQuantityOutputAdapter::QUANTITY_Energy:
            str << ConvertAuToOutputEnergy( Adapter.QuantityInAu );
            if ( Adapter.ShowUnit )
                str << EnergyUnit();
            break;
        case FQuantityOutputAdapter::QUANTITY_Distance:
            str << ConvertAuToOutputDistance( Adapter.QuantityInAu );
            if ( Adapter.ShowUnit )
                str << DistanceUnit();
            break;
        default:
            assert_rt(0);
    };
    return str;
};


// will load an entire file into the stringstream str. Returns false if failed.
bool LoadFileIntoMemory( TArray<char> &pFileContent,
        std::string const &FileName, uint *pFileLength )
{
    // read the entire file into an stringstream object.
    std::ifstream
        File( FileName.c_str() );
    std::size_t
        FileLength;
    if ( false == File.good() )
        return false;
    File.seekg( 0, std::ios::end );
    FileLength = File.tellg();
    if ( 0 != pFileLength )
        *pFileLength = FileLength;
    pFileContent.resize(2 + FileLength);
    memset( &pFileContent[0], 0, 2 + FileLength );
    File.seekg( 0, std::ios::beg );
    File.read( &pFileContent[0], FileLength );
    return true;
};


} // namespace ct

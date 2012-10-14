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

#ifndef _BASISLIBRARY_H
#define _BASISLIBRARY_H

#include <vector>
#include <map>
#include <list>
#include <set>
#include <string>

#include "AicShells.h"

#include "CtCommon.h"
// #include "CtAtomSet.h"
// #include "CtBasisSet.h"

namespace ct {

struct FBasisSet;

typedef std::vector<double>
   FScalarArray;
typedef std::vector<int>
   FIntArray;
typedef std::list<std::string>
   FStringList;

struct FBasisSetLibrary
{
// protected:
   struct FBasisKey {
      std::string const
         *pName; // pointer into this->BasisNames.
      uint
         iElement; // element number
      bool operator < ( FBasisKey const &other ) const {
         if ( pName < other.pName ) return true;
         if ( pName > other.pName ) return false;
         return iElement < other.iElement;
      }
   };
   FBasisKey MakeKey(std::string const &Name, uint iElement, bool AssertExists = true) const;
   void TryMakeBasis(std::string const &BasisDesc, int iElement);

   typedef std::multimap<FBasisKey, aic::FGaussBfnPtr>
      FBasisFnMap;
   FBasisFnMap
      m_BasisFns;

   typedef std::set<std::string>
      FBasisNameSet;
   FBasisNameSet
      m_BasisNames; // list of all names of imported entries (AO sets, Aux sets, Guesses, ECPs...)
public:
   // imports a .libmol file into memory (does not yet write any files)
   void ImportMolproLib( std::string const &FileName );

   // return false if failed. otherwise: _add_s basis functions to the
   // array, does not replace them.
   //   - nAtomIdx: index of atom in atom set allowed to be broken when not used
   //     anyway.
   //   - nElement: nuclear charge of atom.
   void LoadBasisFunctions( std::vector<aic::FGaussShell> &Shells,
            int iElement, std::string const &BasisDesc,
            FVector3 const &vAtomPos, FAtomIndex iAtomIdx ) const;
};

extern FBasisSetLibrary
   g_BasisSetLibrary; // <- no point in having more than one of these things around.

// returns g_BasisSetLibrary. Used for python interface because global vars
// can't really be exported as it seems.
FBasisSetLibrary& GetBasisSetLibrary();

} // namespace ct

#endif // _BASISLIBRARY_H

#ifndef PHF_SCF_H
#define PHF_SCF_H

#include "PhfSolidDef.h"
#include "PhfOperator.h"

struct FScfOptions
{
   double
      // convergence thresholds for SCF orbital gradient and energy
      ThrOrb, ThrDen;

   // (other SCF-option related stuff should also go here)
};




#endif // PHF_SCF_H

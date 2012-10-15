#ifndef DIIS_H
#define DIIS_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"
#include "RSPM.h"

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 */
class DIIS {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * @param output The stream to which you are writing (e.g. cout)
    * @param diis_p de DIIS you want to print
    */
   friend ostream &operator<<(ostream &output,DIIS &diis_p);

   public:
      
      //constructor
      DIIS();

      //copy constructor
      DIIS(const DIIS &);

      //destructor
      virtual ~DIIS();

      void push_comm(const RSPM &);

      void push_F(const RSPM &);

      void construct();

      const vector<RSPM> &gFD() const;

      const RSPM &gFD(int) const;

      const RSPM &gF(int) const;

      int gflag() const;

      int size() const;

      double gb(int) const;

   private:

      //!array containing the commutator [F,D] for all the previous iterations
      vector<RSPM> FD;

      //!array containing the Fock matrices of previous iterations
      vector<RSPM> F;

      double *b;

      Matrix *B;

      int flag;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/

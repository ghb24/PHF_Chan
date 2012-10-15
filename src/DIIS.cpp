#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "include.h"

DIIS::DIIS() {

   flag = 0;
   
}

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
DIIS::DIIS(const DIIS &diis_copy) {

   FD = diis_copy.gFD();
   flag = diis_copy.gflag();
   
}

/**
 * destructor
 */
DIIS::~DIIS(){ 
   
  if(flag == 1){

      delete B;
      delete [] b;

   }
  
}


ostream &operator<<(ostream &output,DIIS &diis_p){

   return output;

}

/**
 * add a new commutator to the DIIS object
 * @param FD_i new commutator to add
 */
void DIIS::push_comm(const RSPM &FD_i){

   FD.push_back(FD_i);

}

/**
 * add a new Fock matrix to the DIIS object
 * @param F_i new Fock matrix to add
 */
void DIIS::push_F(const RSPM &F_i){

   F.push_back(F_i);

}



/** 
 * @return the vector containing the commutator matrices
 */
const vector<RSPM> &DIIS::gFD() const {

   return FD;

}

/**
 * access to the individual commutators
 * @param i index of the iteration
 * @return commutator with index i
 */
const RSPM &DIIS::gFD(int i) const {

   return FD[i];

}

/**
 * access to the individual Fock Matrices
 * @param i index of the iteration
 * @return Fock matrix with index i
 */
const RSPM &DIIS::gF(int i) const {

   return F[i];

}

/**
 * construct the DIIS matrix
 */ 
void DIIS::construct(){

   if(flag == 1){

      delete B;
      delete [] b;

   }

   flag = 1;

   //construct the B matrix
   B = new Matrix(FD.size() + 1);

   for(unsigned int i = 0;i < FD.size();++i){

      for(unsigned int j = i;j < FD.size();++j)
         (*B)(i,j) = FD[i].ddot(FD[j]);

      (*B)(i,FD.size()) = 1.0;

   }

   (*B)(FD.size(),FD.size()) = 0.0;

   B->symmetrize();

   //and the vector on the right hand side
   b = new double [FD.size() + 1];

   for(unsigned int i = 0;i < FD.size();++i)
      b[i] = 0;

   b[FD.size()] = 1;

   //now solve!
   B->solve(b);

}

/**
 * @return the flag that indicates if the B and b vectors have been allocated
 */
int DIIS::gflag() const{

   return flag;

}

/**
 * access to the relaxation coefficients
 * @param i index of the coefficient
 */
double DIIS::gb(int i) const {

   if(flag == 0)
      std::cout << "This is wrong, you should never enter this function when flag = 0. Rubbish will come out." << std::endl;

   return b[i];

}

/**
 * @return the size of the system
 */
int DIIS::size() const{

   return FD.size();

}

/* vim: set ts=3 sw=3 expandtab :*/

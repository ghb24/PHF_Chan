#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "include.h"

int RSPM::N;
int *RSPM::dim;
int RSPM::nr;

/**
 * static function that initializes the static variables 
 * @param N_i Nr of particles
 * @param nr_i input nr of blocks
 * @param dim_i dimensions of the blocks
 */
void RSPM::init(int N_i,int nr_i,int *dim_i){

   nr = nr_i;
   dim = new int [nr];

   for(int B = 0;B < nr_i;++B)
      dim[B] = dim_i[B];

}

/**
 * delete the static stuff
 */
void RSPM::clear(){

   delete [] dim;

}

// takes in the number of blocks and an array of dimensions
RSPM::RSPM() : BlockMatrix(nr,dim) { }

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
RSPM::RSPM(const RSPM &spm_copy) : BlockMatrix(spm_copy) { }

/**
 * destructor
 */
RSPM::~RSPM(){ }

/**
 * @return nr of particles
 */
int RSPM::gN() const
{
   return N;
}

ostream &operator<<(ostream &output,RSPM &spm_p){

   for(int B = 0;B < spm_p.gnr();++B)
      for(int i = 0;i < spm_p.gdim(B);++i)
         for(int j = 0;j < spm_p.gdim(B);++j)
            output << i << "\t" << j << "\t" << spm_p(B,i,j) << endl;

   return output;

}

/**
 * fill the RSPM with the unit matrix on the first N/2 orbitals
 */
void RSPM::unit(){
/*
   *this = 0.0;

   for(int i = 0;i < N/2;++i)
      (*this)(i,i) = 1.0;
*/
}

/**
 * construct the Fock matrix, given the rspm from the previous iteration
 */
void RSPM::construct_Fock(const RSPM &spm){
/*
   for(int a = 0;a < M;++a)
      for(int c = a;c < M;++c){

         //first the SP part
         (*this)(a,c) = CartInt::gT()(a,c);

         for(int core = 0;core < input::gN_Z();++core)
            (*this)(a,c) += CartInt::gU(core)(a,c);

         //the Hartree term
         for(int b = 0;b < M;++b)
            for(int d = 0;d < M;++d)
               (*this)(a,c) += 2.0 *CartInt::gV()(a,b,c,d) * spm(b,d);

         //Fock
         for(int b = 0;b < M;++b)
            for(int d = 0;d < M;++d)
               (*this)(a,c) -= CartInt::gV()(a,b,d,c) * spm(b,d);

      }

   this->symmetrize();
*/
}

/**
 * construct the new rspm guess using the diagonalized Fock matrix
 */
void RSPM::update(const BlockVector &v,const RSPM &F){

   //to create the new rdm, we need the lowest N/2 eigenvectors
   *this = 0.0;

   //list containing the the index in the current block
   int list[nr];

   //start at zero of course
   for(int B = 0;B < nr;++B)
      list[B] = 0;

   double min;

   int B_min,i_min;

   for(int n = 0;n < N/2;++n){

      //make sure the list dimension isn't overflowed!
      int B_start = 0;

      while(list[B_start] >= dim[B_start])
         ++B_start;

      min = v(B_start,list[B_start]);

      B_min = B_start;
      i_min = list[B_start];

      for(int B = B_start + 1;B < nr;++B){

         if(list[B] < dim[B]){

            if(v(B,list[B]) < min){

               //new lowest value
               min = v(B,list[B]);

               //in block B with index list[B]
               B_min = B;
               i_min = list[B];

            }

         }

      }

      //for block B_min eigenvalue list[B_min] has been used:
      list[B_min]++;

      //construct part of rdm
      for(int a = 0;a < dim[B_min];++a)
         for(int b = a;b < dim[B_min];++b)
            (*this)(B_min,a,b) += F(B_min,a,i_min) * F(B_min,b,i_min);

   }

   this->symmetrize();

}

/**
 * construct the Fock matrix, given the rspm from the previous iteration
 */
void RSPM::construct_sp_ham(){
   /*
      for(int a = 0;a < M;++a)
      for(int c = a;c < M;++c){

   //first the SP part
   (*this)(a,c) = CartInt::gT()(a,c);

   for(int core = 0;core < input::gN_Z();++core)
   (*this)(a,c) += CartInt::gU(core)(a,c);

   }

   this->symmetrize();
    */
}

/**
 * take a linear combination of the matrices Fock matrices in the DIIS object using the relaxation coefficients
 */
void RSPM::relax(const DIIS &diis){

   *this = 0.0;

   for(int i = 0;i < diis.size();++i)
      this->daxpy(diis.gb(i),diis.gF(i));

}

/**
 * @return the expectation value of the energy
 */
double RSPM::energy() const{
   /*
      double energy = 0.0;

      for(int a = 0;a < M;++a)
      for(int c = 0;c < M;++c){

      energy += 2.0*(*this)(a,c) * CartInt::gT()(a,c);

      for(int core = 0;core < input::gN_Z();++core)
      energy += 2.0*(*this)(a,c) * CartInt::gU(core)(a,c);

      for(int b = 0;b < M;++b)
      for(int d = 0;d < M;++d)
      energy += (*this)(a,c) * ( 2.0 * CartInt::gV()(a,b,c,d) - CartInt::gV()(a,b,d,c) ) * (*this)(b,d);

      }

      return energy;
    */
   return 0.0;
}

/* vim: set ts=3 sw=3 expandtab :*/

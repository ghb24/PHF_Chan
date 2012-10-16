#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * constructor that takes dimension as input
 * @param nr nr of blocks
 * @param dim array containing the dimensions
 */
BlockVector::BlockVector(int nr,int *dim){

   this->nr = nr;

   this->dim = new int [nr];

   vector = new double * [nr];

   for(int B = 0;B < nr;++B){

      this->dim[B] = dim[B];
      
      vector[B] = new double [dim[B]];

   }

}

/**
 * Construct and initialize the BlockVector object by diagonalizing a Matrix object:
 */
BlockVector::BlockVector(BlockMatrix &matrix){

   //allocate
   this->nr = matrix.gnr();

   this->dim = new int [nr];

   vector = new double * [nr];

   for(int B = 0;B < nr;++B){

      this->dim[B] = matrix.gdim(B);

      vector[B] = new double [dim[B]];

   }

   for(int B = 0;B < nr;++B){

      //initialize
      char jobz = 'V';
      char uplo = 'U';

      int lwork = 3*dim[B] - 1;

      double *work = new double [lwork];

      int info;

      dsyev_(&jobz,&uplo,&dim[B],(matrix.gBlockMatrix())[B][0],&dim[B],vector[B],work,&lwork,&info);

      delete [] work;

   }


}

/**
 * copy constructor 
 * @param vec_copy The vector you want to be copied into the object you are constructing
 */
BlockVector::BlockVector(const BlockVector &vec_copy){

   this->nr = vec_copy.gnr();

   this->dim = new int [nr];

   for(int B = 0;B < nr;++B)
      this->dim[B] = vec_copy.gdim(B);

   vector = new double * [nr];

   for(int B = 0;B < nr;++B){

      vector[B] = new double [dim[B]];

      int inc = 1;

      dcopy_(&dim[B],vec_copy.vector[B],&inc,vector[B],&inc);

   }

}

/**
 * Destructor
 */
BlockVector::~BlockVector(){

   for(int B = 0;B < nr;++B)
      delete [] vector[B];

   delete [] vector;

   delete [] dim;

}

/**
 * overload the equality operator
 * @param vector_copy The vector you want to be copied into this
 */
BlockVector &BlockVector::operator=(const BlockVector &vector_copy){

   for(int B = 0;B < nr;++B){

      int inc = 1;

      dcopy_(&dim[B],vector_copy.vector[B],&inc,vector[B],&inc);

   }

   return *this;

}

/**
 * Make all the number in your vector equal to the number a, e.g. usefull for initialization (BlockVector M = 0)
 * @param a the number
 */
BlockVector &BlockVector::operator=(double a){

   for(int B = 0;B < nr;++B)
      for(int i = 0;i < dim[B];++i)
         vector[B][i] = a;

   return *this;

}

/**
 * overload the += operator for matrices
 * @param vector_pl The vector you want to add to this
 */
BlockVector &BlockVector::operator+=(const BlockVector &vector_pl){

   for(int B = 0;B < nr;++B){

      int inc = 1;
      double alpha = 1.0;

      daxpy_(&dim[B],&alpha,vector_pl.vector[B],&inc,vector[B],&inc);

   }

   return *this;

}

/**
 * overload the -= operator for matrices
 * @param vector_pl The vector you want to add to this
 */
BlockVector &BlockVector::operator-=(const BlockVector &vector_pl){

   for(int B = 0;B < nr;++B){

      int inc = 1;
      double alpha = -1.0;

      daxpy_(&dim[B],&alpha,vector_pl.vector[B],&inc,vector[B],&inc);

   }

   return *this;

}



/**
 * add the vector vector_pl times the constant alpha to this
 * @param alpha the constant to multiply the vector_pl with
 * @param vector_pl the BlockVector to be multiplied by alpha and added to this
 */
BlockVector &BlockVector::daxpy(double alpha,const BlockVector &vector_pl){

   for(int B = 0;B < nr;++B){

      int inc = 1;

      daxpy_(&dim[B],&alpha,vector_pl.vector[B],&inc,vector[B],&inc);

   }

   return *this;

}

/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your vector through
 */
BlockVector &BlockVector::operator/=(double c){

   int inc = 1;

   double alpha = 1.0/c;

   for(int B = 0;B < nr;++B)
      dscal_(&dim[B],&alpha,vector[B],&inc);

   return *this;

}

/**
 * write access to your vector, change the number on index i
 * @param i row number
 * @return the entry on place i
 */
double &BlockVector::operator()(int B,int i){

   return vector[B][i];

}

/**
 * read access to your vector, change the number on index i: const version
 * @param i row number
 * @return the entry on place i
 */
double BlockVector::operator()(int B,int i) const{

   return vector[B][i];

}

/**
 * Diagonalize the Matrix matrix when you have allready allocated the memory of the vector
 * on the correct dimension.
 */
void BlockVector::diagonalize(BlockMatrix &matrix){

   for(int B = 0;B < nr;++B){

      char jobz = 'V';
      char uplo = 'U';

      int lwork = 3*dim[B] - 1;

      double *work = new double [lwork];

      int info;

      dsyev_(&jobz,&uplo,&dim[B],(matrix.gBlockMatrix())[B][0],&dim[B],vector[B],work,&lwork,&info);

      delete [] work;

   }

}

/**
 * @return the underlying pointer to vector, useful for mkl and lapack applications
 */
double **BlockVector::gBlockVector(){

   return vector;

}

/**
 * @return the total nr of blocks
 */
int BlockVector::gnr() const{

   return nr;

}

/**
 * @param B blockindex
 * @return the dimension of block B
 */
int BlockVector::gdim(int B) const{

   return dim[B];

}

/**
 * @return inproduct of (*this) vector with vector_i
 * @param vector_i input vector
 */
double BlockVector::ddot(const BlockVector &vector_i)const{

   int inc = 1;

   double ward = 0.0;

   for(int B = 0;B < nr;++B)
      ward += ddot_(&dim[B],vector[B],&inc,vector_i.vector[B],&inc);

   return ward;

}

/**
 * Scale the vector (*this) with parameter alpha
 * @param alpha scalefactor
 */
void BlockVector::dscal(double alpha){

   int inc = 1;

   for(int B = 0;B < nr;++B)
      dscal_(&dim[B],&alpha,vector[B],&inc);

}

/**
 * Fill the vector with random numbers.
 */
void BlockVector::fill_Random(){

   srand(time(NULL));

   for(int B = 0;B < nr;++B)
      for(int i = 0;i < dim[B];++i)
         vector[B][i] = (double) rand()/RAND_MAX;

}

ostream &operator<<(ostream &output,const BlockVector &vector_p){

   for(int B =0;B < vector_p.gnr();++B){

      output << endl;
      output << "Block\t" << B << endl;
      output << endl;

      for(int i = 0;i < vector_p.gdim(B);++i)
         output << i << "\t" << vector_p(B,i) << endl;

   }

   return output;

}

/**
 * @return the minimal element present in this BlockVector object.
 * watch out, only works when BlockVector is filled with the eigenvalues of a diagonalized Matrix object
 */
double BlockVector::min() const {

   double min = vector[0][0];

   for(int B = 1;B < nr;++B)
      if(vector[B][0] < min)
          min = vector[B][0];

   return min;

}

/**
 * @return the maximal element present in this BlockVector object.
 * watch out, only works when BlockVector is filled with the eigenvalues of a diagonalized Matrix object
 */
double BlockVector::max() const {

   double max = vector[0][dim[0]-1];

   for(int B = 1;B < nr;++B)
      if(vector[B][dim[B]-1] > max)
         vector[B][dim[B]-1] = max;

   return max;

}

/* vim: set ts=3 sw=3 expandtab :*/

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
 * @param n dimension of the vector
 */
Vector::Vector(int n){

   this->n = n;

   vector = new double [n];

}

/**
 * Construct and initialize the Vector object by diagonalizing a Matrix object:
 */
Vector::Vector(Matrix &matrix){

   //allocate
   this->n = matrix.gn();

   vector = new double [n];

   //initialize
   char jobz = 'V';
   char uplo = 'U';

   int lwork = 3*n - 1;

   double *work = new double [lwork];

   int info;

   dsyev_(&jobz,&uplo,&n,(matrix.gMatrix())[0],&n,vector,work,&lwork,&info);

   delete [] work;

}

/**
 * copy constructor 
 * @param vec_copy The vector you want to be copied into the object you are constructing
 */
Vector::Vector(const Vector &vec_copy){

   this->n = vec_copy.n;

   vector = new double [n];

   int inc = 1;

   dcopy_(&n,vec_copy.vector,&inc,vector,&inc);

}


/**
 * Destructor
 */
Vector::~Vector(){

   delete [] vector;

}

/**
 * overload the equality operator
 * @param vector_copy The vector you want to be copied into this
 */
Vector &Vector::operator=(const Vector &vector_copy){

   int inc = 1;

   dcopy_(&n,vector_copy.vector,&inc,vector,&inc);

   return *this;

}

/**
 * Make all the number in your vector equal to the number a, e.g. usefull for initialization (Vector M = 0)
 * @param a the number
 */
Vector &Vector::operator=(double a){

   for(int i = 0;i < n;++i)
      vector[i] = a;

   return *this;

}

/**
 * overload the += operator for matrices
 * @param vector_pl The vector you want to add to this
 */
Vector &Vector::operator+=(const Vector &vector_pl){

   int inc = 1;
   double alpha = 1.0;

   daxpy_(&n,&alpha,vector_pl.vector,&inc,vector,&inc);

   return *this;

}

/**
 * overload the -= operator for matrices
 * @param vector_pl The vector you want to deduct from this
 */
Vector &Vector::operator-=(const Vector &vector_pl){

   int inc = 1;
   double alpha = -1.0;

   daxpy_(&n,&alpha,vector_pl.vector,&inc,vector,&inc);

   return *this;

}

/**
 * add the vector vector_pl times the constant alpha to this
 * @param alpha the constant to multiply the vector_pl with
 * @param vector_pl the Vector to be multiplied by alpha and added to this
 */
Vector &Vector::daxpy(double alpha,const Vector &vector_pl){

   int inc = 1;

   daxpy_(&n,&alpha,vector_pl.vector,&inc,vector,&inc);

   return *this;

}

/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your vector through
 */
Vector &Vector::operator/=(double c){

   int inc = 1;

   double alpha = 1.0/c;

   dscal_(&n,&alpha,vector,&inc);

   return *this;

}

/**
 * write access to your vector, change the number on index i
 * @param i row number
 * @return the entry on place i
 */
double &Vector::operator[](int i){

   return vector[i];

}

/**
 * read access to your vector, change the number on index i: const version
 * @param i row number
 * @return the entry on place i
 */
double Vector::operator[](int i) const{

   return vector[i];

}

/**
 * Diagonalize the Matrix matrix when you have allready allocated the memory of the vector
 * on the correct dimension.
 */
void Vector::diagonalize(Matrix &matrix){

   char jobz = 'V';
   char uplo = 'U';

   int lwork = 3*n - 1;

   double *work = new double [lwork];

   int info;

   dsyev_(&jobz,&uplo,&n,(matrix.gMatrix())[0],&n,vector,work,&lwork,&info);

   delete [] work;

}

/**
 * @return the underlying pointer to vector, useful for mkl and lapack applications
 */
double *Vector::gVector(){

   return vector;

}

/**
 * @return the underlying pointer to vector, useful for mkl and lapack applications: const version
 */
const double *Vector::gVector() const{

   return vector;

}

/**
 * @return the dimension of the vector
 */
int Vector::gn() const{

   return n;

}

/**
 * @return the sum of all the elements in the vector
 */
double Vector::sum() const{

   double ward = 0;

   for(int i = 0;i < n;++i)
      ward += vector[i];

   return ward;

}

/**
 * @return the logarithm of the product of all the elements in the vector (so the sum of all the logarithms)
 */
double Vector::log_product() const{

   double ward = 0;

   for(int i = 0;i < n;++i)
      ward += log(vector[i]);

   return ward;

}

/**
 * @return inproduct of (*this) vector with vector_i
 * @param vector_i input vector
 */
double Vector::ddot(const Vector &vector_i)const{

   int inc = 1;

   return ddot_(&n,vector,&inc,vector_i.vector,&inc);

}

/**
 * Scale the vector (*this) with parameter alpha
 * @param alpha scalefactor
 */
void Vector::dscal(double alpha){

   int inc = 1;

   dscal_(&n,&alpha,vector,&inc);

}

/**
 * Fill the vector with random numbers.
 */
void Vector::fill_Random(){

   srand(time(NULL));

   for(int i = 0;i < n;++i)
      vector[i] = (double) rand()/RAND_MAX;

}

ostream &operator<<(ostream &output,const Vector &vector_p){

   for(int i = 0;i < vector_p.gn();++i)
      output << i << "\t" << vector_p[i] << endl;

   return output;

}

/**
 * @return the minimal element present in this Vector object.
 * watch out, only works when Vector is filled with the eigenvalues of a diagonalized Matrix object
 */
double Vector::min() const {

   return vector[0];

}

/**
 * @return the maximal element present in this Vector object.
 * watch out, only works when Vector is filled with the eigenvalues of a diagonalized Matrix object
 */
double Vector::max() const {

   return vector[n - 1];

}

/* vim: set ts=3 sw=3 expandtab :*/

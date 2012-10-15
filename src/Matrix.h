#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cstdlib>

using std::ostream;

class Vector;

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * This is a class written for symmetric matrices. It is a wrapper around a double pointer and
 * redefines much used lapack and blas routines as memberfunctions
 */
class Matrix{

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << matrix << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << matrix << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param matrix_p de Matrix you want to print
    */
   friend ostream &operator<<(ostream &output,const Matrix &matrix_p);

   public:

      //constructor
      Matrix(int n);

      //copy constructor
      Matrix(const Matrix &);

      //construct with filename
      Matrix(const char *filename);

      //destructor
      virtual ~Matrix();

      //overload equality operator
      Matrix &operator=(const Matrix &);

      Matrix &operator=(double );

      //overload += operator
      Matrix &operator+=(const Matrix &);

      //overload -= operator
      Matrix &operator-=(const Matrix &);

      Matrix &daxpy(double alpha,const Matrix &);

      Matrix &operator/=(double );

      Matrix &mprod(const Matrix &,const Matrix &);

      Matrix &commute(const Matrix &,const Matrix &);

      //easy to change the numbers
      double &operator()(int i,int j);

      //easy to access the numbers
      double operator()(int i,int j) const;

      //get the pointer to the matrix
      double **gMatrix();

      int gn() const;

      double trace() const;

      double ddot(const Matrix &) const;

      void invert();

      void solve(double *);

      void dscal(double alpha);

      void fill_Random();

      //positieve of negatieve vierkantswortel uit de matrix
      void sqrt(int option);

      void mdiag(const Vector &diag);

      void L_map(const Matrix &,const Matrix &);

      void symmetrize();

      void out(const char*) const;

      void sep_pm(Matrix &p,Matrix &m);

   private:

      //!double pointer of doubles, contains the numbers, the matrix
      double **matrix;

      //!dimension of the matrix
      int n;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/

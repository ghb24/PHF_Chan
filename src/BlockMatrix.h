#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H

#include <iostream>
#include <cstdlib>

using std::ostream;

class BlockVector;

/**
 * @author Brecht Verstichel
 * @date 18-02-2010\n\n
 * This is a class written for symmetric matrices. It is a wrapper around a double pointer and
 * redefines much used lapack and blas routines as memberfunctions
 */
class BlockMatrix{

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << matrix << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << matrix << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param matrix_p de BlockMatrix you want to print
    */
   friend ostream &operator<<(ostream &output,const BlockMatrix &matrix_p);

   public:

      //constructor
      BlockMatrix(int nr,int *deg);

      //copy constructor
      BlockMatrix(const BlockMatrix &);

      //destructor
      virtual ~BlockMatrix();

      //overload equality operator
      BlockMatrix &operator=(const BlockMatrix &);

      BlockMatrix &operator=(double );

      //overload += operator
      BlockMatrix &operator+=(const BlockMatrix &);

      //overload -= operator
      BlockMatrix &operator-=(const BlockMatrix &);

      BlockMatrix &daxpy(double alpha,const BlockMatrix &);

      BlockMatrix &operator/=(double );

      BlockMatrix &mprod(const BlockMatrix &,const BlockMatrix &);

      BlockMatrix &commute(const BlockMatrix &,const BlockMatrix &);

      //easy to change the numbers
      double &operator()(int B,int i,int j);

      //easy to access the numbers
      double operator()(int B,int i,int j) const;

      //get the pointer to the matrix
      double ***gBlockMatrix();

      int gnr() const;

      int gdim(int) const;

      const int *gdim() const;

      double trace() const;

      double ddot(const BlockMatrix &) const;

      void invert();

      void solve(double *);

      void dscal(double alpha);

      void fill_Random();

      //positieve of negatieve vierkantswortel uit de matrix
      void sqrt(int option);

      void mdiag(const BlockVector &diag);

      void L_map(const BlockMatrix &,const BlockMatrix &);

      void symmetrize();

      void out(const char*) const;

      void sep_pm(BlockMatrix &p,BlockMatrix &m);

   private:

      //!three-double pointer of doubles, contains the numbers, the matrix
      double ***matrix;

      //!nr of blocks
      int nr;

      //!array of size nr containing the dimensions of the blocks
      int *dim;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/

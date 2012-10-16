#ifndef BLOCKVECTOR_H
#define BLOCKVECTOR_H

#include <iostream>
#include <cstdlib>
#include <cmath>

using std::ostream;

#include "lapack.h"
#include "BlockMatrix.h"

/**
 * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
 * ifstream object and type:\n\n
 * object << vector << endl;\n\n
 * For output onto the screen type: \n\n
 * cout << vector << endl;\n\n
 * @param output The stream to which you are writing (e.g. cout)
 * @param vector_p de BlockVector you want to print
 */
ostream &operator<<(ostream &output,const BlockVector &vector_p);

/**
 * @author Brecht Verstichel
 * @date 15-04-2010\n\n
 * This is a class written for vectors. It will contain the eigenvalues of the TPM, etc. Matrices. It is a template class,
 * corresponding to the different BlockVectorType's that can be put in, it will automatically get the right dimension.
 * It is a wrapper around a pointer and redefines much used lapack and blas routines as memberfunctions
 */
class BlockVector{

   public:

      //construct with as input a dimension
      BlockVector(int ,int*);

      //construct with as input a Matrix
      BlockVector(BlockMatrix &);

      //copy constructor
      BlockVector(const BlockVector &);

      //destructor
      virtual ~BlockVector();

      //overload equality operator
      BlockVector &operator=(const BlockVector &);

      BlockVector &operator=(double );

      //overload += operator
      BlockVector &operator+=(const BlockVector &);

      //overload -= operator
      BlockVector &operator-=(const BlockVector &);

      BlockVector &daxpy(double alpha,const BlockVector &);

      BlockVector &operator/=(double );

      void diagonalize(BlockMatrix &);

      //easy to change the numbers
      double &operator()(int B,int i);

      double operator()(int B,int i) const;

      void fill_Random();

      //get the pointer to the vector
      double **gBlockVector();

      int gn() const;

      double sum() const;

      double log_product() const;

      double ddot(const BlockVector &) const;

      void dscal(double alpha);

      int gnr() const;

      int gdim(int) const;

      double min() const;

      double max() const;

   private:

      //!double pointer of doubles, contains the numbers, the vector
      double **vector;

      //!nr of blocks
      int nr;

      //!dimension of the blocks
      int *dim;


};

#endif

/* vim: set ts=3 sw=3 expandtab :*/

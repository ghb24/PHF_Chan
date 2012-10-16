#include <iostream>
#include <fstream>
#include <cmath>

using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * constructor 
 * @param nr nr of blocks
 * @param dim pointer containing the dimensions of the blocks
 */
BlockMatrix::BlockMatrix(int nr,int *dim){

   this->nr = nr;

   this->dim = new int [nr];

   matrix = new double ** [nr];

   for(int B = 0;B < nr;++B){

      this->dim[B] = dim[B];

      matrix[B] = new double * [dim[B]];
      matrix[B][0] = new double [dim[B]*dim[B]];

      for(int i = 1;i < dim[B];++i)
         matrix[B][i] = matrix[B][i - 1] + dim[B];

   }

}

/**
 * copy constructor 
 * @param mat_copy The matrix you want to be copied into the object you are constructing
 */
BlockMatrix::BlockMatrix(const BlockMatrix &mat_copy){

   this->nr = mat_copy.gnr();

   this->dim = new int [nr];

   for(int B = 0;B < nr;++B)
      dim[B] = mat_copy.gdim(B);

   matrix = new double ** [nr];

   for(int B = 0;B < nr;++B){

      this->dim[B] = dim[B];

      matrix[B] = new double * [dim[B]];
      matrix[B][0] = new double [dim[B]*dim[B]];

      for(int i = 1;i < dim[B];++i)
         matrix[B][i] = matrix[B][i - 1] + dim[B];

      int lapdim = dim[B]*dim[B];

      int incx = 1;
      int incy = 1;

      dcopy_(&lapdim,mat_copy.matrix[B][0],&incx,matrix[B][0],&incy);

   }
}

/**
 * Destructor
 */
BlockMatrix::~BlockMatrix(){

   delete [] dim;

   for(int B = 0;B < nr;++B){

      delete [] matrix[B][0];
      delete [] matrix[B];

   }

   delete [] matrix;

}

/**
 * overload the equality operator
 * @param matrix_copy The matrix you want to be copied into this
 */
BlockMatrix &BlockMatrix::operator=(const BlockMatrix &matrix_copy){

   for(int B = 0;B < nr;++B){

      int lapdim = dim[B]*dim[B];
      int incx = 1;
      int incy = 1;

      dcopy_(&lapdim,matrix_copy.matrix[B][0],&incx,matrix[B][0],&incy);

   }

   return *this;

}

/**
 * Make all the number in your matrix equal to the number a, e.g. usefull for initialization (BlockMatrix M = 0)
 * @param a the number
 */
BlockMatrix &BlockMatrix::operator=(double a){

   for(int B = 0;B < nr;++B)
      for(int i = 0;i < dim[B];++i)
         for(int j = 0;j < dim[B];++j)
            matrix[B][j][i] = a;

   return *this;

}

/**
 * overload the += operator for matrices
 * @param matrix_pl The matrix you want to add to this
 */
BlockMatrix &BlockMatrix::operator+=(const BlockMatrix &matrix_pl){

   for(int B = 0;B < nr;++B){

      int lapdim = dim[B]*dim[B];
      int inc = 1;
      double alpha = 1.0;

      daxpy_(&lapdim,&alpha,matrix_pl.matrix[B][0],&inc,matrix[B][0],&inc);

   }

   return *this;

}

/**
 * overload the -= operator for matrices
 * @param matrix_pl The matrix you want to deduct from this
 */
BlockMatrix &BlockMatrix::operator-=(const BlockMatrix &matrix_pl){

   for(int B = 0;B < nr;++B){

      int lapdim = dim[B]*dim[B];
      int inc = 1;
      double alpha = -1.0;

      daxpy_(&lapdim,&alpha,matrix_pl.matrix[B][0],&inc,matrix[B][0],&inc);

   }

   return *this;

}

/**
 * add the matrix matrix_pl times the constant alpha to this
 * @param alpha the constant to multiply the matrix_pl with
 * @param matrix_pl the BlockMatrix to be multiplied by alpha and added to this
 */
BlockMatrix &BlockMatrix::daxpy(double alpha,const BlockMatrix &matrix_pl){

   for(int B = 0;B < nr;++B){

      int lapdim = dim[B]*dim[B];
      int inc = 1;

      daxpy_(&lapdim,&alpha,matrix_pl.matrix[B][0],&inc,matrix[B][0],&inc);

   }

   return *this;

}
/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your matrix through
 */
BlockMatrix &BlockMatrix::operator/=(double c){

   for(int B = 0;B < nr;++B){

      int lapdim = dim[B]*dim[B];
      int inc = 1;

      double alpha = 1.0/c;

      dscal_(&lapdim,&alpha,matrix[B][0],&inc);

   }

   return *this;

}

/**
 * write access to your matrix, change the number on row i and column j
 * remark that for the conversion to lapack functions the double pointer is transposed!
 * @param B block index
 * @param i row number
 * @param j column number
 * @return the entry on place i,j
 */
double &BlockMatrix::operator()(int B,int i,int j){

   return matrix[B][j][i];

}

/**
 * read access to your matrix, view the number on row i and column j
 * remark that for the conversion to lapack functions the double pointer is transposed!
 * @param B block index
 * @param i row number
 * @param j column number
 * @return the entry on place i,j
 */
double BlockMatrix::operator()(int B,int i,int j) const {

   return matrix[B][j][i];

}

/**
 * @return the underlying pointer to matrix, useful for mkl applications
 */
double ***BlockMatrix::gBlockMatrix(){

   return matrix;

}

/**
 * @return the nr of blocks in the matrix
 */
int BlockMatrix::gnr() const{

   return nr;

}

/**
 * @param B the blockindex
 * @return the dimension of block with index B
 */
int BlockMatrix::gdim(int B) const{

   return dim[B];

}

/**
 * @return the trace of the matrix:
 */
double BlockMatrix::trace() const{

   double ward = 0;

   for(int B = 0;B < nr;++B)
      for(int i = 0;i < dim[B];++i)
         ward += matrix[B][i][i];

   return ward;

}

/**
 * @return inproduct of (*this) matrix with matrix_i, defined as Tr (A B)
 * @param matrix_i input matrix
 */
double BlockMatrix::ddot(const BlockMatrix &matrix_i) const{

   double ward = 0.0;

   for(int B = 0;B < nr;++B){

      int lapdim = dim[B]*dim[B];
      int inc = 1;

      ward += ddot_(&lapdim,matrix[B][0],&inc,matrix_i.matrix[B][0],&inc); 

   }

   return ward;

}

/**
 * Invert positive semidefinite symmetric matrix is stored in (*this), original matrix (*this) is destroyed
 */
void BlockMatrix::invert(){

   char uplo = 'U';

   for(int B = 0;B < nr;++B){

      int INFO;

      dpotrf_(&uplo,&dim[B],matrix[B][0],&dim[B],&INFO);//cholesky decompositie

      dpotri_(&uplo,&dim[B],matrix[B][0],&dim[B],&INFO);//inverse berekenen

   }

   //terug symmetrisch maken:
   this->symmetrize();

}

/**
 * Scale the matrix (*this) with parameter alpha
 * @param alpha scalefactor
 */
void BlockMatrix::dscal(double alpha){

   for(int B = 0;B < nr;++B){

      int lapdim = dim[B]*dim[B];
      int inc = 1;

      dscal_(&lapdim,&alpha,matrix[B][0],&inc);

   }

}

/**
 * Fill the matrix with random numbers.
 */
void BlockMatrix::fill_Random(){

   for(int B = 0;B < nr;++B)
      for(int i = 0;i < dim[B];++i)
         for(int j = i;j < dim[B];++j)
            matrix[B][j][i] = (double) rand()/RAND_MAX;

   for(int B = 0;B < nr;++B)
      for(int i = 0;i < dim[B];++i)
         for(int j = i + 1;j < dim[B];++j)
            matrix[B][i][j] = matrix[B][j][i];

}

/**
 * Take the square root out of the positive semidefinite matrix, destroys original matrix, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void BlockMatrix::sqrt(int option){

   BlockMatrix hulp(*this);

   BlockVector eigen(hulp);

   if(option == 1)
      for(int B = 0;B < nr;++B)
         for(int i = 0;i < dim[B];++i)
            eigen(B,i) = std::sqrt(eigen(B,i));
   else
      for(int B = 0;B < nr;++B)
         for(int i = 0;i < dim[B];++i)
            eigen(B,i) = 1.0/std::sqrt(eigen(B,i));

   //hulp opslaan
   BlockMatrix hulp_c = hulp;

   //vermenigvuldigen met diagonaalmatrix
   hulp_c.mdiag(eigen);

   for(int B = 0;B < nr;++B){

      //en tenslotte de laatste matrixvermenigvuldiging
      char transA = 'N';
      char transB = 'T';

      double alpha = 1.0;
      double beta = 0.0;

      dgemm_(&transA,&transB,&dim[B],&dim[B],&dim[B],&alpha,hulp_c.matrix[B][0],&dim[B],hulp.matrix[B][0],&dim[B],&beta,matrix[B][0],&dim[B]);

   }

}

/**
 * Multiply this matrix with diagonal matrix
 * @param diag Diagonal matrix to multiply with this, has to be allocated on matrix dimension.
 */
void BlockMatrix::mdiag(const BlockVector &diag){

   int inc = 1;

   for(int B = 0;B < nr;++B){
      for(int i = 0;i < dim[B];++i){

         double scal = diag(B,i);

         dscal_(&dim[B],&scal,matrix[B][i],&inc);

      }
   }

}

/**
 * BlockMatrix product of two general matrices A en B, put result in this
 * @param A left matrix
 * @param B right matrix
 */
BlockMatrix &BlockMatrix::mprod(const BlockMatrix &L,const BlockMatrix &R){

   char trans = 'N';

   double alpha = 1.0;
   double beta = 0.0;

   for(int B = 0;B < nr;++B)
      dgemm_(&trans,&trans,&dim[B],&dim[B],&dim[B],&alpha,L.matrix[B][0],&dim[B],R.matrix[B][0],&dim[B],&beta,matrix[B][0],&dim[B]);

   return *this;

}

/**
 * Copy upper triangle into lower triangle.
 */
void BlockMatrix::symmetrize(){

   for(int B = 0;B < nr;++B)
      for(int i = 0;i < dim[B];++i)
         for(int j = i + 1;j < dim[B];++j)
            matrix[B][i][j] = matrix[B][j][i];

}

ostream &operator<<(ostream &output,const BlockMatrix &matrix_p){

   for(int B = 0;B < matrix_p.gnr();++B)
      for(int i = 0;i < matrix_p.gdim(B);++i)
         for(int j = 0;j < matrix_p.gdim(B);++j)
            output << i << "\t" << j << "\t" << matrix_p(B,i,j) << endl;

   return output;

}

/**
 * commute two symmetric matrices
 * @param A left matrix
 * @param B right matrix
 */
BlockMatrix &BlockMatrix::commute(const BlockMatrix &L,const BlockMatrix &R){

   for(int B = 0;B < nr;++B){

      for(int i = 0;i < dim[B];++i)
         for(int k = 0;k < dim[B];++k){

            (*this)(B,i,k) = 0.0;

            for(int j = 0;j < dim[B];++j)
               (*this)(B,i,k) += L(B,i,j) * R(B,j,k)  - R(B,i,j) * L(B,j,k);

         }

   }

   return *this;

}

/* vim: set ts=3 sw=3 expandtab :*/

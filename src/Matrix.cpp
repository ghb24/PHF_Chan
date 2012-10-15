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
 * @param n dimension of the matrix
 */
Matrix::Matrix(int n){

   this->n = n;

   matrix = new double * [n];
   matrix[0] = new double [n*n];

   for(int i = 1;i < n;++i)
      matrix[i] = matrix[i - 1] + n;

}

/**
 * copy constructor 
 * @param mat_copy The matrix you want to be copied into the object you are constructing
 */
Matrix::Matrix(const Matrix &mat_copy){

   this->n = mat_copy.n;

   matrix = new double * [n];
   matrix[0] = new double [n*n];

   for(int i = 1;i < n;++i)
      matrix[i] = matrix[i - 1] + n;

   int dim = n*n;
   int incx = 1;
   int incy = 1;

   dcopy_(&dim,mat_copy.matrix[0],&incx,matrix[0],&incy);

}

/**
 * construct from file: matrix is allocated and is filled with number from the file "filename"
 * @param filename char containing the name of the input file
 */
Matrix::Matrix(const char *filename){

   ifstream input(filename);

   input >> this->n;

   matrix = new double * [n];
   matrix[0] = new double [n*n];

   for(int i = 1;i < n;++i)
      matrix[i] = matrix[i - 1] + n;

   int I,J;

   for(int i = 0;i < n;++i)
      for(int j = 0;j < n;++j)
         input >> I >> J >> matrix[j][i];

}

/**
 * Destructor
 */
Matrix::~Matrix(){

   delete [] matrix[0];
   delete [] matrix;

}

/**
 * overload the equality operator
 * @param matrix_copy The matrix you want to be copied into this
 */
Matrix &Matrix::operator=(const Matrix &matrix_copy){

   int dim = n*n;
   int incx = 1;
   int incy = 1;

   dcopy_(&dim,matrix_copy.matrix[0],&incx,matrix[0],&incy);

   return *this;

}

/**
 * Make all the number in your matrix equal to the number a, e.g. usefull for initialization (Matrix M = 0)
 * @param a the number
 */
Matrix &Matrix::operator=(double a){

   for(int i = 0;i < n;++i)
      for(int j = 0;j < n;++j)
         matrix[j][i] = a;

   return *this;

}

/**
 * overload the += operator for matrices
 * @param matrix_pl The matrix you want to add to this
 */
Matrix &Matrix::operator+=(const Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;
   double alpha = 1.0;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}

/**
 * overload the -= operator for matrices
 * @param matrix_pl The matrix you want to deduct from this
 */
Matrix &Matrix::operator-=(const Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;
   double alpha = -1.0;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}

/**
 * add the matrix matrix_pl times the constant alpha to this
 * @param alpha the constant to multiply the matrix_pl with
 * @param matrix_pl the Matrix to be multiplied by alpha and added to this
 */
Matrix &Matrix::daxpy(double alpha,const Matrix &matrix_pl){

   int dim = n*n;
   int inc = 1;

   daxpy_(&dim,&alpha,matrix_pl.matrix[0],&inc,matrix[0],&inc);

   return *this;

}
/**
 * /= operator overloaded: divide by a constant
 * @param c the number to divide your matrix through
 */
Matrix &Matrix::operator/=(double c){

   int dim = n*n;
   int inc = 1;

   double alpha = 1.0/c;

   dscal_(&dim,&alpha,matrix[0],&inc);

   return *this;

}

/**
 * write access to your matrix, change the number on row i and column j
 * remark that for the conversion to lapack functions the double pointer is transposed!
 * @param i row number
 * @param j column number
 * @return the entry on place i,j
 */
double &Matrix::operator()(int i,int j){

   return matrix[j][i];

}

/**
 * read access to your matrix, view the number on row i and column j
 * remark that for the conversion to lapack functions the double pointer is transposed!
 * @param i row number
 * @param j column number
 * @return the entry on place i,j
 */
double Matrix::operator()(int i,int j) const {

   return matrix[j][i];

}

/**
 * @return the underlying pointer to matrix, useful for mkl applications
 */
double **Matrix::gMatrix(){

   return matrix;

}

/**
 * @return the dimension of the matrix
 */
int Matrix::gn() const{

   return n;

}

/**
 * @return the trace of the matrix:
 */
double Matrix::trace() const{

   double ward = 0;

   for(int i = 0;i < n;++i)
      ward += matrix[i][i];

   return ward;

}

/**
 * @return inproduct of (*this) matrix with matrix_i, defined as Tr (A B)
 * @param matrix_i input matrix
 */
double Matrix::ddot(const Matrix &matrix_i) const{

   int dim = n*n;
   int inc = 1;

   return ddot_(&dim,matrix[0],&inc,matrix_i.matrix[0],&inc);

}

/**
 * Invert positive semidefinite symmetric matrix is stored in (*this), original matrix (*this) is destroyed
 */
void Matrix::invert(){

   char uplo = 'U';

   int INFO;

   dpotrf_(&uplo,&n,matrix[0],&n,&INFO);//cholesky decompositie

   dpotri_(&uplo,&n,matrix[0],&n,&INFO);//inverse berekenen

   //terug symmetrisch maken:
   this->symmetrize();

}

/**
 * Scale the matrix (*this) with parameter alpha
 * @param alpha scalefactor
 */
void Matrix::dscal(double alpha){

   int dim = n*n;
   int inc = 1;

   dscal_(&dim,&alpha,matrix[0],&inc);

}

/**
 * Fill the matrix with random numbers.
 */
void Matrix::fill_Random(){

   for(int i = 0;i < n;++i)
      for(int j = i;j < n;++j)
         matrix[j][i] = (double) rand()/RAND_MAX;

   for(int i = 0;i < n;++i)
      for(int j = i + 1;j < n;++j)
         matrix[i][j] = matrix[j][i];

}

/**
 * Take the square root out of the positive semidefinite matrix, destroys original matrix, square root will be put in (*this)
 * @param option = 1, positive square root, = -1, negative square root.
 */
void Matrix::sqrt(int option){

   Matrix hulp(*this);

   Vector eigen(hulp);

   if(option == 1)
      for(int i = 0;i < n;++i)
         eigen[i] = std::sqrt(eigen[i]);
   else
      for(int i = 0;i < n;++i)
         eigen[i] = 1.0/std::sqrt(eigen[i]);

   //hulp opslaan
   Matrix hulp_c = hulp;

   //vermenigvuldigen met diagonaalmatrix
   hulp_c.mdiag(eigen);

   //en tenslotte de laatste matrixvermenigvuldiging
   char transA = 'N';
   char transB = 'T';

   double alpha = 1.0;
   double beta = 0.0;

   dgemm_(&transA,&transB,&n,&n,&n,&alpha,hulp_c.matrix[0],&n,hulp.matrix[0],&n,&beta,matrix[0],&n);

}

/**
 * Multiply this matrix with diagonal matrix
 * @param diag Diagonal matrix to multiply with this, has to be allocated on matrix dimension.
 */
void Matrix::mdiag(const Vector &diag){

   int inc = 1;

   for(int i = 0;i < n;++i){

      double scal = diag[i];

      dscal_(&n,&scal,matrix[i],&inc);

   }

}

/**
 * Multiply symmetric matrix object left en right with symmetric matrix map to 
 * form another symmetric matrix and put it in (*this): this = map*object*map
 * @param map matrix that will be multiplied to the left en to the right of matrix object
 * @param object central matrix
 */
void Matrix::L_map(const Matrix &map,const Matrix &object){
   
   char side = 'L';
   char uplo = 'U';

   double alpha = 1.0;
   double beta = 0.0;

   double *hulp = new double [n*n];

   dsymm_(&side,&uplo,&n,&n,&alpha,map.matrix[0],&n,object.matrix[0],&n,&beta,hulp,&n);

   side = 'R';

   dsymm_(&side,&uplo,&n,&n,&alpha,map.matrix[0],&n,hulp,&n,&beta,matrix[0],&n);

   delete [] hulp;

   //expliciet symmetriseren van de uit matrix
   this->symmetrize();

}

/**
 * Matrix product of two general matrices A en B, put result in this
 * @param A left matrix
 * @param B right matrix
 */
Matrix &Matrix::mprod(const Matrix &A,const Matrix &B){

   char trans = 'N';

   double alpha = 1.0;
   double beta = 0.0;

   dgemm_(&trans,&trans,&n,&n,&n,&alpha,A.matrix[0],&n,B.matrix[0],&n,&beta,matrix[0],&n);

   return *this;

}

/**
 * Copy upper triangle into lower triangle.
 */
void Matrix::symmetrize(){

   for(int i = 0;i < n;++i)
      for(int j = i + 1;j < n;++j)
         matrix[i][j] = matrix[j][i];

}

ostream &operator<<(ostream &output,const Matrix &matrix_p){

   for(int i = 0;i < matrix_p.gn();++i)
      for(int j = 0;j < matrix_p.gn();++j)
         output << i << "\t" << j << "\t" << matrix_p(i,j) << endl;

   return output;

}

/**
 * print the matrix in a file with name and location filename
 * @param filename char with name and location
 */
void Matrix::out(const char *filename) const{

   ofstream output(filename);
   output.precision(10);

   for(int i = 0;i < n;++i)
      for(int j = 0;j < n;++j)
         output << i << "\t" << j << "\t" << matrix[j][i] << endl;

}

/**
 * Seperate matrix into two matrices, a positive and negative semidefinite part.
 * @param p positive (plus) output part
 * @param m negative (minus) output part
 */
void Matrix::sep_pm(Matrix &p,Matrix &m){

   //init:
#ifdef __INTEL_COMPILER
   p = 0;
   m = *this;
#else
   p = *this;
   m = 0;
#endif

   double *eigenvalues = new double [n];

   //diagonalize orignal matrix:
   char jobz = 'V';
   char uplo = 'U';

   int info;

#ifdef SYEVD
   int lwork = 2*n*n+6*n+1;
   int liwork = 5*n+3;

   double *work = new double [lwork];
   int *iwork = new int [liwork];

   dsyevd_(&jobz,&uplo,&n,matrix[0],&n,eigenvalues,work,&lwork,iwork,&liwork,&info);

   delete [] iwork;
#else
   int lwork = 3*n - 1;

   double *work = new double [lwork];

   dsyev_(&jobz,&uplo,&n,matrix[0],&n,eigenvalues,work,&lwork,&info);

#endif

   delete [] work;

#ifdef __INTEL_COMPILER
   Matrix copy(*this);

   if( eigenvalues[0] >= 0 )
   {
      p = m;
      m = 0;
      delete [] eigenvalues;
      return;
   }

   if( eigenvalues[n-1] < 0)
   {
      delete [] eigenvalues;
      return;
   }

   int i = 0;

   while(i < n && eigenvalues[i] < 0.0)
      eigenvalues[i++] = 0;

   int inc = 1;

   for(i = 0;i < n;++i)
      dscal_(&n,&eigenvalues[i],matrix[i],&inc);

   char transA = 'N';
   char transB = 'T';

   double alpha = 1.0;
   double beta = 0.0;

   dgemm_(&transA,&transB,&n,&n,&n,&alpha,matrix[0],&n,copy.matrix[0],&n,&beta,p.matrix[0],&n);

   m -= p;
#else
   if( eigenvalues[0] >= 0 )
   {
      delete [] eigenvalues;
      return;
   }

   if( eigenvalues[n-1] < 0)
   {
      m = p;
      p = 0;
      delete [] eigenvalues;
      return;
   }

   //fill the plus and minus matrix
   int i = 0;

   while(i < n && eigenvalues[i] < 0.0){

      for(int j = 0;j < n;++j)
         for(int k = j;k < n;++k)
            m(j,k) += eigenvalues[i] * matrix[i][j] * matrix[i][k];

      ++i;

   }

   m.symmetrize();

   p -= m;
#endif

   delete [] eigenvalues;

}

/**
 * commute two symmetric matrices
 * @param A left matrix
 * @param B right matrix
 */
Matrix &Matrix::commute(const Matrix &A,const Matrix &B){

   for(int i = 0;i < n;++i)
      for(int k = 0;k < n;++k){

         (*this)(i,k) = 0.0;

         for(int j = 0;j < n;++j)
            (*this)(i,k) += A(i,j) * B(j,k)  - B(i,j) * A(j,k);

      }

   return *this;

}

/**
 * solve the linear system of equations: Ax = b. (*this) is the A matrix defining the left hand side
 * @param b is the right hand side of the equation
 */
void Matrix::solve(double *b){

   //variabels needed for the lapack call
   char trans = 'T';
   int lda = n;
   int lda2 = n;
   int nrhs = 1;
   int ldb = n;
   int INFO;
   int rows = n;
   int colums = n;

   int *ipiv;
   ipiv = new int[n];

   dgetrf_(&rows,&colums,matrix[0],&lda2,ipiv,&INFO);

   dgetrs_(&trans,&n,&nrhs,matrix[0],&lda,ipiv,b,&ldb,&INFO);

   delete [] ipiv;


}

/* vim: set ts=3 sw=3 expandtab :*/

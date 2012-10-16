#ifndef LAPACK_H
#define LAPACK_H

extern "C" {

   void dgeqrf_(int *m,int *n,double *A,int *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dorgqr_(int *m,int *n, int *k,double *A,int  *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dgelqf_(int *m,int *n,double *A,int *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dorglq_(int *m,int *n, int *k,double *A,int  *LDA,double *tau,double *WORK,int *LWORK,int *INFO);
   void dcopy_(int *n,double *x,int *incx,double *y,int *incy);
   void daxpy_(int *n,double *alpha,double *x,int *incx,double *y,int *incy);
   void dscal_(int *n,double *alpha,double *x,int *incx);
   void dgemm_(char *transA,char *transB,int *m,int *n,int *k,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);
   void dsymm_(char *side,char *uplo,int *m,int *n,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);
   void dgemv_(char *trans,int *m,int *n,double *alpha,double *A,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
   void dsymv_(char *uplo, int *n, double *alpha, double *A, int *lda, double *X, int *incx, double *beta, double *Y, int *incy);
   double ddot_(int *n,double *x,int *incx,double *y,int *incy);
   void dsyev_(char *jobz,char *uplo,int *n,double *A,int *lda,double *W,double *work,int *lwork,int *info);
   void dpotrf_(char *uplo,int *n,double *A,int *lda,int *INFO);
   void dpotri_(char *uplo,int *n,double *A,int *lda,int *INFO);
   void dgesdd_(char* JOBZ, int* M, int* N, double* A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* IWORK, int* INFO);
   void dsygvd_(int *itype, char *jobz, char *uplo, int *N, double *A, int *lda, double *B, int *ldb, double *W, double *work, int *lwork, int *iwork, int *liwork, int *info);
   void dsygv_(int *itype, char *jobz, char *uplo, int *N, double *A, int *lda, double *B, int *ldb, double *W, double *work, int *lwork, int *info);
   void dsytrf_(char *uplo, int *n, double *A, int *lda, int *ipiv, double *work, int *lwork, int *info);
   void dsycon_(char *uplo, int *n, double *A, int *lda, int *ipiv, double *anorm, double *rcond, double *work, int *iwork, int *info);
   double dlansy_(char *norm, char *uplo, int *N, double *A, int *lda, double *work);
   void dsyevd_(char *jobz, char *uplo, int *N, double *G, int *lda, double *eigs, double *work, int *lwork, int *iwork, int *liwork, int *info);
   double dnrm2_(int *N, double *diff, int *inc1);

}

#endif

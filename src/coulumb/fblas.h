/*
 * File: fblas.h
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#if defined __cplusplus
extern "C" {
#endif

void dscal_(const int *n, const double *da, double *dx, const int *incx);
void daxpy_(const int *n, const double *da, const double *dx,
            const int *incx, double *dy, const int *incy);
double ddot_(const int *n, const double *dx, const int *incx,
             const double *dy, const int *incy);
void dgemm_(const char*, const char*, const int*, const int*, const int*,
            const double*, const double*, const int*,
            const double*, const int*, const double*, double*, const int*);
void dgemv_(const char*, const int*, const int*,
            const double*, const double*, const int*,
            const double*, const int*, const double*, double*, const int*);
void dger_(const int *m, const int *n, const double *alpha, const double *x,
           const int *incx, const double *y, const int *incy,
           double *a, const int *lda);
//void dsyrk_
void zgerc_(const int *m, const int *n, const double *alpha, const double *x,
            const int *incx, const double *y, const int *incy,
            double *a, const int *lda);
void zgemv_(const char*, const int*, const int*,
            const double*, const double*, const int*,
            const double*, const int*, const double*, double*, const int*);

#if defined __cplusplus
} // end extern "C"
#endif

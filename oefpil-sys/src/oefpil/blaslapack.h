#ifndef BLASLAPACK_H
#define BLASLAPACK_H

#include <stddef.h> /* size_t */

/* TODO: figure out for other compilers (e.g. ifort + icc/msvc - hidden argument type, position (another hiddden arguments present?) */

/* check LAPACKE/include/lapack.h */
/* check name mangling / length parameter with other compilers, move add respective #defines */

/* note that clang also defines __GNUC__ */


// #if defined(__clang__)

/* TODO */

// #elif defined(__GNUC__) || defined(__GNUG__)
#if defined(__GNUC__) || defined(__GNUG__)

#if (__GNUC__ > 7)
typedef size_t f_char_len;
#else
typedef int f_char_len;
#endif


/* BLAS */

extern double ddot_(const int *n, const double *dx, const int *incx,  const double *dy, const int *incy);

extern void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                   const double *alpha, const double *a, const int *lda,
                   const double *b, const int *ldb,
                   const double *beta, double *c, const int *ldc,
                   f_char_len ltransa, f_char_len ltransb);

extern void dgemv_(const char *transa, const int *m, const int *n,
                   const double *alpha, const double *a, const int *lda,
                   const double *x, const int *incx,
                   const double *beta, double *y, const int *incy,
                   f_char_len ltransa);

extern void dtrmv_(const char *uplo, const char *transa, const char *diag,
                   const int *n, const double *a, const int *lda, 
                   double *x, const int *incx);

extern void dsymm_(const char *side, const char *uplo, const int *m, const int *n,
                   const double *alpha, const double *a, const int *lda,
                   const double *b, const int *ldb,
                   const double *beta, double *c, const int *ldc,
                   f_char_len lside, f_char_len luplo);

extern void dsymv_(const char *uplo, const int *n,
                   const double *alpha, const double *a, const int *lda,
                   const double *x, const int *incx,
                   const double *beta, double *y, const int *incy,
                   f_char_len luplo);

extern void dsyrk_(const char *uplo, const char *trans, const int *n, const int *k,
                   const double *alpha, const double *a, const int *lda,
                   const double *beta, double *c, const int *ldc,
                   f_char_len luplo, f_char_len ltrans);

/* LAPACK */

extern void dgeqr_(const int *m, const int *n, double *a, const int *lda, double *t, const int *tsize,
		   double *work, const int *lwork, int *info);

extern void dgemqr_(const char *side, const char *trans, const int *m, const int *n, const int *k,
		    const double *a, const int *lda, const double *t, const int *tsize,
		    double *c, const int *ldc, double *work, const int *lwork, int *info,
		    f_char_len lside, f_char_len ltrans);

extern void dpotrf_(const char *uplo, const int *n,
                    double *a, const int *lda,
                    int *info,
                    f_char_len luplo);

extern void dpotri_(const char *uplo, const int *n,
                    double *a, const int *lda,
                    int *info,
                    f_char_len luplo);

extern void dtrsm_(const char *side, const char *uplo, const char *transa, const char *diag, const int *m, const int *n,
                   const double *alpha, const double *a, const int *lda,
                   double *b, const int *ldb,
                   f_char_len lside, f_char_len luplo, f_char_len ltransa, f_char_len ldiag);

extern void dtrsv_(const char *uplo, const char *trans, const char *diag,
                   const int *n, const double *a, const int *lda,
                   double *x, const int *incx,
                   f_char_len luplo, f_char_len ltrans, f_char_len ldiag);


#elif defined(_MSC_VER)

typedef size_t f_char_len;

/* BLAS */

extern double DDOT(const int *n, const double *dx, const int *incx,  const double *dy, const int *incy);

extern void DGEMM(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                  const double *alpha, const double *a, const int *lda,
                  const double *b, const int *ldb,
                  const double *beta, double *c, const int *ldc,
                  f_char_len ltransa, f_char_len ltransb);

extern void DGEMV(const char *transa, const int *m, const int *n,
                  const double *alpha, const double *a, const int *lda,
                  const double *x, const int *incx,
                  const double *beta, double *y, const int *incy,
                  f_char_len ltransa);

extern void DTRMV(const char *uplo, const char *transa, const char *diag,
                  const int *n, const double *a, const int *lda, 
                  double *x, const int *incx);

extern void DSYMM(const char *side, const char *uplo, const int *m, const int *n,
                  const double *alpha, const double *a, const int *lda,
                  const double *b, const int *ldb,
                  const double *beta, double *c, const int *ldc,
                  f_char_len lside, f_char_len luplo);

extern void DSYMV(const char *uplo, const int *n,
                  const double *alpha, const double *a, const int *lda,
                  const double *x, const int *incx,
                  const double *beta, double *y, const int *incy,
                  f_char_len luplo);

extern void DSYRK(const char *uplo, const char *trans, const int *n, const int *k,
                  const double *alpha, const double *a, const int *lda,
                  const double *beta, double *c, const int *ldc,
                  f_char_len luplo, f_char_len ltrans);

/* LAPACK */

extern void DGEQR(const int *m, const int *n, double *a, const int *lda, double *t, const int *tsize,
                  double *work, const int *lwork, int *info);

extern void DGEMQR(const char *side, const char *trans, const int *m, const int *n, const int *k,
                   const double *a, const int *lda, const double *t, const int *tsize,
                   double *c, const int *ldc, double *work, const int *lwork, int *info,
                   f_char_len lside, f_char_len ltrans);

extern void DPOTRF(const char *uplo, const int *n,
                   double *a, const int *lda,
                   int *info,
                   f_char_len luplo);

extern void DPOTRI(const char *uplo, const int *n,
                   double *a, const int *lda,
                   int *info,
                   f_char_len luplo);

extern void DTRSM(const char *side, const char *uplo, const char *transa, const char *diag, const int *m, const int *n,
                  const double *alpha, const double *a, const int *lda,
                  double *b, const int *ldb,
                  f_char_len lside, f_char_len luplo, f_char_len ltransa, f_char_len ldiag);

extern void DTRSV(const char *uplo, const char *trans, const char *diag,
                  const int *n, const double *a, const int *lda,
                  double *x, const int *incx,
                  f_char_len luplo, f_char_len ltrans, f_char_len ldiag);

#endif


/* declarations for C */

double ddot(int n, const double *dx, int incx,  const double *dy, int incy);

void dgemm(const char *transa, const char *transb, int m, int n, int k,
           double alpha, const double *a, int lda,
           const double *b, int ldb,
           double beta, double *c, int ldc);

void dgemv(const char *transa, int m, int n,
           double alpha, const double *a, int lda,
           const double *x, int incx,
           double beta, double *y, int incy);

void dtrmv(const char *uplo, const char *transa, const char *diag,
           int n, const double *a, int lda, double *x, int incx);

void dsymm(const char *side, const char *uplo, int m, int n,
           double alpha, const double *a, int lda,
           const double *b, int ldb,
           double beta, double *c, int ldc);

void dsymv(const char *uplo, int n,
           double alpha, const double *a, int lda,
           const double *x, int incx,
           double beta, double *y, int incy);

void dsyrk(const char *uplo, const char *trans, int n, int k,
           double alpha, const double *a, int lda,
           double beta, double *c, int ldc);

void dgesvd(const char *jobu, const char *jobvt, int m, int n,
            double *a, int lda,
            double *s, double *u, int ldu, double *vt, int ldvt,
            double *work, int lwork,
            int *info);

void dgesvdq(const char *joba, const char *jobp, const char *jobr, const char *jobu, const char *jobv,
             int m, int n, double *a, int lda,
             double *s, double *u, int ldu, double *v, int ldv,
             int *numrank,
             int *iwork, int liwork, double *work, int lwork, double *rwork, int lrwork,
             int *info);

void dgeqr(int m, int n, double *a, int lda, double *t, int tsize,
           double *work, int lwork, int *info);

void dgemqr(const char *side, const char *trans, int m, int n, int k,
            const double *a, int lda, const double *t, int tsize,
            double *c, int ldc, double *work, int lwork, int *info);

void dpotrf(const char *uplo, int n, double *a, int lda, int *info);

void dpotri(const char *uplo, int n, double *a, int lda, int *info);

void dtrsm(const char *side, const char *uplo, const char *transa, const char *diag, int m, int n,
           double alpha, const double *a, int lda,
           double *b, int ldb);

void dtrsv(const char *uplo, const char *trans, const char *diag,
           int n, const double *a, int lda,
           double *x, int incx);

#endif


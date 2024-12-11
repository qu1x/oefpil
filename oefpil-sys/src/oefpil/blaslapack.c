#include "blaslapack.h"

#include <string.h>

/* wrapper functions to abstract compiler-dependent name mangling and passing string length parameters */
/* and to allow passing arguments by value when possible  */

/* note that clang also defines __GNUC__ */

///// #if defined(__clang__)

/* TODO */

//// #elif defined(__GNUC__) || defined(__GNUG__)


/* BLAS */

double ddot(int n, const double *dx, int incx,  const double *dy, int incy)
{
#if defined(__GNUC__) || defined(__GNUG__)
    return ddot_(&n, dx, &incx, dy, &incy);
#elif defined(_MSC_VER)
    return DDOT(&n, dx, &incx, dy, &incy);
#endif
}

void dgemm(const char *transa, const char *transb, int m, int n, int k,
           double alpha, const double *a, int lda,
           const double *b, int ldb,
           double beta, double *c, int ldc)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc,
           strlen(transa), strlen(transb));
#elif defined(_MSC_VER)
    DGEMM(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc,
          strlen(transa), strlen(transb));
#endif
}

void dgemv(const char *transa, int m, int n,
           double alpha, const double *a, int lda,
           const double *x, int incx,
           double beta, double *y, int incy)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dgemv_(transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy,
           strlen(transa));
#elif defined(_MSC_VER)
    DGEMV(transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy,
          strlen(transa));
#endif
}

void dtrmv(const char *uplo, const char *transa, const char *diag,
           int n, const double *a, int lda, double *x, int incx)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dtrmv_(uplo, transa, diag, &n, a, &lda, x, &incx);
#elif defined(_MSC_VER)
    DTRMV(uplo, transa, diag, &n, a, &lda, x, &incx);
#endif
}

void dsymm(const char *side, const char *uplo, int m, int n,
           double alpha, const double *a, int lda,
           const double *b, int ldb,
           double beta, double *c, int ldc)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dsymm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc,
           strlen(side), strlen(uplo));
#elif defined(_MSC_VER)
    DSYMM(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc,
          strlen(side), strlen(uplo));
#endif
}

void dsymv(const char *uplo, int n,
           double alpha, const double *a, int lda,
           const double *x, int incx,
           double beta, double *y, int incy)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dsymv_(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy,
           strlen(uplo));
#elif defined(_MSC_VER)
    DSYMV(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy,
          strlen(uplo));
#endif
}

void dsyrk(const char *uplo, const char *trans, int n, int k,
           double alpha, const double *a, int lda,
           double beta, double *c, int ldc)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dsyrk_(uplo, trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc,
           strlen(uplo), strlen(trans));
#elif defined(_MSC_VER)
    DSYRK(uplo, trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc,
          strlen(uplo), strlen(trans));
#endif
}

/* LAPACK */

void dgeqr(int m, int n, double *a, int lda, double *t, int tsize,
           double *work, int lwork, int *info)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dgeqr_(&m, &n, a, &lda, t, &tsize, work, &lwork, info);
#elif defined(_MSC_VER)
    DGEQR(&m, &n, a, &lda, t, &tsize, work, &lwork, info);
#endif
}

void dgemqr(const char *side, const char *trans, int m, int n, int k,
            const double *a, int lda, const double *t, int tsize,
            double *c, int ldc, double *work, int lwork, int *info)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dgemqr_(side, trans, &m, &n, &k, a, &lda, t, &tsize, c, &ldc, work, &lwork, info,
	    strlen(side), strlen(trans));
#elif defined(_MSC_VER)
    DGEMQR(side, trans, &m, &n, &k, a, &lda, t, &tsize, c, &ldc, work, &lwork, info,
           strlen(side), strlen(trans));
#endif
}

void dpotrf(const char *uplo, int n, double *a, int lda, int *info)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dpotrf_(uplo, &n, a, &lda, info,
            strlen(uplo));
#elif defined(_MSC_VER)
    DPOTRF(uplo, &n, a, &lda, info,
           strlen(uplo));
#endif
}

void dpotri(const char *uplo, int n, double *a, int lda, int *info)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dpotri_(uplo, &n, a, &lda, info,
            strlen(uplo));
#elif defined(_MSC_VER)
    DPOTRI(uplo, &n, a, &lda, info,
           strlen(uplo));
#endif
}

void dtrsm(const char *side, const char *uplo, const char *transa, const char *diag, int m, int n,
           double alpha, const double *a, int lda,
           double *b, int ldb)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dtrsm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb,
           strlen(side), strlen(uplo), strlen(transa), strlen(diag));
#elif defined(_MSC_VER)
    DTRSM(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb,
          strlen(side), strlen(uplo), strlen(transa), strlen(diag));
#endif
}

void dtrsv(const char *uplo, const char *trans, const char *diag,
           int n, const double *a, int lda,
           double *x, int incx)
{
#if defined(__GNUC__) || defined(__GNUG__)
    dtrsv_(uplo, trans, diag, &n, a, &lda, x, &incx,
           strlen(uplo), strlen(trans), strlen(diag));
#elif defined(_MSC_VER)
    DTRSV(uplo, trans, diag, &n, a, &lda, x, &incx,
          strlen(uplo), strlen(trans), strlen(diag));
#endif
}

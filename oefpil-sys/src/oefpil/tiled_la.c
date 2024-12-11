#include "tiled_la.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h> /* until debug printfs are removed */

#include "blaslapack.h"


int *tilemap_n_new(int bn)
{
    if (bn > 0)
        return (int*)malloc(bn * sizeof(int));
    else
        return NULL;
}

int *tilemap_mn_new(int bm, int bn)
{
    if ((bm > 0) && (bn > 0))
        return (int*)malloc(bm*bn * sizeof(int));
    else
        return NULL;
}

void tiled_matrix_clear(int bsize, int nblocks, int *tilemap, double *A)
{
    int i;
    
    for (i = 0; i < bsize*nblocks; i++)
        A[i] = 0.0;
        
    for (i = 0; i < nblocks; i++)
        tilemap[i] = TILE_NULL;
}

double* tiled_matrix_diag_new(int n, int bn, int *tilemap)
{
    double *A;
    
    if ((n <= 0) || (bn <= 0))
        return NULL;

    A = (double*)malloc(n*bn * sizeof(double));
    tiled_matrix_clear(n, bn, tilemap, A);
    return A;
}

double* tiled_matrix_blockdiag_new(int n, int bn, int *tilemap)
{
    double *A;
    
    if ((n <= 0) || (bn <= 0))
        return NULL;

    A = (double*)malloc((n*n)*bn * sizeof(double));
    tiled_matrix_clear(n*n, bn, tilemap, A);
    return A;
}

double* tiled_matrix_diags_new(int n, int bm, int bn, int *tilemap)
{
    double *A;
    
    if ((n <= 0) || (bm <= 0) || (bn <= 0))
        return NULL;

    A = (double*)malloc(n*(bm*bn) * sizeof(double));
    tiled_matrix_clear(n, bm*bn, tilemap, A);
    return A;
}

double* tiled_matrix_full_new(int n, int bm, int bn, int *tilemap)
{
    double *A;
    
    if ((n <= 0) ||  (bm <= 0) || (bn <= 0))
        return NULL;

    A = (double*)malloc((n*n)*(bm*bn) * sizeof(double));
    tiled_matrix_clear(n*n, bm*bn, tilemap, A);
    return A;
}

void tiled_matrix_diag_set_tile_diag(int n, int bn, double *A, int *tilemap,
                                      int ii, const double *B)
{
    int i;
    double *A_;

    if (B) {
        A_ = A + ii*n;

        for (i = 0; i < n; i++)
            A_[i] = B[i];

        tilemap[ii] = TILE_DIAG;
    }
    else {
        tilemap[ii] = TILE_NULL;
    }
}

void tiled_matrix_blockdiag_set_tile_diag(int n, int bn, double *A, int *tilemap,
                                           int ii, const double *B)
{
    int i;
    double *A_;

    if (B) {
        A_ = A + ii*n*n;

        for (i = 0; i < n; i++)
            A_[i*n + i] = B[i];

        tilemap[ii] = TILE_DIAG;
    }
    else {
        tilemap[ii] = TILE_NULL;
    }
}

void tiled_matrix_blockdiag_set_tile_full(int n, int bn, double *A, int *tilemap,
                                           int ii, const double *B, const char *trans)
{
    int i, j;
    double *A_;
    bool transpose;

    if (!trans)
        return;

    if ((trans[0] == 'T') || (trans[0] == 't'))
        transpose = true;
    else
        transpose = false;

    if (B) {
        A_ = A + ii*n*n;

        if (!transpose)
            for (i = 0; i < n*n; i++)
                A_[i] = B[i];
        else
            for (j = 0; j < n; j++)
                for (i = 0; i < n; i++)
                    A_[j*n + i] = B[i*n + j];

        tilemap[ii] = TILE_FULL;
    }
    else {
        tilemap[ii] = TILE_NULL;
    }
}

void tiled_matrix_diags_set_tile_diag(int n, int bm, int bn, double *A, int *tilemap,
                                       int ii, int jj, const double *B)
{
    int i;
    double *A_;
    
    if (B) {
        A_ = A + (jj*bm + ii) * n;
        
        for (i = 0; i < n; i++)
            A_[i] = B[i];
        
        tilemap[jj*bm + ii] = TILE_DIAG;
    }
    else {
        tilemap[jj*bm + ii] = TILE_NULL;
    }
}

void tiled_matrix_full_set_tile_diag(int n, int bm, int bn, double *A, int *tilemap,
                                      int ii, int jj, const double *B)
{
    int i;
    double *A_;
    
    if (B) {
        A_ = A + (jj*bm*n*n + ii*n);

        for (i = 0; i < n; i++)
            A_[i*n*bm + i] = B[i];

        tilemap[jj*bm + ii] = TILE_DIAG;
    }
    else {
        tilemap[jj*bm + ii] = TILE_NULL;
    }
}

void tiled_matrix_full_set_tile_full(int n, int bm, int bn, double *A, int *tilemap,
                                 int ii, int jj, const double *B, const char *trans)
{
    int i, j;
    double *A_;
    bool transpose;

    if (!trans)
        return;

    if ((trans[0] == 'T') || (trans[0] == 't'))
        transpose = true;
    else
        transpose = false;

    if (B) {
        A_ = A + (jj*bm*n*n + ii*n);
//        printf("jj %d ii %d \n", jj, ii);

        if (!transpose)
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++){
                    A_[i*n*bn + j] = B[i*n+j];
 //                   printf("i %d n %d bn %d i*n*bn +j %d  B[%d] %g\n", i, n, bn, i*n*bn+j, i*n+j, B[i*n+j]);
                }
        else
            for (j = 0; j < n; j++)
                for (i = 0; i < n; i++)
                    A_[j*n*bn + i] = B[i*n + j];

        tilemap[jj*bm + ii] = TILE_FULL;
    }
    else {
        tilemap[jj*bm + ii] = TILE_NULL;
    }
}


void dgetmv_diag(int n, int bn,
                 double alpha, const double *a, const int *tma,
                 const double *x,
                 double beta, double *y)
{
    int ii, i;
    const double *a_, *x_;
    double *y_;

    for (ii = 0; ii < bn; ii++)
        if (tma[ii] == TILE_DIAG) {
            a_ = a + ii*n;
            x_ = x + ii*n;
            y_ = y + ii*n;
                
            for (i = 0; i < n; i++)
                y_[i] = beta * y_[i] + alpha * a_[i] * x_[i];
        }
}

void dgetmv_blockdiag(int n, int bn,
                      double alpha, const double *a, const int *tma,
                      const double *x,
                      double beta, double *y)
{
    int ii, i;
    const double *a_, *x_;
    double *y_;

    for (ii = 0; ii < bn; ii++)
        switch (tma[ii]) {
        case TILE_DIAG:
            a_ = a + ii*(n*n);
            x_ = x + ii*n;
            y_ = y + ii*n;

            for (i = 0; i < n; i++)
                y_[i] = beta * y_[i] + alpha * a_[i*n + i] * x_[i];
            break;

        case TILE_FULL:
            a_ = a + ii*(n*n);
            x_ = x + ii*n;
            y_ = y + ii*n;

            dgemv("N", n, n, alpha, a_, n, x_, 1, beta, y_, 1);
            /* dsymv("L", n, alpha, a_, n, x_, &i_1, beta, y_, &i_1); */
            break;

        default:
            break;
        }
}

void dgetmv_diags(int n, int bm, int bn,
                  double alpha, const double *a, const int *tma,
                  const double *x,
                  double beta, double *y)
{
    int ii, jj, i;
    const double *a_, *x_;
    double *y_;

    for (ii = 0; ii < bm; ii++) {
        y_ = y + ii*n;

        for (i = 0; i < n; i++)
            y_[i] *= beta;
            
        for (jj = 0; jj < bn; jj++)
            if (tma[jj*bm + ii] == TILE_DIAG) {
                a_ = a + (jj*bm + ii)*n;
                x_ = x + jj*n;
                    
                for (i = 0; i < n; i++)
                    y_[i] += alpha * a_[i] * x_[i];
            }
    }
}

void dgetmv_full(int n, int bm, int bn,
                 double alpha, const double *a, const int *tma,
                 const double *x,
                 double beta, double *y)
{
    int ii, jj, i;
    const double *a_, *x_;
    double *y_;
    double *atmp;

    atmp = (double *)malloc(n*n*sizeof(double));

    for (ii = 0; ii < bm; ii++) {
   //     printf("ii %d \n", ii);
        y_ = y + ii*n;

        for (i = 0; i < n; i++){
            //printf("y_[%d] %g \n", i, y_[i]);
            y_[i] *= beta;
        }

        for (jj = 0; jj < bn; jj++) {
            a_ = a + (jj*bm*n*n + ii*n);
            x_ = x + jj*n;

/*            for (i = 0; i < n; i++)
                printf("x_[%d] %g \n", i, x_[i]); 

            printf("a_ \n");
            */
            for (int ik =0;ik < n; ik ++){
                for (int il=0; il < n; il++){
             //       printf("%g ", a_[ik*n*bn +il]);
                    atmp[ik*n+il] = a_[ik*n*bn + il];
                }
           // printf("\n");
            }


            switch (tma[jj*bm + ii]) {
            case TILE_DIAG:
                for (i = 0; i < n; i++)
                    y_[i] += alpha * a_[i*n + i] * x_[i];
                break;

            case TILE_FULL:
                dgemv("N", n, n, alpha, atmp, n, x_, 1, beta, y_, 1);
                break;

            default:
                break;
            }
        }
    }
    free(atmp);
}

void dpottrf_diag(int n, int bn, double *a, const int *tma, int *info)
{
    int i;

    for (i = 0; i < bn * n; i++)
        a[i] = sqrt(a[i]);

    *info = 0;
}

void dpottrf_blockdiag(int n, int bn, double *a, const int *tma, int *info)
{
    int ii;
    double *a_;

    for (ii = 0; ii < bn; ii++) {
        a_ = a + ii*(n*n);
        dpotrf("L", n, a_, n, info);
    }
}

/* only for bn = 2 !!! */
void dpottrf_diags(int n, int bn, double *a, const int *tma, int *info)
{
    int i;
    double *a_, *L_;

    if (bn != 2)
        printf("WARNING: dpottrf_diags not implemented for bn != 2" "\n");

    /* L11: A11 = L11*L11^T */
    a_ = a + 0;
    
    for (i = 0; i < n; i++)
        a_[i] = sqrt(a_[i]);

    /* L21: solve(X * L11^T = A21) */
    L_ = a + 0;
    a_ = a + 1*n;
    
    for (i = 0; i < n; i++)
        a_[i] = a_[i] / L_[i];

    /* L22: A22 - L21*L21^T = L22*L22^T */
    L_ = a + 1*n;
    a_ = a + 3*n;

    for (i = 0; i < n; i++)
        a_[i] = sqrt(a_[i] - L_[i]*L_[i]);

    *info = 0;
}

/* only for bn = 2 !!! */
void dpottrf_full(int n, int bn, double *a, const int *tma, int *info)
{
    int ii;
    double *a_, *L_;
    double *atmp1, *atmp2;
    int i,j;


    atmp1 = (double *)malloc(n*n*sizeof(double));
    atmp2 = (double *)malloc(n*n*sizeof(double));
    for (i=0;i< n; i++)
       for (j=0;j<n; j++)
           atmp1[i*n+j] = a[i*2*n+j];

    if (bn != 2)
        printf("WARNING: dpottrf_full not implemented for bn != 2" "\n");

    /* L11: A11 = L11*L11^T */
 //   a_ = a + 0; /* A11 -> L11 */
    
    a_ = atmp1;
    dpotrf("L", n, a_, n, info);
    for (i=0;i< n; i++)
       for (j=0;j<n; j++)
           a[i*2*n+j] = atmp1[i*n+j];

    /* L21: solve(X * L11^T = A21) */
    for (i=0;i< n; i++)
       for (j=0;j<n; j++)
           atmp2[i*n+j] = a[i*2*n+j+n];
    L_ = atmp1;
//    L_ = a + 0; /* L11 */
    //a_ = a + 1*(n*n); /* A21 -> L21 */

    a_ = atmp2;

    dtrsm("R", "L", "T", "N", n, n, 1.0 , L_, n, a_, n);
    for (i=0;i< n; i++)
       for (j=0;j<n; j++)
           a[i*2*n+j+n] = atmp2[i*n+j];

    /* L22: A22 - L21*L21^T = L22*L22^T */
 //   L_ = a + 1*(n*n); /* L21 */
  //  a_ = a + 3*(n*n); /* A22 -> L22 */
      L_ = atmp2;

    for (i=0;i< n; i++)
       for (j=0;j<n; j++)
           atmp1[i*n+j] = a[(i+n)*2*n+j+n];

    a_ = atmp1;

    dsyrk("L", "N", n, n, -1.0, L_, n, 1.0, a_, n);
    dpotrf("L", n, a_, n, info);

    for (i=0;i< n; i++)
       for (j=0;j<n; j++)
           a[(i+n)*2*n+j+n] = atmp1[i*n+j];

    free(atmp1);
    free(atmp2);
}

void dtrtsv_diag(int n, int bn, const double *a, const int *tma, double *b, int *info)
{
    int i;

    for (i = 0; i < bn*n; i++)
        b[i] = b[i] / a[i];

    *info = 0;
}

void dtrtsv_blockdiag(int n, int bn, const double *a, const int *tma, double *b, int *info)
{
    int ii;
    const double *a_;
    double *b_;

    /* xi: solve(Aii * xi = bi) */
    for (ii = 0; ii < bn; ii++) {
        a_ = a + ii*(n*n);
        b_ = b + ii*n;
        dtrsv("L", "N", "N", n, a_, n, b_, 1);
    }

    *info = 0;
}

/* only for bn = 2 !!! */
void dtrtsv_diags(int n, int bn, const double *a, const int *tma, double *b, int *info)
{
    int i, ii, jj;
    const double *a_, *x_;
    double *b_;

    if (bn != 2)
        printf("WARNING: dpotrsv_full not implemented for bn != 2" "\n");

    /* x1: solve(A11*x1 = b1) */
    a_ = a + 0;
    b_ = b + 0;
    
    for (i = 0; i < n; i++)
        b_[i] = b_[i] / a_[i];

    /* b2 = b2 - A21*x1 */
    a_ = a + 1*n;
    x_ = b + 0*n; /* b1 = x1 */
    b_ = b + 1*n;

    for (i = 0; i < n; i++)
	b_[i] = b_[i] - a_[i]*x_[i];

    /* x2: solve(A22*x2 = b2) */
    a_ = a + 3*n;
    b_ = b + 1*n;
    
    for (i = 0; i < n; i++)
        b_[i] = b_[i] / a_[i];

    *info = 0;
}

/* only for bn = 2 !!! */
void dtrtsv_full(int n, int bn, const double *a, const int *tma, double *b, int *info)
{
    int i, ii, jj;
    const double *a_, *x_;
    double *atmp, *btmp;
    double *b_;

    if (bn != 2)
        printf("WARNING: dpotrtsv_diags not implemented for bn != 2" "\n");

    atmp = (double *)malloc(n*n*sizeof(double));
    btmp = (double *)malloc(n*sizeof(double));

    /* x1: solve(A11*x1 = b1) */

    //a_ = a + 0;
    //b_ = b + 0;
    b_ = b;
    for (ii =0;ii < n; ii++)
        for (jj =0; jj<n ; jj++)
                atmp[ii*n+jj] = a[ii*2*n +jj];

    a_ = atmp;
    dtrsv("L", "N", "N", n, a_, n, b_, 1);
    
    /* b2 = b2 - A21*x1 */

    //a_ = a + 1*(n*n);
    //x_ = b + 0*n; /* b1 = x1 */
    //b_ = b + 1*n;

    for (ii =0;ii < n; ii++)
        btmp[ii] = b[ii];

    for (ii =0;ii < n; ii++)
        for (jj =0; jj<n ; jj++)
            atmp[ii*n+jj] = a[ii*2*n +jj+n];

    a_ = atmp;   
    x_ = btmp;
    b_ = b +n;

    dgemv("N", n, n, -1.0 , a_, n, x_, 1, 1.0, b_, 1);

    /* x2: solve(A22*x2 = b2) */
    //a_ = a + 3*(n*n);
   // b_ = b + 1*n;
    for (ii =0;ii < n; ii++)
        for (jj =0; jj<n ; jj++)
            atmp[ii*n+jj] = a[(ii+n)*2*n +jj+n];
    a_ = atmp;        
    b_ = b+n;        
    dtrsv("L", "N", "N", n, a_, n, b_, 1);

    *info = 0;
    free(atmp);
    free(btmp);
}

/* symmetric tiled matrix X tiled vector */
void dsytmtv(const char *uplo, int n, int bn,
             double alpha, const double *a, const int *tma,
             const double *x,
             double beta, double *y)
{
    int ii, jj, i;
    const double *a_, *x_;
    double *y_;

    if (!uplo)
        return;

    if ((uplo[0] == 'l') || (uplo[0] == 'L')) /* L */
        for (ii = 0; ii < bn; ii++) {
            y_ = y + n*(bn*n) + ii*n;

            for (i = 0; i < n*n; i++)
                y_[i] = y_[i] * beta;
        
            for (jj = 0; jj < bn; jj++)
		switch (tma[jj*bn + ii]) {
		case TILE_DIAG:
		    /* TODO */
		    break;

		case TILE_FULL:
		    x_ = x + n*(bn*n) + jj*n;

                    if (ii > jj) { /* lower tr. */
                        a_ = a + (jj*n)*(bn*n) + ii*n;
                        dgemm("N", "N", n, n, n, alpha, a_, bn*n, x_, bn*n, 1.0, y_, bn*n); //////// ????????????? ////////////
                    }
                    else if (ii < jj) { /* upper tr. */
                        a_ = a + (ii*n)*(bn*n) + jj*n;
                        dgemm("T", "N", n, n, n, alpha, a_, bn*n, x_, 1, 1.0, y_, bn*n);
                    }
                    else { /* diag */
                        a_ = a + (ii*n)*(bn*n) + ii*n;
                        dsymm("L", "L", bn*n, n, alpha, a_, bn*n, x_, bn*n, 1.0, y_, bn*n); ////// ??????? ///////
                    }
		    break;

		case TILE_NULL:
		default:
		    break;
		}
	}
    else /* U */
        for (ii = 0; ii < bn; ii++) {
            y_ = y + n*(bn*n) + ii*n;

            for (i = 0; i < n*n; i++)
                y_[i] = y_[i] * beta;
            
            for (jj = 0; jj < bn; jj++)
		switch (tma[jj*bn + ii]) {
		case TILE_DIAG:
		    /* TODO */
		    break;

		case TILE_FULL:
		    x_ = x + n*(bn*n) + jj*n;

                    if (ii > jj) { /* lower tr. */
                        a_ = a + (ii*n)*(bn*n) + jj*n; 
                        dgemm("T", "N", n, n, n, alpha, a_, bn*n, x_, 1, 1.0, y_, bn*n);
                    }
                    else if (ii < jj) { /* upper tr. */
                        a_ = a + (jj*n)*(bn*n) + ii*n;
                        dgemm("N", "N", n, n, n, alpha, a_, bn*n, x_, bn*n, 1.0, y_, bn*n); ///// ???????? //////
                    }
                    else { /* diag */
                        a_ = a + (ii*n)*(bn*n) + ii*n;
                        dsymm("L", "L", bn*n, n, alpha, a_, bn*n, x_, bn*n, 1.0, y_, bn*n); ///// ??????? //////
                    }
		    break;

		case TILE_NULL:
		default:
		    break;
		}
	}
}

void dsymmt(const char *side, const char *uplo, const int *m, const int *n, const int *ts,
            const double *alpha, const double *a,
            const double *b,
            const double *beta, double *c)
{
}

void dpotrft(const char *uplo, const int *n, const int *ts,
             double *a,
             int *info)
{
}

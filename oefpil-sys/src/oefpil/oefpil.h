#ifndef OEFPIL_H
#define OEFPIL_H

#include "tiled_la.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

typedef double real;

enum function_mode {
    FUNCTION_EXPLICIT,
    FUNCTION_IMPLICIT
};
enum convergence {
    REL_PAR =0,
    REL_PAR_OR_ABS_PAR,
    REL_X ,
    REL_X_OR_ABS_X,
    REL_PAR_AND_REL_X,
    REL_PAR_OR_ABS_PAR_AND_REL_X_OR_ABS_X,
    CHI2
};



#define PRINT_FORMAT_E3 " % 10.3e"
#define PRINT_FORMAT_E8 " % 10.8e"
#define PRINT_FORMAT_E15 " % 10.15e"
#define PRINT_FORMAT_F8 " % 10.8f"

enum oefpil_status {
    OEFPIL_STATUS_OK,
    OEFPIL_CONVERGED,
    OEFPIL_MAX_ITER_REACHED,
    OEFPIL_ERROR_NAN,
    OEFPIL_ERROR_INF,
    OEFPIL_ERROR_LAPACK,
    OEFPIL_STATUS_UNKNOWN
};

/* mininum sizes: x: n*nv, p: np, fx: n, dfdx: n*nv, dfdp: n*np */
typedef void fitfunc_oefpil(void *data, int n, const double *x, const double *p, double *fx, double *dfdx, double *dfdp);


void printmatrix(FILE *rptfile, int m, int n, const double *a, const char *format, const char *label, int mode);

int* oefpil_tilemap_diagtiles_new(int bn);
int* oefpil_tilemap_alltiles_new(int bn);

/* diag */
double* oefpil_tcm_diag_new(int n, int bn, int *tilemap);
void oefpil_tcm_diag_set_tile_diag(int n, int bn, double *A, int *tilemap, int ii, const double *B);

/* blockdiag */
double* oefpil_tcm_blockdiag_new(int n, int bn, int *tilemap);
void oefpil_tcm_blockdiag_set_tile_diag(int n, int bn, double *A, int *tilemap, int ii, const double *B);
void oefpil_tcm_blockdiag_set_tile_full(int n, int bn, double *A, int *tilemap, int ii, const double *B);

/* diags */
double* oefpil_tcm_diags_new(int n, int bn, int *tilemap);
void oefpil_tcm_diags_set_tile_diag(int n, int bn, double *A, int *tilemap, int ii, int jj, const double *B);

/* full */
double* oefpil_tcm_full_new(int n, int bn, int *tilemap);
void oefpil_tcm_full_set_tile_diag(int n, int bn, double *A, int *tilemap, int ii, int jj, const double *B);
void oefpil_tcm_full_set_tile_full(int n, int bn, double *A, int *tilemap, int ii, int jj, const double *B);



void oefpil(fitfunc_oefpil *fcn, void *data, int fun_mode,
	    int np, double *p, double *pcov,
	    int n, int nv, const double *x, const double *y,
            double *xtrue, double *ytrue,
            const double *K, int tiling, const int *tilemap,
	    int maxit, double tol,
	    int printlevel, FILE *rptfile,
            double *chi2, int conv_crit,
	    int *info, int *niter, double *sigma2, bool usesigma2, double *chi2pval);

#endif

#if 0
/* fun_mode removed - decide from presence of y */
/* yfit = NULL <=> y = NULL */
/* p menit; kdyby se nepovedlo, zkopirovat puvodni? */
/* maxit, tol <= 0 => default */
/* printlevel > 0 & rptfname = NULL => print to stdout */
/* pcov != NULL => calculate it */
/* chi2 != NULL => calculate it */
/* xfit, yfit != NULL => calculate them */

/* zavest vnitrni funkce pro jasnejsi strukturu v ramci hlavni funkce: */
/* po revizi parametru a alokaci poli zavolat oefpil_worker(), v ramci neho oefpil_iter() a oefpil_converged(), na konci pak uklidit */


void oefpil(fitfunc_oefpil *fcn, void *data,
                int np, double *p, const int *pfix, double *pcov, 
                int n, int nv, const double *x, const double *y,
                const double *K, int tiling, const int *tilemap,
                int maxit, double tol,
                int printlevel, const char *rptfname,
                double *xfit, double *yfit, double *chi2,
                int *info, int *niter);

/* nt - no tiling */
/* zavolat oefpil_tile() */

void oefpil_nt(fitfunc_oefpil *fcn, void *data,
               int np, double *p, const int *pfix, double *pcov,
               int n, int nv, const double *x, const double *y, const double *K,
               int maxit, double tol,
               int printlevel, const char *rptfname,
               double *xfit, double *yfit, double *chi2,
               int *info, int *niter);

/* static? */
void oefpil_tile(int n, int nv, const double *A,
                 double *B, int *tiling, int *tilemap);

#endif

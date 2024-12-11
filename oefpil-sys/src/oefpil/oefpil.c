#include "oefpil.h"

#include "blaslapack.h"
#include "tiled_la.h"
#include "../math77_chi2/dcdchi.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


enum matrix_print_mode {
    MATRIX_PRINT_FULL,
    MATRIX_PRINT_LOWER,
    MATRIX_PRINT_UPPER
};

enum {MATRIX_PRINT_PREC = 8};


/* podivat se na ODRPACK - PARTOL jako mozny vychozi threshold */
/* nebo jestli neco nema i minpack */
/* nejake konstanty jsou ve float.h - DBL_EPSILON, ten pak ^2/3 jako ODRPACK? */
/* pro implicitni pak asi ^1/3 */
/* pro sumu ctvercu pak ^1/2 */
/* proverit jeste: Celmins, str. 29 */
static bool convergence_check(int np, const double *p0, const double *pdiff, 
        int nx, const double *x0, const double *xdiff, 
        double chi2, double chi20, double tol, int conv_crit);
static bool conv_param_rel(int np, const double *p0, const double *pdiff, double tol);
static bool conv_param_abs(int np, const double *pdiff, double tol); /* je potreba? odrpack nejspis nema */

static void check_inf_nan(int n, const double *a, const char *where, enum oefpil_status *status);

static void matrix_debug(int iter, const char *name, int m, int n, const double *a, const char *format, int mode);
static void tiled_matrix_debug(int iter, const char *name, int bm, int bn, int n, const double *a,
                               int tiling, const char *format, int mode);


void printmatrix(FILE *rptfile, int m, int n, const double *a, const char *format, const char *label, int mode)
{
    int i, j;
    double x;

    if (m != n)
        mode = MATRIX_PRINT_FULL;

    if (label) {
        fputs(label, rptfile);
        fputs("\n", rptfile);
    }
    
    for (i = 0; i < m; i++) {
	for (j = 0; j < n; j++) {
            switch (mode) {
                
            case MATRIX_PRINT_FULL:
                x = a[j*m + i];
                break;
                
            case MATRIX_PRINT_LOWER:
                if (i >= j)
                    x = a[j*m + i];
                else
                    x = 0.0;
                break;

            case MATRIX_PRINT_UPPER:
                if (j >= i)
                    x = a[j*m + i];
                else
                    x = 0.0;
                break;

            default:
                x = NAN;
                break;
            }
                
            fprintf(rptfile, format, x);
	}
	fputs("\n", rptfile);
    }
}

static void matrix_debug(int iter, const char *name, int m, int n, const double *a, const char *format, int mode)
{
    FILE *f;
    char fname[256];

    sprintf(fname, "%d.%s.txt", iter, name);

    f = fopen(fname, "w");

    if (f) {
        printmatrix(f, m, n, a, format, NULL, mode);
        fclose(f);
    }
    /* else report an error? */
}

static void tiled_matrix_debug(int iter, const char *name, int bm, int bn, int n, const double *a,
                               int tiling, const char *format, int mode)
{
    int i, j;
    FILE *f;
    char fname[256];

    if ((bm <= 0) || (bn <= 0))
        return;

    switch (tiling) {
    case TILING_DIAG:
    case TILING_BLOCKDIAG:
        if (bm != bn)
            return;
        
        for (i = 0; i < bn; i++) {
            sprintf(fname, "%d.%s.%d.txt", iter, name, i);
            f = fopen(fname, "w");

            if (f) {
                if (tiling == TILING_DIAG)
                    printmatrix(f, n, 1, a + i*n, format, NULL, mode);
                else
                    printmatrix(f, n, n, a + i*(n*n), format, NULL, mode);

                fclose(f);
            }
        }
        break;

    case TILING_DIAGS:
    case TILING_FULL:
        for (j = 0; j < bn; j++)
            for (i = 0; i < bm; i++) {
                sprintf(fname, "%d.%s.%d.%d.txt", iter, name, i, j);
                f = fopen(fname, "w");
                
                if (f) {
                    if (tiling == TILING_DIAGS)
                        printmatrix(f, n, 1, a + (j*bm + i)*n, format, NULL, mode);
                    else
                        printmatrix(f, n, n, a + (j*bm + i)*(n*n), format, NULL, mode);
                    
                    fclose(f);
                }
            }
        break;
    } /* switch (tiling) */
}

static void check_inf_nan(int n, const double *a, const char *where, enum oefpil_status *status)
{
    int i;

    i = 0;
    
    do {
	if (isnan(a[i])) {
	    fprintf(stderr, "OEFPIL: NaN encountered in %s evaluation (%d-th value)!\n", where, i);
	    *status = OEFPIL_ERROR_NAN;
	}
	else if (isinf(a[i])) {
	    fprintf(stderr, "OEFPIL: Inf encountered in %s evaluation (%d-th value)!\n", where, i);
	    *status = OEFPIL_ERROR_INF;
	}
	
	i++;
    } while ((status == OEFPIL_STATUS_OK) && (i < n));
}


/* OEFPIL */

/* mozna by stacilo pdiff_hat a p (tj. pocitat relativni zmenu vuci p1) */
static bool conv_param_rel(int np, const double *p0, const double *pdiff, double tol)
{
    int i;
    bool r;

    r = true;

    for (i = 0; i < np; i++)
	r = r && (fabs(pdiff[i] / p0[i]) < tol);

    return r;
}

/* takze pdiff_hat */
static bool conv_param_abs(int np, const double *pdiff, double tol)
{
    int i;
    bool r;

    r = true;

    for (i = 0; i < np; i++)
        r = r && (fabs(pdiff[i]) < tol);

    return r;
}

static bool conv_chi2(double chi20, double chi2, double tol)
{
    bool r ;
    r =   (fabs(chi2-chi20) < tol);
    return r;
}

static bool convergence_check(int np, const double *p0, const double *pdiff, 
        int nx, const double *x0, const double *xdiff, 
        double chi2, double chi20, double tol, int conv_crit)
{

       /* if (conv_param_rel(np, p, w, tol) || conv_param_abs(np, w, DBL_EPSILON))
        if ((conv_param_rel(np,p,w,tol) || conv_param_abs(np, w, DBL_EPSILON)) && 
                    (conv_param_rel(nve*n, x,deltax_hat, tol) || conv_param_abs(nve*n, deltax_hat, tol)))
            */
    bool r;

    switch(conv_crit){
        case REL_PAR:
            r=conv_param_rel(np, p0, pdiff, tol) ;
            break;
        case REL_PAR_OR_ABS_PAR:
            r = (conv_param_rel(np,p0,pdiff, tol) || conv_param_abs(np, pdiff, tol));
            break;
        case REL_X:
            r=conv_param_rel(nx, x0, xdiff, tol) ;
            break;
        case REL_X_OR_ABS_X:
            r = (conv_param_rel(nx,x0,xdiff, tol) || conv_param_abs(nx, xdiff, tol));
            break;
        case REL_PAR_AND_REL_X:
            r=(conv_param_rel(np, p0, pdiff, tol) &&conv_param_rel(nx, x0, xdiff, tol)) ;
            break;
        case REL_PAR_OR_ABS_PAR_AND_REL_X_OR_ABS_X:
            r = ((conv_param_rel(np,p0,pdiff, tol) || conv_param_abs(np, pdiff, tol)) 
                    && (conv_param_rel(nx,x0,xdiff, tol) || conv_param_abs(nx, xdiff, tol)));
            break;
        case CHI2:
            r = conv_chi2(chi20, chi2, tol);
            break;
        default:
            r = true;
            printf( "*** Failure: not a valid convergence criterion %d **** \n", conv_crit);
            break;
    }
    return r;

}


/* after sufficient testing and cleanups, level 3 prints could be removed */
/* printlevel:
   0: no output
   1: print parameter values in each iteration plus stopping information
   2: print also parameter covariance matrix in each iteration
   3: comment on individual steps performed
*/

/* return value set in info:
   1: OK
   2: iteration limit reached
   3: matrix inversion encountered singular matrix
*/

void oefpil(fitfunc_oefpil *fcn, void *data, int fun_mode,
	    int np, double *p, double *pcov, 
	    int n, int nv, const double *x, const double *y,
        double *xtrue, double *ytrue,
        const double *K, int tiling, const int *tilemap,
	    int maxit, double tol,
	    int printlevel, FILE *rptfile,
        double *chi2, int conv_crit,
	    int *info, int *niter, double *sigma2, bool usesigma2,double *chi2pval)
{
    int i, j, ii, jj;
    int nve;
    int iter_num;

    int lpkinfo;
    char *lpkerrsrc = NULL;
    
    int info_; /* local info */
    enum oefpil_status status;

    char *env_printformat;
    char *dblprintformat;

    char *env_maxit;

    double *wx; /* (n*nv) work arrays of x-data */
    double *wf; /* (n) work array for function values */

    double *xdiff; /* (n*nv) */
    double *b; /* (n) */

    double *B1_diags; /* (n*nv) */
    
    const double *K_;
    double *B1_ii, *B1_jj, *xdiff_ii; /* helper sub-blocks */

    double *B2; /* (n, np) */
    double *M11 = NULL, *LM = NULL; /* (n, n) / (n) + (n) upper left block of the large matrix before inverting, its Cholesky factor */

    double *deltax_hat; /* (n*nve) */

    double topt[5];
    double qrworkopt, mqrworkopt;
    double *t, *qrwork;
    int tsize = 0, lqrwork = 0, lmqrwork1 = 0, lmqrwork2 = 0;
    double *mqrwork1 = NULL, *mqrwork2 = NULL;
    bool mqrqueried = false;
    
    double *u, *v, *w, *z;
    double *z_;
    
    double *LK = NULL;
    int *tilemapLK = NULL;
    double chi2_, chi2prev = 0;
    double pchi2, qchi2;
    long ierr;

    bool debug;

    
    /* ENVIRONMENT VARIABLES */
    
    debug = false;
    
    if (getenv("OEFPIL_DEBUG")) {
	if (printlevel > 0)
	    fprintf(rptfile, "### OEFPIL debug mode (printing intermediate matrices) ON ###" "\n\n");
	
        debug = true;
    }

    /* print format */
    dblprintformat = PRINT_FORMAT_E8;
    env_printformat = getenv("OEFPIL_PRINTFORMAT");

    if (env_printformat) {
	if (printlevel > 0)
	    printf("### Format for printing doubles: %% 10.%s ###" "\n\n", env_printformat);

#if defined(__POSIX__)
        dblprintformat = strdup("% 10.");
#elif defined(_MSC_VER)
        dblprintformat = _strdup("% 10.");
#endif
        
        strcat(dblprintformat, env_printformat);
    }

    /* maximum iterations */
    env_maxit = getenv("OEFPIL_MAXIT");

    if (env_maxit)
        maxit = atoi(env_maxit);

    
    info_ = 0;
    status = OEFPIL_STATUS_OK;

    switch (fun_mode) {
    case FUNCTION_IMPLICIT:
        nve = nv;
        break;
    case FUNCTION_EXPLICIT:
        nve = nv + 1;
        break;
    }
    
    if (!rptfile)
	printlevel = 0;
    
    if (printlevel > 0) {
	fprintf(rptfile, "### OEFPIL begin ###" "\n");
        fprintf(rptfile, "\n");
	fprintf(rptfile, "*** Input ***" "\n");
        fprintf(rptfile, "\n");
        fprintf(rptfile, "* function mode: ");

        if (fun_mode == FUNCTION_EXPLICIT)
            fprintf(rptfile, "explicit" "\n");
        else
            fprintf(rptfile, "implicit" "\n");

	fprintf(rptfile, "* covariance matrix tiling mode: %d" "\n", tiling);
	fprintf(rptfile, "* number of data points: %d" "\n", n);
	fprintf(rptfile, "* number of variables: %d" "\n", nv);
	fprintf(rptfile, "* number of parameters: %d" "\n", np);
	fprintf(rptfile, "* max number of iterations: %d" "\n", maxit);
	fprintf(rptfile, "* relative tolerance for convergence: %e" "\n", tol);
	fprintf(rptfile, "* number of decimal digits in printed output: %d" "\n", MATRIX_PRINT_PREC);
        fprintf(rptfile, "\n");
    }

    if (printlevel > 2) {
	printmatrix(rptfile, n, nv, x, dblprintformat, "* x values *", MATRIX_PRINT_FULL);

        if (fun_mode == FUNCTION_EXPLICIT)
            printmatrix(rptfile, n, 1, y, dblprintformat, "* y values *", MATRIX_PRINT_FULL);

	printmatrix(rptfile, n, nv, xtrue, dblprintformat, "* xtrue values *", MATRIX_PRINT_FULL);

        if (fun_mode == FUNCTION_EXPLICIT)
            printmatrix(rptfile, n, 1, ytrue, dblprintformat, "* ytrue values *", MATRIX_PRINT_FULL);
    }

    
    iter_num = 0;

    if (debug) {
        tiled_matrix_debug(iter_num, "K", nve, nve, n, K, tiling, dblprintformat, MATRIX_PRINT_FULL);
        matrix_debug(iter_num, "x", n, nv, x, dblprintformat, MATRIX_PRINT_FULL);

        if (fun_mode == FUNCTION_EXPLICIT)
            matrix_debug(iter_num, "y", n, 1, y, dblprintformat, MATRIX_PRINT_FULL);
    }
    
    if (printlevel > 2)
	fprintf(rptfile, "*** Initialization ***\n\n");

    wf = (double*)malloc(n * sizeof(double));
    b = (double*)malloc(n * sizeof(double));

    B1_diags = (double*)malloc(n*nve * sizeof(double));
    wx = (double*)malloc(n*nve * sizeof(double));
    xdiff = (double*)malloc(n*nve * sizeof(double));
    deltax_hat = (double*)malloc(n*nve * sizeof(double));

    B2 = (double*)malloc(n*np * sizeof(double));

    u = (double*)malloc(n * sizeof(double));
    v = (double*)malloc(n * sizeof(double));
    w = (double*)malloc(np * sizeof(double)); /* deltap_hat = w */
    z = (double*)malloc(nve*n * sizeof(double));
    
    switch (tiling) {
    case TILING_DIAG:
    case TILING_DIAGS:
        /* M11 diagonal: use separate arrays */
        M11 = (double*)malloc(n * sizeof(double));
        LM = (double*)malloc(n * sizeof(double));
        break;

    case TILING_NOTILES:
    case TILING_BLOCKDIAG:
    case TILING_FULL:
        /* M11 nondiagonal: the matrices can replace each other */
        M11 = (double*)malloc(n * n * sizeof(double));
        LM = M11;
        break;
    }

    dgeqr(n, np, B2, n, topt, -1, &qrworkopt, -1, &lpkinfo); /* E = B2 */

    if (lpkinfo != 0) {
        status = OEFPIL_ERROR_LAPACK;
        lpkerrsrc = "DGEQR query";
    }

    tsize = topt[0];
    lqrwork = qrworkopt;

    t = (double*)malloc(tsize * sizeof(double));
    qrwork = (double*)malloc(lqrwork * sizeof(double));

    if (printlevel > 0)
        fprintf(rptfile, "* workspace needed for QR: %zu + %zu bytes *" "\n",
                tsize * sizeof(double), lqrwork * sizeof(double));


    /* for chi2 */
    switch (tiling) {
    case TILING_NOTILES:
        /* TODO - check */
        LK = (double*)malloc((nve*n)*(nve*n) * sizeof(double));

        for (i = 0; i < (nve*n)*(nve*n); i++)
            LK[i] = K[i];

        dpotrf("L", nve*n, LK, nve*n, &lpkinfo);
        break;

    case TILING_DIAG:
        tilemapLK = oefpil_tilemap_diagtiles_new(nve);
        LK = oefpil_tcm_diag_new(n, nve, tilemapLK);

        for (i = 0; i < nve*n; i++)
            LK[i] = K[i];
        
        dpottrf_diag(n, nve, LK, tilemapLK, &lpkinfo);
        break;

    case TILING_BLOCKDIAG:
        tilemapLK = oefpil_tilemap_diagtiles_new(nve);
        LK = oefpil_tcm_blockdiag_new(n, nve, tilemapLK);

        for (i = 0; i < nve*(n*n); i++)
            LK[i] = K[i];
        
        dpottrf_blockdiag(n, nve, LK, tilemapLK, &lpkinfo);
        break;

    case TILING_DIAGS:
        tilemapLK = oefpil_tilemap_alltiles_new(nve);
        LK = oefpil_tcm_diags_new(n, nve, tilemapLK);
        
        for (i = 0; i < (nve*nve)*n; i++)
            LK[i] = K[i];

        dpottrf_diags(n, nve, LK, tilemapLK, &lpkinfo);
        break;

    case TILING_FULL:
        tilemapLK = oefpil_tilemap_alltiles_new(nve);
        LK = oefpil_tcm_full_new(n, nve, tilemapLK);

        for (i = 0; i < (nve*nve)*(n*n); i++)
            LK[i] = K[i];
        
        dpottrf_full(n, nve, LK, tilemapLK, &lpkinfo);
        break;
    }
    

    if (lpkinfo != 0) {
        status = OEFPIL_ERROR_LAPACK;
        lpkerrsrc = "DPOTTRF_xxx";
    }
    
    for (i = 0; i < n*nv; i++)
        wx[i] = xtrue[i];

    if (fun_mode == FUNCTION_EXPLICIT){
        for (i = 0; i < n; i++)
            wx[n*nv + i] = ytrue[i];
    }
    
    if (printlevel > 0) {
        printmatrix(rptfile, 1, np, p, dblprintformat, "* initial parameter estimates *" "\n", MATRIX_PRINT_FULL);
        fprintf(rptfile, "\n");
    }

    if (printlevel > 0)
        fprintf(rptfile, "*** Iterations ***" "\n\n");


    /* ITERATIONS */
    while (status == OEFPIL_STATUS_OK) {
	
        if (printlevel > 2)
            fprintf(rptfile, "** iteration %d begin **" "\n\n", iter_num);

        if (debug)
            matrix_debug(iter_num, "wx", n, nve, wx, dblprintformat, MATRIX_PRINT_FULL);

        fcn(data, n, wx, p, wf, B1_diags, B2);
        check_inf_nan(n, wf, "function values", &status);
        check_inf_nan(n*nv, B1_diags, "function derivatives", &status);
        check_inf_nan(n*np, B2, "function parameter derivatives", &status);
        
        if (fun_mode == FUNCTION_EXPLICIT) {
            /* b(n) = (f(x) - y) */
            /* B1(n,(nv+1)*n) = ( B11 | ... | B1nv | -I_n ) */
            for (i = 0; i < n; i++) {
                b[i] = wf[i] - wx[n*nv + i];
                B1_diags[n*nv + i] = -1.0;
            }
        }
        else {
            /* b(n) = (f(x_i,y_i)) */
            for (i = 0; i < n; i++)
                b[i] = wf[i];
        }

        if (debug) {
            matrix_debug(iter_num, "b", n, 1, b, dblprintformat, MATRIX_PRINT_FULL);
            matrix_debug(iter_num, "B1_diags", n, nve, B1_diags, dblprintformat, MATRIX_PRINT_FULL);
            matrix_debug(iter_num, "B2", n, np, B2, dblprintformat, MATRIX_PRINT_FULL);
        }

        /* M11(n, n) = B1(n, nve*n) * U(nve*n, nve*n) * B1^T(nve*n, n) */
        /* M11 = B11 * Kxx * B11 + B11 * Kxy * B12 + B12 * Kyx * B11 + B12 * Kyy * B12 */
        /* Cholesky M11 = LM LM^T */
        /* E = LM^(-1) B2 (solve LM E = B2), in place */

        switch (tiling) {

	case TILING_NOTILES:
	    /* TODO */
	    break;

        case TILING_DIAG:
        case TILING_DIAGS:
            for (i = 0; i < n; i++)
                M11[i] = 0.0;

            if (tiling == TILING_DIAG) {
                for (ii = 0; ii < nve; ii++)
                    if (tilemap[ii] == TILE_DIAG) {
                        B1_ii = B1_diags + ii*n;
                        K_ = K + ii*n;
                    
                        for (i = 0; i < n; i++)
                            M11[i] += B1_ii[i] * K_[i] * B1_ii[i];
                    }
            }
            else {  /* tiling == TILING_DIAGS */
                for (jj = 0; jj < nve; jj++)
                    for (ii = 0; ii < nve; ii++)
                        if (tilemap[jj*nve + ii] == TILE_DIAG) {
                            B1_ii = B1_diags + ii*n;
                            B1_jj = B1_diags + jj*n;
                            K_ = K + (jj*nve + ii) * n;
                            
                            for (i = 0; i < n; i++)
                                M11[i] += B1_ii[i] * K_[i] * B1_jj[i];
                        }
            }

            /* "diagonal Cholesky" */
            for (i = 0; i < n; i++)
                LM[i] = sqrt(M11[i]);

            /* E = B2 */
            for (i = 0; i < n; i++)
                for (j = 0; j < np; j++)
                    B2[j*n + i] /= LM[i];

            if (debug) {
                matrix_debug(iter_num, "M11_diag", n, 1, M11, dblprintformat, MATRIX_PRINT_FULL);
                matrix_debug(iter_num, "LM_diag", n, 1, LM, dblprintformat, MATRIX_PRINT_FULL);
            }
            break; /* end case TILING_DIAG, TILING_DIAGS */

        case TILING_BLOCKDIAG:
        case TILING_FULL:
            for (i = 0; i < n*n; i++)
                M11[i] = 0.0;

            if (tiling == TILING_BLOCKDIAG)
                for (ii = 0; ii < nve; ii++) {
                    B1_ii = B1_diags + ii*n;
                    K_ = K + ii*(n*n);

                    for (j = 0; j < n; j++)
                        for (i = 0; i < n; i++)
                            M11[j*n + i] += B1_ii[i] * K_[j*n + i] * B1_ii[j];
                }
            else /* tiling == TILING_FULL */
                for (jj = 0; jj < nve; jj++)
                    for (ii = 0; ii < nve; ii++){
                        if (tilemap[jj*nve + ii] != TILE_NULL) {
                            B1_ii = B1_diags + ii*n;
                            B1_jj = B1_diags + jj*n;
                            K_ = K + (jj*nve*n*n + ii*n) ;
                        
                            for (j = 0; j < n; j++)
                                for (i = 0; i < n; i++)
                                    M11[j*n + i] += B1_ii[i] * K_[j*n*nve + i] * B1_jj[j];
                        }
                    }

            if (debug)
                matrix_debug(iter_num, "M11", n, n, M11, dblprintformat, MATRIX_PRINT_LOWER);

            /* Cholesky in-place */
            dpotrf("L", n, M11, n, &lpkinfo);
            /* M11 now contains just the factor LM */
            
            if (lpkinfo != 0) {
                status = OEFPIL_ERROR_LAPACK;
                lpkerrsrc = "DPOTRF(M11)";
                break;
            }

            if (debug)
                matrix_debug(iter_num, "LM_triang", n, n, LM, dblprintformat, MATRIX_PRINT_LOWER);

            /* E */
            dtrsm("L", "L", "N", "N", n, np, 1.0, LM, n, B2, n);
            /* B2 now contains the solution, E */

            if (debug)
                matrix_debug(iter_num, "E", n, np, B2, dblprintformat, MATRIX_PRINT_FULL);
            break; /* end case TILING_BLOCKDIAG, TILING_FULL */
        } /* end switch (tiling) */

        
        /* break the loop when LAPACK error encountered */
        if (status == OEFPIL_ERROR_LAPACK) {
            break;
        }

        /* common part */
        
        /* ### (A1) xdiff ### */
        /* xdiff(nve*n) */
        for (i = 0; i < n*nv; i++)
            xdiff[i] = x[i] - wx[i];

	if (fun_mode == FUNCTION_EXPLICIT)
	    for (i = 0; i < n; i++)
		xdiff[n*nv + i] = y[i] - wx[n*nv + i];


        /* ### (A2) B1 * xdiff + b ### */
        /* u(n) = B1(n, nve*n) * xdiff(nve*n) + b(n) */
        for (i = 0; i < n; i++)
            u[i] = b[i]; // u lze eliminovat, muzeme pouzit b?
        
        for (ii = 0; ii < nve; ii++) {
            B1_ii = B1_diags + ii*n;
            xdiff_ii = xdiff + ii*n;
            
            for (i = 0; i < n; i++)
                u[i] += B1_ii[i] * xdiff_ii[i];
        }


        /* E(n, np) -> QE(n, n) * RE(n, np) */
        dgeqr(n, np, B2, n, t, tsize, qrwork, lqrwork, &lpkinfo);

        if (lpkinfo != 0) {
            status = OEFPIL_ERROR_LAPACK;
            lpkerrsrc = "DGEQR (E)";
            break;
        }

        if (!mqrqueried) {
            dgemqr("L", "T", n, 1, np, B2, n, t, tsize, u, n, &mqrworkopt, -1, &lpkinfo); /* E = B2 */

            if (lpkinfo != 0) {
                status = OEFPIL_ERROR_LAPACK;
                lpkerrsrc = "DGEMQR 1 query";
                break;
            }

            lmqrwork1 = mqrworkopt;
            mqrwork1 = (double*)malloc(lmqrwork1 * sizeof(double));

            if (printlevel > 0)
                fprintf(rptfile, "* workspace needed for MQR1: %zu bytes *" "\n",
                        lmqrwork1 * sizeof(double));

            dgemqr("L", "N", n, 1, np, B2, n, t, tsize, u, n, &mqrworkopt, -1, &lpkinfo); /* E = B2 */

            if (lpkinfo != 0) {
                status = OEFPIL_ERROR_LAPACK;
                lpkerrsrc = "DGEMQR 2 query";
                break;
            }

            lmqrwork2 = mqrworkopt;
            mqrwork2 = (double*)malloc(lmqrwork2 * sizeof(double));

            if (printlevel > 0)
                fprintf(rptfile, "* workspace needed for MQR2: %zu bytes *" "\n\n",
                        lmqrwork2 * sizeof(double));

            mqrqueried = true;
        }

        /* u -> solve(LM, u) */
        switch (tiling) {
        case TILING_DIAG:
        case TILING_DIAGS:
            for (i = 0; i < n; i++)
                u[i] = u[i] / LM[i];
		
            break;

	case TILING_NOTILES:
        case TILING_BLOCKDIAG:
        case TILING_FULL:
            dtrsv("L", "N", "N", n, LM, n, u, 1);
            break;
        }  /* end switch (tiling) */

        /* v(n) = u(n) */
        for (i = 0; i < n; i++)
            v[i] = u[i];
            
        /* u(n) -> QE^T(n, n) * u(n) */
        dgemqr("L", "T", n, 1, np, B2, n, t, tsize, u, n, mqrwork1, lmqrwork1, &lpkinfo); /* E = B2 */

        if (lpkinfo != 0) {
            status = OEFPIL_ERROR_LAPACK;
            lpkerrsrc = "DGEMQR 1";
            break;
        }

        /* QE is square = (Q_R Q_0) */
        /* (Q_R Q_0) (Q_R^T; Q_0^T) u = (Q_R Q_0) (Q_R^T u; Q_0^T u) = Q_R Q_R^T u + Q_0 Q_0^T u */
        /* zero-out Q_0^T u, so that only Q_R Q_R^T u is calculated */
        for (i = np; i < n; i++)
            u[i] = 0.0;

        /* w(np) = -u(1:np); for deltap_hat */
        for (i = 0; i < np; i++)
            w[i] = -u[i];

        /* u(n) -> QE{(n, n) * u(n) */
        dgemqr("L", "N", n, 1, np, B2, n, t, tsize, u, n, mqrwork2, lmqrwork2, &lpkinfo); /* E = B2 */

        if (lpkinfo != 0) {
            status = OEFPIL_ERROR_LAPACK;
            lpkerrsrc = "DGEMQR 2";
            break;
        }

        /* v(n) -> v(n) - u(n) */
        for (i = 0; i < n; i++)
            v[i] = v[i] - u[i];
            
        /* v(n) -> solve(LM^T(n, n), v(n)) */
        switch (tiling) {
        case TILING_DIAG:
        case TILING_DIAGS:
            for (i = 0; i < n; i++)
                v[i] = v[i] / LM[i];
		
            break;

	case TILING_NOTILES:
        case TILING_BLOCKDIAG:
        case TILING_FULL:
            dtrsv("L", "T", "N", n, LM, n, v, 1);
            break;
        } /* end switch (tiling) */
            
        /* deltap_hat(np) = w(np) */
        /* w(np) -> solve(RE(np, np), w(np) */
        dtrsv("U", "N", "N", np, B2, n, w, 1); /* E = B2, deltap_hat = w */

        /* pcov */
        dpotri("U", np, B2, n, &lpkinfo); /* E = B2 */

        for (j = 0; j < np; j++)
            for (i = 0; i <= j; i++) {
                pcov[j*np + i] = B2[j*n + i];
                pcov[i*np + j] = B2[j*n + i];
            }

        /* ### finish ### */

        /* rotate and update p */
        for (i = 0; i < np; i++)
            p[i] = p[i] + w[i];
            
        /* calculate deltax_hat, update wx */

        /* z(nve*n) = B1^T(nve*n, n) * v(n) */
        for (ii = 0; ii < nve; ii++) {
            B1_ii = B1_diags + ii*n;
            z_ = z + ii*n;
            
            for (i = 0; i < n; i++)
                z_[i] = B1_ii[i] * v[i];
        }

            
        /* deltax_hat(nve*n) = xdiff(nve*n) - K(nve*n, nve*n) * z(nve*n) */
        for (i = 0; i < nve*n; i++)
            deltax_hat[i] = xdiff[i];
            
        switch (tiling) {
	case TILING_NOTILES:
	    /* TODO */
	    /* dgemv / dsymv */
	    break;
	    
        case TILING_DIAG:
            dgetmv_diag(n, nve, -1.0, K, tilemap, z, 1.0, deltax_hat);
            break;
                
        case TILING_BLOCKDIAG:
            dgetmv_blockdiag(n, nve, -1.0, K, tilemap, z, 1.0, deltax_hat);
            break;
                
        case TILING_DIAGS:
            dgetmv_diags(n, nve, nve, -1.0, K, tilemap, z, 1.0, deltax_hat);
            break;
                
        case TILING_FULL:
            /*         printf("FULL K \n");
                       for (i=0; i< nve*n; i++){
                       for (j=0; j< nve*n; j++)
                       printf("%g ", K[i*nve*n+j]);
                       printf("\n");
                       }
            */

            dgetmv_full(n, nve, nve, -1.0, K, tilemap, z, 1.0, deltax_hat);
//              dgemv("N", nve, nve, -1.0, K, n, z, 1, 1.0, deltax_hat, 1);            
            break;
        } /* end switch (tiling) */
        
        for (i = 0; i < n*nve; i++)
            wx[i] = wx[i] + deltax_hat[i];

       
        /* chi2 */
        /* update xdiff with new wx */
        for (i = 0; i < n*nv; i++)
            xdiff[i] = x[i] - wx[i];

	if (fun_mode == FUNCTION_EXPLICIT)
	    for (i = 0; i < n; i++)
		xdiff[n*nv + i] = y[i] - wx[n*nv + i];

        /* solve (LK, xdiff) */
        switch (tiling) {
	case TILING_NOTILES:
	    /* TODO */
	    /* dtrsv */
	    break;
	    
        case TILING_DIAG:
            dtrtsv_diag(n, nve, LK, tilemapLK, xdiff, &lpkinfo);
            break;
        case TILING_BLOCKDIAG:
            dtrtsv_blockdiag(n, nve, LK, tilemapLK, xdiff, &lpkinfo);
            break;
        case TILING_DIAGS:
            dtrtsv_diags(n, nve, LK, tilemapLK, xdiff, &lpkinfo);
            break;
        case TILING_FULL:
            dtrtsv_full(n, nve, LK, tilemapLK, xdiff, &lpkinfo);
            break;
        } /* end switch (tiling) */

        if (lpkinfo != 0) {
            status = OEFPIL_ERROR_LAPACK;
            lpkerrsrc = "DTRTSV_xxx query";
            break;
        }


        chi2_ = ddot(nve*n, xdiff, 1, xdiff, 1);

        /* deltap_hat = w */
        if (convergence_check(np, p, w, nve*n, x, deltax_hat, chi2_, chi2prev, tol, conv_crit))
            status = OEFPIL_CONVERGED;
        
 	
        if (printlevel > 0) {
            fprintf(rptfile, "%-4d\t", iter_num);
            printmatrix(rptfile, 1, np, p, dblprintformat, NULL, MATRIX_PRINT_FULL);
        }

        if (printlevel > 1) {
            printmatrix(rptfile, np, np, pcov, dblprintformat, "* pcov *", MATRIX_PRINT_FULL);
            printmatrix(rptfile, 1, 1, &chi2_, PRINT_FORMAT_E15, "* chi2 *", MATRIX_PRINT_FULL);
            fprintf(rptfile, "\n");
        }
         
        iter_num++;
        chi2prev = chi2_;

        if (iter_num >= maxit)
            status = OEFPIL_MAX_ITER_REACHED;

    } /* end of main loop */

    switch (status) {
    case OEFPIL_CONVERGED:
        if (printlevel > 0)
            fprintf(rptfile, "*** Success: parameters converged in %d iterations ***\n", iter_num);

        info_ = 1;
        break;

    case OEFPIL_MAX_ITER_REACHED:
        if (printlevel > 0)
            fprintf(rptfile, "*** Failure: maximum number of iterations (%d) reached ***\n", maxit);

        info_ = 2;
        break;

    case OEFPIL_ERROR_NAN:
        if (printlevel > 0)
            fprintf(rptfile, "*** Failure: NaN encountered in function/derivatives values ***\n");

        info_ = 3;
        break;

    case OEFPIL_ERROR_INF:
        if (printlevel > 0)
            fprintf(rptfile, "*** Failure: Inf encountered in function/derivatives values ***\n");

        info_ = 3;
        break;

    case OEFPIL_ERROR_LAPACK:
        if (printlevel > 0)
            fprintf(rptfile, "*** Failure: error encountered in LAPACK procedure %s, "
                    "info code %d (see LAPACK documentation for details) ***\n", lpkerrsrc, lpkinfo);
        
        info_ = 3;
        break;

    default:
        fprintf(rptfile, "*** Unknown termination status ***\n");
        break;
    }

    if (printlevel > 0)
        fprintf(rptfile, "\n### OEFPIL end ###\n");

    if (info)
        *info = info_;

    if (niter)
        *niter = iter_num;

    if (chi2)
        *chi2 = chi2_;

    if (chi2pval)
    {
	    dcdchi(chi2_, (n-np), &pchi2, &qchi2, &ierr);
	    //if (ierr > 0)
	//	    chi2pval = NULL;
	    *chi2pval = qchi2;
    }

    for (i = 0; i < n*nv; i++)
        xtrue[i] = wx[i];
    
    if (fun_mode == FUNCTION_EXPLICIT)
        for (i = 0; i < n; i++)
            ytrue[i] = wx[n*nv + i];

    if (sigma2){
	    *sigma2 = chi2_/(n-np);
	    if (usesigma2){
		    for (j = 0; j < np; j++)
			    for (i = 0; i < np; i++) {
				    pcov[i*np + j] *= *sigma2;
			    }
	    }
    }


    free(wf);
    free(b);
    free(t);
    free(qrwork);

    if (mqrqueried) {
        free(mqrwork1);
        free(mqrwork2);
    }

    free(B1_diags);
    free(wx);
    free(xdiff);
    free(deltax_hat);

    free(B2);

    switch (tiling) {
    case TILING_DIAG:
    case TILING_DIAGS:
        free(M11);
        free(LM);
        break;

    case TILING_NOTILES:
    case TILING_BLOCKDIAG:
    case TILING_FULL:
        free(M11);
        break;
    }

    free(u);
    free(v);
    free(w); /* deltap_hat */
    free(z);

    free(LK);
    free(tilemapLK);

    return;
}

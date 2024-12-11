#ifndef TILED_LA_H
#define TILED_LA_H

enum tile_type {
    TILE_TRANS_FULL = -2,
    TILE_TRANS_DIAG = -1,
    TILE_NULL = 0,
    TILE_DIAG = 1,
    TILE_FULL = 2
};

enum tiling_mode {
    TILING_NOTILES,
    TILING_DIAG,
    TILING_BLOCKDIAG,
    TILING_DIAGS,
    TILING_FULL
};

int *tilemap_n_new(int bn);
int *tilemap_mn_new(int bm, int bn);

double* tiled_matrix_diag_new(int n, int bn, int *tilemap);
double* tiled_matrix_blockdiag_new(int n, int bn, int *tilemap);
double* tiled_matrix_diags_new(int n, int bm, int bn, int *tilemap);
double* tiled_matrix_full_new(int n, int bm, int bn, int *tilemap);

void tiled_matrix_diag_set_tile_diag(int n, int bn, double *A, int *tilemap, int ii, const double *B);

void tiled_matrix_blockdiag_set_tile_diag(int n, int bn, double *A, int *tilemap, int ii, const double *B);
void tiled_matrix_blockdiag_set_tile_full(int n, int bn, double *A, int *tilemap, int ii, const double *B, const char *trans);

void tiled_matrix_diags_set_tile_diag(int n, int bm, int bn, double *A, int *tilemap, int ii, int jj, const double *B);

void tiled_matrix_full_set_tile_diag(int n, int bm, int bn, double *A, int *tilemap, int ii, int jj, const double *B);
void tiled_matrix_full_set_tile_full(int n, int bm, int bn, double *A, int *tilemap, int ii, int jj, const double *B, const char *trans);

/* no trans at the moment */
void dgetmv_diag(int n, int bn,
                 double alpha, const double *a, const int *tma,
                 const double *x,
                 double beta, double *y);
void dgetmv_blockdiag(int n, int bn,
                      double alpha, const double *a, const int *tma,
                      const double *x,
                      double beta, double *y);
void dgetmv_diags(int n, int bm, int bn,
                  double alpha, const double *a, const int *tma,
                  const double *x,
                  double beta, double *y);
void dgetmv_full(int n, int bm, int bn,
                 double alpha, const double *a, const int *tma,
                 const double *x,
                 double beta, double *y);

void dsytmtv(const char *uplo, int n, int bn,
             double alpha, const double *a, const int *tma,
             const double *x,
             double beta, double *y);

void dpottrf_diag(int n, int bn, double *a, const int *tma, int *info);
void dpottrf_blockdiag(int n, int bn, double *a, const int *tma, int *info);
void dpottrf_diags(int n, int bn, double *a, const int *tma, int *info);
void dpottrf_full(int n, int bn, double *a, const int *tma, int *info);

void dtrtsv_diag(int n, int bn, const double *a, const int *tma, double *b, int *info);
void dtrtsv_blockdiag(int n, int bn, const double *a, const int *tma, double *b, int *info);
void dtrtsv_diags(int n, int bn, const double *a, const int *tma, double *b, int *info);
void dtrtsv_full(int n, int bn, const double *a, const int *tma, double *b, int *info);

#endif

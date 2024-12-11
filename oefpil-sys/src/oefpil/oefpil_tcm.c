#include "oefpil.h"

#include "tiled_la.h"

int* oefpil_tilemap_diagtiles_new(int bn)
{
    return tilemap_n_new(bn);
}

int* oefpil_tilemap_alltiles_new(int bn)
{
    return tilemap_mn_new(bn, bn);
}

double* oefpil_tcm_diag_new(int n, int bn, int *tilemap)
{
    return tiled_matrix_diag_new(n, bn, tilemap);
}

double* oefpil_tcm_blockdiag_new(int n, int bn, int *tilemap)
{
    return tiled_matrix_blockdiag_new(n, bn, tilemap);
}

double* oefpil_tcm_diags_new(int n, int bn, int *tilemap)
{
    return tiled_matrix_diags_new(n, bn, bn, tilemap);
}

double* oefpil_tcm_full_new(int n, int bn, int *tilemap)
{
    return tiled_matrix_full_new(n, bn, bn, tilemap);
}

void oefpil_tcm_diag_set_tile_diag(int n, int bn, double *A, int *tilemap, int ii, const double *B)
{
    tiled_matrix_diag_set_tile_diag(n, bn, A, tilemap, ii, B);
}

void oefpil_tcm_blockdiag_set_tile_diag(int n, int bn, double *A, int *tilemap, int ii, const double *B)
{
    tiled_matrix_blockdiag_set_tile_diag(n, bn, A, tilemap, ii, B);
}

void oefpil_tcm_blockdiag_set_tile_full(int n, int bn, double *A, int *tilemap, int ii, const double *B)
{
    tiled_matrix_blockdiag_set_tile_full(n, bn, A, tilemap, ii, B, "N");
}

void oefpil_tcm_diags_set_tile_diag(int n, int bn, double *A, int *tilemap, int ii, int jj, const double *B)
{
    tiled_matrix_diags_set_tile_diag(n, bn, bn, A, tilemap, ii, jj, B);

    if (ii != jj)
        tiled_matrix_diags_set_tile_diag(n, bn, bn, A, tilemap, jj, ii, B);

}

void oefpil_tcm_full_set_tile_diag(int n, int bn, double *A, int *tilemap, int ii, int jj, const double *B)
{
    tiled_matrix_full_set_tile_diag(n, bn, bn, A, tilemap, ii, jj, B);

    if (ii != jj)
        tiled_matrix_full_set_tile_diag(n, bn, bn, A, tilemap, jj, ii, B);

}

void oefpil_tcm_full_set_tile_full(int n, int bn, double *A, int *tilemap, int ii, int jj, const double *B)
{
    tiled_matrix_full_set_tile_full(n, bn, bn, A, tilemap, ii, jj, B, "N");

    if (ii != jj)
        tiled_matrix_full_set_tile_full(n, bn, bn, A, tilemap, jj, ii, B, "T");
}

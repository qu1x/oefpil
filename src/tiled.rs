//! Tiled linear algebra of the covariance matrix with functions named as in the reference
//! implementation, where `n` is the number of samples per tile and `bn` the number of tiles per row
//! and per column.
//!
//! A diagonal tile stores `n` fields whereas a block tile stores `n * n` fields. Modes
//! [`Mode::BlockDiagonal`] and [`Mode::Full`] allocate block tiles even for diagonal tiles and
//! [`Mode::Full`] is stored as dense column-major matrix of order `bn * n`.

#![allow(
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    clippy::neg_cmp_op_on_partial_ord,
    clippy::similar_names,
    clippy::suboptimal_flops,
    clippy::too_many_arguments
)]

use crate::la::{Diag, Trans, Uplo, dgemv, dpotrf, dtrsv};

/// Tiling mode of covariance matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    /// Diagonal tiles on its diagonal.
    Diagonal = 1,
    /// Diagonal or block tiles on its diagonal.
    BlockDiagonal = 2,
    /// Diagonal tiles all over.
    Diagonals = 3,
    /// Diagonal or block tiles all over.
    Full = 4,
}

/// Mode of tile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tile {
    /// Unset.
    None,
    /// Diagonal.
    Diagonal,
    /// Block.
    Block,
}

pub fn tcm_diag_set_tile_diag(n: usize, a: &mut [f64], tilemap: &mut [Tile], ii: usize, b: &[f64]) {
    a[ii * n..][..n].copy_from_slice(b);
    tilemap[ii] = Tile::Diagonal;
}

pub fn tcm_blockdiag_set_tile_diag(
    n: usize,
    a: &mut [f64],
    tilemap: &mut [Tile],
    ii: usize,
    b: &[f64],
) {
    let a = &mut a[ii * n * n..][..n * n];
    // Clear what a previous block tile left behind.
    a.fill(0.0);
    for i in 0..n {
        a[i * n + i] = b[i];
    }
    tilemap[ii] = Tile::Diagonal;
}

pub fn tcm_blockdiag_set_tile_full(
    n: usize,
    a: &mut [f64],
    tilemap: &mut [Tile],
    ii: usize,
    b: &[f64],
) {
    a[ii * n * n..][..n * n].copy_from_slice(b);
    tilemap[ii] = Tile::Block;
}

pub fn tcm_blockdiag_set_tile_half(
    n: usize,
    a: &mut [f64],
    tilemap: &mut [Tile],
    ii: usize,
    b: &[f64],
) {
    let a = &mut a[ii * n * n..][..n * n];
    let mut b = b.iter();
    for i in 0..n {
        for j in 0..=i {
            let field = *b.next().unwrap();
            a[j * n + i] = field;
            a[i * n + j] = field;
        }
    }
    tilemap[ii] = Tile::Block;
}

pub fn tcm_diags_set_tile_diag(
    n: usize,
    bn: usize,
    a: &mut [f64],
    tilemap: &mut [Tile],
    ii: usize,
    jj: usize,
    b: &[f64],
) {
    for (ii, jj) in [(ii, jj), (jj, ii)] {
        a[(jj * bn + ii) * n..][..n].copy_from_slice(b);
        tilemap[jj * bn + ii] = Tile::Diagonal;
    }
}

pub fn tcm_full_set_tile_diag(
    n: usize,
    bn: usize,
    a: &mut [f64],
    tilemap: &mut [Tile],
    ii: usize,
    jj: usize,
    b: &[f64],
) {
    for (ii, jj) in [(ii, jj), (jj, ii)] {
        let a = &mut a[jj * bn * n * n + ii * n..];
        // Clear what a previous block tile left behind.
        for j in 0..n {
            for i in 0..n {
                a[j * n * bn + i] = 0.0;
            }
        }
        for i in 0..n {
            a[i * n * bn + i] = b[i];
        }
        tilemap[jj * bn + ii] = Tile::Diagonal;
    }
}

pub fn tcm_full_set_tile_full(
    n: usize,
    bn: usize,
    a: &mut [f64],
    tilemap: &mut [Tile],
    ii: usize,
    jj: usize,
    b: &[f64],
) {
    let o = bn * n;
    for i in 0..n {
        for j in 0..n {
            let field = b[i * n + j];
            a[(jj * n + i) * o + ii * n + j] = field;
            a[(ii * n + j) * o + jj * n + i] = field;
        }
    }
    tilemap[jj * bn + ii] = Tile::Block;
    tilemap[ii * bn + jj] = Tile::Block;
}

pub fn tcm_full_set_tile_half(
    n: usize,
    bn: usize,
    a: &mut [f64],
    tilemap: &mut [Tile],
    ii: usize,
    jj: usize,
    b: &[f64],
) {
    let o = bn * n;
    let an = (ii * o + jj) * n;
    let at = (jj * o + ii) * n;
    let mut b = b.iter();
    for i in 0..n {
        for j in 0..=i {
            let field = *b.next().unwrap();
            a[at + j * o + i] = field;
            a[an + i * o + j] = field;
            a[at + i * o + j] = field;
            a[an + j * o + i] = field;
        }
    }
    tilemap[jj * bn + ii] = Tile::Block;
    tilemap[ii * bn + jj] = Tile::Block;
}

pub fn dgetmv_diag(
    n: usize,
    bn: usize,
    alpha: f64,
    a: &[f64],
    tma: &[Tile],
    x: &[f64],
    beta: f64,
    y: &mut [f64],
) {
    for ii in 0..bn {
        if tma[ii] == Tile::Diagonal {
            let a = &a[ii * n..];
            let x = &x[ii * n..];
            let y = &mut y[ii * n..];
            for i in 0..n {
                y[i] = beta * y[i] + alpha * a[i] * x[i];
            }
        }
    }
}

pub fn dgetmv_blockdiag(
    n: usize,
    bn: usize,
    alpha: f64,
    a: &[f64],
    tma: &[Tile],
    x: &[f64],
    beta: f64,
    y: &mut [f64],
) {
    for ii in 0..bn {
        let a = &a[ii * (n * n)..];
        let x = &x[ii * n..];
        let y = &mut y[ii * n..];
        match tma[ii] {
            Tile::Diagonal => {
                for i in 0..n {
                    y[i] = beta * y[i] + alpha * a[i * n + i] * x[i];
                }
            }
            Tile::Block => dgemv(Trans::No, n, n, alpha, a, n, x, 1, beta, y, 1),
            Tile::None => {}
        }
    }
}

pub fn dgetmv_diags(
    n: usize,
    bn: usize,
    alpha: f64,
    a: &[f64],
    tma: &[Tile],
    x: &[f64],
    beta: f64,
    y: &mut [f64],
) {
    for ii in 0..bn {
        let y = &mut y[ii * n..][..n];
        for y in &mut *y {
            *y *= beta;
        }
        for jj in 0..bn {
            if tma[jj * bn + ii] == Tile::Diagonal {
                let a = &a[(jj * bn + ii) * n..];
                let x = &x[jj * n..];
                for i in 0..n {
                    y[i] += alpha * a[i] * x[i];
                }
            }
        }
    }
}

pub fn dgetmv_full(
    n: usize,
    bn: usize,
    alpha: f64,
    a: &[f64],
    tma: &[Tile],
    x: &[f64],
    beta: f64,
    y: &mut [f64],
) {
    for ii in 0..bn {
        let y = &mut y[ii * n..][..n];
        for y in &mut *y {
            *y *= beta;
        }
        for jj in 0..bn {
            let a = &a[jj * bn * n * n + ii * n..];
            let x = &x[jj * n..];
            match tma[jj * bn + ii] {
                Tile::Diagonal => {
                    for i in 0..n {
                        y[i] += alpha * a[i * n * bn + i] * x[i];
                    }
                }
                Tile::Block => dgemv(Trans::No, n, n, alpha, a, bn * n, x, 1, 1.0, y, 1),
                Tile::None => {}
            }
        }
    }
}

/// Returns the one-based index of the first pivot which is not positive or `0` on success.
pub fn dpottrf_diag(n: usize, bn: usize, a: &mut [f64]) -> usize {
    for i in 0..bn * n {
        // Not positive definite, or NaN.
        if !(a[i] > 0.0) {
            return i + 1;
        }
        a[i] = libm::sqrt(a[i]);
    }
    0
}

/// Returns the one-based index of the first pivot which is not positive or `0` on success.
pub fn dpottrf_blockdiag(n: usize, bn: usize, a: &mut [f64]) -> usize {
    for ii in 0..bn {
        let info = dpotrf(Uplo::Lower, n, &mut a[ii * (n * n)..], n);
        if info != 0 {
            return info;
        }
    }
    0
}

/// Tiles are diagonal: One independent `bn` by `bn` Cholesky factorization per sample.
///
/// Returns the one-based index of the first pivot which is not positive or `0` on success.
pub fn dpottrf_diags(n: usize, bn: usize, a: &mut [f64]) -> usize {
    for i in 0..n {
        for jj in 0..bn {
            for ii in jj..bn {
                let mut d = a[(jj * bn + ii) * n + i];
                for kk in 0..jj {
                    d -= a[(kk * bn + ii) * n + i] * a[(kk * bn + jj) * n + i];
                }
                if ii > jj {
                    a[(jj * bn + ii) * n + i] = d / a[(jj * bn + jj) * n + i];
                } else {
                    // Not positive definite, or NaN.
                    if !(d > 0.0) {
                        return jj + 1;
                    }
                    a[(jj * bn + jj) * n + i] = libm::sqrt(d);
                }
            }
        }
    }
    0
}

/// Returns the one-based index of the first pivot which is not positive or `0` on success.
pub fn dpottrf_full(n: usize, bn: usize, a: &mut [f64]) -> usize {
    dpotrf(Uplo::Lower, bn * n, a, bn * n)
}

pub fn dtrtsv_diag(n: usize, bn: usize, a: &[f64], b: &mut [f64]) {
    for i in 0..bn * n {
        b[i] /= a[i];
    }
}

pub fn dtrtsv_blockdiag(n: usize, bn: usize, a: &[f64], b: &mut [f64]) {
    for ii in 0..bn {
        dtrsv(
            Uplo::Lower,
            Trans::No,
            Diag::No,
            n,
            &a[ii * (n * n)..],
            n,
            &mut b[ii * n..],
            1,
        );
    }
}

/// Tiles are diagonal: One independent `bn` by `bn` forward substitution per sample.
pub fn dtrtsv_diags(n: usize, bn: usize, a: &[f64], b: &mut [f64]) {
    for i in 0..n {
        for ii in 0..bn {
            let mut d = b[ii * n + i];
            for jj in 0..ii {
                d -= a[(jj * bn + ii) * n + i] * b[jj * n + i];
            }
            b[ii * n + i] = d / a[(ii * bn + ii) * n + i];
        }
    }
}

pub fn dtrtsv_full(n: usize, bn: usize, a: &[f64], b: &mut [f64]) {
    dtrsv(Uplo::Lower, Trans::No, Diag::No, bn * n, a, bn * n, b, 1);
}

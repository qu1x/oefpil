//! Backend in plain Rust without dependencies besides `libm`.
//!
//! It relies on the compiler to vectorize its kernels for the baseline of the target, whereas the
//! recursive factorization and combining columns keeps memory accesses local. This makes it as fast
//! as the `faer` backend unless the covariance matrix is block-diagonal or full and large.

#![allow(clippy::needless_range_loop)]

use super::{Diag, Side, Trans, Uplo};
use alloc::{vec, vec::Vec};

/// Workspace requirement of [`dgemqr()`], which needs none.
#[derive(Debug, Clone, Copy)]
pub struct WorkspaceReq;

impl WorkspaceReq {
    /// Returns the requirement sufficient for this and `other`.
    pub const fn or(self, _other: Self) -> Self {
        self
    }
}

/// Workspace of [`dgemqr()`], which needs none.
#[derive(Debug)]
pub struct Workspace;

/// Order below which factorizations stop recursing.
const RECURSION_LIMIT: usize = 32;

/// Returns `x^T * y` with independent accumulators as a floating-point sum is not reorderable.
#[inline]
fn dot(x: &[f64], y: &[f64]) -> f64 {
    let (x_chunks, y_chunks) = (x.chunks_exact(4), y.chunks_exact(4));
    let (x_rest, y_rest) = (x_chunks.remainder(), y_chunks.remainder());
    let mut sum = [0f64; 4];
    for (x, y) in x_chunks.zip(y_chunks) {
        for i in 0..4 {
            sum[i] += x[i] * y[i];
        }
    }
    let mut sum = sum[0] + sum[1] + sum[2] + sum[3];
    for (x, y) in x_rest.iter().zip(y_rest) {
        sum += x * y;
    }
    sum
}

/// Computes `y = alpha * x + y`.
#[inline]
fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (y, x) in y.iter_mut().zip(x) {
        *y += alpha * x;
    }
}

/// Computes `y = alpha * y` without reading `y` if `alpha` is zero.
fn scale(alpha: f64, y: &mut [f64]) {
    if alpha != 1.0 {
        for y in y {
            *y = if alpha == 0.0 { 0.0 } else { alpha * *y };
        }
    }
}

/// Copies a vector with increment `inc`.
fn gather(x: &[f64], n: usize, inc: usize) -> Vec<f64> {
    (0..n).map(|i| x[i * inc]).collect()
}

/// Copies a vector back to increment `inc`.
fn scatter(v: &[f64], x: &mut [f64], inc: usize) {
    for (i, v) in v.iter().enumerate() {
        x[i * inc] = *v;
    }
}

/// Subtracts the first `count` columns of `left` from `column`, starting at their rows `start` and
/// scaled by their fields at row `row`.
///
/// Combines four columns per pass to touch memory less often than one [`axpy()`] per column.
#[inline]
fn update(column: &mut [f64], left: &[f64], lda: usize, count: usize, row: usize, start: usize) {
    let len = column.len();
    let mut k = 0;
    while k + 4 <= count {
        let (c0, c1) = (&left[k * lda..], &left[(k + 1) * lda..]);
        let (c2, c3) = (&left[(k + 2) * lda..], &left[(k + 3) * lda..]);
        let (f0, f1, f2, f3) = (c0[row], c1[row], c2[row], c3[row]);
        let (c0, c1) = (&c0[start..][..len], &c1[start..][..len]);
        let (c2, c3) = (&c2[start..][..len], &c3[start..][..len]);
        for i in 0..len {
            column[i] -= f0 * c0[i] + f1 * c1[i] + f2 * c2[i] + f3 * c3[i];
        }
        k += 4;
    }
    while k < count {
        let c = &left[k * lda..];
        axpy(-c[row], &c[start..][..len], column);
        k += 1;
    }
}

/// `DDOT`: Returns `x^T * y`.
pub fn ddot(n: usize, x: &[f64], incx: usize, y: &[f64], incy: usize) -> f64 {
    if n == 0 {
        0.0
    } else if incx == 1 && incy == 1 {
        dot(&x[..n], &y[..n])
    } else {
        dot(&gather(x, n, incx), &gather(y, n, incy))
    }
}

/// `DGEMV`: Computes `y = alpha * op(A) * x + beta * y`.
pub fn dgemv(
    trans: Trans,
    m: usize,
    n: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    x: &[f64],
    incx: usize,
    beta: f64,
    y: &mut [f64],
    incy: usize,
) {
    let (rows, columns) = match trans {
        Trans::No => (m, n),
        Trans::Yes => (n, m),
    };
    if rows == 0 {
        return;
    }
    let mut copy = if incy == 1 {
        Vec::new()
    } else {
        gather(y, rows, incy)
    };
    let out: &mut [f64] = if incy == 1 { &mut y[..rows] } else { &mut copy };
    scale(beta, out);
    if alpha != 0.0 && columns > 0 {
        let owned;
        let x: &[f64] = if incx == 1 {
            &x[..columns]
        } else {
            owned = gather(x, columns, incx);
            &owned
        };
        match trans {
            Trans::No => {
                for j in 0..n {
                    axpy(alpha * x[j], &a[j * lda..][..m], out);
                }
            }
            Trans::Yes => {
                for j in 0..n {
                    out[j] += alpha * dot(&a[j * lda..][..m], x);
                }
            }
        }
    }
    if incy != 1 {
        scatter(&copy, y, incy);
    }
}

/// `DTRMV`: Computes `x = op(A) * x`.
pub fn dtrmv(
    uplo: Uplo,
    trans: Trans,
    diag: Diag,
    n: usize,
    a: &[f64],
    lda: usize,
    x: &mut [f64],
    incx: usize,
) {
    let mut v = gather(x, n, incx);
    let pivot = |j: usize| {
        if diag == Diag::Unit {
            1.0
        } else {
            a[j * lda + j]
        }
    };
    match (uplo, trans) {
        (Uplo::Lower, Trans::No) => {
            for j in (0..n).rev() {
                let field = v[j];
                axpy(field, &a[j * lda + j + 1..][..n - j - 1], &mut v[j + 1..]);
                v[j] = pivot(j) * field;
            }
        }
        (Uplo::Lower, Trans::Yes) => {
            for j in 0..n {
                v[j] = pivot(j) * v[j] + dot(&a[j * lda + j + 1..][..n - j - 1], &v[j + 1..]);
            }
        }
        (Uplo::Upper, Trans::No) => {
            for j in 0..n {
                let field = v[j];
                axpy(field, &a[j * lda..][..j], &mut v[..j]);
                v[j] = pivot(j) * field;
            }
        }
        (Uplo::Upper, Trans::Yes) => {
            for j in (0..n).rev() {
                v[j] = pivot(j) * v[j] + dot(&a[j * lda..][..j], &v[..j]);
            }
        }
    }
    scatter(&v, x, incx);
}

/// Solves `op(A) * x = b` in place with dense `b`.
fn solve(uplo: Uplo, trans: Trans, diag: Diag, n: usize, a: &[f64], lda: usize, b: &mut [f64]) {
    let unit = diag == Diag::Unit;
    match (uplo, trans) {
        (Uplo::Lower, Trans::No) => {
            for j in 0..n {
                if !unit {
                    b[j] /= a[j * lda + j];
                }
                let field = b[j];
                axpy(-field, &a[j * lda + j + 1..][..n - j - 1], &mut b[j + 1..n]);
            }
        }
        (Uplo::Lower, Trans::Yes) => {
            for j in (0..n).rev() {
                b[j] -= dot(&a[j * lda + j + 1..][..n - j - 1], &b[j + 1..n]);
                if !unit {
                    b[j] /= a[j * lda + j];
                }
            }
        }
        (Uplo::Upper, Trans::No) => {
            for j in (0..n).rev() {
                if !unit {
                    b[j] /= a[j * lda + j];
                }
                let field = b[j];
                axpy(-field, &a[j * lda..][..j], &mut b[..j]);
            }
        }
        (Uplo::Upper, Trans::Yes) => {
            for j in 0..n {
                b[j] -= dot(&a[j * lda..][..j], &b[..j]);
                if !unit {
                    b[j] /= a[j * lda + j];
                }
            }
        }
    }
}

/// `DTRSM`: Solves `op(A) * X = alpha * B` or `X * op(A) = alpha * B` in place.
pub fn dtrsm(
    side: Side,
    uplo: Uplo,
    trans: Trans,
    diag: Diag,
    m: usize,
    n: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
) {
    for j in 0..n {
        scale(alpha, &mut b[j * ldb..][..m]);
    }
    if alpha == 0.0 {
        return;
    }
    match side {
        Side::Left => {
            for j in 0..n {
                solve(uplo, trans, diag, m, a, lda, &mut b[j * ldb..][..m]);
            }
        }
        Side::Right => {
            // `X * op(A) = B` is `op(A)^T * X^T = B^T`.
            let trans = match trans {
                Trans::No => Trans::Yes,
                Trans::Yes => Trans::No,
            };
            let mut row = vec![0f64; n];
            for i in 0..m {
                for j in 0..n {
                    row[j] = b[j * ldb + i];
                }
                solve(uplo, trans, diag, n, a, lda, &mut row);
                for j in 0..n {
                    b[j * ldb + i] = row[j];
                }
            }
        }
    }
}

/// `DTRSV`: Solves `op(A) * x = b` in place.
pub fn dtrsv(
    uplo: Uplo,
    trans: Trans,
    diag: Diag,
    n: usize,
    a: &[f64],
    lda: usize,
    x: &mut [f64],
    incx: usize,
) {
    if incx == 1 {
        solve(uplo, trans, diag, n, a, lda, &mut x[..n]);
    } else {
        let mut v = gather(x, n, incx);
        solve(uplo, trans, diag, n, a, lda, &mut v);
        scatter(&v, x, incx);
    }
}

/// Recursively computes the lower Cholesky factorization.
fn dpotrf_lower(n: usize, a: &mut [f64], lda: usize) -> usize {
    if n <= RECURSION_LIMIT {
        for j in 0..n {
            let (left, right) = a.split_at_mut(j * lda);
            let column = &mut right[j..n];
            update(column, left, lda, j, j, j);
            // Not positive definite, or NaN.
            #[allow(clippy::neg_cmp_op_on_partial_ord)]
            if !(column[0] > 0.0) {
                return j + 1;
            }
            let pivot = libm::sqrt(column[0]);
            column[0] = pivot;
            let inverse = 1.0 / pivot;
            for field in &mut column[1..] {
                *field *= inverse;
            }
        }
        return 0;
    }
    let (n1, n2) = (n / 2, n - n / 2);
    let info = dpotrf_lower(n1, a, lda);
    if info != 0 {
        return info;
    }
    // A21 = A21 * L11^-T
    for j in 0..n1 {
        let (left, right) = a.split_at_mut(j * lda);
        update(&mut right[n1..n], left, lda, j, j, n1);
        let inverse = 1.0 / right[j];
        for field in &mut right[n1..n] {
            *field *= inverse;
        }
    }
    // A22 = A22 - A21 * A21^T
    let (left, right) = a.split_at_mut(n1 * lda);
    for j in 0..n2 {
        update(
            &mut right[j * lda + n1 + j..][..n2 - j],
            left,
            lda,
            n1,
            n1 + j,
            n1 + j,
        );
    }
    match dpotrf_lower(n2, &mut a[n1 * lda + n1..], lda) {
        0 => 0,
        info => n1 + info,
    }
}

/// `DPOTRF`: Computes the Cholesky factorization in place.
///
/// The opposite triangle is not referenced. Returns the one-based index of the first pivot which is
/// not positive or `0` on success.
pub fn dpotrf(uplo: Uplo, n: usize, a: &mut [f64], lda: usize) -> usize {
    match uplo {
        Uplo::Lower => dpotrf_lower(n, a, lda),
        Uplo::Upper => {
            for j in 0..n {
                let (left, right) = a.split_at_mut(j * lda);
                for i in 0..j {
                    let (head, tail) = right.split_at_mut(i);
                    tail[0] = (tail[0] - dot(&left[i * lda..][..i], head)) / left[i * lda + i];
                }
                let pivot = right[j] - dot(&right[..j], &right[..j]);
                // Not positive definite, or NaN.
                #[allow(clippy::neg_cmp_op_on_partial_ord)]
                if !(pivot > 0.0) {
                    return j + 1;
                }
                right[j] = libm::sqrt(pivot);
            }
            0
        }
    }
}

/// `DPOTRI`: Computes the inverse from the Cholesky factorization in place.
///
/// The opposite triangle is not referenced. Returns the one-based index of the first pivot which is
/// zero and leaves the matrix unchanged or `0` on success.
pub fn dpotri(uplo: Uplo, n: usize, a: &mut [f64], lda: usize) -> usize {
    if let Some(pivot) = (0..n).find(|&pivot| a[pivot * lda + pivot] == 0.0) {
        return pivot + 1;
    }
    // The matrix is `U^T * U` for both triangles and its inverse `U^-1 * U^-T`.
    let mut u = vec![0f64; n * n];
    for j in 0..n {
        for i in 0..=j {
            u[j * n + i] = match uplo {
                Uplo::Upper => a[j * lda + i],
                Uplo::Lower => a[i * lda + j],
            };
        }
    }
    // U = U^-1 like `DTRTI2`.
    for j in 0..n {
        let pivot = 1.0 / u[j * n + j];
        u[j * n + j] = pivot;
        let (left, right) = u.split_at_mut(j * n);
        let column = &mut right[..j];
        for i in 0..j {
            let field = column[i];
            axpy(field, &left[i * n..][..i], &mut column[..i]);
            column[i] = field * left[i * n + i];
        }
        for field in column {
            *field *= -pivot;
        }
    }
    // A = U * U^T like `DLAUUM`.
    for j in 0..n {
        for i in 0..=j {
            let mut sum = 0.0;
            for k in j..n {
                sum += u[k * n + i] * u[k * n + j];
            }
            match uplo {
                Uplo::Upper => a[j * lda + i] = sum,
                Uplo::Lower => a[i * lda + j] = sum,
            }
        }
    }
    0
}

/// `DNRM2`: Returns the Euclidean norm of `x` without unnecessary overflow or underflow.
fn dnrm2(x: &[f64]) -> f64 {
    const TINY: f64 = f64::MIN_POSITIVE / f64::EPSILON;
    let square = dot(x, x);
    if square.is_nan() || (square.is_finite() && square >= TINY) {
        return libm::sqrt(square);
    }
    let scale = x.iter().fold(0f64, |scale, x| scale.max(x.abs()));
    if scale == 0.0 || scale.is_infinite() {
        return scale;
    }
    let mut square = 0.0;
    for x in x {
        let x = x / scale;
        square += x * x;
    }
    scale * libm::sqrt(square)
}

/// Generates the Householder reflector of `(alpha, x)`, overwriting `alpha` with the resulting field
/// and `x` with the reflector, and returns its factor.
///
/// Adds the squares directly unless this would not be accurate, see [`reflector_scaled()`].
fn reflector(alpha: &mut f64, x: &mut [f64]) -> f64 {
    const TINY: f64 = f64::MIN_POSITIVE / f64::EPSILON;
    let below = dot(x, x);
    let square = *alpha * *alpha + below;
    if below >= TINY && square.is_finite() {
        let beta = -libm::sqrt(square).copysign(*alpha);
        let inverse = 1.0 / (*alpha - beta);
        for x in &mut *x {
            *x *= inverse;
        }
        let tau = (beta - *alpha) / beta;
        *alpha = beta;
        return tau;
    }
    reflector_scaled(alpha, x)
}

/// Generates the Householder reflector like [`reflector()`] at extreme scales.
///
/// Scales by the largest magnitude as in Golub and Van Loan, as both the sums of squares and the
/// fields of the reflector would not be accurate at extreme scales otherwise.
fn reflector_scaled(alpha: &mut f64, x: &mut [f64]) -> f64 {
    // Written such that NaN is not considered zero and hence propagates.
    if x.iter().all(|x| *x == 0.0) {
        return 0.0;
    }
    let largest = x
        .iter()
        .fold(alpha.abs(), |largest, x| largest.max(x.abs()));
    let scaled = *alpha / largest;
    for x in &mut *x {
        *x /= largest;
    }
    let beta = -libm::hypot(scaled, dnrm2(x)).copysign(scaled);
    let inverse = 1.0 / (scaled - beta);
    for x in &mut *x {
        *x *= inverse;
    }
    *alpha = beta * largest;
    (beta - scaled) / beta
}

/// `DGEQR`: Computes the QR factorization in place with `min(m, n)` reflectors stored in `t`.
pub fn dgeqr(m: usize, n: usize, a: &mut [f64], lda: usize, t: &mut [f64]) {
    for j in 0..m.min(n) {
        let (left, right) = a.split_at_mut((j + 1) * lda);
        let (head, below) = left[j * lda + j..j * lda + m].split_at_mut(1);
        let tau = reflector(&mut head[0], below);
        t[j] = tau;
        if tau != 0.0 {
            for k in 0..n - j - 1 {
                let (first, rest) = right[k * lda + j..k * lda + m].split_at_mut(1);
                let factor = -tau * (first[0] + dot(below, rest));
                first[0] += factor;
                axpy(factor, below, rest);
            }
        }
    }
}

/// `DGEMQR` workspace query: Returns the workspace requirement of [`dgemqr()`].
pub const fn dgemqr_work(_side: Side, _trans: Trans, _m: usize, _n: usize) -> WorkspaceReq {
    WorkspaceReq
}

/// `DGEMQR`: Computes `C = op(Q) * C` or `C = C * op(Q)` with `k` reflectors of [`dgeqr()`].
pub fn dgemqr(
    side: Side,
    trans: Trans,
    m: usize,
    n: usize,
    k: usize,
    a: &[f64],
    lda: usize,
    t: &[f64],
    c: &mut [f64],
    ldc: usize,
    _work: &mut Workspace,
) {
    let order = match side {
        Side::Left => m,
        Side::Right => n,
    };
    // `Q = H_1 * ... * H_k`
    let forward = matches!(
        (side, trans),
        (Side::Left, Trans::Yes) | (Side::Right, Trans::No)
    );
    for step in 0..k {
        let j = if forward { step } else { k - 1 - step };
        let v = &a[j * lda + j + 1..j * lda + order];
        let tau = t[j];
        match side {
            Side::Left => {
                for column in 0..n {
                    let (first, rest) = c[column * ldc + j..column * ldc + m].split_at_mut(1);
                    let factor = -tau * (first[0] + dot(v, rest));
                    first[0] += factor;
                    axpy(factor, v, rest);
                }
            }
            Side::Right => {
                for row in 0..m {
                    let mut sum = c[j * ldc + row];
                    for (i, v) in v.iter().enumerate() {
                        sum += c[(j + 1 + i) * ldc + row] * v;
                    }
                    let factor = -tau * sum;
                    c[j * ldc + row] += factor;
                    for (i, v) in v.iter().enumerate() {
                        c[(j + 1 + i) * ldc + row] += factor * v;
                    }
                }
            }
        }
    }
}

/// Allocates the workspace of [`dgemqr()`] for `req`.
pub const fn workspace(_req: WorkspaceReq) -> Workspace {
    Workspace
}

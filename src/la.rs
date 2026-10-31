//! BLAS and LAPACK routines required by the algorithm, implemented on [`faer`].
//!
//! Each function mirrors the routine of the same name for comparison with the reference
//! implementation: It takes the same arguments in the same order and expects column-major matrices
//! with leading dimensions and vectors with increments. It differs as follows:
//!
//!   * Character arguments are enums.
//!   * `INFO` is returned if the routine has failure modes.
//!   * The workspace of [`dgemqr()`] is queried via [`dgemqr_work()`] instead of `LWORK = -1`.
//!   * Increments must be positive.
//!   * The factor `T` of [`dgeqr()`] and [`dgemqr()`] is their private business and hence has its own
//!     layout. Its size is given by the number of reflectors instead of being queried.
//!
//! Only this module knows about [`faer`] which must not leak into the public API.

#![allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::suboptimal_flops,
    clippy::too_many_arguments
)]

use alloc::{vec, vec::Vec};
use faer::{
    Accum, Conj, Mat, MatMut, MatRef, Par, Spec,
    dyn_stack::{MemBuffer, MemStack, StackReq},
    linalg::{
        cholesky::llt::factor::{
            LltError, LltRegularization, cholesky_in_place, cholesky_in_place_scratch,
        },
        householder::{
            apply_block_householder_sequence_on_the_left_in_place_scratch,
            apply_block_householder_sequence_on_the_left_in_place_with_conj,
            apply_block_householder_sequence_on_the_right_in_place_scratch,
            apply_block_householder_sequence_on_the_right_in_place_with_conj,
            apply_block_householder_sequence_transpose_on_the_left_in_place_scratch,
            apply_block_householder_sequence_transpose_on_the_left_in_place_with_conj,
            apply_block_householder_sequence_transpose_on_the_right_in_place_scratch,
            apply_block_householder_sequence_transpose_on_the_right_in_place_with_conj,
            make_householder_in_place,
        },
        matmul::{
            dot::inner_prod,
            matmul,
            triangular::{self, BlockStructure},
        },
        triangular_inverse::invert_upper_triangular,
        triangular_solve::{
            solve_lower_triangular_in_place, solve_unit_lower_triangular_in_place,
            solve_unit_upper_triangular_in_place, solve_upper_triangular_in_place,
        },
    },
    prelude::{Reborrow, ReborrowMut},
};

/// Argument `SIDE`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(not(test), allow(dead_code))]
pub enum Side {
    /// `'L'`
    Left,
    /// `'R'`
    Right,
}

/// Argument `UPLO`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Uplo {
    /// `'L'`
    Lower,
    /// `'U'`
    Upper,
}

/// Argument `TRANS`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trans {
    /// `'N'`
    No,
    /// `'T'`
    Yes,
}

/// Argument `DIAG`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(not(test), allow(dead_code))]
pub enum Diag {
    /// `'N'`
    No,
    /// `'U'`
    Unit,
}

impl Uplo {
    const fn transpose(self) -> Self {
        match self {
            Self::Lower => Self::Upper,
            Self::Upper => Self::Lower,
        }
    }
}

impl Trans {
    const fn transpose(self) -> Self {
        match self {
            Self::No => Self::Yes,
            Self::Yes => Self::No,
        }
    }
}

fn matrix(a: &[f64], m: usize, n: usize, lda: usize) -> MatRef<'_, f64> {
    MatRef::from_column_major_slice_with_stride(a, m, n, lda)
}

fn matrix_mut(a: &mut [f64], m: usize, n: usize, lda: usize) -> MatMut<'_, f64> {
    MatMut::from_column_major_slice_with_stride_mut(a, m, n, lda)
}

/// Views a vector as a single column.
fn vector(x: &[f64], n: usize, incx: usize) -> MatRef<'_, f64> {
    MatRef::from_column_major_slice_with_stride(x, 1, n, incx).transpose()
}

/// Views a vector as a single column.
fn vector_mut(x: &mut [f64], n: usize, incx: usize) -> MatMut<'_, f64> {
    MatMut::from_column_major_slice_with_stride_mut(x, 1, n, incx).transpose_mut()
}

fn scale(alpha: f64, mut a: MatMut<'_, f64>) {
    if alpha != 1.0 {
        for column in a.rb_mut().col_iter_mut() {
            for field in column.iter_mut() {
                *field = if alpha == 0.0 { 0.0 } else { alpha * *field };
            }
        }
    }
}

/// Solves `op(A) * X = B` in place.
fn solve_on_the_left(uplo: Uplo, trans: Trans, diag: Diag, a: MatRef<'_, f64>, b: MatMut<'_, f64>) {
    let (a, uplo) = match trans {
        Trans::No => (a, uplo),
        Trans::Yes => (a.transpose(), uplo.transpose()),
    };
    match (uplo, diag) {
        (Uplo::Lower, Diag::No) => solve_lower_triangular_in_place(a, b, Par::Seq),
        (Uplo::Lower, Diag::Unit) => solve_unit_lower_triangular_in_place(a, b, Par::Seq),
        (Uplo::Upper, Diag::No) => solve_upper_triangular_in_place(a, b, Par::Seq),
        (Uplo::Upper, Diag::Unit) => solve_unit_upper_triangular_in_place(a, b, Par::Seq),
    }
}

/// `DDOT`: Returns `x^T * y`.
pub fn ddot(n: usize, x: &[f64], incx: usize, y: &[f64], incy: usize) -> f64 {
    inner_prod(
        vector(x, n, incx).transpose().row(0),
        Conj::No,
        vector(y, n, incy).col(0),
        Conj::No,
    )
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
    let a = matrix(a, m, n, lda);
    let a = match trans {
        Trans::No => a,
        Trans::Yes => a.transpose(),
    };
    let x = vector(x, a.ncols(), incx);
    let mut y = vector_mut(y, a.nrows(), incy);
    scale(beta, y.rb_mut());
    if alpha != 0.0 {
        matmul(y, Accum::Add, a, x, alpha, Par::Seq);
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
    let a = matrix(a, n, n, lda);
    let (a, uplo) = match trans {
        Trans::No => (a, uplo),
        Trans::Yes => (a.transpose(), uplo.transpose()),
    };
    let structure = match (uplo, diag) {
        (Uplo::Lower, Diag::No) => BlockStructure::TriangularLower,
        (Uplo::Lower, Diag::Unit) => BlockStructure::UnitTriangularLower,
        (Uplo::Upper, Diag::No) => BlockStructure::TriangularUpper,
        (Uplo::Upper, Diag::Unit) => BlockStructure::UnitTriangularUpper,
    };
    let mut x = vector_mut(x, n, incx);
    let copy = x.rb().to_owned();
    triangular::matmul(
        x.rb_mut(),
        BlockStructure::Rectangular,
        Accum::Replace,
        a,
        structure,
        copy.as_ref(),
        BlockStructure::Rectangular,
        1.0,
        Par::Seq,
    );
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
    let order = match side {
        Side::Left => m,
        Side::Right => n,
    };
    let a = matrix(a, order, order, lda);
    let mut b = matrix_mut(b, m, n, ldb);
    scale(alpha, b.rb_mut());
    if alpha != 0.0 {
        match side {
            Side::Left => solve_on_the_left(uplo, trans, diag, a, b),
            Side::Right => solve_on_the_left(uplo, trans.transpose(), diag, a, b.transpose_mut()),
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
    solve_on_the_left(
        uplo,
        trans,
        diag,
        matrix(a, n, n, lda),
        vector_mut(x, n, incx),
    );
}

/// `DPOTRF`: Computes the Cholesky factorization in place.
///
/// The opposite triangle is not referenced. Returns the one-based index of the first pivot which is
/// not positive or `0` on success.
pub fn dpotrf(uplo: Uplo, n: usize, a: &mut [f64], lda: usize) -> usize {
    let a = matrix_mut(a, n, n, lda);
    let a = match uplo {
        Uplo::Lower => a,
        Uplo::Upper => a.transpose_mut(),
    };
    let mut work = MemBuffer::new(cholesky_in_place_scratch::<f64>(
        n,
        Par::Seq,
        Spec::default(),
    ));
    cholesky_in_place(
        a,
        LltRegularization::default(),
        Par::Seq,
        MemStack::new(&mut work),
        Spec::default(),
    )
    .map_or_else(|LltError::NonPositivePivot { index }| index + 1, |_info| 0)
}

/// `DPOTRI`: Computes the inverse from the Cholesky factorization in place.
///
/// The opposite triangle is not referenced. Returns the one-based index of the first pivot which is
/// zero and leaves the matrix unchanged or `0` on success.
pub fn dpotri(uplo: Uplo, n: usize, a: &mut [f64], lda: usize) -> usize {
    let a = matrix_mut(a, n, n, lda);
    // The inverse of a lower factor is the one of its transposed upper factor.
    let mut a = match uplo {
        Uplo::Lower => a.transpose_mut(),
        Uplo::Upper => a,
    };
    if let Some(pivot) = (0..n).find(|&pivot| a[(pivot, pivot)] == 0.0) {
        return pivot + 1;
    }
    let mut inverse = Mat::<f64>::zeros(n, n);
    invert_upper_triangular(inverse.as_mut(), a.rb(), Par::Seq);
    triangular::matmul(
        a.rb_mut(),
        BlockStructure::TriangularUpper,
        Accum::Replace,
        inverse.as_ref(),
        BlockStructure::TriangularUpper,
        inverse.transpose(),
        BlockStructure::TriangularLower,
        1.0,
        Par::Seq,
    );
    0
}

/// `DGEQR`: Computes the QR factorization in place with `min(m, n)` reflectors stored in `t`.
///
/// Deliberately not on `faer::linalg::qr::no_pivoting::factor::qr_in_place()` which estimates the
/// rank and skips columns it considers dependent. This yields a row-echelon factor for
/// ill-conditioned matrices whereas the algorithm relies on an upper triangular factor of any
/// condition.
pub fn dgeqr(m: usize, n: usize, a: &mut [f64], lda: usize, t: &mut [f64]) {
    let mut a = matrix_mut(a, m, n, lda);
    for (reflector, t) in t.iter_mut().enumerate().take(m.min(n)) {
        let (mut head, right, below, trailing) =
            a.rb_mut().split_at_mut(reflector + 1, reflector + 1);
        let mut below = below.col_mut(reflector);
        let info = make_householder_in_place(&mut head[(reflector, reflector)], below.rb_mut());
        *t = info.tau;
        let tau_inv = info.tau.recip();
        for (head, mut column) in right
            .row_mut(reflector)
            .iter_mut()
            .zip(trailing.col_iter_mut())
        {
            let dot = *head + inner_prod(below.rb().transpose(), Conj::No, column.rb(), Conj::No);
            let factor = -dot * tau_inv;
            *head += factor;
            for (field, below) in column.rb_mut().iter_mut().zip(below.rb().iter()) {
                *field += factor * below;
            }
        }
    }
}

/// `DGEMQR` workspace query: Returns the workspace requirement of [`dgemqr()`].
pub fn dgemqr_work(side: Side, trans: Trans, m: usize, n: usize) -> StackReq {
    match (side, trans) {
        (Side::Left, Trans::No) => {
            apply_block_householder_sequence_on_the_left_in_place_scratch::<f64>(m, 1, n)
        }
        (Side::Left, Trans::Yes) => {
            apply_block_householder_sequence_transpose_on_the_left_in_place_scratch::<f64>(m, 1, n)
        }
        (Side::Right, Trans::No) => {
            apply_block_householder_sequence_on_the_right_in_place_scratch::<f64>(n, 1, m)
        }
        (Side::Right, Trans::Yes) => {
            apply_block_householder_sequence_transpose_on_the_right_in_place_scratch::<f64>(n, 1, m)
        }
    }
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
    work: &mut MemStack,
) {
    let order = match side {
        Side::Left => m,
        Side::Right => n,
    };
    let a = matrix(a, order, k, lda);
    let t = MatRef::from_column_major_slice(&t[..k], 1, k);
    let c = matrix_mut(c, m, n, ldc);
    match (side, trans) {
        (Side::Left, Trans::No) => apply_block_householder_sequence_on_the_left_in_place_with_conj(
            a,
            t,
            Conj::No,
            c,
            Par::Seq,
            work,
        ),
        (Side::Left, Trans::Yes) => {
            apply_block_householder_sequence_transpose_on_the_left_in_place_with_conj(
                a,
                t,
                Conj::No,
                c,
                Par::Seq,
                work,
            );
        }
        (Side::Right, Trans::No) => {
            apply_block_householder_sequence_on_the_right_in_place_with_conj(
                a,
                t,
                Conj::No,
                c,
                Par::Seq,
                work,
            );
        }
        (Side::Right, Trans::Yes) => {
            apply_block_householder_sequence_transpose_on_the_right_in_place_with_conj(
                a,
                t,
                Conj::No,
                c,
                Par::Seq,
                work,
            );
        }
    }
}

/// Allocates the workspace of [`dgemqr()`] for `req`.
pub fn workspace(req: StackReq) -> MemBuffer {
    MemBuffer::new(req)
}

/// Allocates the factor `T` of [`dgeqr()`].
pub fn dgeqr_t(m: usize, n: usize) -> Vec<f64> {
    vec![0.0; m.min(n)]
}

#[cfg(test)]
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::needless_pass_by_value,
    clippy::too_many_lines,
    clippy::tuple_array_conversions
)]
mod tests {
    use super::*;
    use alloc::format;
    use core::ffi::{c_char, c_int};
    use oefpil_sys as sys;

    /// Deterministic pseudorandom numbers in `-1.0..1.0`.
    struct Random(u64);

    impl Random {
        fn field(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            #[allow(clippy::cast_precision_loss)]
            let unit = (self.0 >> 11) as f64 / (1u64 << 53) as f64;
            2.0 * unit - 1.0
        }
        fn matrix(&mut self, fields: usize) -> Vec<f64> {
            (0..fields).map(|_| self.field()).collect()
        }
    }

    const fn int(value: usize) -> c_int {
        #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
        {
            value as c_int
        }
    }

    fn side_c(side: Side) -> *const c_char {
        match side {
            Side::Left => c"L".as_ptr(),
            Side::Right => c"R".as_ptr(),
        }
    }

    fn uplo_c(uplo: Uplo) -> *const c_char {
        match uplo {
            Uplo::Lower => c"L".as_ptr(),
            Uplo::Upper => c"U".as_ptr(),
        }
    }

    fn trans_c(trans: Trans) -> *const c_char {
        match trans {
            Trans::No => c"N".as_ptr(),
            Trans::Yes => c"T".as_ptr(),
        }
    }

    fn diag_c(diag: Diag) -> *const c_char {
        match diag {
            Diag::No => c"N".as_ptr(),
            Diag::Unit => c"U".as_ptr(),
        }
    }

    const UPLOS: [Uplo; 2] = [Uplo::Lower, Uplo::Upper];
    const TRANSES: [Trans; 2] = [Trans::No, Trans::Yes];
    const DIAGS: [Diag; 2] = [Diag::No, Diag::Unit];
    const SIDES: [Side; 2] = [Side::Left, Side::Right];
    /// Orders across thresholds of blocking and recursion in the backend.
    const ORDERS: [usize; 8] = [0, 1, 2, 3, 7, 17, 64, 130];

    #[track_caller]
    fn assert_close(actual: &[f64], expected: &[f64], tolerance: f64, what: &str) {
        assert_eq!(expected.len(), actual.len(), "{what}");
        let scale = expected.iter().fold(1f64, |scale, e| scale.max(e.abs()));
        for (field, (a, e)) in actual.iter().zip(expected).enumerate() {
            assert!(
                (a - e).abs() <= tolerance * scale,
                "{what}: field {field}: {a} != {e}",
            );
        }
    }

    /// Well-conditioned triangle in the full matrix, garbage in the opposite triangle.
    fn triangular(random: &mut Random, n: usize, lda: usize, uplo: Uplo) -> Vec<f64> {
        let mut a = random.matrix(lda * n);
        for j in 0..n {
            a[j * lda + j] = 2.0 * n as f64 + 1.0;
            for i in 0..n {
                let referenced = match uplo {
                    Uplo::Lower => i >= j,
                    Uplo::Upper => i <= j,
                };
                if !referenced {
                    a[j * lda + i] = f64::NAN;
                }
            }
        }
        a
    }

    /// Symmetric positive definite matrix with NaN in the opposite triangle.
    fn definite(random: &mut Random, n: usize, lda: usize, uplo: Uplo) -> Vec<f64> {
        let mut a = vec![0f64; lda * n];
        let b = random.matrix(n * n);
        for j in 0..n {
            for i in 0..n {
                let field = (0..n).map(|k| b[k * n + i] * b[k * n + j]).sum::<f64>();
                let field = if i == j { field + 1.0 } else { field };
                let referenced = match uplo {
                    Uplo::Lower => i >= j,
                    Uplo::Upper => i <= j,
                };
                a[j * lda + i] = if referenced { field } else { f64::NAN };
            }
        }
        a
    }

    #[test]
    fn ddot_agrees() {
        let mut random = Random(1);
        for n in ORDERS {
            for (incx, incy) in [(1, 1), (2, 3)] {
                let x = random.matrix((n * incx).max(1));
                let y = random.matrix((n * incy).max(1));
                let expected =
                    unsafe { sys::ddot(int(n), x.as_ptr(), int(incx), y.as_ptr(), int(incy)) };
                assert_close(&[ddot(n, &x, incx, &y, incy)], &[expected], 1e-13, "ddot");
            }
        }
    }

    #[test]
    fn dgemv_agrees() {
        let mut random = Random(2);
        for (m, n) in [(0, 0), (1, 1), (3, 2), (2, 3), (17, 7), (7, 17), (130, 64)] {
            for trans in TRANSES {
                for (alpha, beta) in [(1.0, 1.0), (-1.0, 1.0), (0.5, 0.0), (0.0, 2.0), (2.0, -1.5)]
                {
                    for (incx, incy) in [(1, 1), (2, 3)] {
                        let lda = m + 2;
                        let (rows, columns) = match trans {
                            Trans::No => (m, n),
                            Trans::Yes => (n, m),
                        };
                        let a = random.matrix(lda * n);
                        let x = random.matrix((columns * incx).max(1));
                        let y = random.matrix((rows * incy).max(1));
                        let mut actual = y.clone();
                        let mut expected = y;
                        dgemv(
                            trans,
                            m,
                            n,
                            alpha,
                            &a,
                            lda,
                            &x,
                            incx,
                            beta,
                            &mut actual,
                            incy,
                        );
                        unsafe {
                            sys::dgemv(
                                trans_c(trans),
                                int(m),
                                int(n),
                                alpha,
                                a.as_ptr(),
                                int(lda),
                                x.as_ptr(),
                                int(incx),
                                beta,
                                expected.as_mut_ptr(),
                                int(incy),
                            );
                        }
                        assert_close(&actual, &expected, 1e-13, "dgemv");
                    }
                }
            }
        }
    }

    #[test]
    fn dtrmv_agrees() {
        let mut random = Random(3);
        for n in ORDERS {
            for (uplo, trans, diag) in flags() {
                for incx in [1, 2] {
                    let lda = n + 2;
                    let a = triangular(&mut random, n, lda, uplo);
                    let x = random.matrix((n * incx).max(1));
                    let mut actual = x.clone();
                    let mut expected = x;
                    dtrmv(uplo, trans, diag, n, &a, lda, &mut actual, incx);
                    unsafe {
                        sys::dtrmv(
                            uplo_c(uplo),
                            trans_c(trans),
                            diag_c(diag),
                            int(n),
                            a.as_ptr(),
                            int(lda),
                            expected.as_mut_ptr(),
                            int(incx),
                        );
                    }
                    assert_close(&actual, &expected, 1e-13, "dtrmv");
                }
            }
        }
    }

    fn flags() -> impl Iterator<Item = (Uplo, Trans, Diag)> {
        UPLOS.into_iter().flat_map(|uplo| {
            TRANSES
                .into_iter()
                .flat_map(move |trans| DIAGS.into_iter().map(move |diag| (uplo, trans, diag)))
        })
    }

    #[test]
    fn dtrsv_agrees() {
        let mut random = Random(4);
        for n in ORDERS {
            for (uplo, trans, diag) in flags() {
                for incx in [1, 2] {
                    let lda = n + 2;
                    let a = triangular(&mut random, n, lda, uplo);
                    let x = random.matrix((n * incx).max(1));
                    let mut actual = x.clone();
                    let mut expected = x;
                    dtrsv(uplo, trans, diag, n, &a, lda, &mut actual, incx);
                    unsafe {
                        sys::dtrsv(
                            uplo_c(uplo),
                            trans_c(trans),
                            diag_c(diag),
                            int(n),
                            a.as_ptr(),
                            int(lda),
                            expected.as_mut_ptr(),
                            int(incx),
                        );
                    }
                    assert_close(&actual, &expected, 1e-11, "dtrsv");
                }
            }
        }
    }

    #[test]
    fn dtrsm_agrees() {
        let mut random = Random(5);
        for (m, n) in [(0, 0), (1, 1), (3, 2), (2, 3), (17, 7), (7, 17), (130, 64)] {
            for side in SIDES {
                for (uplo, trans, diag) in flags() {
                    for alpha in [1.0, 0.5, 0.0] {
                        let order = if side == Side::Left { m } else { n };
                        let lda = order + 2;
                        let ldb = m + 3;
                        let a = triangular(&mut random, order, lda, uplo);
                        let b = random.matrix(ldb * n);
                        let mut actual = b.clone();
                        let mut expected = b;
                        dtrsm(
                            side,
                            uplo,
                            trans,
                            diag,
                            m,
                            n,
                            alpha,
                            &a,
                            lda,
                            &mut actual,
                            ldb,
                        );
                        unsafe {
                            sys::dtrsm(
                                side_c(side),
                                uplo_c(uplo),
                                trans_c(trans),
                                diag_c(diag),
                                int(m),
                                int(n),
                                alpha,
                                a.as_ptr(),
                                int(lda),
                                expected.as_mut_ptr(),
                                int(ldb),
                            );
                        }
                        assert_close(&actual, &expected, 1e-11, "dtrsm");
                    }
                }
            }
        }
    }

    #[test]
    fn dpotrf_agrees_and_does_not_reference_opposite_triangle() {
        let mut random = Random(6);
        for n in ORDERS {
            for uplo in UPLOS {
                let lda = n + 2;
                let a = definite(&mut random, n, lda, uplo);
                let mut actual = a.clone();
                let mut expected = a;
                let info = dpotrf(uplo, n, &mut actual, lda);
                let mut expected_info = 0;
                unsafe {
                    sys::dpotrf(
                        uplo_c(uplo),
                        int(n),
                        expected.as_mut_ptr(),
                        int(lda),
                        &raw mut expected_info,
                    );
                }
                assert_eq!(0, expected_info);
                assert_eq!(0, info);
                // Opposite triangle stays NaN, so the comparison of the triangle is exact in shape.
                for j in 0..n {
                    for i in 0..n {
                        let (a, e) = (actual[j * lda + i], expected[j * lda + i]);
                        assert_eq!(e.is_nan(), a.is_nan(), "n {n} {uplo:?} ({i}, {j})");
                        if !e.is_nan() {
                            assert!((a - e).abs() <= 1e-11 * e.abs().max(1.0), "{a} != {e}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn dpotrf_reports_first_failing_pivot() {
        let mut random = Random(7);
        for n in ORDERS {
            for uplo in UPLOS {
                for pivot in [0, n / 2, n.saturating_sub(1)] {
                    for bad in [-1.0, 0.0, f64::NAN] {
                        if n == 0 {
                            continue;
                        }
                        let mut actual = definite(&mut random, n, n, uplo);
                        // Below every leading minor of smaller order remains positive definite.
                        actual[pivot * n + pivot] = bad;
                        let mut expected = actual.clone();
                        let info = dpotrf(uplo, n, &mut actual, n);
                        let mut expected_info = 0;
                        unsafe {
                            sys::dpotrf(
                                uplo_c(uplo),
                                int(n),
                                expected.as_mut_ptr(),
                                int(n),
                                &raw mut expected_info,
                            );
                        }
                        assert_eq!(
                            usize::try_from(expected_info).unwrap(),
                            info,
                            "n {n} {uplo:?} pivot {pivot} bad {bad}",
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn dpotri_agrees_and_does_not_reference_opposite_triangle() {
        let mut random = Random(8);
        for n in ORDERS {
            for uplo in UPLOS {
                let lda = n + 2;
                let mut factor = definite(&mut random, n, lda, uplo);
                assert_eq!(0, dpotrf(uplo, n, &mut factor, lda));
                let mut actual = factor.clone();
                let mut expected = factor;
                let info = dpotri(uplo, n, &mut actual, lda);
                let mut expected_info = 0;
                unsafe {
                    sys::dpotri(
                        uplo_c(uplo),
                        int(n),
                        expected.as_mut_ptr(),
                        int(lda),
                        &raw mut expected_info,
                    );
                }
                assert_eq!(0, expected_info);
                assert_eq!(0, info);
                for j in 0..n {
                    for i in 0..n {
                        let (a, e) = (actual[j * lda + i], expected[j * lda + i]);
                        assert_eq!(e.is_nan(), a.is_nan(), "n {n} {uplo:?} ({i}, {j})");
                        if !e.is_nan() {
                            assert!((a - e).abs() <= 1e-10 * e.abs().max(1.0), "{a} != {e}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn dpotri_reports_zero_pivot_and_leaves_matrix_unchanged() {
        let mut a = [2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 3.0];
        let unchanged = a;
        assert_eq!(2, dpotri(Uplo::Upper, 3, &mut a, 3));
        assert_eq!(unchanged, a);
    }

    /// Fortran factorization `(a, t)` of LAPACK with its workspace.
    fn dgeqr_fortran(m: usize, n: usize, a: &mut [f64], lda: usize) -> Vec<f64> {
        let mut query = [0f64; 5];
        let mut work = [0f64; 1];
        let mut info = 0;
        unsafe {
            sys::dgeqr(
                int(m),
                int(n),
                a.as_mut_ptr(),
                int(lda),
                query.as_mut_ptr(),
                -1,
                work.as_mut_ptr(),
                -1,
                &raw mut info,
            );
        }
        assert_eq!(0, info);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let (tsize, lwork) = (query[0] as usize, work[0] as usize);
        let mut t = vec![0f64; tsize];
        let mut work = vec![0f64; lwork];
        unsafe {
            sys::dgeqr(
                int(m),
                int(n),
                a.as_mut_ptr(),
                int(lda),
                t.as_mut_ptr(),
                int(tsize),
                work.as_mut_ptr(),
                int(lwork),
                &raw mut info,
            );
        }
        assert_eq!(0, info);
        t
    }

    /// Applies Fortran `op(Q)` of `(a, t)` on the left of `c`.
    fn dgemqr_fortran(
        trans: Trans,
        m: usize,
        k: usize,
        a: &[f64],
        lda: usize,
        t: &[f64],
        c: &mut [f64],
    ) {
        let (mut query, mut info) = (0f64, 0);
        unsafe {
            sys::dgemqr(
                c"L".as_ptr(),
                trans_c(trans),
                int(m),
                1,
                int(k),
                a.as_ptr(),
                int(lda),
                t.as_ptr(),
                int(t.len()),
                c.as_mut_ptr(),
                int(m),
                &raw mut query,
                -1,
                &raw mut info,
            );
        }
        assert_eq!(0, info);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let lwork = query as usize;
        let mut work = vec![0f64; lwork];
        unsafe {
            sys::dgemqr(
                c"L".as_ptr(),
                trans_c(trans),
                int(m),
                1,
                int(k),
                a.as_ptr(),
                int(lda),
                t.as_ptr(),
                int(t.len()),
                c.as_mut_ptr(),
                int(m),
                work.as_mut_ptr(),
                int(lwork),
                &raw mut info,
            );
        }
        assert_eq!(0, info);
    }

    /// Fits to `rhs` along the steps of the algorithm, compared up to the signs of the reflectors.
    fn assert_qr_agrees(m: usize, n: usize, a: &[f64], lda: usize, rhs: &[f64], tolerance: f64) {
        let mut factor = a.to_vec();
        let mut t = dgeqr_t(m, n);
        dgeqr(m, n, &mut factor, lda, &mut t);
        let mut fortran = a.to_vec();
        let fortran_t = dgeqr_fortran(m, n, &mut fortran, lda);
        // Same triangular factor: `R^T R` is unique and so are the signs by the same convention.
        for j in 0..n {
            for i in 0..=j.min(m.saturating_sub(1)) {
                let (a, e) = (factor[j * lda + i], fortran[j * lda + i]);
                assert!(
                    (a.abs() - e.abs()).abs() <= tolerance * e.abs().max(1.0),
                    "R ({i}, {j}) of {m} x {n}: {a} != {e}",
                );
            }
        }
        // Same projection onto the range `Q_R Q_R^T` used by the algorithm.
        let k = m.min(n);
        let project = |transposed: &mut dyn FnMut(&mut [f64]), back: &mut dyn FnMut(&mut [f64])| {
            let mut c = rhs.to_vec();
            transposed(&mut c);
            c[k..].fill(0.0);
            back(&mut c);
            c
        };
        let mut work = workspace(dgemqr_work(Side::Left, Trans::Yes, m, 1).or(dgemqr_work(
            Side::Left,
            Trans::No,
            m,
            1,
        )));
        let mut work_t = workspace(dgemqr_work(Side::Left, Trans::Yes, m, 1));
        let actual = project(
            &mut |c| {
                dgemqr(
                    Side::Left,
                    Trans::Yes,
                    m,
                    1,
                    k,
                    &factor,
                    lda,
                    &t,
                    c,
                    m,
                    MemStack::new(&mut work_t),
                );
            },
            &mut |c| {
                dgemqr(
                    Side::Left,
                    Trans::No,
                    m,
                    1,
                    k,
                    &factor,
                    lda,
                    &t,
                    c,
                    m,
                    MemStack::new(&mut work),
                );
            },
        );
        let expected = project(
            &mut |c| dgemqr_fortran(Trans::Yes, m, k, &fortran, lda, &fortran_t, c),
            &mut |c| dgemqr_fortran(Trans::No, m, k, &fortran, lda, &fortran_t, c),
        );
        assert_close(
            &actual,
            &expected,
            tolerance,
            &format!("projection of {m} x {n}"),
        );
    }

    #[test]
    fn qr_agrees() {
        let mut random = Random(9);
        for (m, n) in [
            (1, 1),
            (2, 1),
            (3, 2),
            (5, 5),
            (17, 7),
            (64, 16),
            (130, 40),
            (400, 3),
        ] {
            for lda in [m, m + 2] {
                let a = random.matrix(lda * n);
                let rhs = random.matrix(m);
                assert_qr_agrees(m, n, &a, lda, &rhs, 1e-11);
            }
        }
    }

    /// `faer::linalg::qr::no_pivoting::factor::qr_in_place()` would fail this.
    #[test]
    fn qr_keeps_upper_triangular_factor_of_ill_conditioned_matrix() {
        // Vandermonde matrix whose condition exceeds the rank threshold of the backend.
        let (m, n) = (40, 14);
        let mut a = vec![0f64; m * n];
        for j in 0..n {
            for i in 0..m {
                a[j * m + i] = (1.0 + i as f64 / m as f64).powi(i32::try_from(j).unwrap());
            }
        }
        let mut factor = a.clone();
        let mut t = dgeqr_t(m, n);
        dgeqr(m, n, &mut factor, m, &mut t);
        let mut fortran = a;
        dgeqr_fortran(m, n, &mut fortran, m);
        for pivot in 0..n {
            let (a, e) = (factor[pivot * m + pivot], fortran[pivot * m + pivot]);
            assert!(a != 0.0, "pivot {pivot}");
            // Only the leading pivots are well-determined.
            if pivot < 6 {
                assert!((a.abs() - e.abs()).abs() <= 1e-9 * e.abs(), "{a} != {e}");
            }
        }
    }

    #[test]
    fn dgemqr_of_both_sides_is_orthogonal() {
        let mut random = Random(10);
        let (m, n, k) = (17, 11, 7);
        for side in SIDES {
            let order = if side == Side::Left { m } else { n };
            let lda = order + 2;
            let mut factor = random.matrix(lda * k);
            let mut t = dgeqr_t(order, k);
            dgeqr(order, k, &mut factor, lda, &mut t);
            let c = random.matrix(m * n);
            let mut actual = c.clone();
            let mut work = workspace(dgemqr_work(side, Trans::No, m, n).or(dgemqr_work(
                side,
                Trans::Yes,
                m,
                n,
            )));
            let (first, second) = match side {
                Side::Left => (Trans::Yes, Trans::No),
                Side::Right => (Trans::No, Trans::Yes),
            };
            for trans in [first, second] {
                dgemqr(
                    side,
                    trans,
                    m,
                    n,
                    k,
                    &factor,
                    lda,
                    &t,
                    &mut actual,
                    m,
                    MemStack::new(&mut work),
                );
            }
            assert_close(&actual, &c, 1e-12, "round trip");
            // Agrees with explicit `Q^T Q` being compared against dense `C` through Fortran.
            let mut dense = c.clone();
            let mut fortran = random.matrix(lda * k);
            let mut factor = fortran.clone();
            let fortran_t = dgeqr_fortran(order, k, &mut fortran, lda);
            dgeqr(order, k, &mut factor, lda, &mut t);
            let mut actual = c;
            dgemqr(
                side,
                Trans::Yes,
                m,
                n,
                k,
                &factor,
                lda,
                &t,
                &mut actual,
                m,
                MemStack::new(&mut work),
            );
            let (mut query, mut info) = (0f64, 0);
            unsafe {
                sys::dgemqr(
                    side_c(side),
                    c"T".as_ptr(),
                    int(m),
                    int(n),
                    int(k),
                    fortran.as_ptr(),
                    int(lda),
                    fortran_t.as_ptr(),
                    int(fortran_t.len()),
                    dense.as_mut_ptr(),
                    int(m),
                    &raw mut query,
                    -1,
                    &raw mut info,
                );
            }
            assert_eq!(0, info);
            // The workspace this revision of LAPACK reports for the right side corrupts the heap.
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let mut fortran_work = vec![0f64; query as usize + m * m * n];
            unsafe {
                sys::dgemqr(
                    side_c(side),
                    c"T".as_ptr(),
                    int(m),
                    int(n),
                    int(k),
                    fortran.as_ptr(),
                    int(lda),
                    fortran_t.as_ptr(),
                    int(fortran_t.len()),
                    dense.as_mut_ptr(),
                    int(m),
                    fortran_work.as_mut_ptr(),
                    int(fortran_work.len()),
                    &raw mut info,
                );
            }
            assert_eq!(0, info);
            // Same convention for the signs of the reflectors, so even `Q` agrees.
            assert_close(&actual, &dense, 1e-11, "op(Q) C");
        }
    }
}

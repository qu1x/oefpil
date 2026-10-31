//! Backend on [`faer`](::faer), which is the fastest if the covariance matrix is block-diagonal or
//! full and large but takes a long time to compile.

use super::{Diag, Side, Trans, Uplo};
use ::faer::{
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

/// Workspace requirement of [`dgemqr()`].
#[derive(Debug, Clone, Copy)]
pub struct WorkspaceReq(StackReq);

impl WorkspaceReq {
    /// Returns the requirement sufficient for this and `other`.
    pub const fn or(self, other: Self) -> Self {
        Self(self.0.or(other.0))
    }
}

/// Workspace of [`dgemqr()`].
pub struct Workspace(MemBuffer);

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
pub fn dgemqr_work(side: Side, trans: Trans, m: usize, n: usize) -> WorkspaceReq {
    WorkspaceReq(match (side, trans) {
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
    })
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
    work: &mut Workspace,
) {
    let work = MemStack::new(&mut work.0);
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
pub fn workspace(req: WorkspaceReq) -> Workspace {
    Workspace(MemBuffer::new(req.0))
}

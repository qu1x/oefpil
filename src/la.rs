//! BLAS and LAPACK routines required by the algorithm.
//!
//! Each function mirrors the routine of the same name for comparison with the reference
//! implementation: It takes the same arguments in the same order and expects column-major matrices
//! with leading dimensions and vectors with increments. It differs as follows:
//!
//!   * Character arguments are enums.
//!   * `INFO` is returned if the routine has failure modes.
//!   * The workspace of `dgemqr()` is queried via `dgemqr_work()` instead of `LWORK = -1`.
//!   * Increments must be positive.
//!   * The factor `T` of `dgeqr()` and `dgemqr()` is their private business and hence has its own
//!     layout. Its size is given by the number of reflectors instead of being queried.
//!
//! The implementation is either the plain backend or the one on `faer` if its feature is enabled.
//! Only this module knows about `faer` which must not leak into the public API.

#![allow(
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::suboptimal_flops,
    clippy::too_many_arguments
)]

#[cfg(feature = "faer")]
mod faer;
#[cfg(not(feature = "faer"))]
mod plain;

#[cfg(feature = "faer")]
pub use self::faer::*;
#[cfg(not(feature = "faer"))]
pub use self::plain::*;

use alloc::{vec, vec::Vec};

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

    const fn side_c(side: Side) -> *const c_char {
        match side {
            Side::Left => c"L".as_ptr(),
            Side::Right => c"R".as_ptr(),
        }
    }

    const fn uplo_c(uplo: Uplo) -> *const c_char {
        match uplo {
            Uplo::Lower => c"L".as_ptr(),
            Uplo::Upper => c"U".as_ptr(),
        }
    }

    const fn trans_c(trans: Trans) -> *const c_char {
        match trans {
            Trans::No => c"N".as_ptr(),
            Trans::Yes => c"T".as_ptr(),
        }
    }

    const fn diag_c(diag: Diag) -> *const c_char {
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
                        // Reported by the reference but not by every implementation (e.g., OpenBLAS).
                        if bad.is_nan() {
                            assert_eq!(pivot + 1, info, "n {n} {uplo:?} pivot {pivot} bad {bad}");
                            continue;
                        }
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
                    &mut work_t,
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
                    &mut work,
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

    /// Sums of squares overflow or underflow at these scales unless the norm guards against it and
    /// LAPACK even rescales at the smallest one.
    #[test]
    fn qr_scales_with_matrix_at_extreme_scales() {
        let mut random = Random(11);
        let (m, n) = (9, 4);
        let a = random.matrix(m * n);
        let mut factor = a.clone();
        let mut t = dgeqr_t(m, n);
        dgeqr(m, n, &mut factor, m, &mut t);
        for scale in [1e200, 1e-200, 1e-300] {
            let mut scaled = a.iter().map(|a| a * scale).collect::<Vec<f64>>();
            let mut fortran = scaled.clone();
            dgeqr_fortran(m, n, &mut fortran, m);
            let mut t = dgeqr_t(m, n);
            dgeqr(m, n, &mut scaled, m, &mut t);
            for j in 0..n {
                for i in 0..=j {
                    let (actual, expected) = (scaled[j * m + i], factor[j * m + i]);
                    assert!(
                        (actual / scale - expected).abs() <= 1e-12 * expected.abs().max(1.0),
                        "R ({i}, {j}) at {scale}: {} != {expected}",
                        actual / scale,
                    );
                    let expected = fortran[j * m + i];
                    assert!(
                        (actual.abs() - expected.abs()).abs() <= 1e-12 * scale,
                        "R ({i}, {j}) at {scale}: {actual} != {expected}",
                    );
                }
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
                    &mut work,
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
                &mut work,
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

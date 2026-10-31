//! The algorithm with variable names as in the reference implementation.

#![allow(
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    clippy::similar_names,
    clippy::suboptimal_flops,
    clippy::too_many_arguments
)]

use crate::{
    Criterion, Format, Model, Notation, Verbosity,
    la::{
        Diag, Side, Trans, Uplo, ddot, dgemqr, dgemqr_work, dgeqr, dgeqr_t, dpotrf, dpotri, dtrsm,
        dtrsv, workspace,
    },
    log::{Log, c_float},
    tiled::{
        Mode, Tile, dgetmv_blockdiag, dgetmv_diag, dgetmv_diags, dgetmv_full, dgetmv_grouped,
        dpottrf_blockdiag, dpottrf_diag, dpottrf_diags, dpottrf_full, dpottrf_grouped,
        dtrtsv_blockdiag, dtrtsv_diag, dtrtsv_diags, dtrtsv_full, dtrtsv_grouped, spans,
    },
};
use alloc::{vec, vec::Vec};

/// Why the iterations ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// Iterating.
    Ok,
    /// Satisfied convergence criterion.
    Converged,
    /// Reached iteration limit.
    MaxIterReached,
    /// Model evaluated to NaN.
    ErrorNan,
    /// Model evaluated to infinity.
    ErrorInf,
    /// Failed factorization.
    ErrorLapack,
}

/// What the iterations yielded.
#[derive(Debug, Clone, Copy)]
pub struct Outcome {
    pub status: Status,
    pub iterations: usize,
    pub chi2: f64,
    pub sigma2: f64,
}

fn conv_param_rel(p0: &[f64], pdiff: &[f64], tol: f64) -> bool {
    p0.iter()
        .zip(pdiff)
        .all(|(p0, pdiff)| (pdiff / p0).abs() < tol)
}

fn conv_param_abs(pdiff: &[f64], tol: f64) -> bool {
    pdiff.iter().all(|pdiff| pdiff.abs() < tol)
}

fn conv_chi2(chi20: f64, chi2: f64, tol: f64) -> bool {
    (chi2 - chi20).abs() < tol
}

fn convergence_check(
    p0: &[f64],
    pdiff: &[f64],
    x0: &[f64],
    xdiff: &[f64],
    chi2: f64,
    chi20: f64,
    tol: f64,
    conv_crit: Criterion,
) -> bool {
    match conv_crit {
        Criterion::RelP => conv_param_rel(p0, pdiff, tol),
        Criterion::RelPOrAbsP => conv_param_rel(p0, pdiff, tol) || conv_param_abs(pdiff, tol),
        Criterion::RelX => conv_param_rel(x0, xdiff, tol),
        Criterion::RelXOrAbsX => conv_param_rel(x0, xdiff, tol) || conv_param_abs(xdiff, tol),
        Criterion::RelPAndRelX => conv_param_rel(p0, pdiff, tol) && conv_param_rel(x0, xdiff, tol),
        Criterion::RelPOrAbsPAndRelXOrAbsX => {
            (conv_param_rel(p0, pdiff, tol) || conv_param_abs(pdiff, tol))
                && (conv_param_rel(x0, xdiff, tol) || conv_param_abs(xdiff, tol))
        }
        Criterion::ChiSquared => conv_chi2(chi20, chi2, tol),
    }
}

/// Evaluates `model` at `n` samples `x` of `xi` and parameters `p`.
fn fcn(
    model: Model,
    n: usize,
    x: &[f64],
    p: &[f64],
    fx: &mut [f64],
    dfdx: &mut [f64],
    dfdp: &mut [f64],
    xi: &mut [f64],
) {
    let nv = model.dfdx.len();
    let np = model.dfdp.len();
    for i in 0..n {
        for (ii, xi) in xi.iter_mut().enumerate() {
            *xi = x[ii * n + i];
        }
        fx[i] = (model.fx)(xi, p);
        for ii in 0..nv {
            dfdx[ii * n + i] = (model.dfdx[ii])(xi, p);
        }
        for ii in 0..np {
            dfdp[ii * n + i] = (model.dfdp[ii])(xi, p);
        }
    }
}

/// Returns the status of `a` not being finite and logs the first such field.
fn check_inf_nan(a: &[f64], location: &str, log: &mut Log) -> Option<Status> {
    let (i, field) = a.iter().enumerate().find(|(_i, a)| !a.is_finite())?;
    let (name, status) = if field.is_nan() {
        ("NaN", Status::ErrorNan)
    } else {
        ("Inf", Status::ErrorInf)
    };
    if log.at(Verbosity::ParameterMean) {
        log.print(format_args!(
            "OEFPIL: {name} encountered in {location} evaluation ({i}-th value)!\n"
        ));
    }
    Some(status)
}

/// Fits parameter `p` with covariance `pcov` to samples `x` where `xtrue` serves as initial estimate
/// and yields the optimum estimate.
///
/// Both `x` and `xtrue` comprise the dependent variable if `model` is explicit.
#[allow(clippy::too_many_lines)]
pub fn oefpil(
    model: Model,
    p: &mut [f64],
    pcov: &mut [f64],
    n: usize,
    x: &[f64],
    xtrue: &mut [f64],
    k: &[f64],
    tiling: Mode,
    tilemap: &[Tile],
    groups: &[usize],
    maxit: usize,
    tol: f64,
    conv_crit: Criterion,
    usesigma2: bool,
    log: &mut Log,
) -> Outcome {
    let nv = model.dfdx.len();
    let nve = model.variables();
    let np = p.len();

    let mut lpkerrsrc = "";
    let mut status = Status::Ok;

    if log.at(Verbosity::ParameterMean) {
        log.print(format_args!("### OEFPIL begin ###\n\n*** Input ***\n\n"));
        log.print(format_args!(
            "* function mode: {}\n",
            if model.implicit {
                "implicit"
            } else {
                "explicit"
            }
        ));
        log.print(format_args!(
            "* covariance matrix tiling mode: {}\n",
            tiling as i32
        ));
        log.print(format_args!("* number of data points: {n}\n"));
        log.print(format_args!("* number of variables: {nv}\n"));
        log.print(format_args!("* number of parameters: {np}\n"));
        log.print(format_args!("* max number of iterations: {maxit}\n"));
        log.print(format_args!(
            "* relative tolerance for convergence: {}\n",
            c_float(tol, Notation::Exponential, 6, false, 0)
        ));
        log.print(format_args!(
            "* number of decimal digits in printed output: {}\n\n",
            log.format().precision
        ));
    }
    if log.at(Verbosity::IndividualSteps) {
        log.matrix(n, nv, x, Some("* x values *"));
        if !model.implicit {
            log.matrix(n, 1, &x[n * nv..], Some("* y values *"));
        }
        log.matrix(n, nv, xtrue, Some("* xtrue values *"));
        if !model.implicit {
            log.matrix(n, 1, &xtrue[n * nv..], Some("* ytrue values *"));
        }
        log.print(format_args!("*** Initialization ***\n\n"));
    }

    let mut wf = vec![0f64; n];
    let mut b = vec![0f64; n];
    let mut b1_diags = vec![0f64; n * nve];
    let mut xdiff = vec![0f64; n * nve];
    let mut deltax_hat = vec![0f64; n * nve];
    let mut b2 = vec![0f64; n * np];
    let mut u = vec![0f64; n];
    let mut v = vec![0f64; n];
    let mut w = vec![0f64; np];
    let mut z = vec![0f64; nve * n];
    let mut xi = vec![0f64; nv];

    // M11 is diagonal for some modes and hence has an array of its own, otherwise the factor LM
    // replaces it. It has one block per group of samples, which is a single one unless grouped.
    let diagonal = matches!(tiling, Mode::Diagonal | Mode::Diagonals);
    let whole = [n];
    let blocks = if tiling == Mode::Grouped {
        groups
    } else {
        &whole
    };
    let m11_fields = blocks.iter().map(|&block| block * block).sum();
    let mut m11 = vec![0f64; if diagonal { n } else { m11_fields }];
    let mut lm = vec![0f64; if diagonal { n } else { 0 }];
    let largest = groups.iter().copied().max().unwrap_or(0);
    let mut lk_work = vec![0f64; nve * largest];

    let mut t = dgeqr_t(n, np);
    let mut work = workspace(dgemqr_work(Side::Left, Trans::Yes, n, 1).or(dgemqr_work(
        Side::Left,
        Trans::No,
        n,
        1,
    )));

    // For chi-squared.
    let mut lk = k.to_vec();
    let mut lpkinfo = match tiling {
        Mode::Diagonal => dpottrf_diag(n, nve, &mut lk),
        Mode::BlockDiagonal => dpottrf_blockdiag(n, nve, &mut lk),
        Mode::Diagonals => dpottrf_diags(n, nve, &mut lk),
        Mode::Full => dpottrf_full(n, nve, &mut lk),
        Mode::Grouped => dpottrf_grouped(groups, nve, &mut lk),
    };
    if lpkinfo != 0 {
        status = Status::ErrorLapack;
        lpkerrsrc = "DPOTTRF_xxx";
    }

    let mut wx: Vec<f64> = xtrue[..n * nve].to_vec();
    let mut chi2 = f64::NAN;
    let mut chi2prev = 0f64;
    let mut iter_num = 0;

    if log.at(Verbosity::ParameterMean) {
        log.matrix(1, np, p, Some("* initial parameter estimates *\n"));
        log.print(format_args!("\n*** Iterations ***\n\n"));
    }

    while status == Status::Ok {
        if log.at(Verbosity::IndividualSteps) {
            log.print(format_args!("** iteration {iter_num} begin **\n\n"));
        }

        fcn(model, n, &wx, p, &mut wf, &mut b1_diags, &mut b2, &mut xi);
        if let Some(error) = check_inf_nan(&wf, "function values", log)
            .or_else(|| check_inf_nan(&b1_diags[..n * nv], "function derivatives", log))
            .or_else(|| check_inf_nan(&b2, "function parameter derivatives", log))
        {
            status = error;
            break;
        }

        if model.implicit {
            // b(n) = f(x_i, y_i)
            b.copy_from_slice(&wf);
        } else {
            // b(n) = f(x) - y
            // B1(n, (nv + 1) * n) = (B11 | ... | B1nv | -I_n)
            for i in 0..n {
                b[i] = wf[i] - wx[n * nv + i];
                b1_diags[n * nv + i] = -1.0;
            }
        }

        // M11(n, n) = B1(n, nve * n) * U(nve * n, nve * n) * B1^T(nve * n, n)
        // M11 = B11 * Kxx * B11 + B11 * Kxy * B12 + B12 * Kyx * B11 + B12 * Kyy * B12
        // Cholesky M11 = LM * LM^T
        // E = LM^(-1) * B2 (solve LM * E = B2), in place
        m11.fill(0.0);
        match tiling {
            Mode::Diagonal => {
                for ii in 0..nve {
                    if tilemap[ii] == Tile::Diagonal {
                        let b1_ii = &b1_diags[ii * n..];
                        let k_ = &k[ii * n..];
                        for i in 0..n {
                            m11[i] += b1_ii[i] * k_[i] * b1_ii[i];
                        }
                    }
                }
            }
            Mode::Diagonals => {
                for jj in 0..nve {
                    for ii in 0..nve {
                        if tilemap[jj * nve + ii] == Tile::Diagonal {
                            let b1_ii = &b1_diags[ii * n..];
                            let b1_jj = &b1_diags[jj * n..];
                            let k_ = &k[(jj * nve + ii) * n..];
                            for i in 0..n {
                                m11[i] += b1_ii[i] * k_[i] * b1_jj[i];
                            }
                        }
                    }
                }
            }
            Mode::BlockDiagonal => {
                for ii in 0..nve {
                    let b1_ii = &b1_diags[ii * n..];
                    let k_ = &k[ii * (n * n)..];
                    for j in 0..n {
                        for i in 0..n {
                            m11[j * n + i] += b1_ii[i] * k_[j * n + i] * b1_ii[j];
                        }
                    }
                }
            }
            Mode::Full => {
                for jj in 0..nve {
                    for ii in 0..nve {
                        if tilemap[jj * nve + ii] != Tile::None {
                            let b1_ii = &b1_diags[ii * n..];
                            let b1_jj = &b1_diags[jj * n..];
                            let k_ = &k[jj * nve * n * n + ii * n..];
                            for j in 0..n {
                                for i in 0..n {
                                    m11[j * n + i] += b1_ii[i] * k_[j * n * nve + i] * b1_jj[j];
                                }
                            }
                        }
                    }
                }
            }
            Mode::Grouped => {
                for ((m, start, offset), (_, _, m11_offset)) in
                    spans(groups, nve).zip(spans(groups, 1))
                {
                    let m11 = &mut m11[m11_offset..];
                    for jj in 0..nve {
                        for ii in 0..nve {
                            if tilemap[jj * nve + ii] != Tile::None {
                                let b1_ii = &b1_diags[ii * n + start..];
                                let b1_jj = &b1_diags[jj * n + start..];
                                let k_ = &k[offset + jj * nve * m * m + ii * m..];
                                for j in 0..m {
                                    for i in 0..m {
                                        m11[j * m + i] += b1_ii[i] * k_[j * m * nve + i] * b1_jj[j];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if diagonal {
            // Diagonal Cholesky.
            for i in 0..n {
                lm[i] = libm::sqrt(m11[i]);
            }
            // E = B2
            for j in 0..np {
                for i in 0..n {
                    b2[j * n + i] /= lm[i];
                }
            }
        } else {
            // Cholesky in place, M11 now contains just the factor LM.
            lpkinfo = 0;
            for (block, _, offset) in spans(blocks, 1) {
                lpkinfo = dpotrf(Uplo::Lower, block, &mut m11[offset..], block);
                if lpkinfo != 0 {
                    break;
                }
            }
            if lpkinfo != 0 {
                status = Status::ErrorLapack;
                lpkerrsrc = "DPOTRF(M11)";
                break;
            }
            // E, B2 now contains the solution.
            for (block, start, offset) in spans(blocks, 1) {
                dtrsm(
                    Side::Left,
                    Uplo::Lower,
                    Trans::No,
                    Diag::No,
                    block,
                    np,
                    1.0,
                    &m11[offset..],
                    block,
                    &mut b2[start..],
                    n,
                );
            }
        }

        // xdiff(nve * n)
        for i in 0..n * nve {
            xdiff[i] = x[i] - wx[i];
        }

        // u(n) = B1(n, nve * n) * xdiff(nve * n) + b(n)
        u.copy_from_slice(&b);
        for ii in 0..nve {
            let b1_ii = &b1_diags[ii * n..];
            let xdiff_ii = &xdiff[ii * n..];
            for i in 0..n {
                u[i] += b1_ii[i] * xdiff_ii[i];
            }
        }

        // E(n, np) -> QE(n, n) * RE(n, np)
        dgeqr(n, np, &mut b2, n, &mut t);

        // u -> solve(LM, u)
        if diagonal {
            for i in 0..n {
                u[i] /= lm[i];
            }
        } else {
            for (block, start, offset) in spans(blocks, 1) {
                dtrsv(
                    Uplo::Lower,
                    Trans::No,
                    Diag::No,
                    block,
                    &m11[offset..],
                    block,
                    &mut u[start..],
                    1,
                );
            }
        }

        // v(n) = u(n)
        v.copy_from_slice(&u);

        // u(n) -> QE^T(n, n) * u(n)
        dgemqr(
            Side::Left,
            Trans::Yes,
            n,
            1,
            np,
            &b2,
            n,
            &t,
            &mut u,
            n,
            &mut work,
        );

        // QE is square = (Q_R Q_0)
        // (Q_R Q_0) (Q_R^T; Q_0^T) u = (Q_R Q_0) (Q_R^T u; Q_0^T u) = Q_R Q_R^T u + Q_0 Q_0^T u
        // Zero out Q_0^T u, so that only Q_R Q_R^T u is calculated.
        u[np..].fill(0.0);

        // w(np) = -u(1:np), for deltap_hat
        for i in 0..np {
            w[i] = -u[i];
        }

        // u(n) -> QE(n, n) * u(n)
        dgemqr(
            Side::Left,
            Trans::No,
            n,
            1,
            np,
            &b2,
            n,
            &t,
            &mut u,
            n,
            &mut work,
        );

        // v(n) -> v(n) - u(n)
        for i in 0..n {
            v[i] -= u[i];
        }

        // v(n) -> solve(LM^T(n, n), v(n))
        if diagonal {
            for i in 0..n {
                v[i] /= lm[i];
            }
        } else {
            for (block, start, offset) in spans(blocks, 1) {
                dtrsv(
                    Uplo::Lower,
                    Trans::Yes,
                    Diag::No,
                    block,
                    &m11[offset..],
                    block,
                    &mut v[start..],
                    1,
                );
            }
        }

        // deltap_hat(np) = w(np)
        // w(np) -> solve(RE(np, np), w(np))
        dtrsv(Uplo::Upper, Trans::No, Diag::No, np, &b2, n, &mut w, 1);

        // pcov, like the reference this does not check for an exactly singular factor which yields
        // non-finite estimates failing the next evaluation of the model.
        dpotri(Uplo::Upper, np, &mut b2, n);
        for j in 0..np {
            for i in 0..=j {
                pcov[j * np + i] = b2[j * n + i];
                pcov[i * np + j] = b2[j * n + i];
            }
        }

        // Update p.
        for i in 0..np {
            p[i] += w[i];
        }

        // z(nve * n) = B1^T(nve * n, n) * v(n)
        for ii in 0..nve {
            let b1_ii = &b1_diags[ii * n..];
            let z_ = &mut z[ii * n..];
            for i in 0..n {
                z_[i] = b1_ii[i] * v[i];
            }
        }

        // deltax_hat(nve * n) = xdiff(nve * n) - K(nve * n, nve * n) * z(nve * n)
        deltax_hat.copy_from_slice(&xdiff);
        match tiling {
            Mode::Diagonal => dgetmv_diag(n, nve, -1.0, k, tilemap, &z, 1.0, &mut deltax_hat),
            Mode::BlockDiagonal => {
                dgetmv_blockdiag(n, nve, -1.0, k, tilemap, &z, 1.0, &mut deltax_hat);
            }
            Mode::Diagonals => dgetmv_diags(n, nve, -1.0, k, tilemap, &z, 1.0, &mut deltax_hat),
            Mode::Full => dgetmv_full(n, nve, -1.0, k, tilemap, &z, 1.0, &mut deltax_hat),
            Mode::Grouped => {
                dgetmv_grouped(groups, n, nve, -1.0, k, tilemap, &z, 1.0, &mut deltax_hat);
            }
        }

        for i in 0..n * nve {
            wx[i] += deltax_hat[i];
        }

        // chi2, update xdiff with new wx
        for i in 0..n * nve {
            xdiff[i] = x[i] - wx[i];
        }

        // solve(LK, xdiff)
        match tiling {
            Mode::Diagonal => dtrtsv_diag(n, nve, &lk, &mut xdiff),
            Mode::BlockDiagonal => dtrtsv_blockdiag(n, nve, &lk, &mut xdiff),
            Mode::Diagonals => dtrtsv_diags(n, nve, &lk, &mut xdiff),
            Mode::Full => dtrtsv_full(n, nve, &lk, &mut xdiff),
            Mode::Grouped => dtrtsv_grouped(groups, n, nve, &lk, &mut xdiff, &mut lk_work),
        }

        chi2 = ddot(nve * n, &xdiff, 1, &xdiff, 1);

        // deltap_hat = w
        if convergence_check(p, &w, x, &deltax_hat, chi2, chi2prev, tol, conv_crit) {
            status = Status::Converged;
        }

        if log.at(Verbosity::ParameterMean) {
            log.print(format_args!("{iter_num:<4}\t"));
            log.matrix(1, np, p, None);
        }
        if log.at(Verbosity::ParameterCovariance) {
            log.matrix(np, np, pcov, Some("* pcov *"));
            log.matrix_with_format(
                1,
                1,
                core::slice::from_ref(&chi2),
                Some("* chi2 *"),
                Format {
                    notation: Notation::Exponential,
                    precision: 15,
                },
            );
            log.print(format_args!("\n"));
        }

        iter_num += 1;
        chi2prev = chi2;

        if status == Status::Ok && iter_num >= maxit {
            status = Status::MaxIterReached;
        }
    }

    if log.at(Verbosity::ParameterMean) {
        match status {
            Status::Converged => log.print(format_args!(
                "*** Success: parameters converged in {iter_num} iterations ***\n"
            )),
            Status::MaxIterReached => log.print(format_args!(
                "*** Failure: maximum number of iterations ({maxit}) reached ***\n"
            )),
            Status::ErrorNan => log.print(format_args!(
                "*** Failure: NaN encountered in function/derivatives values ***\n"
            )),
            Status::ErrorInf => log.print(format_args!(
                "*** Failure: Inf encountered in function/derivatives values ***\n"
            )),
            Status::ErrorLapack => log.print(format_args!(
                "*** Failure: error encountered in LAPACK procedure {lpkerrsrc}, info code \
                 {lpkinfo} (see LAPACK documentation for details) ***\n"
            )),
            Status::Ok => unreachable!("still iterating"),
        }
        log.print(format_args!("\n### OEFPIL end ###\n"));
    }

    xtrue[..n * nve].copy_from_slice(&wx);

    #[allow(clippy::cast_precision_loss)]
    let sigma2 = chi2 / (n - np) as f64;
    if usesigma2 {
        for field in &mut *pcov {
            *field *= sigma2;
        }
    }

    Outcome {
        status,
        iterations: iter_num,
        chi2,
        sigma2,
    }
}

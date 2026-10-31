//! Whole fits compared against the C/Fortran library via its `oefpil-sys` bindings.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::needless_pass_by_value,
    clippy::similar_names,
    clippy::suboptimal_flops,
    clippy::too_many_lines
)]

use core::ffi::{c_int, c_void};
use oefpil::{
    Algorithm, AlgorithmError, Covariance, Criterion, Model, OefpilError, Parameter, Report,
    Variable, Verbosity,
};
use oefpil_sys as sys;
use std::{fs, process, ptr, slice};

/// Tiling mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    Diagonal,
    BlockDiagonal,
    Diagonals,
    Full,
}

/// Problem being fitted by either implementation.
#[derive(Clone)]
struct Problem {
    model: Model<'static>,
    samples: usize,
    mode: Mode,
    /// Tiles at row and column with fields as passed to [`Covariance::set_tile()`].
    tiles: Vec<(usize, usize, Vec<f64>)>,
    sample: Vec<f64>,
    parameter: Vec<f64>,
    unknown_scale: bool,
    algorithm: Algorithm,
}

/// What a fit yielded, whether successful or not.
#[derive(Debug, Clone)]
struct Fit {
    status: Result<(), AlgorithmError>,
    iterations: usize,
    chi_squared: f64,
    chi_squared_reduced: f64,
    parameter: Vec<f64>,
    parameter_covariance: Vec<f64>,
    mean: Vec<f64>,
}

impl Problem {
    fn covariance(&self) -> Covariance {
        let variables = self.model.variables();
        let mut covariance = match self.mode {
            Mode::Diagonal => Covariance::new_diagonal(self.samples, variables),
            Mode::BlockDiagonal => Covariance::new_block_diagonal(self.samples, variables),
            Mode::Diagonals => Covariance::new_diagonals(self.samples, variables),
            Mode::Full => Covariance::new_full(self.samples, variables),
        };
        for (row, column, fields) in &self.tiles {
            covariance.set_tile(*row, *column, fields).unwrap();
        }
        if self.unknown_scale {
            covariance = covariance.with_unknown_scale();
        }
        covariance
    }

    fn fit_rust(&self, log: &mut dyn core::fmt::Write) -> Fit {
        let covariance = self.covariance();
        let mut fit = Fit {
            status: Ok(()),
            iterations: 0,
            chi_squared: f64::NAN,
            chi_squared_reduced: f64::NAN,
            parameter: self.parameter.clone(),
            parameter_covariance: vec![0f64; self.parameter.len().pow(2)],
            mean: self.sample.clone(),
        };
        let result = self.algorithm.fit_with_log(
            self.model,
            &mut Variable {
                sample: &self.sample,
                mean: &mut fit.mean,
                covariance: &covariance,
            },
            &mut Parameter {
                mean: &mut fit.parameter,
                covariance: &mut fit.parameter_covariance,
                ..Parameter::default()
            },
            log,
        );
        match result {
            Ok(Report {
                chi_squared,
                chi_squared_reduced,
                iterations,
                ..
            }) => {
                fit.chi_squared = chi_squared;
                fit.chi_squared_reduced = chi_squared_reduced;
                fit.iterations = iterations;
            }
            Err(OefpilError::Algorithm(error)) => fit.status = Err(error),
            Err(error) => panic!("{error}"),
        }
        fit
    }

    /// Fits via C library with `log` path if not silent.
    fn fit_c(&self, log: Option<&std::path::Path>) -> Fit {
        let variables = self.model.variables();
        let x_variables = self.model.dfdx.len();
        let (n, nv) = (self.samples, variables);
        let (fields, tilemap, mode) = match self.mode {
            Mode::Diagonal => (n * nv, nv, 1),
            Mode::BlockDiagonal => (n * n * nv, nv, 2),
            Mode::Diagonals => (n * nv * nv, nv * nv, 3),
            Mode::Full => (n * n * nv * nv, nv * nv, 4),
        };
        let mut tcm = vec![0f64; fields];
        let mut map = vec![0 as c_int; tilemap];
        let (cn, cnv) = (n as c_int, nv as c_int);
        for (row, column, b) in &self.tiles {
            let (r, c) = (*row as c_int, *column as c_int);
            let (a, m, b) = (tcm.as_mut_ptr(), map.as_mut_ptr(), b.as_ptr());
            let length = self
                .tiles
                .iter()
                .find(|tile| tile.0 == *row && tile.1 == *column)
                .unwrap()
                .2
                .len();
            unsafe {
                match self.mode {
                    Mode::Diagonal => sys::oefpil_tcm_diag_set_tile_diag(cn, cnv, a, m, r, b),
                    Mode::BlockDiagonal if length == n => {
                        sys::oefpil_tcm_blockdiag_set_tile_diag(cn, cnv, a, m, r, b);
                    }
                    Mode::BlockDiagonal if length == n * (n + 1) / 2 => {
                        sys::oefpil_tcm_blockdiag_set_tile_half(cn, cnv, a, m, r, b);
                    }
                    Mode::BlockDiagonal => {
                        sys::oefpil_tcm_blockdiag_set_tile_full(cn, cnv, a, m, r, b);
                    }
                    Mode::Diagonals => sys::oefpil_tcm_diags_set_tile_diag(cn, cnv, a, m, r, c, b),
                    Mode::Full if length == n => {
                        sys::oefpil_tcm_full_set_tile_diag(cn, cnv, a, m, c, r, b);
                    }
                    Mode::Full if length == n * (n + 1) / 2 => {
                        sys::oefpil_tcm_full_set_tile_half(cn, cnv, a, m, c, r, b);
                    }
                    Mode::Full => sys::oefpil_tcm_full_set_tile_full(cn, cnv, a, m, c, r, b),
                }
            }
        }
        let x_fields = n * x_variables;
        let mut parameter = self.parameter.clone();
        let mut parameter_covariance = vec![0f64; parameter.len().pow(2)];
        let mut mean = self.sample.clone();
        let (x_mean, y_mean) = mean.split_at_mut(x_fields);
        let mut data = Data {
            model: self.model,
            x: vec![0f64; x_variables],
        };
        let (mut chi_squared, mut chi_squared_reduced) = (f64::NAN, f64::NAN);
        let (mut info, mut iterations) = (0, 0);
        let file = log.map_or(ptr::null_mut(), |path| {
            let path = std::ffi::CString::new(path.as_os_str().as_encoded_bytes()).unwrap();
            let file = unsafe { sys::libc::fopen(path.as_ptr(), c"w".as_ptr()) };
            assert!(!file.is_null());
            file
        });
        unsafe {
            sys::oefpil(
                Some(evaluate),
                ptr::from_mut(&mut data).cast(),
                self.model.implicit.into(),
                parameter.len() as c_int,
                parameter.as_mut_ptr(),
                parameter_covariance.as_mut_ptr(),
                cn,
                x_variables as c_int,
                self.sample.as_ptr(),
                self.sample[x_fields..].as_ptr(),
                x_mean.as_mut_ptr(),
                y_mean.as_mut_ptr(),
                tcm.as_ptr(),
                mode,
                map.as_ptr(),
                self.algorithm.iteration_limit as c_int,
                self.algorithm.tolerance,
                self.algorithm.verbosity as c_int,
                file,
                &raw mut chi_squared,
                self.algorithm.criterion as c_int,
                &raw mut info,
                &raw mut iterations,
                &raw mut chi_squared_reduced,
                self.unknown_scale,
                ptr::null_mut(),
            );
            if !file.is_null() {
                sys::libc::fclose(file);
            }
        }
        Fit {
            status: match info {
                1 => Ok(()),
                2 => Err(AlgorithmError::IterationLimit),
                3 => Err(AlgorithmError::NumericalError),
                _ => unreachable!(),
            },
            iterations: iterations as usize,
            chi_squared,
            chi_squared_reduced,
            parameter,
            parameter_covariance,
            mean,
        }
    }

    /// Asserts both implementations agree within `tolerance`, except for hard problems where they
    /// may need a different number of iterations.
    #[track_caller]
    fn assert_agreement(&self, tolerance: f64, what: &str) -> Fit {
        let rust = self.fit_rust(&mut String::new());
        let c = self.fit_c(None);
        assert_eq!(c.status, rust.status, "{what}: status");
        if c.status.is_ok() {
            assert_eq!(c.iterations, rust.iterations, "{what}: iterations");
            assert_close(c.chi_squared, rust.chi_squared, tolerance, what);
            // Division by zero degrees of freedom depends on whether rounding yields zero.
            if self.samples > self.parameter.len() {
                assert_close(
                    c.chi_squared_reduced,
                    rust.chi_squared_reduced,
                    tolerance,
                    what,
                );
            }
            assert_all_close(&rust.parameter, &c.parameter, tolerance, what, "parameter");
            assert_all_close(
                &rust.parameter_covariance,
                &c.parameter_covariance,
                tolerance,
                what,
                "parameter covariance",
            );
            assert_all_close(&rust.mean, &c.mean, tolerance, what, "mean");
        }
        rust
    }
}

#[track_caller]
fn assert_close(expected: f64, actual: f64, tolerance: f64, what: &str) {
    assert!(
        actual == expected
            || (actual.is_nan() && expected.is_nan())
            || (actual - expected).abs() <= tolerance * expected.abs().max(1e-6),
        "{what}: {actual} != {expected}",
    );
}

#[track_caller]
fn assert_all_close(actual: &[f64], expected: &[f64], tolerance: f64, what: &str, field: &str) {
    let scale = expected.iter().fold(0f64, |scale, e| scale.max(e.abs()));
    for (index, (a, e)) in actual.iter().zip(expected).enumerate() {
        assert!(
            (a.is_nan() && e.is_nan()) || (a - e).abs() <= tolerance * scale,
            "{what}: {field} {index}: {a} != {e}",
        );
    }
}

struct Data {
    model: Model<'static>,
    x: Vec<f64>,
}

unsafe extern "C" fn evaluate(
    data: *mut c_void,
    samples: c_int,
    x_sample: *const f64,
    parameter: *const f64,
    fx: *mut f64,
    dfdx: *mut f64,
    dfdp: *mut f64,
) {
    let Data { model, x } = unsafe { &mut *data.cast::<Data>() };
    let samples = samples as usize;
    let x_variables = model.dfdx.len();
    let parameters = model.dfdp.len();
    let x_sample = unsafe { slice::from_raw_parts(x_sample, samples * x_variables) };
    let p = unsafe { slice::from_raw_parts(parameter, parameters) };
    let fx = unsafe { slice::from_raw_parts_mut(fx, samples) };
    let dfdx = unsafe { slice::from_raw_parts_mut(dfdx, samples * x_variables) };
    let dfdp = unsafe { slice::from_raw_parts_mut(dfdp, samples * parameters) };
    for sample in 0..samples {
        for x_variable in 0..x_variables {
            x[x_variable] = x_sample[x_variable * samples + sample];
        }
        fx[sample] = (model.fx)(x, p);
        for x_variable in 0..x_variables {
            dfdx[x_variable * samples + sample] = (model.dfdx[x_variable])(x, p);
        }
        for parameter in 0..parameters {
            dfdp[parameter * samples + sample] = (model.dfdp[parameter])(x, p);
        }
    }
}

const LINE: Model = Model {
    fx: &|x, p| p[0] * x[0] + p[1],
    dfdx: &[&|_, p| p[0]],
    dfdp: &[&|x, _| x[0], &|_, _| 1.0],
    implicit: false,
};

const CUBIC: Model = Model {
    fx: &|x, p| p[0] + x[0] * (p[1] + x[0] * (p[2] + x[0] * p[3])),
    dfdx: &[&|x, p| p[1] + x[0] * (2.0 * p[2] + x[0] * 3.0 * p[3])],
    dfdp: &[&|_, _| 1.0, &|x, _| x[0], &|x, _| x[0] * x[0], &|x, _| {
        x[0].powi(3)
    }],
    implicit: false,
};

const QUINTIC: Model = Model {
    fx: &|x, p| p[0] + x[0] * (p[1] + x[0] * (p[2] + x[0] * (p[3] + x[0] * (p[4] + x[0] * p[5])))),
    dfdx: &[&|x, p| {
        p[1] + x[0] * (2.0 * p[2] + x[0] * (3.0 * p[3] + x[0] * (4.0 * p[4] + x[0] * 5.0 * p[5])))
    }],
    dfdp: &[
        &|_, _| 1.0,
        &|x, _| x[0],
        &|x, _| x[0].powi(2),
        &|x, _| x[0].powi(3),
        &|x, _| x[0].powi(4),
        &|x, _| x[0].powi(5),
    ],
    implicit: false,
};

const BOXBOD: Model = Model {
    fx: &|x, p| p[0] * (1.0 - (-p[1] * x[0]).exp()),
    dfdx: &[&|x, p| p[0] * (-p[1] * x[0]).exp() * p[1]],
    dfdp: &[&|x, p| 1.0 - (-p[1] * x[0]).exp(), &|x, p| {
        p[0] * (-p[1] * x[0]).exp() * x[0]
    }],
    implicit: false,
};

/// Explicit bivariate model, so three variables.
const PLANE: Model = Model {
    fx: &|x, p| p[0] + p[1] * x[0] + p[2] * x[1],
    dfdx: &[&|_, p| p[1], &|_, p| p[2]],
    dfdp: &[&|_, _| 1.0, &|x, _| x[0], &|x, _| x[1]],
    implicit: false,
};

/// Implicit nonlinear model of three variables.
const SURFACE: Model = Model {
    fx: &|x, p| x[0] + p[0] * x[1] + p[1] * x[2] + p[2] * x[0] * x[2] - 1.0,
    dfdx: &[&|x, p| 1.0 + p[2] * x[2], &|_, p| p[0], &|x, p| {
        p[1] + p[2] * x[0]
    }],
    dfdp: &[&|x, _| x[1], &|x, _| x[2], &|x, _| x[0] * x[2]],
    implicit: true,
};

/// Implicit circle of two variables.
const CIRCLE: Model = Model {
    fx: &|x, p| (x[0] - p[0]).powi(2) + (x[1] - p[1]).powi(2) - p[2] * p[2],
    dfdx: &[&|x, p| 2.0 * (x[0] - p[0]), &|x, p| 2.0 * (x[1] - p[1])],
    dfdp: &[
        &|x, p| -2.0 * (x[0] - p[0]),
        &|x, p| -2.0 * (x[1] - p[1]),
        &|_, p| -2.0 * p[2],
    ],
    implicit: true,
};

/// Well-conditioned noisy sample of `model` at `samples` points with `exact` being applied on them.
fn sample(
    model: Model,
    samples: usize,
    p: &[f64],
    exact: impl Fn(&mut [f64], usize, &[f64]),
) -> Vec<f64> {
    let nv = model.dfdx.len();
    let mut v = vec![0f64; samples * model.variables()];
    let mut x = vec![0f64; nv];
    for sample in 0..samples {
        let s = sample as f64;
        for (variable, x) in x.iter_mut().enumerate() {
            *x =
                1.0 + 0.3 * (s * (1.1 + variable as f64)).sin() + 0.1 * s / (1.0 + variable as f64);
        }
        exact(&mut x, sample, p);
        for (variable, x) in x.iter().enumerate() {
            v[variable * samples + sample] =
                x + 0.02 * ((variable * 7 + sample * 3) as f64 * 1.7).sin();
        }
        if !model.implicit {
            v[nv * samples + sample] =
                (model.fx)(&x, p) + 0.02 * ((nv * 7 + sample * 3) as f64 * 1.7).sin();
        }
    }
    v
}

fn among(r: usize, c: usize) -> f64 {
    if r == c {
        1.0
    } else {
        0.4 / (1.0 + r.abs_diff(c) as f64)
    }
}

/// Full tile at `r` and `c` of `n` samples with correlation `within` variables.
fn block(n: usize, r: usize, c: usize, within: f64) -> Vec<f64> {
    let mut tile = vec![0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let correlation = within.powi(i.abs_diff(j) as i32);
            tile[i * n + j] = among(r, c) * correlation * 1e-3 * (1.0 + 0.1 * (r + c) as f64);
        }
    }
    tile
}

/// Lower triangle of [`block()`] as row-major triangular slice.
fn half(n: usize, r: usize, c: usize, within: f64) -> Vec<f64> {
    let tile = block(n, r, c, within);
    let mut half = Vec::new();
    for i in 0..n {
        for j in 0..=i {
            half.push(tile[i * n + j]);
        }
    }
    half
}

fn diagonal(n: usize, r: usize, c: usize) -> Vec<f64> {
    vec![among(r, c) * 1e-3 * (1.0 + 0.1 * (r + c) as f64); n]
}

/// Tiles of all upper positions by `tile`.
fn upper(
    variables: usize,
    tile: impl Fn(usize, usize) -> Vec<f64>,
) -> Vec<(usize, usize, Vec<f64>)> {
    let mut tiles = Vec::new();
    for r in 0..variables {
        for c in r..variables {
            tiles.push((r, c, tile(r, c)));
        }
    }
    tiles
}

/// Problems of `model` of each tiling mode inclusive any kind of tile.
fn problems(
    model: Model<'static>,
    samples: usize,
    sampled: Vec<f64>,
    p: Vec<f64>,
) -> Vec<(&'static str, Problem)> {
    let n = samples;
    let nv = model.variables();
    let problem = |mode, tiles| Problem {
        model,
        samples,
        mode,
        tiles,
        sample: sampled.clone(),
        parameter: p.clone(),
        unknown_scale: false,
        algorithm: Algorithm::default(),
    };
    let diagonals_on_diagonal = (0..nv)
        .map(|r| (r, r, diagonal(n, r, r)))
        .collect::<Vec<_>>();
    let blocks_on_diagonal = (0..nv)
        .map(|r| (r, r, block(n, r, r, 0.3)))
        .collect::<Vec<_>>();
    let halves_on_diagonal = (0..nv)
        .map(|r| (r, r, half(n, r, r, 0.3)))
        .collect::<Vec<_>>();
    let mixed_on_diagonal = (0..nv)
        .map(|r| {
            let tile = match r % 3 {
                0 => diagonal(n, r, r),
                1 => block(n, r, r, 0.3),
                _ => half(n, r, r, 0.3),
            };
            (r, r, tile)
        })
        .collect::<Vec<_>>();
    vec![
        (
            "diagonal",
            problem(Mode::Diagonal, diagonals_on_diagonal.clone()),
        ),
        (
            "block-diagonal of diagonal tiles",
            problem(Mode::BlockDiagonal, diagonals_on_diagonal.clone()),
        ),
        (
            "block-diagonal of block tiles",
            problem(Mode::BlockDiagonal, blocks_on_diagonal.clone()),
        ),
        (
            "block-diagonal of triangular slices",
            problem(Mode::BlockDiagonal, halves_on_diagonal),
        ),
        (
            "block-diagonal of mixed tiles",
            problem(Mode::BlockDiagonal, mixed_on_diagonal),
        ),
        (
            "diagonals on diagonal",
            problem(Mode::Diagonals, diagonals_on_diagonal.clone()),
        ),
        (
            "diagonals all over",
            problem(Mode::Diagonals, upper(nv, |r, c| diagonal(n, r, c))),
        ),
        (
            "full of diagonal tiles on diagonal",
            problem(Mode::Full, diagonals_on_diagonal),
        ),
        (
            "full of block tiles on diagonal",
            problem(Mode::Full, blocks_on_diagonal),
        ),
        (
            "full of diagonal tiles",
            problem(Mode::Full, upper(nv, |r, c| diagonal(n, r, c))),
        ),
        (
            "full of block tiles",
            problem(Mode::Full, upper(nv, |r, c| block(n, r, c, 0.3))),
        ),
        (
            "full of triangular slices",
            problem(Mode::Full, upper(nv, |r, c| half(n, r, c, 0.3))),
        ),
        (
            "full of mixed tiles",
            problem(
                Mode::Full,
                upper(nv, |r, c| match (r + c) % 3 {
                    0 => block(n, r, c, 0.3),
                    1 => diagonal(n, r, c),
                    _ => half(n, r, c, 0.2),
                }),
            ),
        ),
    ]
}

fn suite() -> Vec<(&'static str, Problem)> {
    let line = |n: usize| {
        let p = vec![1.0, 2.0];
        (sample(LINE, n, &p, |_, _, _| {}), vec![1.1, 1.9])
    };
    let plane = |n: usize| {
        let p = vec![1.0, 0.5, -0.5];
        (sample(PLANE, n, &p, |_, _, _| {}), vec![0.9, 0.6, -0.4])
    };
    let surface = |n: usize| {
        let p = vec![0.8, 0.5, -0.2];
        let exact = |x: &mut [f64], _: usize, p: &[f64]| {
            x[1] = (1.0 - x[0] - p[1] * x[2] - p[2] * x[0] * x[2]) / p[0];
        };
        (sample(SURFACE, n, &p, exact), vec![0.85, 0.45, -0.15])
    };
    let circle = |n: usize| {
        let p = vec![1.0, 2.0, 3.0];
        let exact = |x: &mut [f64], sample: usize, p: &[f64]| {
            let angle = sample as f64 * 0.7;
            x[0] = p[0] + p[2] * angle.cos();
            x[1] = p[1] + p[2] * angle.sin();
        };
        (sample(CIRCLE, n, &p, exact), vec![1.1, 1.9, 3.2])
    };
    let mut suite = Vec::new();
    for n in [2, 5, 13, 40] {
        let (v, p) = line(n);
        suite.extend(problems(LINE, n, v, p));
    }
    for n in [3, 7, 21] {
        let (v, p) = plane(n);
        suite.extend(problems(PLANE, n, v, p));
        let (v, p) = surface(n);
        suite.extend(problems(SURFACE, n, v, p));
    }
    for n in [3, 8, 30] {
        let (v, p) = circle(n);
        suite.extend(problems(CIRCLE, n, v, p));
    }
    suite
}

#[test]
fn every_tiling_mode_agrees() {
    let suite = suite();
    assert!(suite.len() > 100);
    for (name, problem) in &suite {
        problem.assert_agreement(
            1e-7,
            &format!(
                "{name} of {} samples and {} variables",
                problem.samples,
                problem.model.variables()
            ),
        );
    }
}

#[test]
fn most_problems_converge() {
    // Guards against agreement in merely failing together.
    let suite = suite();
    let converged = suite
        .iter()
        .filter(|(_, problem)| problem.fit_rust(&mut String::new()).status.is_ok())
        .count();
    assert!(
        converged * 10 >= suite.len() * 9,
        "{converged} of {} converged",
        suite.len()
    );
}

#[test]
fn every_criterion_agrees() {
    for (name, problem) in suite() {
        for criterion in [
            Criterion::RelP,
            Criterion::RelPOrAbsP,
            Criterion::RelX,
            Criterion::RelXOrAbsX,
            Criterion::RelPAndRelX,
            Criterion::RelPOrAbsPAndRelXOrAbsX,
            Criterion::ChiSquared,
        ] {
            let mut problem = problem.clone();
            problem.algorithm.criterion = criterion;
            problem.assert_agreement(1e-7, &format!("{name} with {criterion:?}"));
        }
    }
}

#[test]
fn unknown_scale_agrees() {
    for (name, mut problem) in suite() {
        if problem.samples > problem.parameter.len() {
            problem.unknown_scale = true;
            problem.assert_agreement(1e-7, name);
        }
    }
}

#[test]
fn iteration_limit_agrees() {
    let (_name, problem) = &suite()[13 * 3 + 10];
    let iterations = problem.fit_rust(&mut String::new()).iterations;
    assert!((2..12).contains(&iterations));
    for limit in 0..12 {
        let mut problem = problem.clone();
        problem.algorithm.iteration_limit = limit;
        let fit = problem.assert_agreement(1e-7, &format!("limit {limit}"));
        assert_eq!(fit.status.is_ok(), limit >= iterations, "limit {limit}");
    }
}

#[test]
fn not_positive_definite_agrees() {
    let n = 6;
    let p = vec![1.0, 2.0];
    let v = sample(LINE, n, &p, |_, _, _| {});
    for (mode, tiles) in [
        (Mode::Diagonal, vec![(0, 0, vec![1e-3; n])]),
        (
            Mode::BlockDiagonal,
            vec![(0, 0, vec![1e-3; n]), (1, 1, vec![1.0; n * n])],
        ),
        (
            Mode::Diagonals,
            vec![
                (0, 0, vec![1e-3; n]),
                (1, 1, vec![1e-3; n]),
                (0, 1, vec![1.0; n]),
            ],
        ),
        (
            Mode::Full,
            vec![
                (0, 0, vec![1e-3; n]),
                (1, 1, vec![1e-3; n]),
                (0, 1, vec![1.0; n]),
            ],
        ),
    ] {
        let problem = Problem {
            model: LINE,
            samples: n,
            mode,
            tiles,
            sample: v.clone(),
            parameter: p.clone(),
            unknown_scale: false,
            algorithm: Algorithm::default(),
        };
        let fit = problem.assert_agreement(1e-7, "not positive definite");
        assert_eq!(Err(AlgorithmError::NumericalError), fit.status);
    }
}

fn not_finite(problem: Model<'static>) {
    let n = 6;
    let p = vec![1.0, 2.0];
    let v = sample(LINE, n, &p, |_, _, _| {});
    let problem = Problem {
        model: problem,
        samples: n,
        mode: Mode::Diagonal,
        tiles: vec![(0, 0, vec![1e-3; n]), (1, 1, vec![1e-3; n])],
        sample: v,
        parameter: p,
        unknown_scale: false,
        algorithm: Algorithm::default(),
    };
    let fit = problem.assert_agreement(1e-7, "not finite");
    assert_eq!(Err(AlgorithmError::NumericalError), fit.status);
}

#[test]
fn not_finite_agrees() {
    not_finite(Model {
        fx: &|_, _| f64::NAN,
        ..LINE
    });
    not_finite(Model {
        fx: &|_, _| f64::INFINITY,
        ..LINE
    });
    not_finite(Model {
        dfdx: &[&|_, _| f64::NAN],
        ..LINE
    });
    not_finite(Model {
        dfdp: &[&|_, _| f64::INFINITY, &|_, _| 1.0],
        ..LINE
    });
}

/// Parameter without effect on the model yielding an exactly singular factor.
#[test]
fn exactly_singular_factor_agrees() {
    const MODEL: Model = Model {
        fx: &|x, p| p[0] * x[0],
        dfdx: &[&|_, p| p[0]],
        dfdp: &[&|x, _| x[0], &|_, _| 0.0],
        implicit: false,
    };
    let n = 6;
    let p = vec![2.0, 1.0];
    let problem = Problem {
        model: MODEL,
        samples: n,
        mode: Mode::Diagonal,
        tiles: vec![(0, 0, vec![1e-3; n]), (1, 1, vec![1e-3; n])],
        sample: sample(MODEL, n, &p, |_, _, _| {}),
        parameter: p,
        unknown_scale: false,
        algorithm: Algorithm::default(),
    };
    let rust = problem.fit_rust(&mut String::new());
    let c = problem.fit_c(None);
    assert_eq!(c.status, rust.status);
    assert!(rust.status.is_err());
}

const PEARSON_YORK_X: [f64; 10] = [0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4];
const PEARSON_YORK_Y: [f64; 10] = [5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5];
const PEARSON_YORK_XCOV: [f64; 10] = [
    0.001,
    0.001,
    0.002,
    0.00125,
    0.005,
    0.0125,
    0.016_666_666_666_7,
    0.05,
    0.555_555_555_556,
    1.0,
];
const PEARSON_YORK_YCOV: [f64; 10] = [
    1.0,
    0.555_555_555_556,
    0.25,
    0.125,
    0.05,
    0.05,
    0.014_285_714_285_7,
    0.014_285_714_285_7,
    0.01,
    0.002,
];

fn pearson_york(model: Model<'static>, parameter: Vec<f64>) -> Problem {
    Problem {
        model,
        samples: 10,
        mode: Mode::Diagonal,
        tiles: vec![
            (0, 0, PEARSON_YORK_XCOV.to_vec()),
            (1, 1, PEARSON_YORK_YCOV.to_vec()),
        ],
        sample: PEARSON_YORK_X
            .iter()
            .chain(&PEARSON_YORK_Y)
            .copied()
            .collect(),
        parameter,
        unknown_scale: false,
        algorithm: Algorithm::default(),
    }
}

/// `examples/straight-line/pearson-york` of the C/Fortran library.
#[test]
fn pearson_york_line_agrees_with_literature() {
    let fit = pearson_york(LINE, vec![0.0, 0.0]).assert_agreement(1e-10, "line");
    fit.status.unwrap();
    // York's least-squares-cubic solution, Lybanon (1984).
    assert!(
        (fit.parameter[0] + 0.480_533_4).abs() < 1e-6,
        "{}",
        fit.parameter[0]
    );
    assert!(
        (fit.parameter[1] - 5.479_910_2).abs() < 1e-6,
        "{}",
        fit.parameter[1]
    );
}

/// Polynomials whose matrix of derivatives in parameters is badly conditioned.
#[test]
fn pearson_york_polynomials_agree() {
    let cubic = pearson_york(CUBIC, vec![5.0, -0.5, 0.0, 0.0]);
    cubic.assert_agreement(1e-8, "cubic").status.unwrap();
    let quintic = pearson_york(QUINTIC, vec![5.0, -0.5, 0.0, 0.0, 0.0, 0.0]);
    quintic.assert_agreement(1e-8, "quintic").status.unwrap();
}

/// `examples/straight-line/ISO_TS_28037` of the C/Fortran library.
#[test]
fn iso_ts_28037_agrees() {
    let x = [50.4, 99.0, 149.9, 200.4, 248.5, 299.7, 349.1];
    let y = [52.3, 97.8, 149.7, 200.1, 250.4, 300.9, 349.2];
    #[rustfmt::skip]
    let xcov = vec![
        0.50, 0.00, 0.25, 0.00, 0.25, 0.00, 0.25,
        0.00, 1.25, 1.00, 0.00, 0.00, 1.00, 1.00,
        0.25, 1.00, 1.50, 0.00, 0.25, 1.00, 1.25,
        0.00, 0.00, 0.00, 1.25, 1.00, 1.00, 1.00,
        0.25, 0.00, 0.25, 1.00, 1.50, 1.00, 1.25,
        0.00, 1.00, 1.00, 1.00, 1.00, 2.25, 2.00,
        0.25, 1.00, 1.25, 1.00, 1.25, 2.00, 2.50,
    ];
    let mut ycov = vec![1.0; 49];
    for i in 0..7 {
        ycov[i * 7 + i] = 5.0;
    }
    for mode in [Mode::BlockDiagonal, Mode::Full] {
        let problem = Problem {
            model: LINE,
            samples: 7,
            mode,
            tiles: vec![(0, 0, xcov.clone()), (1, 1, ycov.clone())],
            sample: x.iter().chain(&y).copied().collect(),
            parameter: vec![1.0, 0.0],
            unknown_scale: false,
            algorithm: Algorithm::default(),
        };
        let fit = problem.assert_agreement(1e-10, "ISO/TS 28037");
        fit.status.unwrap();
    }
}

/// `examples/nonlinear/nist_nls/diff-higher/BoxBOD` of the C/Fortran library from all its starts.
#[test]
fn boxbod_agrees_from_every_start() {
    let x = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0];
    let y = [109.0, 149.0, 149.0, 191.0, 213.0, 224.0];
    let mut converged = 0;
    for start in [[1.0, 1.0], [100.0, 0.75], [213.8, 0.5472]] {
        let problem = Problem {
            model: BOXBOD,
            samples: 6,
            mode: Mode::Diagonal,
            tiles: vec![(0, 0, vec![0.01; 6]), (1, 1, vec![1.0; 6])],
            sample: x.iter().chain(&y).copied().collect(),
            parameter: start.to_vec(),
            unknown_scale: false,
            algorithm: Algorithm::default(),
        };
        let fit = problem.assert_agreement(1e-7, "BoxBOD");
        converged += usize::from(fit.status.is_ok());
    }
    assert!(converged >= 2);
}

/// Text with numbers being compared as such since they may differ by rounding.
#[track_caller]
fn assert_logs_agree(actual: &str, expected: &str, what: &str) {
    // Workspace of the C library, whose last line is followed by an empty line.
    let mut lines = expected.lines();
    let mut expected = Vec::new();
    while let Some(line) = lines.next() {
        if line.starts_with("* workspace needed") {
            if line.contains("MQR2") {
                lines.next();
            }
        } else {
            expected.push(line);
        }
    }
    let actual = actual.lines().collect::<Vec<_>>();
    assert_eq!(
        expected.len(),
        actual.len(),
        "{what}\n{actual:#?}\n{expected:#?}"
    );
    for (a, e) in actual.iter().zip(&expected) {
        let (a_words, e_words) = (
            a.split_whitespace().collect::<Vec<_>>(),
            e.split_whitespace().collect::<Vec<_>>(),
        );
        assert_eq!(e_words.len(), a_words.len(), "{what}: {a:?} != {e:?}");
        for (a_word, e_word) in a_words.iter().zip(&e_words) {
            if a_word != e_word {
                let (a_number, e_number) = (
                    a_word.parse::<f64>().unwrap(),
                    e_word.parse::<f64>().unwrap(),
                );
                assert!(
                    (a_number.is_nan() && e_number.is_nan())
                        || (a_number - e_number).abs() <= 1e-6 * e_number.abs().max(1e-6),
                    "{what}: {a:?} != {e:?}",
                );
            }
        }
        // Same layout unless numbers differ.
        if a_words == e_words {
            assert_eq!(e, a, "{what}");
        }
    }
}

#[test]
fn log_agrees_at_every_verbosity() {
    let path = std::env::temp_dir().join(format!("oefpil-{}-differential", process::id()));
    let mut suite = suite();
    suite.retain(|(_, problem)| problem.samples <= 8);
    for (name, mut problem) in suite {
        for verbosity in [
            Verbosity::Silent,
            Verbosity::ParameterMean,
            Verbosity::ParameterCovariance,
            Verbosity::IndividualSteps,
        ] {
            problem.algorithm.verbosity = verbosity;
            let mut actual = String::new();
            problem.fit_rust(&mut actual);
            let is_silent = verbosity == Verbosity::Silent;
            problem.fit_c((!is_silent).then_some(path.as_path()));
            let expected = if is_silent {
                String::new()
            } else {
                fs::read_to_string(&path).unwrap()
            };
            assert_logs_agree(&actual, &expected, &format!("{name} at {verbosity:?}"));
        }
    }
    fs::remove_file(&path).unwrap();
}

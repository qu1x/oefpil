//! Correlations among more than two variables, checked against a dense reference.

#![allow(
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::many_single_char_names,
    clippy::needless_range_loop,
    clippy::suboptimal_flops
)]

use oefpil::{
    Algorithm, AlgorithmError, Covariance, CovarianceError, Criterion, Distribution, Model,
    OefpilError, Parameter, Report, Variable, Verbosity,
};
use rand::RngExt;
use rand_chacha::{ChaChaRng, rand_core::SeedableRng};
use rand_distr::StandardNormal;

const SAMPLES: usize = 12;

/// Explicit bivariate model $`Y=P^1+P^2X^1+P^3X^2`$, so $`|V|=3`$.
const EXPLICIT: Model = Model {
    fx: &|x, p| p[0] + p[1] * x[0] + p[2] * x[1],
    dfdx: &[&|_, p| p[1], &|_, p| p[2]],
    dfdp: &[&|_, _| 1.0, &|x, _| x[0], &|x, _| x[1]],
    implicit: false,
};

/// Implicit nonlinear model $`0=X^1+P^1X^2+P^2X^3+P^3X^1X^3-1`$, so $`|V|=3`$.
const IMPLICIT: Model = Model {
    fx: &|x, p| x[0] + p[0] * x[1] + p[1] * x[2] + p[2] * x[0] * x[2] - 1.0,
    dfdx: &[&|x, p| 1.0 + p[2] * x[2], &|_, p| p[0], &|x, p| {
        p[1] + p[2] * x[0]
    }],
    dfdp: &[&|x, _| x[1], &|x, _| x[2], &|x, _| x[0] * x[2]],
    implicit: true,
};

/// Explicit trivariate model $`Y=P^1+P^2X^1+P^3X^2+P^4X^3`$, so $`|V|=4`$.
const QUADRIVARIATE: Model = Model {
    fx: &|x, p| p[0] + p[1] * x[0] + p[2] * x[1] + p[3] * x[2],
    dfdx: &[&|_, p| p[1], &|_, p| p[2], &|_, p| p[3]],
    dfdp: &[&|_, _| 1.0, &|x, _| x[0], &|x, _| x[1], &|x, _| x[2]],
    implicit: false,
};

/// Deterministic observational noise.
fn noise(variable: usize, sample: usize) -> f64 {
    0.02 * ((variable * 7 + sample * 3) as f64 * 1.7).sin()
}

/// Noisy sample matrix of `model` around well-conditioned points.
fn sample(model: Model, p: &[f64]) -> Vec<f64> {
    let x_variables = model.dfdx.len();
    let mut v = vec![0f64; SAMPLES * model.variables()];
    let mut x = vec![0f64; x_variables];
    for sample in 0..SAMPLES {
        let s = sample as f64;
        for (variable, x) in x.iter_mut().enumerate() {
            *x =
                1.0 + 0.3 * (s * (1.1 + variable as f64)).sin() + 0.1 * s / (1.0 + variable as f64);
        }
        if model.implicit {
            // Solve for the last independent variable (the model is linear in `x[1]`).
            x[1] = (1.0 - x[0] - p[1] * x[2] - p[2] * x[0] * x[2]) / p[0];
        }
        for (variable, x) in x.iter().enumerate() {
            v[variable * SAMPLES + sample] = x + noise(variable, sample);
        }
        if !model.implicit {
            v[x_variables * SAMPLES + sample] = (model.fx)(&x, p) + noise(x_variables, sample);
        }
    }
    v
}

/// Correlation coefficient among variables `r` and `c`.
fn among(r: usize, c: usize) -> f64 {
    if r == c {
        1.0
    } else {
        0.4 / (1.0 + r.abs_diff(c) as f64)
    }
}

/// Fields of diagonal tile at `r` and `c`.
fn diagonal_tile(r: usize, c: usize) -> Vec<f64> {
    vec![among(r, c) * 1e-3; SAMPLES]
}

/// Fields of block tile at `r` and `c` being correlated `within` variables.
fn block_tile(r: usize, c: usize, within: f64) -> Vec<f64> {
    let mut tile = vec![0f64; SAMPLES * SAMPLES];
    for row in 0..SAMPLES {
        for column in 0..SAMPLES {
            tile[row * SAMPLES + column] =
                among(r, c) * 1e-3 * within.powi(row.abs_diff(column).try_into().unwrap());
        }
    }
    tile
}

/// Covariance matrix with upper tiles from `tile`.
fn tiled(
    mut covariance: Covariance,
    tile: impl Fn(usize, usize) -> Vec<f64>,
) -> Result<Covariance, OefpilError> {
    for r in 0..covariance.variables() {
        for c in r..covariance.variables() {
            covariance.set_tile(r, c, &tile(r, c))?;
        }
    }
    Ok(covariance)
}

struct Fit {
    report: Report,
    p: Vec<f64>,
    p_covariance: Vec<f64>,
    mean: Vec<f64>,
}

fn fit(model: Model, v: &[f64], covariance: &Covariance, p: &[f64]) -> Result<Fit, OefpilError> {
    let mut fit = Fit {
        report: Report::default(),
        p: p.to_vec(),
        p_covariance: vec![0f64; p.len() * p.len()],
        mean: v.to_vec(),
    };
    fit.report = Algorithm {
        tolerance: 1e-13,
        criterion: Criterion::RelPOrAbsPAndRelXOrAbsX,
        ..Algorithm::default()
    }
    .fit(
        model,
        &mut Variable {
            sample: v,
            mean: &mut fit.mean,
            covariance,
        },
        &mut Parameter {
            mean: &mut fit.p,
            covariance: &mut fit.p_covariance,
            ..Parameter::default()
        },
    )?;
    Ok(fit)
}

/// Asserts `a` and `b` agree to `tolerance` relative to `scale`.
#[track_caller]
fn assert_close(a: &[f64], b: &[f64], scale: f64, tolerance: f64) {
    assert_eq!(a.len(), b.len());
    for (index, (a, b)) in a.iter().zip(b).enumerate() {
        assert!(
            (a - b).abs() <= tolerance * scale,
            "field {index}: {a} != {b}"
        );
    }
}

/// Solves symmetric positive definite `a` of `order` for `b` by dense Cholesky decomposition.
fn solve(a: &[f64], order: usize, b: &mut [f64]) {
    let mut l = a.to_vec();
    for j in 0..order {
        for i in j..order {
            let d = (0..j).fold(l[i * order + j], |d, k| {
                d - l[i * order + k] * l[j * order + k]
            });
            l[i * order + j] = if i == j {
                assert!(d > 0.0);
                d.sqrt()
            } else {
                d / l[j * order + j]
            };
        }
    }
    for i in 0..order {
        let d = (0..i).fold(b[i], |d, k| d - l[i * order + k] * b[k]);
        b[i] = d / l[i * order + i];
    }
    for i in (0..order).rev() {
        let d = (i + 1..order).fold(b[i], |d, k| d - l[k * order + i] * b[k]);
        b[i] = d / l[i * order + i];
    }
}

/// Verifies the optimality conditions of the generalized least squares problem
///
/// ```text
/// minimize (v - mu)^T K^-1 (v - mu) subject to f(mu; p) = 0
/// ```
///
/// with dense linear algebra independent of the tiled routines, that is, the constraint, stationary
/// in `mu` and `p`, the sum of squares, and the parameter covariance matrix.
#[track_caller]
#[allow(clippy::too_many_lines)]
fn assert_optimal(model: Model, v: &[f64], covariance: &Covariance, fit: &Fit) {
    let variables = model.variables();
    let x_variables = model.dfdx.len();
    let parameters = fit.p.len();
    let order = SAMPLES * variables;
    let full = covariance.to_full();
    let at = |sample: usize| -> Vec<f64> {
        (0..x_variables)
            .map(|variable| fit.mean[variable * SAMPLES + sample])
            .collect()
    };
    // Constraint and Jacobians at optimum `B1` of `SAMPLES * order` and `B2` of `SAMPLES * parameters`.
    let mut b1 = vec![0f64; SAMPLES * order];
    let mut b2 = vec![0f64; SAMPLES * parameters];
    for sample in 0..SAMPLES {
        let x = at(sample);
        let mut constraint = (model.fx)(&x, &fit.p);
        if !model.implicit {
            constraint -= fit.mean[x_variables * SAMPLES + sample];
            b1[sample * order + x_variables * SAMPLES + sample] = -1.0;
        }
        assert!(constraint.abs() < 1e-9, "constraint {constraint}");
        for (variable, dfdx) in model.dfdx.iter().enumerate() {
            b1[sample * order + variable * SAMPLES + sample] = dfdx(&x, &fit.p);
        }
        for (parameter, dfdp) in model.dfdp.iter().enumerate() {
            b2[sample * parameters + parameter] = dfdp(&x, &fit.p);
        }
    }
    // Weighted residuals `w = K^-1 (v - mu)` yielding `chi^2 = (v - mu)^T w`.
    let residual = v
        .iter()
        .zip(&fit.mean)
        .map(|(v, m)| v - m)
        .collect::<Vec<_>>();
    let mut w = residual.clone();
    solve(&full, order, &mut w);
    let chi_squared = residual.iter().zip(&w).map(|(r, w)| r * w).sum::<f64>();
    assert_close(&[fit.report.chi_squared], &[chi_squared], chi_squared, 1e-8);
    // Schur complement `M = B1 K B1^T` of `SAMPLES * SAMPLES`.
    let mut m = vec![0f64; SAMPLES * SAMPLES];
    for r in 0..SAMPLES {
        for c in 0..SAMPLES {
            for i in 0..order {
                for j in 0..order {
                    m[r * SAMPLES + c] +=
                        b1[r * order + i] * full[i * order + j] * b1[c * order + j];
                }
            }
        }
    }
    // Lagrange multipliers from least squares of `w = B1^T l`.
    let mut l = vec![0f64; SAMPLES];
    let k_w = (0..order)
        .map(|i| (0..order).map(|j| full[i * order + j] * w[j]).sum::<f64>())
        .collect::<Vec<_>>();
    for r in 0..SAMPLES {
        l[r] = (0..order).map(|i| b1[r * order + i] * k_w[i]).sum();
    }
    solve(&m, SAMPLES, &mut l);
    // Stationarity in `mu`: `w = B1^T l` and in `p`: `B2^T l = 0`.
    let b1_l = (0..order)
        .map(|i| (0..SAMPLES).map(|r| b1[r * order + i] * l[r]).sum::<f64>())
        .collect::<Vec<_>>();
    let scale = w.iter().fold(1f64, |scale, w| scale.max(w.abs()));
    assert_close(&w, &b1_l, scale, 1e-7);
    for parameter in 0..parameters {
        let gradient = (0..SAMPLES)
            .map(|r| b2[r * parameters + parameter] * l[r])
            .sum::<f64>();
        assert!(gradient.abs() < 1e-7 * scale, "gradient {gradient}");
    }
    // Parameter covariance matrix `(B2^T M^-1 B2)^-1`.
    let mut information = vec![0f64; parameters * parameters];
    for c in 0..parameters {
        let mut column = (0..SAMPLES)
            .map(|r| b2[r * parameters + c])
            .collect::<Vec<_>>();
        solve(&m, SAMPLES, &mut column);
        for r in 0..parameters {
            information[r * parameters + c] = (0..SAMPLES)
                .map(|s| b2[s * parameters + r] * column[s])
                .sum();
        }
    }
    for c in 0..parameters {
        let mut column = vec![0f64; parameters];
        column[c] = 1.0;
        solve(&information, parameters, &mut column);
        for r in 0..parameters {
            let expected = column[r];
            let actual = fit.p_covariance[r * parameters + c];
            assert!(
                (actual - expected).abs() <= 1e-6 * expected.abs().max(1e-12),
                "covariance ({r}, {c}): {actual} != {expected}"
            );
        }
    }
}

#[test]
fn diagonals_and_full_are_optimal_for_any_number_of_variables() -> Result<(), OefpilError> {
    for (model, p, start) in [
        (EXPLICIT, &[1.0, 2.0, -1.0][..], &[0.9, 2.1, -0.9][..]),
        (IMPLICIT, &[0.8, -0.5, 0.3][..], &[0.9, -0.4, 0.2][..]),
        (
            QUADRIVARIATE,
            &[1.0, 2.0, -1.0, 0.5][..],
            &[0.9, 2.1, -0.9, 0.6][..],
        ),
    ] {
        let variables = model.variables();
        let v = sample(model, p);
        let diagonals = tiled(Covariance::new_diagonals(SAMPLES, variables), diagonal_tile)?;
        let full_diagonal = tiled(Covariance::new_full(SAMPLES, variables), diagonal_tile)?;
        let full_block = tiled(Covariance::new_full(SAMPLES, variables), |r, c| {
            block_tile(r, c, 0.0)
        })?;
        let correlated = tiled(Covariance::new_full(SAMPLES, variables), |r, c| {
            block_tile(r, c, 0.3)
        })?;
        // All tiling modes of the same matrix yield the same fit which is optimal.
        assert_eq!(diagonals.to_full(), full_diagonal.to_full());
        assert_eq!(diagonals.to_full(), full_block.to_full());
        let reference = fit(model, &v, &diagonals, start)?;
        assert_optimal(model, &v, &diagonals, &reference);
        for covariance in [&full_diagonal, &full_block] {
            let fit = fit(model, &v, covariance, start)?;
            assert_optimal(model, &v, covariance, &fit);
            assert_close(&fit.p, &reference.p, 1.0, 1e-9);
            assert_close(&fit.mean, &reference.mean, 1.0, 1e-9);
        }
        // Correlations within and among variables at once.
        let fit = fit(model, &v, &correlated, start)?;
        assert_optimal(model, &v, &correlated, &fit);
        assert!((fit.p[0] - reference.p[0]).abs() > 1e-9);
    }
    Ok(())
}

#[test]
fn null_tiles_among_variables_coincide_with_diagonal() -> Result<(), OefpilError> {
    let p = [1.0, 2.0, -1.0];
    let v = sample(EXPLICIT, &p);
    let tile = |r, c| {
        if r == c {
            diagonal_tile(r, c)
        } else {
            vec![0f64; SAMPLES]
        }
    };
    let mut diagonal = Covariance::new_diagonal(SAMPLES, 3);
    let mut block_diagonal = Covariance::new_block_diagonal(SAMPLES, 3);
    for variable in 0..3 {
        diagonal.set_tile(variable, variable, &diagonal_tile(variable, variable))?;
        block_diagonal.set_tile(variable, variable, &diagonal_tile(variable, variable))?;
    }
    let reference = fit(EXPLICIT, &v, &diagonal, &p)?;
    assert_close(
        &fit(EXPLICIT, &v, &block_diagonal, &p)?.p,
        &reference.p,
        1.0,
        1e-9,
    );
    let diagonals = tiled(Covariance::new_diagonals(SAMPLES, 3), tile)?;
    assert_close(
        &fit(EXPLICIT, &v, &diagonals, &p)?.p,
        &reference.p,
        1.0,
        1e-9,
    );
    let full = tiled(Covariance::new_full(SAMPLES, 3), tile)?;
    assert_close(&fit(EXPLICIT, &v, &full, &p)?.p, &reference.p, 1.0, 1e-9);
    Ok(())
}

#[test]
fn not_positive_definite_is_numerical_error() -> Result<(), OefpilError> {
    // Pairwise valid correlation coefficients which are jointly invalid.
    let tile = |r: usize, c: usize| vec![if r == c { 1e-3 } else { -0.9e-3 }; SAMPLES];
    let p = [1.0, 2.0, -1.0];
    let v = sample(EXPLICIT, &p);
    for covariance in [
        tiled(Covariance::new_diagonals(SAMPLES, 3), tile)?,
        tiled(Covariance::new_full(SAMPLES, 3), tile)?,
    ] {
        assert!(matches!(
            fit(EXPLICIT, &v, &covariance, &p),
            Err(OefpilError::Algorithm(AlgorithmError::NumericalError))
        ));
    }
    Ok(())
}

#[test]
fn decomposition_of_three_variables() -> Result<(), OefpilError> {
    let covariance = tiled(Covariance::new_diagonals(SAMPLES, 3), diagonal_tile)?;
    covariance.with_decomposition()?;
    let covariance = tiled(Covariance::new_full(SAMPLES, 3), |r, c| {
        block_tile(r, c, 0.3)
    })?;
    covariance.with_decomposition()?;
    Ok(())
}

#[test]
fn asymmetric_tile_among_variables_is_row_major() -> Result<(), OefpilError> {
    let covariance = Covariance::new_full(2, 2)
        .with_tile(0, 0, &[1.0, 0.0, 0.0, 1.0])?
        .with_tile(0, 1, &[0.1, 0.2, 0.3, 0.4])?
        .with_tile(1, 1, &[1.0, 0.0, 0.0, 1.0])?;
    let full = [
        1.0, 0.0, 0.1, 0.2, //
        0.0, 1.0, 0.3, 0.4, //
        0.1, 0.3, 1.0, 0.0, //
        0.2, 0.4, 0.0, 1.0, //
    ];
    assert_eq!(full, covariance.to_full().as_slice());
    Ok(())
}

#[test]
fn asymmetric_tile_among_variables_is_optimal() -> Result<(), OefpilError> {
    let p = [1.0, 2.0, -1.0];
    let v = sample(EXPLICIT, &p);
    let covariance = tiled(Covariance::new_full(SAMPLES, 3), |r, c| {
        let mut tile = block_tile(r, c, 0.3);
        if r != c {
            // Upper triangle of tile only, hence asymmetric.
            for row in 0..SAMPLES {
                for column in row + 1..SAMPLES {
                    tile[row * SAMPLES + column] *= 0.5;
                }
            }
        }
        tile
    })?;
    let fit = fit(EXPLICIT, &v, &covariance, &[0.9, 2.1, -0.9])?;
    assert_optimal(EXPLICIT, &v, &covariance, &fit);
    Ok(())
}

#[test]
fn not_finite_beyond_first_sample_is_numerical_error() -> Result<(), OefpilError> {
    const THRESHOLD: f64 = 1.2;
    const FX: Model = Model {
        fx: &|x, p| {
            if x[0] > THRESHOLD {
                f64::NAN
            } else {
                p[0] + p[1] * x[0] + p[2] * x[1]
            }
        },
        ..EXPLICIT
    };
    const DFDX: Model = Model {
        dfdx: &[
            &|x, p| {
                if x[0] > THRESHOLD {
                    f64::INFINITY
                } else {
                    p[1]
                }
            },
            &|_, p| p[2],
        ],
        ..EXPLICIT
    };
    const DFDP: Model = Model {
        dfdp: &[
            &|x, _| if x[0] > THRESHOLD { f64::NAN } else { 1.0 },
            &|x, _| x[0],
            &|x, _| x[1],
        ],
        ..EXPLICIT
    };
    let p = [1.0, 2.0, -1.0];
    let v = sample(EXPLICIT, &p);
    // First sample is finite, a later one is not.
    assert!(v[0] <= THRESHOLD && v[..SAMPLES].iter().any(|&x| x > THRESHOLD));
    let covariance = tiled(Covariance::new_diagonals(SAMPLES, 3), diagonal_tile)?;
    // Stops at once rather than discovering it via the parameters in the next iteration, if any.
    let oefpil = Algorithm {
        iteration_limit: 1,
        ..Algorithm::default()
    };
    for model in [FX, DFDX, DFDP] {
        let result = oefpil.fit(
            model,
            &mut Variable {
                sample: &v,
                mean: &mut [],
                covariance: &covariance,
            },
            &mut Parameter {
                mean: &mut p.clone(),
                covariance: &mut [0f64; 9],
                ..Parameter::default()
            },
        );
        assert!(matches!(
            result,
            Err(OefpilError::Algorithm(AlgorithmError::NumericalError))
        ));
    }
    Ok(())
}

#[test]
fn not_positive_definite_tile_before_last_is_numerical_error() -> Result<(), OefpilError> {
    let p = [1.0, 2.0, -1.0];
    let v = sample(EXPLICIT, &p);
    let mut invalid = block_tile(0, 0, 0.3);
    // Tile is indefinite on its own whereas `B1 K B1^T` remains positive definite.
    invalid[1] = 1.2e-3;
    invalid[SAMPLES] = 1.2e-3;
    for variable in 0..2 {
        let mut covariance = Covariance::new_block_diagonal(SAMPLES, 3);
        for tile in 0..3 {
            let fields = if tile == variable {
                invalid.clone()
            } else {
                block_tile(tile, tile, 0.3)
            };
            covariance.set_tile(tile, tile, &fields)?;
        }
        assert!(matches!(
            fit(EXPLICIT, &v, &covariance, &p),
            Err(OefpilError::Algorithm(AlgorithmError::NumericalError))
        ));
    }
    Ok(())
}

#[test]
fn not_positive_variance_of_diagonal_is_numerical_error() -> Result<(), OefpilError> {
    let p = [1.0, 2.0, -1.0];
    let v = sample(EXPLICIT, &p);
    for variance in [0.0, -1e-3, f64::NAN] {
        for variable in 0..3 {
            let mut covariance = Covariance::new_diagonal(SAMPLES, 3);
            for tile in 0..3 {
                let mut fields = diagonal_tile(tile, tile);
                if tile == variable {
                    fields[SAMPLES / 2] = variance;
                }
                covariance.set_tile(tile, tile, &fields)?;
            }
            assert!(matches!(
                fit(EXPLICIT, &v, &covariance, &p),
                Err(OefpilError::Algorithm(AlgorithmError::NumericalError))
            ));
        }
    }
    Ok(())
}

/// Number of samples of each group adding up to [`SAMPLES`].
const GROUPS: [usize; 4] = [3, 1, 5, 3];

/// Fields of the blocks of [`GROUPS`] of a block tile at `r` and `c` being asymmetric among
/// variables.
fn grouped_tile(r: usize, c: usize) -> Vec<f64> {
    let mut tile = Vec::new();
    let mut start = 0;
    for group in GROUPS {
        for row in 0..group {
            for column in 0..group {
                let asymmetry = if r != c && row < column { 0.5 } else { 1.0 };
                tile.push(
                    asymmetry
                        * among(r, c)
                        * (1.0 + 0.1 * (start + row.min(column)) as f64)
                        * 1e-3
                        * 0.3f64.powi(row.abs_diff(column).try_into().unwrap()),
                );
            }
        }
        start += group;
    }
    tile
}

/// Fields of a block tile of all samples of [`grouped_tile()`].
fn ungrouped_tile(r: usize, c: usize) -> Vec<f64> {
    let grouped = grouped_tile(r, c);
    let mut tile = vec![0f64; SAMPLES * SAMPLES];
    let (mut start, mut offset) = (0, 0);
    for group in GROUPS {
        for row in 0..group {
            for column in 0..group {
                tile[(start + row) * SAMPLES + start + column] =
                    grouped[offset + row * group + column];
            }
        }
        start += group;
        offset += group * group;
    }
    tile
}

#[test]
fn grouped_is_optimal_and_agrees_with_full() -> Result<(), OefpilError> {
    for (model, p, start) in [
        (EXPLICIT, &[1.0, 2.0, -1.0][..], &[0.9, 2.1, -0.9][..]),
        (IMPLICIT, &[0.8, -0.5, 0.3][..], &[0.9, -0.4, 0.2][..]),
        (
            QUADRIVARIATE,
            &[1.0, 2.0, -1.0, 0.5][..],
            &[0.9, 2.1, -0.9, 0.6][..],
        ),
    ] {
        let variables = model.variables();
        let v = sample(model, p);
        let grouped = tiled(Covariance::new_grouped(&GROUPS, variables), grouped_tile)?;
        let full = tiled(Covariance::new_full(SAMPLES, variables), ungrouped_tile)?;
        assert_eq!(full.to_full(), grouped.to_full());
        let fit_grouped = fit(model, &v, &grouped, start)?;
        assert_optimal(model, &v, &grouped, &fit_grouped);
        let fit_full = fit(model, &v, &full, start)?;
        assert_close(&fit_grouped.p, &fit_full.p, 1.0, 1e-11);
        assert_close(&fit_grouped.mean, &fit_full.mean, 1.0, 1e-11);
        assert_close(
            &fit_grouped.p_covariance,
            &fit_full.p_covariance,
            1e-3,
            1e-9,
        );
        assert_close(
            &[fit_grouped.report.chi_squared],
            &[fit_full.report.chi_squared],
            fit_full.report.chi_squared,
            1e-9,
        );
    }
    Ok(())
}

#[test]
fn grouped_coincides_with_other_modes_at_its_limits() -> Result<(), OefpilError> {
    let p = [1.0, 2.0, -1.0];
    let start = [0.9, 2.1, -0.9];
    let v = sample(EXPLICIT, &p);
    // One group of all samples is a full matrix.
    let tile = |r, c| block_tile(r, c, 0.3);
    let one = tiled(Covariance::new_grouped(&[SAMPLES], 3), tile)?;
    let full = tiled(Covariance::new_full(SAMPLES, 3), tile)?;
    assert_eq!(full.to_full(), one.to_full());
    let reference = fit(EXPLICIT, &v, &full, &start)?;
    let fitted = fit(EXPLICIT, &v, &one, &start)?;
    assert_close(&fitted.p, &reference.p, 1.0, 1e-11);
    assert_close(&fitted.mean, &reference.mean, 1.0, 1e-11);
    // Groups of one sample each are diagonal tiles all over.
    let singles = tiled(Covariance::new_grouped(&[1; SAMPLES], 3), diagonal_tile)?;
    let diagonals = tiled(Covariance::new_diagonals(SAMPLES, 3), diagonal_tile)?;
    assert_eq!(diagonals.to_full(), singles.to_full());
    let reference = fit(EXPLICIT, &v, &diagonals, &start)?;
    let fitted = fit(EXPLICIT, &v, &singles, &start)?;
    assert_close(&fitted.p, &reference.p, 1.0, 1e-11);
    assert_close(&fitted.mean, &reference.mean, 1.0, 1e-11);
    Ok(())
}

#[test]
fn grouped_tile_forms_agree() -> Result<(), OefpilError> {
    let block = [4.0, 1.0, 1.0, 3.0, 2.0];
    let half = [4.0, 1.0, 3.0, 2.0];
    let diagonal = [4.0, 3.0, 2.0];
    let first = Covariance::new_grouped(&[2, 1], 1).with_tile(0, 0, &block)?;
    let second = Covariance::new_grouped(&[2, 1], 1).with_tile(0, 0, &half)?;
    assert_eq!(first.to_full(), second.to_full());
    assert_eq!(
        [4.0, 1.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 2.0],
        first.to_full().as_slice()
    );
    // A diagonal tile clears what a previous block tile left behind.
    let third = first.with_tile(0, 0, &diagonal)?;
    assert_eq!(
        [4.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 2.0],
        third.to_full().as_slice()
    );
    Ok(())
}

#[test]
fn grouped_rejects_wrong_tiles() {
    let covariance = || Covariance::new_grouped(&[2, 1], 2);
    assert!(matches!(
        covariance().with_tile(0, 0, &[1.0; 2]),
        Err(OefpilError::Covariance(CovarianceError::MismatchingTile))
    ));
    assert!(matches!(
        covariance().with_tile(0, 0, &[1.0, 2.0, 3.0, 4.0, 5.0]),
        Err(OefpilError::Covariance(CovarianceError::AsymmetricTile))
    ));
    assert!(matches!(
        covariance().with_tile(0, 2, &[1.0; 3]),
        Err(OefpilError::Covariance(
            CovarianceError::OutOfBoundTileIndex
        ))
    ));
    // Asymmetry among variables is valid.
    assert!(
        covariance()
            .with_tile(0, 1, &[1.0, 2.0, 3.0, 4.0, 5.0])
            .is_ok()
    );
}

#[test]
#[should_panic(expected = "assertion failed")]
fn grouped_rejects_empty_group() {
    let _ = Covariance::new_grouped(&[2, 0, 1], 2);
}

#[test]
#[should_panic(expected = "assertion")]
fn grouped_rejects_no_groups() {
    let _ = Covariance::new_grouped(&[], 2);
}

#[test]
fn grouped_not_positive_definite_is_numerical_error() -> Result<(), OefpilError> {
    let p = [1.0, 2.0, -1.0];
    let v = sample(EXPLICIT, &p);
    // Only the last group fails.
    let tile = |r: usize, c: usize| {
        let mut tile = vec![if r == c { 1e-3 } else { 0.0 }; SAMPLES];
        if r != c {
            tile[SAMPLES - 1] = -0.9e-3;
        }
        tile
    };
    let covariance = tiled(Covariance::new_grouped(&GROUPS, 3), tile)?;
    assert!(matches!(
        fit(EXPLICIT, &v, &covariance, &p),
        Err(OefpilError::Algorithm(AlgorithmError::NumericalError))
    ));
    assert!(covariance.with_decomposition().is_err());
    Ok(())
}

#[test]
fn grouped_samples_like_full() -> Result<(), OefpilError> {
    let grouped = tiled(Covariance::new_grouped(&GROUPS, 3), grouped_tile)?.with_decomposition()?;
    let full = tiled(Covariance::new_full(SAMPLES, 3), ungrouped_tile)?.with_decomposition()?;
    let mean = (0..3 * SAMPLES)
        .map(|field| field as f64)
        .collect::<Vec<f64>>();
    let draw = |covariance| -> Result<Vec<f64>, OefpilError> {
        let mut sample = vec![0f64; 3 * SAMPLES];
        let mut field = 0;
        Distribution {
            sample: &mut sample,
            mean: &mean,
            covariance,
        }
        .sample(|| {
            field += 1;
            (f64::from(field) * 1.7).sin()
        })?;
        Ok(sample)
    };
    assert_close(&draw(&grouped)?, &draw(&full)?, 1.0, 1e-13);
    Ok(())
}

/// Draws a sample with reproducible standard normal `values`.
fn draw(covariance: &Covariance, mean: &[f64], values: &[f64]) -> Result<Vec<f64>, OefpilError> {
    let mut sample = vec![0f64; mean.len()];
    let mut values = values.iter();
    Distribution {
        sample: &mut sample,
        mean,
        covariance,
    }
    .sample(|| *values.next().unwrap())?;
    Ok(sample)
}

#[test]
#[allow(clippy::too_many_lines)]
fn grouped_agrees_with_full_for_random_groups_and_matrices() -> Result<(), OefpilError> {
    const TRIALS: usize = 300;
    let problems = [
        (EXPLICIT, &[1.0, 2.0, -1.0][..], &[0.9, 2.1, -0.9][..]),
        (IMPLICIT, &[0.8, -0.5, 0.3][..], &[0.9, -0.4, 0.2][..]),
        (
            QUADRIVARIATE,
            &[1.0, 2.0, -1.0, 0.5][..],
            &[0.9, 2.1, -0.9, 0.6][..],
        ),
    ];
    let mut random = ChaChaRng::seed_from_u64(0);
    let mut converged = 0;
    for trial in 0..TRIALS {
        let (model, p, start) = problems[trial % problems.len()];
        let variables = model.variables();
        // A random partition of the samples into groups.
        let mut groups = Vec::new();
        let mut left = SAMPLES;
        while left > 0 {
            let group = random.random_range(1..=left.min(6));
            groups.push(group);
            left -= group;
        }
        // A random positive definite block per group, which is asymmetric among variables.
        let mut grouped_tiles = vec![Vec::new(); variables * variables];
        let mut full_tiles = vec![vec![0f64; SAMPLES * SAMPLES]; variables * variables];
        let mut first = 0;
        for &group in &groups {
            let order = variables * group;
            let root = (0..order * order)
                .map(|_| random.sample(StandardNormal))
                .collect::<Vec<f64>>();
            for r in 0..variables {
                for c in 0..variables {
                    for i in 0..group {
                        for j in 0..group {
                            let (row, column) = (r * group + i, c * group + j);
                            let product = (0..order)
                                .map(|k| root[row * order + k] * root[column * order + k])
                                .sum::<f64>();
                            let variance = if row == column { 1e-3 } else { 0.0 };
                            let field = 1e-4 * product + variance;
                            grouped_tiles[r * variables + c].push(field);
                            full_tiles[r * variables + c][(first + i) * SAMPLES + first + j] =
                                field;
                        }
                    }
                }
            }
            first += group;
        }
        let grouped = tiled(Covariance::new_grouped(&groups, variables), |r, c| {
            grouped_tiles[r * variables + c].clone()
        })?
        .with_decomposition()?;
        let full = tiled(Covariance::new_full(SAMPLES, variables), |r, c| {
            full_tiles[r * variables + c].clone()
        })?
        .with_decomposition()?;
        assert_eq!(full.to_full(), grouped.to_full(), "groups {groups:?}");
        // Same draws from both.
        let order = SAMPLES * variables;
        let mean = sample(model, p);
        let values = (0..order)
            .map(|_| random.sample(StandardNormal))
            .collect::<Vec<f64>>();
        let v = draw(&grouped, &mean, &values)?;
        assert_close(&v, &draw(&full, &mean, &values)?, 1.0, 1e-11);
        // Same fits of which every seventh is of unknown scale.
        let (grouped, full, unknown_scale) = if trial % 7 == 0 {
            (
                grouped.with_unknown_scale(),
                full.with_unknown_scale(),
                true,
            )
        } else {
            (grouped, full, false)
        };
        match (
            fit(model, &v, &grouped, start),
            fit(model, &v, &full, start),
        ) {
            (Ok(fit_grouped), Ok(fit_full)) => {
                converged += 1;
                assert_close(&fit_grouped.p, &fit_full.p, 1.0, 1e-9);
                assert_close(&fit_grouped.mean, &fit_full.mean, 1.0, 1e-9);
                assert_close(
                    &fit_grouped.p_covariance,
                    &fit_full.p_covariance,
                    fit_full.p_covariance[0].abs(),
                    1e-6,
                );
                assert_close(
                    &[fit_grouped.report.chi_squared],
                    &[fit_full.report.chi_squared],
                    fit_full.report.chi_squared,
                    1e-6,
                );
                if !unknown_scale {
                    assert_optimal(model, &v, &grouped, &fit_grouped);
                }
            }
            (Err(_), Err(_)) => {}
            _ => panic!("only one fit converged with groups {groups:?}"),
        }
    }
    assert!(converged > TRIALS / 2, "converged {converged}");
    Ok(())
}

#[test]
fn grouped_logs_its_mode() -> Result<(), OefpilError> {
    let p = [1.0, 2.0, -1.0];
    let v = sample(EXPLICIT, &p);
    let covariance = tiled(Covariance::new_grouped(&GROUPS, 3), grouped_tile)?;
    let mut log = String::new();
    Algorithm {
        verbosity: Verbosity::IndividualSteps,
        ..Algorithm::default()
    }
    .fit_with_log(
        EXPLICIT,
        &mut Variable {
            sample: &v,
            mean: &mut [],
            covariance: &covariance,
        },
        &mut Parameter {
            mean: &mut p.clone(),
            covariance: &mut [0f64; 9],
            ..Parameter::default()
        },
        &mut log,
    )?;
    assert!(log.contains("tiling mode: 5"));
    Ok(())
}

//! Fits halving steps which do not improve them.

#![allow(
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::suboptimal_flops
)]

use core::cell::RefCell;
use oefpil::{
    Algorithm, AlgorithmError, Covariance, Model, OefpilError, Parameter, Report, Variable,
    Verbosity,
};

const GAUSSIAN: Model = Model {
    fx: &|x, p| p[0] * (-(x[0] - p[1]).powi(2) / (2.0 * p[2] * p[2])).exp(),
    dfdx: &[],
    dfdp: &[],
    implicit: false,
};

const LOGISTIC: Model = Model {
    fx: &|x, p| p[0] / (1.0 + (-p[1] * (x[0] - p[2])).exp()),
    dfdx: &[],
    dfdp: &[],
    implicit: false,
};

const EXPONENTIAL: Model = Model {
    fx: &|x, p| p[0] * (-p[1] * x[0]).exp() + p[2],
    dfdx: &[],
    dfdp: &[],
    implicit: false,
};

/// Samples along `model` at `truth` for evenly spaced `x` with deterministic noise.
fn sample(model: Model, truth: &[f64], end: f64, samples: usize) -> Vec<f64> {
    let x = (0..samples).map(|i| end * i as f64 / (samples - 1) as f64);
    let noise = |i: usize| (1.7 * i as f64 + 0.3).sin();
    x.clone()
        .enumerate()
        .map(|(i, x)| x + 1e-3 * noise(i))
        .chain(
            x.enumerate()
                .map(|(i, x)| (model.fx)(&[x], truth) + 1e-2 * noise(i + samples)),
        )
        .collect()
}

fn diagonal(samples: usize) -> Result<Covariance, OefpilError> {
    Covariance::new_diagonal(samples, 2)
        .with_tile(0, 0, &vec![1e-6; samples])?
        .with_tile(1, 1, &vec![1e-4; samples])
}

/// Fits `model` from `start` returning the report and parameter mean.
fn fit(
    algorithm: Algorithm,
    model: Model,
    sample: &[f64],
    covariance: &Covariance,
    start: &[f64],
) -> Result<(Report, Vec<f64>), OefpilError> {
    let mut mean = start.to_vec();
    let report = algorithm.fit(
        model,
        &mut Variable {
            sample,
            mean: &mut [],
            covariance,
        },
        &mut Parameter {
            mean: &mut mean,
            covariance: &mut vec![0f64; start.len().pow(2)],
            ..Parameter::default()
        },
    )?;
    Ok((report, mean))
}

fn damped() -> Algorithm {
    Algorithm {
        damping: true,
        ..Algorithm::default()
    }
}

#[test]
fn is_disabled_by_default() {
    assert!(!Algorithm::default().damping);
}

#[test]
fn poor_starts_converge_about_twice_as_often_to_the_same_optimum() -> Result<(), OefpilError> {
    let factors = [0.1, 0.3, 0.5, 0.7, 1.5, 2.0, 3.0, 10.0];
    let (mut plain_count, mut damped_count) = (0, 0);
    for (model, truth, end) in [
        (GAUSSIAN, [5.0, 2.0, 0.7], 4.0),
        (LOGISTIC, [10.0, 1.5, 3.0], 6.0),
        (EXPONENTIAL, [4.0, 1.2, 1.0], 5.0),
    ] {
        let sample = sample(model, &truth, end, 25);
        let covariance = diagonal(25)?;
        let optimum = fit(Algorithm::default(), model, &sample, &covariance, &truth)?;
        for a in factors {
            for b in factors {
                for c in factors {
                    let start = [a * truth[0], b * truth[1], c * truth[2]];
                    let reached = |result: Result<(Report, Vec<f64>), OefpilError>| {
                        result.is_ok_and(|(report, _)| {
                            (report.chi_squared - optimum.0.chi_squared).abs()
                                <= 1e-6 * optimum.0.chi_squared
                        })
                    };
                    let plain = Algorithm::default();
                    plain_count +=
                        u32::from(reached(fit(plain, model, &sample, &covariance, &start)));
                    damped_count +=
                        u32::from(reached(fit(damped(), model, &sample, &covariance, &start)));
                }
            }
        }
    }
    assert!(
        f64::from(damped_count) > 1.5 * f64::from(plain_count),
        "damped {damped_count}, undamped {plain_count}"
    );
    Ok(())
}

#[test]
fn rescues_step_leaving_domain_of_model() -> Result<(), OefpilError> {
    // A full step from this start is beyond the domain of the logarithm of the parameter.
    const LOGARITHM: Model = Model {
        fx: &|x, p| p[0].ln() * x[0],
        dfdx: &[],
        dfdp: &[],
        implicit: false,
    };
    let x = [1.0, 2.0, 3.0, 4.0];
    let v = x
        .iter()
        .copied()
        .chain(x.iter().map(|x| x * 2f64.ln()))
        .collect::<Vec<f64>>();
    let covariance = Covariance::new_diagonal(4, 2)
        .with_tile(0, 0, &[1e-6; 4])?
        .with_tile(1, 1, &[1e-6; 4])?;
    let start = [20.0];
    assert!(matches!(
        fit(Algorithm::default(), LOGARITHM, &v, &covariance, &start),
        Err(OefpilError::Algorithm(AlgorithmError::NumericalError))
    ));
    let (_, mean) = fit(damped(), LOGARITHM, &v, &covariance, &start)?;
    assert!((mean[0] - 2.0).abs() < 1e-6, "{mean:?}");
    Ok(())
}

#[test]
fn halves_step_of_parameters_and_variables_alike() -> Result<(), OefpilError> {
    let truth = [10.0, 1.5, 3.0];
    let v = sample(LOGISTIC, &truth, 6.0, 25);
    let covariance = diagonal(25)?;
    // The estimates of the first sample and parameter per iteration.
    let estimates = RefCell::new(Vec::new());
    let exp = |x: &[f64], p: &[f64]| (-p[1] * (x[0] - p[2])).exp();
    let fx = |x: &[f64], p: &[f64]| {
        estimates.borrow_mut().push((x[0], p[0]));
        p[0] / (1.0 + exp(x, p))
    };
    let dfdx = |x: &[f64], p: &[f64]| p[0] * p[1] * exp(x, p) / (1.0 + exp(x, p)).powi(2);
    let dfdp0 = |x: &[f64], p: &[f64]| 1.0 / (1.0 + exp(x, p));
    let dfdp1 = |x: &[f64], p: &[f64]| p[0] * (x[0] - p[2]) * exp(x, p) / (1.0 + exp(x, p)).powi(2);
    let dfdp2 = |x: &[f64], p: &[f64]| -p[0] * p[1] * exp(x, p) / (1.0 + exp(x, p)).powi(2);
    let model = Model {
        fx: &fx,
        dfdx: &[&dfdx],
        dfdp: &[&dfdp0, &dfdp1, &dfdp2],
        implicit: false,
    };
    fit(damped(), model, &v, &covariance, &[5.0, 0.75, 1.5])?;
    let estimates = estimates.into_inner();
    let first = estimates
        .iter()
        .step_by(25)
        .copied()
        .collect::<Vec<(f64, f64)>>();
    let midway =
        |a: f64, b: f64, c: f64| (c - f64::midpoint(a, b)).abs() <= 1e-12 * a.abs().max(b.abs());
    let mut halved = 0;
    for window in first.windows(3) {
        let [(x_0, p_0), (x_1, p_1), (x_2, p_2)] = window else {
            unreachable!();
        };
        // The start of a step, its rejected end, and the next trial.
        if (p_1 - p_0).abs() > 1e-6 && midway(*p_0, *p_1, *p_2) {
            halved += 1;
            assert!((x_1 - x_0).abs() > 1e-9, "{x_0} {x_1}");
            assert!(midway(*x_0, *x_1, *x_2), "{x_0} {x_1} {x_2}");
        }
    }
    assert!(halved > 0);
    Ok(())
}

#[test]
fn does_not_stall_near_solution() -> Result<(), OefpilError> {
    // Tiny correct steps can raise the effective chi-squared of a tiny residual.
    const LOGARITHM: Model = Model {
        fx: &|x, p| p[0] * x[0].ln() + p[1],
        dfdx: &[&|x, p| p[0] / x[0]],
        dfdp: &[&|x, _| x[0].ln(), &|_, _| 1.0],
        implicit: false,
    };
    let x = [0.5f64, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = x
        .iter()
        .enumerate()
        .map(|(i, x)| 2.0 * x.ln() + 1.0 + 0.01 * (1.7 * i as f64).sin());
    let v = x.iter().copied().chain(y).collect::<Vec<f64>>();
    let covariance = Covariance::new_diagonal(6, 2)
        .with_tile(0, 0, &[0.03 * 0.03; 6])?
        .with_tile(1, 1, &[1.0; 6])?;
    let start = [-2.0, -20.0];
    let (plain, plain_mean) = fit(Algorithm::default(), LOGARITHM, &v, &covariance, &start)?;
    let (report, mean) = fit(damped(), LOGARITHM, &v, &covariance, &start)?;
    assert_eq!(plain.iterations, report.iterations);
    assert_eq!(plain_mean, mean);
    Ok(())
}

#[test]
fn counts_trials_as_iterations() -> Result<(), OefpilError> {
    let truth = [10.0, 1.5, 3.0];
    let v = sample(LOGISTIC, &truth, 6.0, 25);
    let covariance = diagonal(25)?;
    let mut log = String::new();
    let report = Algorithm {
        verbosity: Verbosity::ParameterMean,
        ..damped()
    }
    .fit_with_log(
        LOGISTIC,
        &mut Variable {
            sample: &v,
            mean: &mut [],
            covariance: &covariance,
        },
        &mut Parameter {
            mean: &mut [5.0, 0.75, 1.5],
            covariance: &mut [0f64; 9],
            ..Parameter::default()
        },
        &mut log,
    )?;
    let trials = log.matches("damping: trying").count();
    let steps = log
        .lines()
        .filter(|line| line.starts_with(|c: char| c.is_ascii_digit()) && line.contains('\t'))
        .count();
    assert!(trials > 0);
    assert_eq!(steps + trials, report.iterations);
    Ok(())
}

#[test]
fn agrees_across_tiling_modes() -> Result<(), OefpilError> {
    let truth = [10.0, 1.5, 3.0];
    let v = sample(LOGISTIC, &truth, 6.0, 25);
    let block = Covariance::new_block_diagonal(25, 2)
        .with_tile(0, 0, &[1e-6; 25])?
        .with_tile(1, 1, &[1e-4; 25])?;
    let grouped = Covariance::new_grouped(&[5, 1, 12, 7], 2)
        .with_tile(0, 0, &[1e-6; 25])?
        .with_tile(1, 1, &[1e-4; 25])?;
    let start = [5.0, 0.75, 1.5];
    let (reference, mean) = fit(damped(), LOGISTIC, &v, &diagonal(25)?, &start)?;
    // Requires halving, so these cover its paths of a factorized dense matrix.
    assert!(fit(Algorithm::default(), LOGISTIC, &v, &diagonal(25)?, &start).is_err());
    for covariance in [&block, &grouped] {
        let (report, other) = fit(damped(), LOGISTIC, &v, covariance, &start)?;
        assert_eq!(reference.iterations, report.iterations);
        for (a, b) in mean.iter().zip(&other) {
            assert!((a - b).abs() < 1e-9, "{a} != {b}");
        }
    }
    Ok(())
}

#[test]
fn respects_iteration_limit() -> Result<(), OefpilError> {
    let truth = [10.0, 1.5, 3.0];
    let v = sample(LOGISTIC, &truth, 6.0, 25);
    let covariance = diagonal(25)?;
    let start = [5.0, 0.75, 1.5];
    let (report, _) = fit(damped(), LOGISTIC, &v, &covariance, &start)?;
    assert!(report.iterations > 3);
    let limited = Algorithm {
        iteration_limit: report.iterations - 1,
        ..damped()
    };
    assert!(matches!(
        fit(limited, LOGISTIC, &v, &covariance, &start),
        Err(OefpilError::Algorithm(AlgorithmError::IterationLimit))
    ));
    let exact = Algorithm {
        iteration_limit: report.iterations,
        ..damped()
    };
    assert_eq!(
        report.iterations,
        fit(exact, LOGISTIC, &v, &covariance, &start)?.0.iterations
    );
    Ok(())
}

#[test]
fn hopeless_start_still_fails() -> Result<(), OefpilError> {
    // Not defined at the initial estimate, so there is no step to halve.
    let truth = [5.0, 2.0, 0.7];
    let v = sample(GAUSSIAN, &truth, 4.0, 25);
    assert!(matches!(
        fit(damped(), GAUSSIAN, &v, &diagonal(25)?, &[5.0, 2.0, 0.0]),
        Err(OefpilError::Algorithm(AlgorithmError::NumericalError))
    ));
    Ok(())
}

#[test]
fn logs_halved_steps() -> Result<(), OefpilError> {
    let truth = [10.0, 1.5, 3.0];
    let v = sample(LOGISTIC, &truth, 6.0, 25);
    let covariance = diagonal(25)?;
    let run = |damping| -> Result<String, OefpilError> {
        let mut log = String::new();
        let _ = Algorithm {
            damping,
            verbosity: Verbosity::ParameterMean,
            ..Algorithm::default()
        }
        .fit_with_log(
            LOGISTIC,
            &mut Variable {
                sample: &v,
                mean: &mut [],
                covariance: &covariance,
            },
            &mut Parameter {
                mean: &mut [5.0, 0.75, 1.5],
                covariance: &mut [0f64; 9],
                ..Parameter::default()
            },
            &mut log,
        );
        Ok(log)
    };
    assert!(run(true)?.contains("damping: trying fraction 5.0"));
    assert!(!run(false)?.contains("damping"));
    Ok(())
}

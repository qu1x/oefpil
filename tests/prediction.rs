//! Predictions of the dependent variable at new samples after a fit.

#![allow(
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::suboptimal_flops
)]

use oefpil::{
    Algorithm, Covariance, Distribution, Model, OefpilError, Parameter, ParameterError,
    PredictionError, Variable,
};
use rand::RngExt;
use rand_chacha::{ChaChaRng, rand_core::SeedableRng};
use rand_distr::StandardNormal;

const LINE: Model = Model {
    fx: &|x, p| p[0] + p[1] * x[0],
    dfdx: &[&|_, p| p[1]],
    dfdp: &[&|_, _| 1.0, &|x, _| x[0]],
    implicit: false,
};

const GROWTH: Model = Model {
    fx: &|x, p| p[0] * (p[1] * x[0]).exp(),
    dfdx: &[&|x, p| p[0] * p[1] * (p[1] * x[0]).exp()],
    dfdp: &[&|x, p| (p[1] * x[0]).exp(), &|x, p| {
        p[0] * x[0] * (p[1] * x[0]).exp()
    }],
    implicit: false,
};

/// Samples along `model` at `p` for `x` from zero to two with deterministic noise.
fn sample(model: Model, p: &[f64], samples: usize) -> Vec<f64> {
    let x = (0..samples).map(|i| 2.0 * i as f64 / (samples - 1) as f64);
    let noise = |i: usize| 1e-2 * (1.7 * i as f64).sin();
    x.clone()
        .enumerate()
        .map(|(i, x)| x + 0.1 * noise(i))
        .chain(
            x.enumerate()
                .map(|(i, x)| (model.fx)(&[x], p) + noise(i + samples)),
        )
        .collect()
}

fn covariance(samples: usize) -> Result<Covariance, OefpilError> {
    Covariance::new_diagonal(samples, 2)
        .with_tile(0, 0, &vec![1e-6; samples])?
        .with_tile(1, 1, &vec![1e-4; samples])
}

/// Fits `model` to `sample` from `start` yielding the parameter mean and covariance.
fn fit(
    model: Model,
    sample: &[f64],
    covariance: &Covariance,
    start: &[f64],
) -> Result<(Vec<f64>, Vec<f64>), OefpilError> {
    let mut mean = start.to_vec();
    let mut parameter_covariance = vec![0f64; start.len().pow(2)];
    Algorithm::default().fit(
        model,
        &mut Variable {
            sample,
            mean: &mut [],
            covariance,
        },
        &mut Parameter {
            mean: &mut mean,
            covariance: &mut parameter_covariance,
            ..Parameter::default()
        },
    )?;
    Ok((mean, parameter_covariance))
}

/// Parameter of fitted `mean` and `covariance`.
fn parameter<'a>(mean: &'a mut [f64], covariance: &'a mut [f64]) -> Parameter<'a> {
    Parameter {
        mean,
        covariance,
        ..Parameter::default()
    }
}

#[test]
fn line_agrees_with_exact_propagation() -> Result<(), OefpilError> {
    let sample = sample(LINE, &[1.0, 2.0], 12);
    let (mut mean, mut c) = fit(LINE, &sample, &covariance(12)?, &[0.9, 2.1])?;
    // The parameters of a line starting at the origin of its data are correlated.
    assert!(c[1] / (c[0] * c[3]).sqrt() < -0.5);
    let parameter = parameter(&mut mean, &mut c);
    for x in [-3.0, 0.0, 1.0, 2.0, 30.0] {
        let prediction = LINE.predict(&[x], &parameter)?;
        let (m, c) = (&*parameter.mean, &*parameter.covariance);
        let variance = c[0] + 2.0 * x * c[1] + x * x * c[3];
        assert!((prediction.mean - (m[0] + m[1] * x)).abs() < 1e-14);
        assert!(
            (prediction.deviation - variance.sqrt()).abs() < 1e-14,
            "{x}: {} != {}",
            prediction.deviation,
            variance.sqrt(),
        );
    }
    Ok(())
}

#[test]
fn approximated_derivative_agrees_with_analytic_one() -> Result<(), OefpilError> {
    const APPROXIMATED: Model = Model {
        dfdp: &[],
        ..GROWTH
    };
    let sample = sample(GROWTH, &[2.0, 0.5], 12);
    let (mut mean, mut c) = fit(GROWTH, &sample, &covariance(12)?, &[1.8, 0.6])?;
    let parameter = parameter(&mut mean, &mut c);
    for x in [-1.0, 0.0, 1.0, 5.0] {
        let analytic = GROWTH.predict(&[x], &parameter)?;
        let approximated = APPROXIMATED.predict(&[x], &parameter)?;
        assert_eq!(analytic.mean.to_bits(), approximated.mean.to_bits());
        assert!(
            (analytic.deviation - approximated.deviation).abs() <= 1e-9 * analytic.deviation,
            "{x}: {} != {}",
            analytic.deviation,
            approximated.deviation,
        );
    }
    Ok(())
}

#[test]
fn step_follows_deviation_of_parameter_at_zero() -> Result<(), OefpilError> {
    const ONE: Model = Model {
        fx: &|x, p| p[0] * x[0],
        dfdx: &[],
        dfdp: &[],
        implicit: false,
    };
    const TINY: Model = Model {
        fx: &|x, p| (p[0] * 1e9).sin() * x[0],
        dfdx: &[],
        dfdp: &[],
        implicit: false,
    };
    // A step of magnitude one would wrap around the sine many times.
    let (mut mean, mut c) = ([0.0], [1e-20]);
    let prediction = TINY.predict(&[2.0], &parameter(&mut mean, &mut c))?;
    assert!((prediction.deviation - 2.0 * 1e9 * 1e-10).abs() < 1e-9);
    // Without any magnitude, it is assumed to be of magnitude one.
    let (mut mean, mut c) = ([0.0], [0.0]);
    assert_eq!(
        0.0,
        ONE.predict(&[2.0], &parameter(&mut mean, &mut c))?
            .deviation
    );
    Ok(())
}

#[test]
fn unknown_scale_scales_deviation() -> Result<(), OefpilError> {
    let sample = sample(LINE, &[1.0, 2.0], 12);
    let (mut mean, mut c) = fit(LINE, &sample, &covariance(12)?, &[0.9, 2.1])?;
    let known = LINE.predict(&[1.0], &parameter(&mut mean, &mut c))?;
    let scaled = covariance(12)?.with_unknown_scale();
    let (mut mean, mut scaled_c) = fit(LINE, &sample, &scaled, &[0.9, 2.1])?;
    let ratio = (scaled_c[0] / c[0]).sqrt();
    assert!((ratio - 1.0).abs() > 1e-3);
    let unknown = LINE.predict(&[1.0], &parameter(&mut mean, &mut scaled_c))?;
    assert!((unknown.deviation - ratio * known.deviation).abs() < 1e-12);
    Ok(())
}

#[test]
fn deviation_agrees_with_monte_carlo_simulation() -> Result<(), OefpilError> {
    const SAMPLES: usize = 20;
    const SIMULATIONS: u32 = 2000;
    let truth = [2.0, 0.5];
    let x = (0..SAMPLES).map(|i| 2.0 * i as f64 / (SAMPLES - 1) as f64);
    let true_mean = x
        .clone()
        .chain(x.map(|x| (GROWTH.fx)(&[x], &truth)))
        .collect::<Vec<f64>>();
    // Correlated within and among variables.
    let block = |deviation: f64| {
        let mut block = vec![0.5 * deviation * deviation; SAMPLES * SAMPLES];
        for sample in 0..SAMPLES {
            block[sample * SAMPLES + sample] = deviation * deviation;
        }
        block
    };
    let covariance = Covariance::new_full(SAMPLES, 2)
        .with_tile(0, 0, &block(1e-3))?
        .with_tile(0, 1, &[2e-6; SAMPLES])?
        .with_tile(1, 1, &block(1e-2))?
        .with_decomposition()?;
    let mut random = ChaChaRng::seed_from_u64(0);
    let mut sample = vec![0f64; 2 * SAMPLES];
    // Within the data and extrapolating.
    let points = [0.5, 1.5, 4.0];
    let mut squared_error = [0f64; 3];
    let mut reported = [0f64; 3];
    let mut covered = [0u32; 3];
    for _ in 0..SIMULATIONS {
        Distribution {
            sample: &mut sample,
            mean: &true_mean,
            covariance: &covariance,
        }
        .sample(|| random.sample(StandardNormal))?;
        let (mut mean, mut c) = fit(GROWTH, &sample, &covariance, &[1.8, 0.6])?;
        let parameter = parameter(&mut mean, &mut c);
        for (index, x) in points.into_iter().enumerate() {
            let prediction = GROWTH.predict(&[x], &parameter)?;
            let error = prediction.mean - (GROWTH.fx)(&[x], &truth);
            squared_error[index] += error * error;
            reported[index] += prediction.deviation;
            covered[index] += u32::from(error.abs() <= 1.96 * prediction.deviation);
        }
    }
    let count = f64::from(SIMULATIONS);
    for index in 0..points.len() {
        let actual = (squared_error[index] / count).sqrt();
        let reported = reported[index] / count;
        let coverage = f64::from(covered[index]) / count;
        assert!(
            (actual / reported - 1.0).abs() < 0.05,
            "{}: actual {actual}, reported {reported}",
            points[index],
        );
        assert!(
            (0.93..0.97).contains(&coverage),
            "{}: coverage {coverage}",
            points[index],
        );
    }
    Ok(())
}

#[test]
fn not_finite_model_yields_not_finite_prediction() -> Result<(), OefpilError> {
    const ROOT: Model = Model {
        fx: &|x, p| p[0] * x[0].sqrt(),
        dfdx: &[&|x, p| 0.5 * p[0] / x[0].sqrt()],
        dfdp: &[&|x, _| x[0].sqrt()],
        implicit: false,
    };
    let (mut mean, mut c) = ([1.0], [1e-4]);
    let prediction = ROOT.predict(&[-1.0], &parameter(&mut mean, &mut c))?;
    assert!(prediction.mean.is_nan() && prediction.deviation.is_nan());
    Ok(())
}

#[test]
fn invalid_requests_are_errors() {
    const IMPLICIT: Model = Model {
        fx: &|x, p| x[0] - p[0],
        dfdx: &[&|_, _| 1.0],
        dfdp: &[&|_, _| -1.0],
        implicit: true,
    };
    let (mut mean, mut c) = ([1.0], [1e-4]);
    assert!(matches!(
        IMPLICIT.predict(&[1.0], &parameter(&mut mean, &mut c)),
        Err(OefpilError::Prediction(PredictionError::ImplicitModel))
    ));
    let (mut mean, mut c) = ([1.0, 2.0], [1e-4; 4]);
    for x in [&[][..], &[1.0, 2.0]] {
        assert!(matches!(
            LINE.predict(x, &parameter(&mut mean, &mut c)),
            Err(OefpilError::Prediction(
                PredictionError::MismatchingVariable
            ))
        ));
    }
    let (mut mean, mut c) = ([1.0], [1e-4]);
    assert!(matches!(
        LINE.predict(&[1.0], &parameter(&mut mean, &mut c)),
        Err(OefpilError::Parameter(ParameterError::MismatchingMean))
    ));
    let (mut mean, mut c) = ([1.0, 2.0], [1e-4; 3]);
    assert!(matches!(
        LINE.predict(&[1.0], &parameter(&mut mean, &mut c)),
        Err(OefpilError::Parameter(
            ParameterError::MismatchingCovariance
        ))
    ));
}

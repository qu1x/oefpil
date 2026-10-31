//! Fits approximating derivatives via central differences.

#![allow(
    clippy::cast_precision_loss,
    clippy::similar_names,
    clippy::suboptimal_flops
)]

use oefpil::{
    Algorithm, AlgorithmError, Covariance, Model, OefpilError, Parameter, ParameterError, Report,
    Variable, VariableError,
};

/// Fits `model` to `sample` of `samples` with independent errors of `deviation` per variable.
fn fit(
    model: Model,
    sample: &[f64],
    deviation: &[f64],
    p: &[f64],
) -> Result<(Report, Vec<f64>, Vec<f64>), OefpilError> {
    let samples = sample.len() / deviation.len();
    let mut covariance = Covariance::new_diagonal(samples, deviation.len());
    for (variable, deviation) in deviation.iter().enumerate() {
        covariance.set_tile(variable, variable, &vec![deviation * deviation; samples])?;
    }
    let mut mean = p.to_vec();
    let mut parameter_covariance = vec![0f64; p.len().pow(2)];
    let report = Algorithm::default().fit(
        model,
        &mut Variable {
            sample,
            mean: &mut [],
            covariance: &covariance,
        },
        &mut Parameter {
            mean: &mut mean,
            covariance: &mut parameter_covariance,
            ..Parameter::default()
        },
    )?;
    Ok((report, mean, parameter_covariance))
}

/// Noisy samples of `y = a * sin(x / s) + b` along `x` around `s`.
fn sine(a: f64, s: f64, b: f64, shift: f64) -> Vec<f64> {
    let x = (0..16).map(|i| s * (shift + 0.25 * f64::from(i)));
    let noise = |i: usize, scale: f64| scale * 1e-3 * (1.7 * i as f64).sin();
    let xs = x.clone().enumerate().map(|(i, x)| x + noise(i, s));
    let ys = x
        .enumerate()
        .map(|(i, x)| a * (x / s).sin() + b + noise(i + 16, a));
    xs.chain(ys).collect()
}

#[test]
fn any_magnitude_needs_no_tuning() {
    // Scales of the independent variable, the amplitude, and the offset.
    for (s, a, b, shift) in [
        (1.0, 1.0, 0.5, 0.25),
        (1e-9, 1.0, 0.5, 0.25),
        (1e9, 1.0, 0.5, 0.25),
        (1.0, 1e-12, 5e-13, 0.25),
        (1.0, 1e12, 5e11, 0.25),
        (1e-9, 1e12, 5e11, 0.25),
        (1e9, 1e-12, 5e-13, 0.25),
        // Samples at and around zero.
        (1.0, 1.0, 0.5, 0.0),
        (1e-9, 1.0, 0.5, -2.0),
    ] {
        let sample = sine(a, s, b, shift);
        let fx = move |x: &[f64], p: &[f64]| p[0] * (x[0] / s).sin() + p[1];
        let dfdx = move |x: &[f64], p: &[f64]| p[0] * (x[0] / s).cos() / s;
        let dfdp = [
            &(move |x: &[f64], _: &[f64]| (x[0] / s).sin()) as &dyn Fn(&[f64], &[f64]) -> f64,
            &|_, _| 1.0,
        ];
        let analytic = Model {
            fx: &fx,
            dfdx: &[&dfdx],
            dfdp: &dfdp,
            implicit: false,
        };
        let numerical = Model {
            fx: &fx,
            dfdx: &[],
            dfdp: &[],
            implicit: false,
        };
        let (deviation, start) = ([s * 1e-3, a * 1e-3], [a * 1.1, b * 0.9]);
        let (_, p_a, c_a) = fit(analytic, &sample, &deviation, &start).unwrap();
        let (_, p_n, c_n) = fit(numerical, &sample, &deviation, &start).unwrap();
        for (i, (p_a, p_n)) in p_a.iter().zip(&p_n).enumerate() {
            assert!(
                (p_a - p_n).abs() <= 1e-8 * p_a.abs().max(a * 1e-3),
                "{s} {a} {b} parameter {i}: {p_a} != {p_n}"
            );
        }
        let scale = c_a.iter().fold(0f64, |scale, c| scale.max(c.abs()));
        for (c_a, c_n) in c_a.iter().zip(&c_n) {
            assert!(
                (c_a - c_n).abs() <= 1e-6 * scale,
                "{s} {a} {b} covariance: {c_a} != {c_n}"
            );
        }
    }
}

#[test]
fn initial_parameter_of_zero_assumes_magnitude_one() {
    const MODEL: Model = Model {
        fx: &|x, p| p[0] + p[1] * x[0],
        dfdx: &[],
        dfdp: &[],
        implicit: false,
    };
    let sample = [0.0, 1.0, 2.0, 3.0, 4.0, 1.01, 2.99, 5.02, 6.98, 9.0];
    let (_, p, _) = fit(MODEL, &sample, &[1e-2, 1e-2], &[0.0, 0.0]).unwrap();
    assert!(
        (p[0] - 1.0).abs() < 0.05 && (p[1] - 2.0).abs() < 0.05,
        "{p:?}"
    );
}

#[test]
fn implicit_model_takes_variables_from_covariance() {
    const CIRCLE: Model = Model {
        fx: &|x, p| x[0].hypot(x[1]) - p[0],
        dfdx: &[],
        dfdp: &[],
        implicit: true,
    };
    let sample = [
        3.01, 0.0, -2.99, 0.0, 2.12, -2.12, //
        0.0, 3.0, 0.01, -3.01, 2.12, 2.13,
    ];
    let (report, p, _) = fit(CIRCLE, &sample, &[1e-2, 1e-2], &[2.0]).unwrap();
    assert_eq!(5, report.degrees_of_freedom);
    assert!((p[0] - 3.0).abs() < 0.01, "{p:?}");
}

#[test]
fn not_finite_within_a_step_is_numerical_error() {
    const MODEL: Model = Model {
        fx: &|x, p| p[0] * x[0].ln(),
        dfdx: &[],
        dfdp: &[],
        implicit: false,
    };
    // The first sample is below any step.
    let sample = [1e-300, 1.0, 2.0, 3.0, 0.0, 0.0, 0.7, 1.1];
    assert!(matches!(
        fit(MODEL, &sample, &[1e-2, 1e-2], &[1.0]),
        Err(OefpilError::Algorithm(AlgorithmError::NumericalError))
    ));
}

#[test]
fn counts_are_still_validated() {
    const NUMERICAL: Model = Model {
        fx: &|x, p| p[0] * x[0],
        dfdx: &[],
        dfdp: &[],
        implicit: false,
    };
    const ANALYTIC: Model = Model {
        fx: &|x, p| p[0] * x[0],
        dfdx: &[&|_, p| p[0]],
        dfdp: &[&|x, _| x[0]],
        implicit: false,
    };
    // An explicit model requires an independent variable besides the dependent one.
    assert!(matches!(
        fit(NUMERICAL, &[1.0, 2.0], &[1e-2], &[1.0]),
        Err(OefpilError::Variable(VariableError::MismatchingCovariance))
    ));
    // At least one parameter is required.
    assert!(matches!(
        fit(NUMERICAL, &[1.0, 2.0, 1.0, 2.0], &[1e-2, 1e-2], &[]),
        Err(OefpilError::Parameter(ParameterError::EmptyMean))
    ));
    // Non-empty derivatives still have to match.
    assert!(matches!(
        fit(ANALYTIC, &[1.0, 2.0, 1.0, 2.0], &[1e-2, 1e-2], &[1.0, 2.0]),
        Err(OefpilError::Parameter(ParameterError::MismatchingMean))
    ));
    assert!(matches!(
        fit(ANALYTIC, &[1.0; 6], &[1e-2, 1e-2, 1e-2], &[1.0]),
        Err(OefpilError::Variable(VariableError::MismatchingCovariance))
    ));
}

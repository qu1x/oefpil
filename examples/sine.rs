//! Explicit univariate biparametric sinusoidal model.

use oefpil::{
    Algorithm, Covariance, Distribution, Logfile, Model, OefpilError, Parameter, Variable,
};
use rand_chacha::{ChaChaRng, rand_core::SeedableRng};

#[allow(clippy::float_cmp)]
fn main() -> Result<(), OefpilError> {
    // Define explicit model of one independent variable and two parameters.
    const MODEL: Model = Model {
        fx: |x, p| p[1] * (x[0] + p[0]).to_radians().sin(),
        dfdx: &[|x, p| p[1] * (x[0] + p[0]).to_radians().cos().to_radians()],
        dfdp: &[
            |x, p| p[1] * (x[0] + p[0]).to_radians().cos().to_radians(),
            |x, p| (x[0] + p[0]).to_radians().sin(),
        ],
        implicit: false,
    };
    // Define true values.
    let p_mean = [8f64, 10.0];
    let x_mean = [-60.0, -40.0, -20.0, -10.0, 10.0, 20.0, 40.0, 60.0];
    let y_mean = x_mean.map(|x_mean| (MODEL.fx)(&[x_mean], &p_mean));
    // Define covariance matrix with constant Pearson `correlation` `coefficient` (PCC).
    let x_deviation = x_mean.map(|x| x.abs() * 2e-2);
    let y_deviation = y_mean.map(|y| y.abs() * 8e-2);
    let coefficient = 0.5;
    let correlation = x_deviation
        .iter()
        .zip(&y_deviation)
        .map(|(x, y)| x * y * coefficient)
        .collect::<Vec<f64>>();
    let covariance = &Covariance::new_diagonals(8, 2)
        .with_tile(0, 0, &x_deviation.map(|x| x * x))?
        .with_tile(1, 1, &y_deviation.map(|y| y * y))?
        .with_tile(0, 1, &correlation)?
        .with_decomposition()?;
    // Seed of 256 bits for deterministic high-quality randomness.
    let seed = b"85a4d500a53d81b65255dc718d7df63a3cd3d6bb15e6f724af04a23d74cc02d3";
    // Cryptographically secure pseudorandom number generator (CSPRNG).
    let mut csprng = ChaChaRng::from_seed(const_hex::decode_to_array(seed).unwrap());
    // Define correlated normal distributions of variables.
    let mut distribution = Distribution {
        sample: &mut vec![0f64; x_mean.len() + y_mean.len()],
        mean: &mut x_mean.iter().chain(&y_mean).copied().collect::<Vec<f64>>(),
        covariance,
    };
    // Simulate measurement by sampling observations from distributions of variables.
    distribution.sample(&mut csprng)?;
    // Pretend mediocre estimates of variable means (empty slices default to samples).
    let mut variable = Variable {
        sample: distribution.sample,
        mean: &mut [],
        covariance,
    };
    // Pretend mediocre estimates of parameter means and assign output slices.
    let p_mean_initial = [9.5f64, 12.0];
    let mut parameter = Parameter {
        mean: &mut p_mean_initial.clone(),
        deviation: &mut [0f64; 2],
        covariance: &mut [0f64; 4],
        correlation: &mut [0f64; 4],
    };
    // Perform fitting with default settings.
    let oefpil = Algorithm::default();
    let report = oefpil.fit(MODEL, &mut variable, &mut parameter, Logfile::default())?;
    // Assert expected results.
    let round = |value: f64| (value * 1e2).round() / 1e2;
    assert_eq!(7, report.iterations);
    assert_eq!(0.73, report.chi_squared_p_value().map(round)?);
    // `parameter.mean` deviates from `p_mean` in agreement with `parameter.deviation`.
    assert_eq!(8.16, round(parameter.mean[0]));
    assert_eq!(0.19, round(parameter.deviation[0]));
    assert_eq!(9.95, round(parameter.mean[1]));
    assert_eq!(0.28, round(parameter.deviation[1]));
    assert_eq!(0.30, round(parameter.correlation[1]));
    // Confirm `parameter.deviation` via Monte Carlo simulation.
    let samples = 10_000;
    let mut p_sample = vec![0f64; samples * p_mean.len()];
    let mut p_correlation = vec![0f64; samples];
    for sample in 0..samples {
        distribution.sample(&mut csprng)?;
        let mut variable = Variable {
            sample: distribution.sample,
            mean: &mut [],
            covariance,
        };
        let mut parameter = Parameter {
            mean: &mut p_mean_initial.clone(),
            covariance: &mut [0f64; 4],
            deviation: &mut [0f64; 2],
            correlation: &mut [0f64; 4],
        };
        oefpil.fit(MODEL, &mut variable, &mut parameter, Logfile::default())?;
        p_sample[sample] = parameter.mean[0];
        p_sample[sample + samples] = parameter.mean[1];
        p_correlation[sample] = parameter.correlation[1];
    }
    let mean_deviation = |s: &[f64]| {
        let n = f64::from(u32::try_from(s.len()).unwrap());
        let m = s.iter().sum::<f64>() / n;
        let v = s.iter().map(|s| s - m).map(|d| d * d).sum::<f64>() / (n - 1.0);
        [m, v.sqrt()]
    };
    let [p_mean_0, p_deviation_0] = mean_deviation(&p_sample[..samples]);
    let [p_mean_1, p_deviation_1] = mean_deviation(&p_sample[samples..]);
    let [p_correlation_mean, p_correlation_deviation] = mean_deviation(&p_correlation);
    // Simulation agrees well with `p_mean` and `parameter.deviation`.
    assert_eq!(8.00, round(p_mean_0));
    assert_eq!(0.19, round(p_deviation_0));
    assert_eq!(10.0, round(p_mean_1));
    assert_eq!(0.28, round(p_deviation_1));
    assert_eq!(0.30, round(p_correlation_mean));
    assert_eq!(0.02, round(p_correlation_deviation));
    Ok(())
}

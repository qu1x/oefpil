//! Calibration of a barometric altimeter against GPS along a climb.
//!
//! The barometric formula relates pressure to altitude via the sea level pressure, which is what a
//! pilot sets on an altimeter. Both variables are noisy but the error of GPS drifts slowly and is
//! hence largely shared by all samples of a flight. More samples do not average out what they share,
//! which is why ignoring this correlation yields about the same estimate but claims a precision the
//! flight cannot deliver. This example confirms both via Monte Carlo simulation.

#![allow(clippy::suboptimal_flops)]

use oefpil::{Algorithm, Covariance, Distribution, Model, OefpilError, Parameter, Variable};
use rand::Rng;
use rand_chacha::{ChaChaRng, rand_core::SeedableRng};
use rand_distr::StandardNormal;

/// Exponent of the barometric formula in the international standard atmosphere.
const EXPONENT: f64 = 0.190_263;

/// Barometric formula with sea level pressure `p[0]` in hPa and scale height `p[1]` in m.
const MODEL: Model = Model {
    fx: &|x, p| p[1] * (1.0 - (x[0] / p[0]).powf(EXPONENT)),
    dfdx: &[&|x, p| -p[1] * EXPONENT * (x[0] / p[0]).powf(EXPONENT) / x[0]],
    dfdp: &[
        &|x, p| p[1] * EXPONENT * (x[0] / p[0]).powf(EXPONENT) / p[0],
        &|x, p| 1.0 - (x[0] / p[0]).powf(EXPONENT),
    ],
    implicit: false,
};

/// Fits [`MODEL`] returning the parameter mean and deviation vector.
fn fit(sample: &[f64], covariance: &Covariance) -> Result<([f64; 2], [f64; 2]), OefpilError> {
    // Mediocre estimates of parameter means as they would be guessed.
    let mut mean = [1000.0, 40_000.0];
    let mut deviation = [0f64; 2];
    let mut variable = Variable {
        sample,
        mean: &mut [],
        covariance,
    };
    let mut parameter = Parameter {
        mean: &mut mean,
        deviation: &mut deviation,
        covariance: &mut [0f64; 4],
        correlation: &mut [0f64; 4],
    };
    Algorithm::default().fit(MODEL, &mut variable, &mut parameter)?;
    Ok((mean, deviation))
}

#[allow(clippy::cast_precision_loss)]
fn main() -> Result<(), OefpilError> {
    // Define true values of a climb to 3 km.
    const SAMPLES: usize = 40;
    let p_mean = [1013.25, 44_330.8f64];
    let y_mean: [f64; SAMPLES] = core::array::from_fn(|i| 3000.0 * i as f64 / 39.0);
    let x_mean = y_mean.map(|y| p_mean[0] * (1.0 - y / p_mean[1]).powf(1.0 / EXPONENT));
    // Define deviations of pressure in hPa and of altitude in m where 80 % of the variance of the
    // latter is shared by all samples.
    let (x_deviation, y_deviation, coefficient) = (0.15f64, 8.0f64, 0.8);
    let x_variance = [x_deviation * x_deviation; SAMPLES];
    let y_variance = [y_deviation * y_deviation; SAMPLES];
    let mut y_block = vec![coefficient * y_deviation * y_deviation; SAMPLES * SAMPLES];
    for sample in 0..SAMPLES {
        y_block[sample * SAMPLES + sample] = y_variance[sample];
    }
    // Define covariance matrix with its correlation within the altitude variable.
    let covariance = &Covariance::new_block_diagonal(SAMPLES, 2)
        .with_tile(0, 0, &x_variance)?
        .with_tile(1, 1, &y_block)?
        .with_decomposition()?;
    // Define covariance matrix which ignores this correlation.
    let ignorant = &Covariance::new_diagonal(SAMPLES, 2)
        .with_tile(0, 0, &x_variance)?
        .with_tile(1, 1, &y_variance)?;
    // Seed of 256 bits for deterministic high-quality randomness.
    let seed = b"2f6a1c0e9b7d4853a1c6e0f27b94d3581e7ac24f60d9b3175c8e2a4d06f1b937";
    // Cryptographically secure pseudorandom number generator (CSPRNG).
    let mut csprng = ChaChaRng::from_seed(const_hex::decode_to_array(seed).unwrap());
    // Define correlated normal distributions of variables.
    let mut distribution = Distribution {
        sample: &mut [0f64; 2 * SAMPLES],
        mean: &x_mean.iter().chain(&y_mean).copied().collect::<Vec<f64>>(),
        covariance,
    };
    // Fit many simulated flights both ways.
    let flights = 2000;
    let mut squared_error = [0f64; 2];
    let mut reported = [0f64; 2];
    let mut covered = [0u32; 2];
    for _ in 0..flights {
        distribution.sample(|| csprng.sample(StandardNormal))?;
        for (way, covariance) in [covariance, ignorant].into_iter().enumerate() {
            let (mean, deviation) = fit(distribution.sample, covariance)?;
            let error = mean[0] - p_mean[0];
            squared_error[way] += error * error;
            reported[way] += deviation[0];
            covered[way] += u32::from(error.abs() <= 1.96 * deviation[0]);
        }
    }
    let flights = f64::from(flights);
    let actual = squared_error.map(|squared_error| (squared_error / flights).sqrt());
    let reported = reported.map(|reported| reported / flights);
    let coverage = covered.map(|covered| f64::from(covered) / flights);
    // A difference in sea level pressure shifts the indicated altitude at sea level as follows.
    let meters_per_hectopascal = p_mean[1] * EXPONENT / p_mean[0];
    println!("sea level pressure of {flights} simulated flights:");
    for (way, name) in ["modeled", "ignored"].into_iter().enumerate() {
        println!(
            "  correlation {name}: actual error {:.2} hPa, reported error {:.2} hPa ({:.1} m), \
             95 % interval holds true value in {:.1} %",
            actual[way],
            reported[way],
            reported[way] * meters_per_hectopascal,
            coverage[way] * 100.0,
        );
    }
    // Both ways yield about the same estimate as their actual errors agree.
    assert!((actual[0] - actual[1]).abs() < 0.01);
    // Modeling the correlation reports the actual error and its interval holds what it claims.
    assert!((reported[0] - actual[0]).abs() < 0.03);
    assert!((0.93..0.97).contains(&coverage[0]));
    // Ignoring it claims about three times the precision and its interval is wrong half the time.
    assert!(actual[1] / reported[1] > 2.5);
    assert!((0.45..0.55).contains(&coverage[1]));
    Ok(())
}

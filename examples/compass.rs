//! Calibration of the magnetometer of a drone.
//!
//! Turning a drone traces the magnetic field of the earth on a sphere, which the magnetic parts of
//! the drone shift by the hard iron offset and distort to an ellipsoid by the soft iron scales. A
//! point on an ellipsoid has no dependent variable, so this implicit model treats all three axes
//! alike. Their noise is correlated as interference of the drone affects them together. As the
//! drone is only tilted a little during calibration, the vertical axis is less well determined,
//! which shows in its parameter deviation. This example confirms them via Monte Carlo simulation.

#![allow(clippy::suboptimal_flops)]

use oefpil::{Algorithm, Covariance, Distribution, Model, OefpilError, Parameter, Variable};
use rand::Rng;
use rand_chacha::{ChaChaRng, rand_core::SeedableRng};
use rand_distr::StandardNormal;

/// Axis-aligned ellipsoid with center `p[0..3]` and radii `p[3..6]`.
const MODEL: Model = Model {
    fx: &|x, p| {
        let axis = |i: usize| ((x[i] - p[i]) / p[3 + i]).powi(2);
        axis(0) + axis(1) + axis(2) - 1.0
    },
    dfdx: &[
        &|x, p| 2.0 * (x[0] - p[0]) / p[3].powi(2),
        &|x, p| 2.0 * (x[1] - p[1]) / p[4].powi(2),
        &|x, p| 2.0 * (x[2] - p[2]) / p[5].powi(2),
    ],
    dfdp: &[
        &|x, p| -2.0 * (x[0] - p[0]) / p[3].powi(2),
        &|x, p| -2.0 * (x[1] - p[1]) / p[4].powi(2),
        &|x, p| -2.0 * (x[2] - p[2]) / p[5].powi(2),
        &|x, p| -2.0 * (x[0] - p[0]).powi(2) / p[3].powi(3),
        &|x, p| -2.0 * (x[1] - p[1]).powi(2) / p[4].powi(3),
        &|x, p| -2.0 * (x[2] - p[2]).powi(2) / p[5].powi(3),
    ],
    implicit: true,
};

/// Fits [`MODEL`] returning the parameter mean and deviation vector.
fn fit(sample: &[f64], covariance: &Covariance) -> Result<([f64; 6], [f64; 6]), OefpilError> {
    // Mediocre estimates of parameter means: No offset and the field strength of the earth.
    let mut mean = [0.0, 0.0, 0.0, 50.0, 50.0, 50.0];
    let mut deviation = [0f64; 6];
    let mut variable = Variable {
        sample,
        mean: &mut [],
        covariance,
    };
    let mut parameter = Parameter {
        mean: &mut mean,
        deviation: &mut deviation,
        covariance: &mut [0f64; 36],
        correlation: &mut [0f64; 36],
    };
    Algorithm::default().fit(MODEL, &mut variable, &mut parameter)?;
    Ok((mean, deviation))
}

/// Heading in degrees of a horizontal field after removing `center` and dividing by `radius`.
fn heading(x: f64, y: f64, center: [f64; 2], radius: [f64; 2]) -> f64 {
    ((y - center[1]) / radius[1])
        .atan2((x - center[0]) / radius[0])
        .to_degrees()
}

fn main() -> Result<(), OefpilError> {
    // Define true values in uT: The offset and the radii of the ellipsoid.
    const SAMPLES: usize = 60;
    let p_mean = [12.0f64, -7.0, 20.0, 48.0, 52.0, 45.0];
    // Turn the drone evenly while tilting it by up to 35 degrees.
    let mut v_mean = [0f64; 3 * SAMPLES];
    for sample in 0..SAMPLES {
        let turn = f64::from(u32::try_from(sample).unwrap());
        let (yaw, tilt) = (turn * 2.399_963, 35f64.to_radians() * (turn * 0.7).sin());
        v_mean[sample] = p_mean[0] + p_mean[3] * tilt.cos() * yaw.cos();
        v_mean[SAMPLES + sample] = p_mean[1] + p_mean[4] * tilt.cos() * yaw.sin();
        v_mean[2 * SAMPLES + sample] = p_mean[2] + p_mean[5] * tilt.sin();
    }
    // Define covariance matrix with noise of 0.5 uT correlated among the axes by 0.6.
    let (variance, covariance_among) = ([0.25; SAMPLES], [0.6 * 0.25; SAMPLES]);
    let covariance = &Covariance::new_diagonals(SAMPLES, 3)
        .with_tile(0, 0, &variance)?
        .with_tile(1, 1, &variance)?
        .with_tile(2, 2, &variance)?
        .with_tile(0, 1, &covariance_among)?
        .with_tile(0, 2, &covariance_among)?
        .with_tile(1, 2, &covariance_among)?
        .with_decomposition()?;
    // Seed of 256 bits for deterministic high-quality randomness.
    let seed = b"c41d7e2096ab35f8d07e16b93c5a48f2e6071bd9a2c43f58d1e7096bb3a2c5f4";
    // Cryptographically secure pseudorandom number generator (CSPRNG).
    let mut csprng = ChaChaRng::from_seed(const_hex::decode_to_array(seed).unwrap());
    // Define correlated normal distributions of variables.
    let mut distribution = Distribution {
        sample: &mut [0f64; 3 * SAMPLES],
        mean: &v_mean,
        covariance,
    };
    // Calibrate on a single simulated flight.
    distribution.sample(|| csprng.sample(StandardNormal))?;
    let (mean, deviation) = fit(distribution.sample, covariance)?;
    println!("calibration of a single flight in uT:");
    for (index, name) in ["x", "y", "z"].into_iter().enumerate() {
        println!(
            "  {name}: offset {:6.2} +- {:.2} (true {:5.1}), radius {:6.2} +- {:.2} (true {:5.1})",
            mean[index],
            deviation[index],
            p_mean[index],
            mean[3 + index],
            deviation[3 + index],
            p_mean[3 + index],
        );
    }
    // Correct the reading of a level drone facing north-east.
    let (x, y) = (
        p_mean[0] + p_mean[3] * 45f64.to_radians().cos(),
        p_mean[1] + p_mean[4] * 45f64.to_radians().sin(),
    );
    let raw = heading(x, y, [0.0, 0.0], [1.0, 1.0]);
    let calibrated = heading(x, y, [mean[0], mean[1]], [mean[3], mean[4]]);
    println!("heading of north-east (45 deg): raw {raw:.1} deg, calibrated {calibrated:.1} deg");
    // The calibration removes an error of more than ten degrees.
    assert!((raw - 45.0).abs() > 10.0);
    assert!((calibrated - 45.0).abs() < 0.5);
    // Confirm `deviation` via Monte Carlo simulation.
    let flights = 2000;
    let mut squared_error = [0f64; 6];
    let mut reported = [0f64; 6];
    let mut covered = [0u32; 6];
    for _ in 0..flights {
        distribution.sample(|| csprng.sample(StandardNormal))?;
        let (mean, deviation) = fit(distribution.sample, covariance)?;
        for index in 0..6 {
            let error = mean[index] - p_mean[index];
            squared_error[index] += error * error;
            reported[index] += deviation[index];
            covered[index] += u32::from(error.abs() <= 1.96 * deviation[index]);
        }
    }
    let flights = f64::from(flights);
    let actual = squared_error.map(|squared_error| (squared_error / flights).sqrt());
    let reported = reported.map(|reported| reported / flights);
    let coverage = covered.map(|covered| f64::from(covered) / flights);
    println!("{flights} simulated flights in uT:");
    for (index, name) in [
        "offset x", "offset y", "offset z", "radius x", "radius y", "radius z",
    ]
    .into_iter()
    .enumerate()
    {
        println!(
            "  {name}: actual error {:.3}, reported error {:.3}, \
             95 % interval holds true value in {:.1} %",
            actual[index],
            reported[index],
            coverage[index] * 100.0,
        );
        // Reports the actual error and its interval holds what it claims.
        assert!((reported[index] / actual[index] - 1.0).abs() < 0.1);
        assert!((0.92..0.98).contains(&coverage[index]));
    }
    // The vertical axis is less well determined as the drone is only tilted a little.
    assert!(reported[2] > reported[0] && reported[2] > reported[1]);
    assert!(reported[5] > reported[3] && reported[5] > reported[4]);
    Ok(())
}

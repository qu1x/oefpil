//! Matter density of the universe from real Type Ia supernovae of the Pantheon+ data release.
//!
//! Supernovae of this type are of about the same brightness, so how faint they appear measures how
//! far away they are and hence how the universe has expanded. Their systematic errors (e.g., of the
//! calibration of a telescope) are shared among many of them, which the release provides as a full
//! covariance matrix. This example fits a flat universe to the supernovae beyond a redshift of 0.01
//! and reproduces the published matter density of `0.334 +- 0.018` of [Brout et al.
//! (2022)](https://doi.org/10.3847/1538-4357/ac8e04). Ignoring the correlation shifts this estimate
//! and claims too much precision.
//!
//! Run it via:
//!
//! ```text
//! cargo run --features faer --example supernovae
//! ```
//!
//! Feature `faer` is optional and about three times faster here.
//!
//! The data is not part of this crate as it is large and its release does not state a license. The
//! first run downloads it from <https://github.com/PantheonPlusSH0ES/DataRelease> via `curl` into
//! `target/pantheon` or the directory given as first argument.
//!
//! The redshifts are treated as exact as their errors already are part of the covariance matrix.
//! With their variances being tiny, this reduces to generalized least squares.

#![allow(clippy::suboptimal_flops)]

use oefpil::{Algorithm, Covariance, Model, Parameter, Variable};
use std::{error::Error, f64::consts::LN_10, fs, path::Path, process::Command};

/// Location of the data.
const RELEASE: &str = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/\
                       Pantheon%2B_Data/4_DISTANCES_AND_COVAR";

/// Positive nodes of the 16-point Gauss-Legendre quadrature on `[-1, 1]`.
const NODES: [f64; 8] = [
    0.095_012_509_837_637_4,
    0.281_603_550_779_258_9,
    0.458_016_777_657_227_4,
    0.617_876_244_402_643_8,
    0.755_404_408_355_003,
    0.865_631_202_387_831_8,
    0.944_575_023_073_232_6,
    0.989_400_934_991_649_9,
];

/// Weights of [`NODES`].
const WEIGHTS: [f64; 8] = [
    0.189_450_610_455_068_5,
    0.182_603_415_044_923_6,
    0.169_156_519_395_002_5,
    0.149_595_988_816_576_7,
    0.124_628_971_255_533_9,
    0.095_158_511_682_492_8,
    0.062_253_523_938_647_9,
    0.027_152_459_411_754_1,
];

/// Integrates `f` from zero to `z`.
fn integrate(z: f64, f: impl Fn(f64) -> f64) -> f64 {
    let half = 0.5 * z;
    let sum = NODES
        .iter()
        .zip(WEIGHTS)
        .map(|(node, weight)| weight * (f(half * (1.0 + node)) + f(half * (1.0 - node))))
        .sum::<f64>();
    half * sum
}

/// Expansion rate of a flat universe of matter density `density` at redshift `z` relative to now.
fn expansion(z: f64, density: f64) -> f64 {
    (density * (1.0 + z).powi(3) + 1.0 - density).sqrt()
}

/// Comoving distance in units of the Hubble distance.
fn distance(z: f64, density: f64) -> f64 {
    integrate(z, |z| 1.0 / expansion(z, density))
}

/// Derivative of [`distance()`] in `density`.
fn distance_in_density(z: f64, density: f64) -> f64 {
    integrate(z, |z| {
        -0.5 * ((1.0 + z).powi(3) - 1.0) / expansion(z, density).powi(3)
    })
}

/// Magnitudes per unit of the natural logarithm of a distance.
const SCALE: f64 = 5.0 / LN_10;

/// Apparent magnitude of redshift `x[0]` of the expansion and of redshift `x[1]` of the observer
/// with matter density `p[0]` and magnitude offset `p[1]`, which comprises the unknown Hubble
/// constant and the unknown absolute magnitude.
const MODEL: Model = Model {
    fx: &|x, p| p[1] + SCALE * ((1.0 + x[1]) * distance(x[0], p[0])).ln(),
    dfdx: &[
        &|x, p| SCALE / expansion(x[0], p[0]) / distance(x[0], p[0]),
        &|x, _| SCALE / (1.0 + x[1]),
    ],
    dfdp: &[
        &|x, p| SCALE * distance_in_density(x[0], p[0]) / distance(x[0], p[0]),
        &|_, _| 1.0,
    ],
    implicit: false,
};

/// Fits [`MODEL`] returning the matter density and its deviation.
fn fit(sample: &[f64], covariance: &Covariance) -> Result<(f64, f64), Box<dyn Error>> {
    // Mediocre estimates of parameter means.
    let mut mean = [0.3, 24.0];
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
    Ok((mean[0], deviation[0]))
}

/// Runs `command` and requires it to succeed.
fn run(command: &mut Command) -> Result<(), Box<dyn Error>> {
    let program = command.get_program().to_string_lossy().into_owned();
    let status = command
        .status()
        .map_err(|error| format!("cannot run `{program}`: {error}"))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("`{program}` failed with {status}").into())
    }
}

/// Downloads `url` to `path` via `curl`, keeping `path` absent unless complete.
fn download(url: &str, path: &Path) -> Result<(), Box<dyn Error>> {
    println!("Downloading {url}");
    let part = path.with_extension("part");
    run(Command::new("curl")
        .args(["-fsSL", "-o"])
        .arg(&part)
        .arg(url))?;
    fs::rename(part, path)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let directory = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "target/pantheon".into());
    let directory = Path::new(&directory);
    let (table, matrix) = (
        directory.join("Pantheon+SH0ES.dat"),
        directory.join("Pantheon+SH0ES_STAT+SYS.cov"),
    );
    fs::create_dir_all(directory)?;
    if !table.exists() {
        download(&format!("{RELEASE}/Pantheon%2BSH0ES.dat"), &table)?;
    }
    if !matrix.exists() {
        download(
            &format!("{RELEASE}/Pantheon%2BSH0ES_STAT%2BSYS.cov"),
            &matrix,
        )?;
    }
    // Read the table of observations.
    let table = fs::read_to_string(table)?;
    let mut lines = table.lines();
    let header = lines.next().ok_or("empty table")?;
    let header = header.split_whitespace().collect::<Vec<&str>>();
    let column = |name: &str| {
        header
            .iter()
            .position(|column| *column == name)
            .ok_or_else(|| format!("missing column {name}"))
    };
    let columns = [column("zHD")?, column("zHEL")?, column("m_b_corr")?];
    let mut observations = Vec::new();
    for line in lines {
        let fields = line.split_whitespace().collect::<Vec<&str>>();
        let mut observation = [0f64; 3];
        for (observation, column) in observation.iter_mut().zip(columns) {
            *observation = fields[column].parse()?;
        }
        observations.push(observation);
    }
    // Read the full covariance matrix which is preceded by its order.
    let matrix = fs::read_to_string(matrix)?;
    let mut fields = matrix.split_ascii_whitespace();
    let order = fields.next().ok_or("empty matrix")?.parse::<usize>()?;
    assert_eq!(order, observations.len());
    let matrix = fields.map(str::parse).collect::<Result<Vec<f64>, _>>()?;
    assert_eq!(order * order, matrix.len());
    // Keep the supernovae beyond the redshift where local motions dominate.
    let keep = (0..order)
        .filter(|&index| observations[index][0] > 0.01)
        .collect::<Vec<usize>>();
    let samples = keep.len();
    let mut sample = vec![0f64; 3 * samples];
    let mut block = vec![0f64; samples * samples];
    for (row, &index) in keep.iter().enumerate() {
        for variable in 0..3 {
            sample[variable * samples + row] = observations[index][variable];
        }
        for (column, &other) in keep.iter().enumerate() {
            // Symmetric up to the rounding of the file.
            block[row * samples + column] =
                f64::midpoint(matrix[index * order + other], matrix[other * order + index]);
        }
    }
    let variance = (0..samples)
        .map(|row| block[row * samples + row])
        .collect::<Vec<f64>>();
    // Tiny positive variances of the redshifts which are treated as exact.
    let exact = vec![1e-14; samples];
    // Define covariance matrix with its correlation within the magnitude variable.
    let covariance = Covariance::new_block_diagonal(samples, 3)
        .with_tile(0, 0, &exact)?
        .with_tile(1, 1, &exact)?
        .with_tile(2, 2, &block)?;
    // Define covariance matrix which ignores this correlation.
    let ignorant = Covariance::new_diagonal(samples, 3)
        .with_tile(0, 0, &exact)?
        .with_tile(1, 1, &exact)?
        .with_tile(2, 2, &variance)?;
    let (density, deviation) = fit(&sample, &covariance)?;
    let (ignorant_density, ignorant_deviation) = fit(&sample, &ignorant)?;
    println!("matter density of {samples} supernovae:");
    println!("  published:           0.334 +- 0.018");
    println!("  correlation modeled: {density:.3} +- {deviation:.3}");
    println!("  correlation ignored: {ignorant_density:.3} +- {ignorant_deviation:.3}");
    // Reproduces the published matter density and its deviation.
    assert_eq!(1590, samples);
    assert!((density - 0.334).abs() < 0.01);
    assert!((deviation - 0.018).abs() < 0.001);
    // Ignoring the correlation shifts the estimate and claims too much precision.
    assert!((ignorant_density - density).abs() > 0.5 * deviation);
    assert!(ignorant_deviation < 0.8 * deviation);
    Ok(())
}

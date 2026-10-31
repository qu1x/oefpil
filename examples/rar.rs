//! Acceleration scale of the radial acceleration relation from real galaxies of the SPARC database.
//!
//! The takeaway: All points of a galaxy share errors in both variables. A fit which ignores this
//! reports an error about four times smaller than its actual one. This fit models these errors,
//! reports its actual error, and agrees with a Bayesian sampler of the same model within about a
//! third of that error.
//!
//! # The problem
//!
//! Stars and gas alone cannot explain how fast galaxies rotate, but the observed acceleration still
//! follows the acceleration expected from them along a tight curve with a single parameter, the
//! acceleration scale `g` in units of `1e-10 m/s^2`. This relation is at the heart of the debate
//! about dark matter and modified gravity.
//!
//! A galaxy has one distance, one inclination, and one mass-to-light ratio, so their errors move
//! all of its points together. This example assigns the published errors of the former two to the
//! observed acceleration and an error of `25 %` of the latter to the expected acceleration, which
//! also has an independent error of `0.01` dex. Orthogonal distance regression handles errors in
//! both variables but cannot express this correlation.
//!
//! # The result
//!
//! ```text
//!                      real data       300 simulations with g = 1.2
//!                      g               reported error  actual error  95 % interval holds g
//! correlated errors    1.078 +- 0.032  0.035           0.035         95 %
//! independent errors   1.117 +- 0.015  0.016           0.068         34 %
//! ```
//!
//! Both fits find the true value on average, but only the fit with correlated errors knows how far
//! off it can be. The difference of `0.040` between both fits on the real data is therefore no
//! surprise, as it is typically `0.057` in the simulations. These assume the errors stated above
//! and are only as good as these.
//!
//! Two further results support this:
//!
//!   * A NUTS sampler in the `bayes` directory of the repository, which is not part of this crate,
//!     yields `1.089 +- 0.032` on the same model. On 40 simulations, this fit was on average
//!     `0.002` below the truth and the sampler `0.010` above it.
//!   * Scaling all mass-to-light ratios by 20 % moves the fit with independent errors by about
//!     `0.26` but the fit with correlated errors only by about `0.03`, presumably as the
//!     mass-to-light ratio of each galaxy absorbs most of it.
//!
//! # What this does not show
//!
//!   * Whether the curve is right. It is assumed and an intrinsic scatter of `0.05` dex is chosen
//!     such that chi-squared matches its degrees of freedom. So this example neither tests the
//!     curve nor decides between dark matter and modified gravity.
//!   * That covariance matrices suffice in general. [Desmond (2023)] jointly infers all parameters
//!     of each galaxy with a closely related curve and finds `1.070 +- 0.032` without external
//!     field effect, while the `1.19` in his abstract averages over models which include that
//!     effect. With his curve and his scatter of `0.035` dex, this fit yields `1.110 +- 0.030` with
//!     a chi-squared per degree of freedom of `1.40`. The above sampler with his main ingredients
//!     yields `1.077 +- 0.030`. Most of the gap to this fit is the exactly log-normal mass-to-light
//!     ratio, which lowers the sampled value by about `0.04` and which this algorithm cannot
//!     express.
//!   * Why [McGaugh et al. (2016)] report `1.20 +- 0.02 (random) +- 0.24 (systematic)`. Fitting
//!     their published table with their independent errors yields `1.11`, so correlation is not
//!     the cause. [Zhou et al. (2017)] report `1.02` from the same kind of fit. All of these values
//!     agree within the systematic uncertainty, but the difference remains unexplained.
//!
//! # How to run it
//!
//! Run the fits, optionally followed by the simulations which take minutes:
//!
//! ```text
//! cargo run --features faer --example rar
//! cargo run --features faer --example rar -- target/sparc 300
//! ```
//!
//! Feature `faer` is optional and about three times faster here.
//!
//! The data is not part of this crate as its release does not state a license. If you use it, cite
//! [Lelli et al. (2016)] as requested by <https://astroweb.case.edu/SPARC>. The first run downloads
//! it via `curl` and `unzip` into `target/sparc` or the directory given as first argument.
//!
//! [McGaugh et al. (2016)]: https://doi.org/10.1103/PhysRevLett.117.201101
//! [Lelli et al. (2016)]: https://doi.org/10.3847/0004-6256/152/6/157
//! [Zhou et al. (2017)]: https://arxiv.org/abs/1707.00417
//! [Desmond (2023)]: https://arxiv.org/abs/2303.11314

#![allow(
    clippy::cast_precision_loss,
    clippy::suboptimal_flops,
    clippy::too_many_lines
)]

use oefpil::{Algorithm, Covariance, Model, Parameter, Variable};
use rand::Rng;
use rand_chacha::{ChaChaRng, rand_core::SeedableRng};
use rand_distr::StandardNormal;
use std::{
    collections::HashMap, error::Error, f64::consts::LN_10, fs, path::Path, process::Command,
    thread,
};

/// Mass-to-light ratios of disks and bulges in solar units adopted by the reference.
const MASS_TO_LIGHT: [f64; 2] = [0.5, 0.7];
/// Relative uncertainty of the mass-to-light ratio of a galaxy.
const MASS_TO_LIGHT_ERROR: f64 = 0.25;
/// Independent error of the logarithm of the expected acceleration in dex.
const EXPECTED_ERROR: f64 = 0.01;
/// Intrinsic scatter of the logarithm of the observed acceleration in dex.
const SCATTER: f64 = 0.05;
/// Square kilometers per square second per kiloparsec in meters per square second.
const UNIT: f64 = 1e6 / 3.085_677_581e19;
/// Location of the data.
const RELEASE: &str = "https://astroweb.case.edu/SPARC";

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

/// Downloads what is missing of the data to `directory`.
fn fetch(directory: &Path) -> Result<(), Box<dyn Error>> {
    fs::create_dir_all(directory)?;
    for file in ["SPARC_Lelli2016c.mrt", "RAR.mrt"] {
        let path = directory.join(file);
        if !path.exists() {
            download(&format!("{RELEASE}/{file}"), &path)?;
        }
    }
    let curves = directory.join("Rotmod_LTG");
    if !curves.exists() {
        let archive = directory.join("Rotmod_LTG.zip");
        download(&format!("{RELEASE}/Rotmod_LTG.zip"), &archive)?;
        let part = directory.join("Rotmod_LTG.unzip");
        run(Command::new("unzip")
            .arg("-qo")
            .arg(&archive)
            .arg("-d")
            .arg(&part))?;
        fs::rename(part, curves)?;
        fs::remove_file(archive)?;
    }
    Ok(())
}

/// Square root of the expected acceleration of logarithm `x` over the acceleration scale `scale`.
fn ratio(x: f64, scale: f64) -> f64 {
    (10f64.powf(x + 10.0) / scale).sqrt()
}

/// Derivative of the logarithm of `1 - exp(-ratio)` in the logarithm of the expected acceleration.
fn slope(ratio: f64) -> f64 {
    0.5 * ratio * (-ratio).exp() / -(-ratio).exp_m1()
}

/// Logarithm of the observed acceleration as a function of the logarithm `x[0]` of the expected
/// acceleration with acceleration scale `p[0]` in units of `1e-10 m/s^2`.
const MODEL: Model = Model {
    fx: &|x, p| x[0] - (-(-ratio(x[0], p[0])).exp_m1()).log10(),
    dfdx: &[&|x, p| 1.0 - slope(ratio(x[0], p[0]))],
    dfdp: &[&|x, p| slope(ratio(x[0], p[0])) / (p[0] * LN_10)],
    implicit: false,
};

/// Terms of `g_obs = g_bar / 2 + sqrt(g_bar^2 / 4 + g_bar * a)` in logarithms with `a = p[0]`: The
/// value, its derivative in the logarithm of the expected acceleration, and its derivative in `a`.
fn simple(x: f64, a: f64) -> [f64; 3] {
    let ratio = a / 10f64.powf(x + 10.0);
    let root = (0.25 + ratio).sqrt();
    let derivative = 1.0 / (2.0 * root * (0.5 + root) * LN_10);
    [
        x + (0.5 + root).log10(),
        1.0 - LN_10 * ratio * derivative,
        derivative * ratio / a,
    ]
}

/// The curve of [`simple()`] used by [Desmond (2023)].
///
/// [Desmond (2023)]: https://arxiv.org/abs/2303.11314
const SIMPLE: Model = Model {
    fx: &|x, p| simple(x[0], p[0])[0],
    dfdx: &[&|x, p| simple(x[0], p[0])[1]],
    dfdp: &[&|x, p| simple(x[0], p[0])[2]],
    implicit: false,
};

/// What a fit yields.
struct Fit {
    /// Acceleration scale.
    scale: f64,
    /// Deviation of [`Self::scale`].
    deviation: f64,
    /// Chi-squared per degree of freedom.
    reduced: f64,
}

/// Fits [`MODEL`].
fn fit(sample: &[f64], covariance: &Covariance) -> Result<Fit, Box<dyn Error>> {
    fit_model(MODEL, sample, covariance)
}

/// Fits `model`.
fn fit_model(model: Model, sample: &[f64], covariance: &Covariance) -> Result<Fit, Box<dyn Error>> {
    let (mut mean, mut deviation) = ([1.0], [0.0]);
    let mut variable = Variable {
        sample,
        mean: &mut [],
        covariance,
    };
    let mut parameter = Parameter {
        mean: &mut mean,
        deviation: &mut deviation,
        covariance: &mut [0.0],
        correlation: &mut [0.0],
    };
    // The default tolerance needs about twice the iterations for digits which are not shown.
    let algorithm = Algorithm {
        tolerance: 1e-6,
        ..Algorithm::default()
    };
    let report = algorithm.fit(model, &mut variable, &mut parameter)?;
    Ok(Fit {
        scale: mean[0],
        deviation: deviation[0],
        reduced: report.chi_squared_reduced,
    })
}

/// A point on a rotation curve.
#[derive(Clone, Copy)]
struct Point {
    /// Index of its galaxy.
    galaxy: usize,
    /// Logarithm of the expected acceleration.
    x: f64,
    /// Logarithm of the observed acceleration.
    y: f64,
    /// Fraction of the expected acceleration due to stars.
    stars: f64,
    /// Independent error of [`Self::y`].
    y_error: f64,
    /// Error of [`Self::y`] shared by its galaxy due to its distance and inclination.
    y_shared: f64,
}

/// Reads the points of the galaxies passing the quality criteria of the reference.
fn read(directory: &Path, mass_to_light: [f64; 2]) -> Result<Vec<Point>, Box<dyn Error>> {
    let table = fs::read_to_string(directory.join("SPARC_Lelli2016c.mrt"))?;
    let mut galaxies = HashMap::new();
    for line in table.lines() {
        let fields = line.split_whitespace().collect::<Vec<&str>>();
        if fields.len() >= 18 && fields[1].parse::<u32>().is_ok() {
            let number = |index: usize| fields[index].parse::<f64>();
            let galaxy = (number(2)?, number(3)?, number(5)?, number(6)?, number(17)?);
            galaxies.insert(fields[0].to_string(), galaxy);
        }
    }
    let mut names = fs::read_dir(directory.join("Rotmod_LTG"))?
        .map(|entry| Ok(entry?.file_name().to_string_lossy().into_owned()))
        .collect::<Result<Vec<String>, std::io::Error>>()?;
    names.sort();
    let mut points = Vec::new();
    let mut galaxy = 0;
    for name in names {
        let Some(&(distance, distance_error, inclination, inclination_error, quality)) =
            galaxies.get(name.trim_end_matches("_rotmod.dat"))
        else {
            continue;
        };
        // Reject asymmetric rotation curves and face-on galaxies.
        #[allow(clippy::float_cmp)]
        if quality == 3.0 || inclination < 30.0 {
            continue;
        }
        let inclination_error =
            2.0 * inclination_error.to_radians() / inclination.to_radians().tan();
        let y_shared = (distance_error / distance).hypot(inclination_error) / LN_10;
        let before = points.len();
        let curve = fs::read_to_string(directory.join("Rotmod_LTG").join(&name))?;
        for line in curve.lines().filter(|line| !line.starts_with('#')) {
            let fields = line
                .split_whitespace()
                .map(str::parse)
                .collect::<Result<Vec<f64>, _>>()?;
            let [radius, speed, speed_error, gas, disk, bulge, ..] = fields[..] else {
                continue;
            };
            let stars =
                mass_to_light[0] * disk.abs() * disk + mass_to_light[1] * bulge.abs() * bulge;
            let baryons = gas.abs() * gas + stars;
            // Require a minimum precision of 10 % in velocity.
            if speed <= 0.0 || speed_error / speed >= 0.1 || baryons <= 0.0 {
                continue;
            }
            points.push(Point {
                galaxy,
                x: (baryons * UNIT / radius).log10(),
                y: (speed * speed * UNIT / radius).log10(),
                stars: (stars / baryons).clamp(0.0, 1.0),
                y_error: 2.0 * speed_error / speed / LN_10,
                y_shared,
            });
        }
        if points.len() > before {
            galaxy += 1;
        }
    }
    Ok(points)
}

/// Variances of the published table of the relation in the same logarithms and the sample.
fn published(table: &str) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let rows = table
        .lines()
        .filter_map(|line| {
            let fields = line
                .split_whitespace()
                .filter_map(|field| field.parse::<f64>().ok())
                .collect::<Vec<f64>>();
            (fields.len() == 4 && fields[0] < 0.0).then_some(fields)
        })
        .collect::<Vec<Vec<f64>>>();
    let column = |index: usize| rows.iter().map(move |row| row[index]);
    let sample = column(0).chain(column(2)).collect();
    let square = |index: usize| column(index).map(|field| field * field).collect();
    (sample, square(1), square(3))
}

/// Covariance matrix tiles of the expected and observed accelerations.
fn tiles(points: &[Point], correlated: bool, scatter: f64) -> (Vec<f64>, Vec<f64>) {
    let samples = points.len();
    let stars = MASS_TO_LIGHT_ERROR / LN_10;
    let (mut x_block, mut y_block) = (vec![0f64; samples * samples], vec![0f64; samples * samples]);
    for (row, a) in points.iter().enumerate() {
        for (column, b) in points.iter().enumerate() {
            if a.galaxy == b.galaxy && (correlated || row == column) {
                x_block[row * samples + column] = a.stars * b.stars * stars * stars;
                y_block[row * samples + column] = a.y_shared * b.y_shared;
            }
        }
        x_block[row * samples + row] += EXPECTED_ERROR * EXPECTED_ERROR;
        y_block[row * samples + row] += a.y_error * a.y_error + scatter * scatter;
    }
    (x_block, y_block)
}

/// Prints a row of results named `name`.
fn row(name: &str, result: &str) {
    println!("  {name:<32}{result}");
}

/// Prints a row of `fit` named `name`.
fn fit_row(name: &str, fit: &Fit) {
    row(
        name,
        &format!(
            "{:.3} +- {:.3} (chi-squared per degree {:.3})",
            fit.scale, fit.deviation, fit.reduced,
        ),
    );
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut arguments = std::env::args().skip(1);
    let directory = arguments.next().unwrap_or_else(|| "target/sparc".into());
    let simulations: u32 = arguments.next().map_or(Ok(0), |number| number.parse())?;
    let directory = Path::new(&directory);
    fetch(directory)?;
    // Fit the published table the way the reference describes.
    let table = directory.join("RAR.mrt");
    let (sample, x_variance, y_variance) = published(&fs::read_to_string(table)?);
    let rows = x_variance.len();
    let covariance = Covariance::new_diagonal(rows, 2)
        .with_tile(0, 0, &x_variance)?
        .with_tile(1, 1, &y_variance)?;
    let table_fit = fit(&sample, &covariance)?;
    println!("acceleration scale in 1e-10 m/s^2:");
    row("published:", "1.20  +- 0.02 (random) +- 0.24 (systematic)");
    fit_row("published table, independent:", &table_fit);
    // Build the points from the rotation curves and fit them with and without their correlation.
    let points = read(directory, MASS_TO_LIGHT)?;
    let samples = points.len();
    let sample = points
        .iter()
        .map(|point| point.x)
        .chain(points.iter().map(|point| point.y))
        .collect::<Vec<f64>>();
    let (x_block, y_block) = tiles(&points, true, SCATTER);
    let correlated = Covariance::new_block_diagonal(samples, 2)
        .with_tile(0, 0, &x_block)?
        .with_tile(1, 1, &y_block)?;
    let diagonal = |block: &[f64]| {
        (0..samples)
            .map(|row| block[row * samples + row])
            .collect::<Vec<f64>>()
    };
    let (x_diagonal, y_diagonal) = (diagonal(&x_block), diagonal(&y_block));
    let independent = Covariance::new_diagonal(samples, 2)
        .with_tile(0, 0, &x_diagonal)?
        .with_tile(1, 1, &y_diagonal)?;
    let correlated_fit = fit(&sample, &correlated)?;
    let independent_fit = fit(&sample, &independent)?;
    fit_row("rotation curves, correlated:", &correlated_fit);
    fit_row("rotation curves, independent:", &independent_fit);
    // Compare with the published value of the same curve without external field effect.
    let (x_block, y_block) = tiles(&points, true, 0.035);
    let covariance = Covariance::new_block_diagonal(samples, 2)
        .with_tile(0, 0, &x_block)?
        .with_tile(1, 1, &y_block)?;
    let simple_fit = fit_model(SIMPLE, &sample, &covariance)?;
    row("Desmond (2023), same curve:", "1.070 +- 0.032");
    fit_row("rotation curves, same curve:", &simple_fit);
    // The published systematic uncertainty represents a normalization of the mass-to-light ratio.
    let mut moved = [[0f64; 2]; 2];
    for (index, factor) in [0.8, 1.2].into_iter().enumerate() {
        let points = read(directory, MASS_TO_LIGHT.map(|ratio| ratio * factor))?;
        let (x_block, y_block) = tiles(&points, true, SCATTER);
        let sample = points
            .iter()
            .map(|point| point.x)
            .chain(points.iter().map(|point| point.y))
            .collect::<Vec<f64>>();
        let covariance = Covariance::new_block_diagonal(points.len(), 2)
            .with_tile(0, 0, &x_block)?
            .with_tile(1, 1, &y_block)?;
        let correlated_fit = fit(&sample, &covariance)?;
        let (x_block, y_block) = tiles(&points, false, SCATTER);
        let diagonal = |block: &[f64]| {
            (0..points.len())
                .map(|row| block[row * points.len() + row])
                .collect::<Vec<f64>>()
        };
        let covariance = Covariance::new_diagonal(points.len(), 2)
            .with_tile(0, 0, &diagonal(&x_block))?
            .with_tile(1, 1, &diagonal(&y_block))?;
        let independent_fit = fit(&sample, &covariance)?;
        row(
            &format!("mass-to-light ratios times {factor}:"),
            &format!(
                "correlated {:.3}, independent {:.3}",
                correlated_fit.scale, independent_fit.scale,
            ),
        );
        moved[index] = [correlated_fit.scale, independent_fit.scale];
    }
    // It agrees with the comparable joint inference within about one and a half deviations.
    assert!((simple_fit.scale - 1.07).abs() < 1.5 * simple_fit.deviation);
    // The fits only agree with the published value within its systematic uncertainty.
    assert!((1.09..1.14).contains(&table_fit.scale));
    assert!((1.05..1.11).contains(&correlated_fit.scale));
    assert!((correlated_fit.reduced - 1.0).abs() < 0.05);
    // Scaling the mass-to-light ratios moves a fit with independent errors like the published
    // systematic uncertainty but hardly moves a fit with correlated errors.
    assert!((0.2..0.3).contains(&(0.5 * (moved[0][1] - moved[1][1]))));
    assert!(0.5 * (moved[0][0] - moved[1][0]) < 0.05);
    // Their correlation doubles the deviation.
    assert!(correlated_fit.deviation > 1.5 * independent_fit.deviation);
    if simulations == 0 {
        return Ok(());
    }
    // Simulate these errors around the curve of the published acceleration scale.
    let truth = 1.2;
    println!("running {simulations} simulations with an acceleration scale of {truth}:");
    // Seed of 256 bits for deterministic high-quality randomness.
    let seed = b"5a1e8c37d2f04b96a7c31e58d09f6b24c83a7d15e62f09b4d17c8a3e5f60b2d9";
    // Cryptographically secure pseudorandom number generator (CSPRNG).
    let mut csprng = ChaChaRng::from_seed(const_hex::decode_to_array(seed).unwrap());
    let galaxies = points
        .iter()
        .map(|point| point.galaxy + 1)
        .max()
        .unwrap_or(0);
    let stars = MASS_TO_LIGHT_ERROR / LN_10;
    let simulated = (0..simulations)
        .map(|_| {
            let mut sample = vec![0f64; 2 * samples];
            let shared = (0..2 * galaxies)
                .map(|_| csprng.sample(StandardNormal))
                .collect::<Vec<f64>>();
            for (index, point) in points.iter().enumerate() {
                let x_noise: f64 = csprng.sample(StandardNormal);
                let y_noise: f64 = csprng.sample(StandardNormal);
                sample[index] =
                    point.x + point.stars * stars * shared[point.galaxy] + EXPECTED_ERROR * x_noise;
                sample[samples + index] = (MODEL.fx)(&[point.x], &[truth])
                    + point.y_shared * shared[galaxies + point.galaxy]
                    + point.y_error.hypot(SCATTER) * y_noise;
            }
            sample
        })
        .collect::<Vec<Vec<f64>>>();
    // Fit them in parallel as the fits with correlated errors take long. More than four threads
    // gain little or even slow down, presumably as dense matrices saturate the memory bandwidth.
    let threads = thread::available_parallelism()?.get().min(4);
    let fits = thread::scope(|scope| {
        let mut workers = Vec::new();
        for chunk in simulated.chunks(simulated.len().div_ceil(threads)) {
            let covariances = [&correlated, &independent];
            workers.push(scope.spawn(move || {
                chunk
                    .iter()
                    .map(|sample| covariances.map(|covariance| fit(sample, covariance).unwrap()))
                    .collect::<Vec<[Fit; 2]>>()
            }));
        }
        workers
            .into_iter()
            .flat_map(|worker| worker.join().unwrap())
            .collect::<Vec<[Fit; 2]>>()
    });
    let mut estimates = [Vec::new(), Vec::new()];
    let mut reported = [0f64; 2];
    let mut covered = [0u32; 2];
    for fits in &fits {
        for (way, fit) in fits.iter().enumerate() {
            estimates[way].push(fit.scale);
            reported[way] += fit.deviation;
            covered[way] += u32::from((fit.scale - truth).abs() <= 1.96 * fit.deviation);
        }
    }
    let count = f64::from(simulations);
    let mut actual = [0f64; 2];
    let mut coverage = [0f64; 2];
    for (way, name) in ["correlated:", "independent:"].into_iter().enumerate() {
        let mean = estimates[way].iter().sum::<f64>() / count;
        let variance = estimates[way]
            .iter()
            .map(|estimate| (estimate - mean).powi(2))
            .sum::<f64>();
        actual[way] = (variance / (count - 1.0)).sqrt();
        coverage[way] = f64::from(covered[way]) / count;
        println!(
            "  {name:<13}mean estimate {mean:.3}, actual error {:.3}, reported error {:.3}, \
             95 % interval holds true value in {:.1} %",
            actual[way],
            reported[way] / count,
            coverage[way] * 100.0,
        );
    }
    let differences = estimates[0].iter().zip(&estimates[1]).map(|(a, b)| a - b);
    let difference = (differences
        .map(|difference| difference * difference)
        .sum::<f64>()
        / count)
        .sqrt();
    println!(
        "  {:<13}typical {difference:.3}, on the real data {:.3}",
        "difference:",
        (correlated_fit.scale - independent_fit.scale).abs(),
    );
    // Fewer simulations are too noisy to assert anything.
    if simulations >= 100 {
        // A fit with independent errors understates its error and its interval misses.
        assert!(reported[1] / count < 0.5 * actual[1]);
        assert!(coverage[1] < 0.7);
        // A fit with correlated errors reports its actual error and its interval holds.
        assert!((reported[0] / count / actual[0] - 1.0).abs() < 0.25);
        assert!(coverage[0] > 0.88);
    }
    Ok(())
}

//! NUTS sampler of the models of `examples/rar.rs`, which confirms what that example states.
//!
//! It samples what the fit assumes with the galaxy terms and the true expected accelerations as
//! parameters, so about 3000 of them, and the variations of it that the fit cannot express. It is a
//! one-off validation, not part of the published crate or of its CI, so it may need adapting to
//! later versions of the crate. It has its own workspace as `nuts-rs` is not to become a dependency
//! of the crate.
//!
//! Run `cargo run --example rar` in the parent directory once to download the data, and then here:
//!
//! ```text
//! cargo run --release
//! cargo run --release -- ../target/sparc 40
//! ```
//!
//! The first samples the real data in about two minutes. The second takes about a minute per ten
//! simulations with an acceleration scale of `1.2`. It reproduces these statements of the example:
//!
//! ```text
//! 2696 points of 147 galaxies
//! OEFPIL McGaugh, correlated errors                   1.0777 +- 0.0316
//! OEFPIL McGaugh, independent errors                  1.1172 +- 0.0152
//! OEFPIL Simple, correlated errors, Desmond's scatter 1.1103 +- 0.0299
//! NUTS   McGaugh, correlated errors                   1.0887 +- 0.0319
//! NUTS   McGaugh, independent errors                  1.1302 +- 0.0152
//! NUTS   McGaugh, exactly log-normal masses           1.0500 +- 0.0308
//! NUTS   Simple, correlated errors, free scatter      1.1205 +- 0.0306, scatter 0.038
//! NUTS   Simple, Desmond's main ingredients           1.0766 +- 0.0301, scatter 0.037
//! 40 simulations: mean error of OEFPIL -0.0019, of NUTS +0.0097
//! ```
//!
//! Each sampled value has a Monte Carlo error of about `0.001`, so its third decimal is uncertain.
//! The seeds are fixed, so these values are reproducible on a platform as long as the versions in
//! `Cargo.lock` are kept.

use nuts_rs::{
    Chain, CpuLogpFunc, CpuMath, CpuMathError, DiagNutsSettings, HasDims, LogpError, Settings,
    rand::{Rng, SeedableRng, rngs::ChaCha8Rng},
};
use oefpil::{Algorithm, Covariance, Model, Parameter, Variable};
use std::{
    collections::HashMap,
    f64::consts::{LN_10, TAU},
    fs,
    path::Path,
    sync::Arc,
    thread,
    time::Instant,
};

/// Mass-to-light ratios of disks and bulges in solar units adopted by the reference.
const MASS_TO_LIGHT: [f64; 2] = [0.5, 0.7];
/// Relative uncertainty of the mass-to-light ratio of a galaxy.
const MASS_TO_LIGHT_ERROR: f64 = 0.25;
/// Error of the logarithm of the stellar mass of a galaxy in dex.
const STARS_ERROR: f64 = MASS_TO_LIGHT_ERROR / LN_10;
/// Relative uncertainty of the gas mass of a galaxy, which the fit does not comprise.
const GAS_ERROR: f64 = 0.1;
/// Independent error of the logarithm of the expected acceleration in dex.
const EXPECTED_ERROR: f64 = 0.01;
/// Intrinsic scatter of the logarithm of the observed acceleration in dex.
const SCATTER: f64 = 0.05;
/// Square kilometers per square second per kiloparsec in meters per square second.
const UNIT: f64 = 1e6 / 3.085_677_581e19;
/// Number of chains per sampling.
const CHAINS: u64 = 4;
/// Acceleration scale of the simulations.
const TRUTH: f64 = 1.2;

/// Logarithm of the observed acceleration of McGaugh et al. as a function of the logarithm of the
/// expected acceleration `x` and the acceleration scale `g`: Its value and its derivatives in `x`
/// and in `g`.
fn mcgaugh(x: f64, g: f64) -> (f64, f64, f64) {
    let ratio = (10f64.powf(x + 10.0) / g).sqrt();
    let slope = 0.5 * ratio * (-ratio).exp() / -(-ratio).exp_m1();
    (
        x - (-(-ratio).exp_m1()).log10(),
        1.0 - slope,
        slope / (g * LN_10),
    )
}

/// Terms of `g_obs = g_bar / 2 + sqrt(g_bar^2 / 4 + g_bar * a)` in logarithms as a function of the
/// logarithm of the expected acceleration `x` and the acceleration scale `a`: Its value and its
/// derivatives in `x` and in `a`.
fn simple(x: f64, a: f64) -> (f64, f64, f64) {
    let ratio = a / 10f64.powf(x + 10.0);
    let root = (0.25 + ratio).sqrt();
    let sum = 0.5 + root;
    let derivative = 1.0 / (2.0 * root * sum * LN_10);
    (
        x + sum.log10(),
        1.0 - LN_10 * ratio * derivative,
        derivative * ratio / a,
    )
}

/// The curve of [`mcgaugh()`] as a model.
const MCGAUGH: Model = Model {
    fx: &|x, p| mcgaugh(x[0], p[0]).0,
    dfdx: &[&|x, p| mcgaugh(x[0], p[0]).1],
    dfdp: &[&|x, p| mcgaugh(x[0], p[0]).2],
    implicit: false,
};

/// The curve of [`simple()`] as a model.
const SIMPLE: Model = Model {
    fx: &|x, p| simple(x[0], p[0]).0,
    dfdx: &[&|x, p| simple(x[0], p[0]).1],
    dfdp: &[&|x, p| simple(x[0], p[0]).2],
    implicit: false,
};

/// A point on a rotation curve.
#[derive(Clone, Copy)]
struct Point {
    /// Index of its galaxy.
    galaxy: usize,
    /// Logarithm of the expected acceleration.
    x: f64,
    /// Logarithm of the observed acceleration.
    y: f64,
    /// Fraction of the expected acceleration due to stars within zero and one.
    stars: f64,
    /// Fraction of the expected acceleration due to stars, which can be beyond these bounds.
    stars_unclamped: f64,
    /// Independent error of [`Self::y`].
    y_error: f64,
    /// Error of [`Self::y`] shared by its galaxy due to its distance and inclination.
    y_shared: f64,
}

/// Reads the points of the galaxies passing the quality criteria of the reference.
fn read(directory: &Path) -> Vec<Point> {
    let table = fs::read_to_string(directory.join("SPARC_Lelli2016c.mrt")).unwrap();
    let mut galaxies = HashMap::new();
    for line in table.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 18 && fields[1].parse::<u32>().is_ok() {
            let number = |index: usize| fields[index].parse::<f64>().unwrap();
            let galaxy = (number(2), number(3), number(5), number(6), number(17));
            galaxies.insert(fields[0].to_string(), galaxy);
        }
    }
    let mut names: Vec<String> = fs::read_dir(directory.join("Rotmod_LTG"))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    let (mut points, mut galaxy) = (Vec::new(), 0);
    for name in names {
        let Some(&(distance, distance_error, inclination, inclination_error, quality)) =
            galaxies.get(name.trim_end_matches("_rotmod.dat"))
        else {
            continue;
        };
        // Reject asymmetric rotation curves and face-on galaxies.
        if quality == 3.0 || inclination < 30.0 {
            continue;
        }
        let inclination_error =
            2.0 * inclination_error.to_radians() / inclination.to_radians().tan();
        let y_shared = (distance_error / distance).hypot(inclination_error) / LN_10;
        let before = points.len();
        let curve = fs::read_to_string(directory.join("Rotmod_LTG").join(&name)).unwrap();
        for line in curve.lines().filter(|line| !line.starts_with('#')) {
            let fields: Vec<f64> = line
                .split_whitespace()
                .map(|field| field.parse().unwrap())
                .collect();
            let [radius, speed, speed_error, gas, disk, bulge, ..] = fields[..] else {
                continue;
            };
            let stars =
                MASS_TO_LIGHT[0] * disk.abs() * disk + MASS_TO_LIGHT[1] * bulge.abs() * bulge;
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
                stars_unclamped: stars / baryons,
                y_error: 2.0 * speed_error / speed / LN_10,
                y_shared,
            });
        }
        if points.len() > before {
            galaxy += 1;
        }
    }
    points
}

/// Logarithms of the expected and observed accelerations as sample matrix of a fit.
fn sample_matrix(points: &[Point]) -> Vec<f64> {
    points
        .iter()
        .map(|point| point.x)
        .chain(points.iter().map(|point| point.y))
        .collect()
}

/// Covariance matrix of errors being shared within each galaxy or, if not `correlated`, of the same
/// variances with independent errors.
fn covariance(points: &[Point], correlated: bool, scatter: f64) -> Covariance {
    let x_own = EXPECTED_ERROR * EXPECTED_ERROR;
    let y_own = |point: &Point| point.y_error * point.y_error + scatter * scatter;
    if !correlated {
        let x_tile: Vec<f64> = points
            .iter()
            .map(|point| (point.stars * STARS_ERROR).powi(2) + x_own)
            .collect();
        let y_tile: Vec<f64> = points
            .iter()
            .map(|point| point.y_shared.powi(2) + y_own(point))
            .collect();
        return Covariance::new_diagonal(points.len(), 2)
            .with_tile(0, 0, &x_tile)
            .unwrap()
            .with_tile(1, 1, &y_tile)
            .unwrap();
    }
    let galaxies = || points.chunk_by(|a, b| a.galaxy == b.galaxy);
    let groups: Vec<usize> = galaxies().map(<[Point]>::len).collect();
    let (mut x_tile, mut y_tile) = (Vec::new(), Vec::new());
    for galaxy in galaxies() {
        for (row, a) in galaxy.iter().enumerate() {
            for (column, b) in galaxy.iter().enumerate() {
                let (x_diagonal, y_diagonal) = if row == column {
                    (x_own, y_own(a))
                } else {
                    (0.0, 0.0)
                };
                x_tile.push(a.stars * b.stars * STARS_ERROR * STARS_ERROR + x_diagonal);
                y_tile.push(a.y_shared * b.y_shared + y_diagonal);
            }
        }
    }
    Covariance::new_grouped(&groups, 2)
        .with_tile(0, 0, &x_tile)
        .unwrap()
        .with_tile(1, 1, &y_tile)
        .unwrap()
}

/// Fits `model` returning the acceleration scale, its deviation, and the seconds it took.
fn fit(model: Model, sample: &[f64], covariance: &Covariance) -> (f64, f64, f64) {
    let start = Instant::now();
    let (mut mean, mut deviation) = ([1.0], [0.0]);
    Algorithm::default()
        .fit(
            model,
            &mut Variable {
                sample,
                mean: &mut [],
                covariance,
            },
            &mut Parameter {
                mean: &mut mean,
                deviation: &mut deviation,
                covariance: &mut [0.0],
                correlation: &mut [0.0],
            },
        )
        .unwrap();
    (mean[0], deviation[0], start.elapsed().as_secs_f64())
}

/// A position outside of what the density is defined for.
#[derive(Debug, thiserror::Error)]
#[error("outside of the support")]
struct Outside;

impl LogpError for Outside {
    fn is_recoverable(&self) -> bool {
        true
    }
}

/// What the model comprises beyond what the fit assumes.
#[derive(Clone, Copy, Default)]
struct Variant {
    /// Log-normal stellar mass instead of its linearization.
    exact_mass: bool,
    /// Log-normal gas mass per galaxy, requires [`Self::exact_mass`].
    gas: bool,
    /// Intrinsic scatter as a parameter instead of [`SCATTER`].
    free_scatter: bool,
    /// The curve of [`simple()`] instead of [`mcgaugh()`].
    simple: bool,
}

/// Posterior density of the acceleration scale and all parameters a fit eliminates.
///
/// Its position is the acceleration scale, optionally the logarithm of the intrinsic scatter, and
/// unless independent the terms per galaxy of mass, distance, and optionally gas, followed by the
/// term per point of its expected acceleration. These terms have standard normal priors, while the
/// acceleration scale and the intrinsic scatter have flat ones.
#[derive(Clone)]
struct Density {
    /// Points of all galaxies.
    points: Arc<Vec<Point>>,
    /// Number of galaxies.
    galaxies: usize,
    /// Whether errors are shared within each galaxy, otherwise they are independent.
    correlated: bool,
    /// What the model comprises beyond what the fit assumes.
    variant: Variant,
}

impl Density {
    /// Number of parameters.
    fn size(&self) -> usize {
        let per_galaxy = if !self.correlated {
            0
        } else if self.variant.gas {
            3
        } else {
            2
        };
        1 + usize::from(self.variant.free_scatter) + per_galaxy * self.galaxies + self.points.len()
    }
}

impl HasDims for Density {
    fn dim_sizes(&self) -> HashMap<String, u64> {
        let size = self.size() as u64;
        HashMap::from([
            ("unconstrained_parameter".to_string(), size),
            ("dim".to_string(), size),
        ])
    }
}

impl CpuLogpFunc for Density {
    type LogpError = Outside;
    type FlowParameters = ();
    type ExpandedVector = Vec<f64>;

    fn dim(&self) -> usize {
        self.size()
    }

    fn logp(&mut self, position: &[f64], gradient: &mut [f64]) -> Result<f64, Outside> {
        let scale = position[0];
        if !(scale > 0.0) {
            return Err(Outside);
        }
        let variant = self.variant;
        let scatter = if variant.free_scatter {
            position[1].exp()
        } else {
            SCATTER
        };
        // Indices of the terms per galaxy and per point.
        let mass_index = 1 + usize::from(variant.free_scatter);
        let galaxies = if self.correlated { self.galaxies } else { 0 };
        let distance_index = mass_index + galaxies;
        let gas_index = distance_index + galaxies;
        let expected_index = gas_index + if variant.gas { galaxies } else { 0 };
        let mass = &position[mass_index..distance_index];
        let distance = &position[distance_index..gas_index];
        let gas = &position[gas_index..expected_index];
        let expected = &position[expected_index..];
        gradient.fill(0.0);
        let mut logp = 0.0;
        // Standard normal priors of the terms per galaxy.
        for (index, term) in mass.iter().enumerate() {
            logp -= 0.5 * term * term;
            gradient[mass_index + index] = -term;
        }
        for (index, term) in distance.iter().enumerate() {
            logp -= 0.5 * term * term;
            gradient[distance_index + index] = -term;
        }
        for (index, term) in gas.iter().enumerate() {
            logp -= 0.5 * term * term;
            gradient[gas_index + index] = -term;
        }
        // Flat prior of the intrinsic scatter rather than of its logarithm.
        if variant.free_scatter {
            logp += position[1];
            gradient[1] += 1.0;
        }
        for (index, point) in self.points.iter().enumerate() {
            let term = expected[index];
            let shared = if self.correlated {
                0.0
            } else {
                point.y_shared * point.y_shared
            };
            let variance = point.y_error * point.y_error + scatter * scatter + shared;
            let (x_error, y_shift) = if self.correlated {
                (EXPECTED_ERROR, point.y_shared * distance[point.galaxy])
            } else {
                ((point.stars * STARS_ERROR).hypot(EXPECTED_ERROR), 0.0)
            };
            // Shift of the expected acceleration and its derivatives in the mass and gas term.
            let (x_shift, x_in_mass, x_in_gas) = if !self.correlated {
                (0.0, 0.0, 0.0)
            } else if variant.exact_mass {
                let stars =
                    point.stars_unclamped * (LN_10 * STARS_ERROR * mass[point.galaxy]).exp();
                let gas_factor = if variant.gas {
                    (GAS_ERROR * gas[point.galaxy]).exp()
                } else {
                    1.0
                };
                let other = (1.0 - point.stars_unclamped) * gas_factor;
                let sum = stars + other;
                if !(sum > 1e-6) {
                    return Err(Outside);
                }
                let x_in_gas = if variant.gas {
                    other * GAS_ERROR / sum / LN_10
                } else {
                    0.0
                };
                (sum.log10(), stars * STARS_ERROR / sum, x_in_gas)
            } else {
                (
                    -point.stars * STARS_ERROR * mass[point.galaxy],
                    -point.stars * STARS_ERROR,
                    0.0,
                )
            };
            let x = point.x + x_shift - x_error * term;
            let (value, y_in_x, y_in_scale) = if variant.simple {
                simple(x, scale)
            } else {
                mcgaugh(x, scale)
            };
            if !value.is_finite() {
                return Err(Outside);
            }
            let residual = point.y - value - y_shift;
            let weighted = residual / variance;
            logp -= 0.5 * residual * weighted + 0.5 * term * term;
            gradient[0] += weighted * y_in_scale;
            if variant.free_scatter {
                logp -= 0.5 * variance.ln();
                gradient[1] +=
                    (residual * residual / variance - 1.0) * scatter * scatter / variance;
            }
            if self.correlated {
                gradient[mass_index + point.galaxy] += weighted * y_in_x * x_in_mass;
                gradient[distance_index + point.galaxy] += weighted * point.y_shared;
                if variant.gas {
                    gradient[gas_index + point.galaxy] += weighted * y_in_x * x_in_gas;
                }
            }
            gradient[expected_index + index] = -term - weighted * y_in_x * x_error;
        }
        Ok(logp)
    }

    fn expand_vector<R: Rng + ?Sized>(
        &mut self,
        _rng: &mut R,
        position: &[f64],
    ) -> Result<Vec<f64>, CpuMathError> {
        Ok(position.to_vec())
    }
}

/// Reports how well the gradient of `density` agrees with finite differences.
fn check_gradient(mut density: Density) {
    let mut random = Random(9);
    let size = density.size();
    let position: Vec<f64> = (0..size)
        .map(|index| {
            if index == 0 {
                1.1
            } else {
                0.5 * random.normal()
            }
        })
        .collect();
    let (mut gradient, mut unused) = (vec![0.0; size], vec![0.0; size]);
    density.logp(&position, &mut gradient).unwrap();
    let mut worst = 0f64;
    for index in [
        0,
        1,
        2,
        5,
        size / 3,
        size / 2,
        2 * size / 3,
        size - 2,
        size - 1,
    ] {
        let step = 1e-5;
        let mut moved = position.clone();
        moved[index] += step;
        let above = density.logp(&moved, &mut unused).unwrap();
        moved[index] -= 2.0 * step;
        let below = density.logp(&moved, &mut unused).unwrap();
        let difference = (above - below) / (2.0 * step);
        let scale = gradient[index].abs().max(1.0);
        worst = worst.max((difference - gradient[index]).abs() / scale);
    }
    println!(
        "  gradient against finite differences (correlated {}, exact mass {}, gas {}, free \
         scatter {}): worst relative error {worst:.1e}",
        density.correlated,
        density.variant.exact_mass,
        density.variant.gas,
        density.variant.free_scatter,
    );
}

/// Linear congruential generator for reproducible simulations and initial positions.
struct Random(u64);

impl Random {
    /// Uniform value from zero up to one.
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Standard normal value.
    fn normal(&mut self) -> f64 {
        (-2.0 * (1.0 - self.unit()).ln()).sqrt() * (TAU * self.unit()).cos()
    }
}

/// Split potential scale reduction of `chains`, which is about one if they have converged.
fn split_rhat(chains: &[Vec<f64>]) -> f64 {
    let parts: Vec<&[f64]> = chains
        .iter()
        .flat_map(|chain| {
            let (first, second) = chain.split_at(chain.len() / 2);
            [first, second]
        })
        .collect();
    let length = parts[0].len() as f64;
    let count = parts.len() as f64;
    let means: Vec<f64> = parts
        .iter()
        .map(|part| part.iter().sum::<f64>() / part.len() as f64)
        .collect();
    let mean = means.iter().sum::<f64>() / count;
    let between = length / (count - 1.0)
        * means
            .iter()
            .map(|part_mean| (part_mean - mean).powi(2))
            .sum::<f64>();
    let within = parts
        .iter()
        .zip(&means)
        .map(|(part, part_mean)| {
            part.iter()
                .map(|draw| (draw - part_mean).powi(2))
                .sum::<f64>()
                / (length - 1.0)
        })
        .sum::<f64>()
        / count;
    (((length - 1.0) / length * within + between / length) / within).sqrt()
}

/// Effective sample size of `chains` via Geyer's initial positive sequence.
fn effective_sample_size(chains: &[Vec<f64>]) -> f64 {
    let mut total = 0.0;
    for chain in chains {
        let length = chain.len();
        let mean = chain.iter().sum::<f64>() / length as f64;
        let variance = chain.iter().map(|draw| (draw - mean).powi(2)).sum::<f64>() / length as f64;
        let autocorrelation = |lag: usize| {
            (0..length - lag)
                .map(|index| (chain[index] - mean) * (chain[index + lag] - mean))
                .sum::<f64>()
                / (length as f64 * variance)
        };
        let mut time = -1.0;
        let mut lag = 0;
        while lag + 1 < length / 2 {
            let pair = autocorrelation(lag) + autocorrelation(lag + 1);
            if pair < 0.0 {
                break;
            }
            time += 2.0 * pair;
            lag += 2;
        }
        total += length as f64 / time.max(1.0);
    }
    total
}

/// What a sampling yields.
struct Sampling {
    /// Posterior mean of the acceleration scale.
    mean: f64,
    /// Posterior deviation of the acceleration scale.
    deviation: f64,
    /// Convergence diagnostic of the acceleration scale, see [`split_rhat()`].
    rhat: f64,
    /// Effective sample size of the acceleration scale.
    effective_sample_size: f64,
    /// Seconds it took.
    seconds: f64,
    /// Number of divergent transitions.
    divergences: u64,
    /// Posterior mean of the free intrinsic scatter, otherwise zero.
    scatter: f64,
}

/// Samples [`Density`] with [`CHAINS`] chains of `tune` and `draws` draws from `seed`.
fn sample(
    points: &[Point],
    galaxies: usize,
    correlated: bool,
    variant: Variant,
    [tune, draws, seed]: [u64; 3],
) -> Sampling {
    let shared = Arc::new(points.to_vec());
    let begin = Instant::now();
    let mut workers = Vec::new();
    for chain in 0..CHAINS {
        let density = Density {
            points: shared.clone(),
            galaxies,
            correlated,
            variant,
        };
        let size = density.size();
        workers.push(thread::spawn(move || {
            let mut settings = DiagNutsSettings::default();
            settings.num_tune = tune;
            settings.num_draws = draws;
            settings.maxdepth = 10;
            let mut rng = ChaCha8Rng::seed_from_u64(seed + chain);
            let mut sampler = settings
                .new_chain(chain, CpuMath::new(density), &mut rng)
                .unwrap();
            let mut random = Random(seed * 7 + chain + 1);
            let start: Vec<f64> = (0..size)
                .map(|index| {
                    if index == 0 {
                        1.0 + 0.2 * random.unit()
                    } else {
                        0.1 * random.normal()
                    }
                })
                .collect();
            sampler.set_position(&start).unwrap();
            let (mut scales, mut divergences, mut scatter) = (Vec::new(), 0u64, 0f64);
            for index in 0..tune + draws {
                let (draw, progress) = sampler.draw().unwrap();
                if index >= tune {
                    scales.push(draw[0]);
                    if variant.free_scatter {
                        scatter += draw[1].exp();
                    }
                    if progress.diverging {
                        divergences += 1;
                    }
                }
            }
            (scales, divergences, scatter)
        }));
    }
    let results: Vec<_> = workers
        .into_iter()
        .map(|worker| worker.join().unwrap())
        .collect();
    let chains: Vec<Vec<f64>> = results.iter().map(|result| result.0.clone()).collect();
    let all: Vec<f64> = chains.iter().flatten().copied().collect();
    let mean = all.iter().sum::<f64>() / all.len() as f64;
    let variance =
        all.iter().map(|draw| (draw - mean).powi(2)).sum::<f64>() / (all.len() - 1) as f64;
    Sampling {
        mean,
        deviation: variance.sqrt(),
        rhat: split_rhat(&chains),
        effective_sample_size: effective_sample_size(&chains),
        seconds: begin.elapsed().as_secs_f64(),
        divergences: results.iter().map(|result| result.1).sum(),
        scatter: results.iter().map(|result| result.2).sum::<f64>() / all.len() as f64,
    }
}

/// Compares fits and samplings of the real data.
fn compare(points: &[Point], galaxies: usize, desmond: Variant) {
    let sample_matrix = sample_matrix(points);
    for (name, model, correlated, scatter) in [
        ("McGaugh, correlated errors", MCGAUGH, true, SCATTER),
        ("McGaugh, independent errors", MCGAUGH, false, SCATTER),
        (
            "Simple, correlated errors, Desmond's scatter",
            SIMPLE,
            true,
            0.035,
        ),
    ] {
        let covariance = covariance(points, correlated, scatter);
        let (mean, deviation, seconds) = fit(model, &sample_matrix, &covariance);
        println!("OEFPIL {name:49} {mean:.4} +- {deviation:.4} in {seconds:.3} s");
    }
    for (name, correlated, variant) in [
        ("McGaugh, correlated errors", true, Variant::default()),
        ("McGaugh, independent errors", false, Variant::default()),
        (
            "McGaugh, exactly log-normal masses",
            true,
            Variant {
                exact_mass: true,
                ..Variant::default()
            },
        ),
        (
            "Simple, correlated errors, free scatter",
            true,
            Variant {
                simple: true,
                free_scatter: true,
                ..Variant::default()
            },
        ),
        ("Simple, Desmond's main ingredients", true, desmond),
    ] {
        let run = sample(points, galaxies, correlated, variant, [1000, 3000, 3]);
        println!(
            "NUTS   {name:49} {:.4} +- {:.4} rhat {:.3} ess {:5.0} divergences {} scatter {:.3} \
             in {:.0} s",
            run.mean,
            run.deviation,
            run.rhat,
            run.effective_sample_size,
            run.divergences,
            run.scatter,
            run.seconds,
        );
    }
    println!(
        "Desmond (2023), Simple, no external field effect:       1.070 +- 0.032 scatter 0.035"
    );
}

/// Compares fits and samplings of `simulations` around the curve of [`TRUTH`].
fn simulate(points: &[Point], galaxies: usize, simulations: u64) {
    let mut random = Random(77);
    let (mut fit_error, mut sampling_error) = (0.0, 0.0);
    let (mut difference, mut squared_difference) = (0.0, 0.0);
    let (mut worst_rhat, mut worst_size) = (0f64, f64::INFINITY);
    println!("simulation  OEFPIL +-      NUTS +-        rhat  ess");
    for simulation in 0..simulations {
        let x_shared: Vec<f64> = (0..galaxies).map(|_| random.normal()).collect();
        let y_shared: Vec<f64> = (0..galaxies).map(|_| random.normal()).collect();
        let simulated: Vec<Point> = points
            .iter()
            .map(|point| Point {
                x: point.x
                    + point.stars * STARS_ERROR * x_shared[point.galaxy]
                    + EXPECTED_ERROR * random.normal(),
                y: mcgaugh(point.x, TRUTH).0
                    + point.y_shared * y_shared[point.galaxy]
                    + point.y_error.hypot(SCATTER) * random.normal(),
                ..*point
            })
            .collect();
        let covariance = covariance(&simulated, true, SCATTER);
        let (mean, deviation, _) = fit(MCGAUGH, &sample_matrix(&simulated), &covariance);
        let run = sample(
            &simulated,
            galaxies,
            true,
            Variant::default(),
            [700, 500, 100 + simulation],
        );
        println!(
            "{simulation:10}  {mean:.4} {deviation:.4}  {:.4} {:.4}  {:.3} {:5.0}",
            run.mean, run.deviation, run.rhat, run.effective_sample_size,
        );
        fit_error += mean - TRUTH;
        sampling_error += run.mean - TRUTH;
        difference += run.mean - mean;
        squared_difference += (run.mean - mean).powi(2);
        worst_rhat = worst_rhat.max(run.rhat);
        worst_size = worst_size.min(run.effective_sample_size);
    }
    let count = simulations as f64;
    println!(
        "mean error from {TRUTH}: OEFPIL {:+.4}, NUTS {:+.4}",
        fit_error / count,
        sampling_error / count,
    );
    println!(
        "NUTS minus OEFPIL: mean {:+.4}, root mean square {:.4}",
        difference / count,
        (squared_difference / count).sqrt(),
    );
    println!("worst rhat {worst_rhat:.3}, worst ess {worst_size:.0}");
}

fn main() {
    let mut arguments = std::env::args().skip(1);
    let directory = arguments.next().unwrap_or_else(|| "../target/sparc".into());
    let simulations: Option<u64> = arguments.next().map(|number| number.parse().unwrap());
    assert!(
        Path::new(&directory).join("Rotmod_LTG").exists(),
        "no data in `{directory}`, run `cargo run --example rar` in the parent directory first",
    );
    let points = read(Path::new(&directory));
    let galaxies = points
        .iter()
        .map(|point| point.galaxy + 1)
        .max()
        .unwrap_or(0);
    println!(
        "{} points of {galaxies} galaxies; {} and {} parameters",
        points.len(),
        1 + 2 * galaxies + points.len(),
        1 + points.len(),
    );
    let desmond = Variant {
        exact_mass: true,
        gas: true,
        free_scatter: true,
        simple: true,
    };
    for (correlated, variant) in [
        (true, Variant::default()),
        (false, Variant::default()),
        (true, desmond),
    ] {
        check_gradient(Density {
            points: Arc::new(points.clone()),
            galaxies,
            correlated,
            variant,
        });
    }
    match simulations {
        Some(simulations) => simulate(&points, galaxies, simulations),
        None => compare(&points, galaxies, desmond),
    }
}

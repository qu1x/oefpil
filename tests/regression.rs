//! Regressions of input validation, input parsing, tile replacement, and statistics.

#![allow(clippy::float_cmp, clippy::suboptimal_flops)]

use oefpil::{
    Algorithm, AlgorithmError, Covariance, Logfile, Model, OefpilError, Parameter, Report,
    Variable, VariableError, read_matrix,
};
use std::{
    fs,
    io::ErrorKind,
    path::{Path, PathBuf},
    process,
};

const LINE: Model = Model {
    fx: |x, p| p[0] + p[1] * x[0],
    dfdx: &[|_, p| p[1]],
    dfdp: &[|_, _| 1.0, |x, _| x[0]],
    implicit: false,
};

/// Fits [`LINE`] to `samples` observations.
fn fit(
    covariance: &Covariance,
    samples: usize,
) -> Result<(Report, [f64; 2], [f64; 4]), OefpilError> {
    let x = (0..samples).map(|i| 1.0 + f64::from(u32::try_from(i).unwrap()));
    let y = x
        .clone()
        .enumerate()
        .map(|(i, x)| 1.0 + 2.0 * x + 0.05 * (2.3 * f64::from(u32::try_from(i).unwrap())).sin());
    let sample = x.chain(y).collect::<Vec<f64>>();
    let mut mean = [0.9, 2.1];
    let mut parameter_covariance = [0f64; 4];
    let report = Algorithm::default().fit(
        LINE,
        &mut Variable {
            sample: &sample,
            mean: &mut [],
            covariance,
        },
        &mut Parameter {
            mean: &mut mean,
            covariance: &mut parameter_covariance,
            ..Parameter::default()
        },
        Logfile::default(),
    )?;
    Ok((report, mean, parameter_covariance))
}

fn diagonal(samples: usize) -> Result<Covariance, OefpilError> {
    Covariance::new_diagonal(samples, 2)
        .with_tile(0, 0, &vec![1e-3; samples])?
        .with_tile(1, 1, &vec![1e-3; samples])
}

#[test]
fn known_covariance_requires_no_degree_of_freedom() -> Result<(), OefpilError> {
    let (report, mean, covariance) = fit(&diagonal(2)?, 2)?;
    assert_eq!(0, report.degrees_of_freedom);
    assert!(
        mean.iter()
            .chain(&covariance)
            .all(|field| field.is_finite())
    );
    Ok(())
}

#[test]
fn unknown_scale_requires_one_degree_of_freedom() -> Result<(), OefpilError> {
    assert!(matches!(
        fit(&diagonal(2)?.with_unknown_scale(), 2),
        Err(OefpilError::Variable(
            VariableError::InsufficientDegreesOfFreedom
        ))
    ));
    let (report, _mean, covariance) = fit(&diagonal(3)?.with_unknown_scale(), 3)?;
    assert_eq!(1, report.degrees_of_freedom);
    assert!(covariance.iter().all(|field| field.is_finite()));
    Ok(())
}

/// Text file removed when dropped.
struct TextFile(PathBuf);

impl TextFile {
    fn new(name: &str, text: &str) -> Self {
        let path = std::env::temp_dir().join(format!("oefpil-{}-{name}", process::id()));
        fs::write(&path, text).unwrap();
        Self(path)
    }
}

impl Drop for TextFile {
    fn drop(&mut self) {
        fs::remove_file(&self.0).unwrap();
    }
}

#[test]
fn read_matrix_of_any_shape() {
    let file = TextFile::new("wide", "1 2 3\n4 5 6\n");
    assert_eq!(
        vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        read_matrix(&file.0, 2, 3, false).unwrap()
    );
    assert_eq!(
        vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0],
        read_matrix(&file.0, 2, 3, true).unwrap()
    );
    let file = TextFile::new("tall", "1\n2\n3\n");
    assert_eq!(
        vec![1.0, 2.0, 3.0],
        read_matrix(&file.0, 3, 1, false).unwrap()
    );
    assert_eq!(
        vec![1.0, 2.0, 3.0],
        read_matrix(&file.0, 3, 1, true).unwrap()
    );
}

#[test]
fn read_matrix_reports_found_dimensions() {
    let file = TextFile::new("short", "1 2\n");
    let error = read_matrix(&file.0, 2, 2, false).unwrap_err();
    assert_eq!(ErrorKind::InvalidData, error.kind());
    assert_eq!("expected 2 row(s), found 1 row(s)", error.to_string());
    let error = read_matrix(&file.0, 1, 3, false).unwrap_err();
    assert_eq!(ErrorKind::InvalidData, error.kind());
    assert_eq!("expected 3 column(s), found 2 column(s)", error.to_string());
}

#[test]
fn diagonal_tile_replaces_block_tile() -> Result<(), OefpilError> {
    const SAMPLES: usize = 8;
    let variance = [1e-3; SAMPLES];
    let mut block = [0.5e-3; SAMPLES * SAMPLES];
    for sample in 0..SAMPLES {
        block[sample * SAMPLES + sample] = 1e-3;
    }
    let block_diagonal = || -> Result<_, OefpilError> {
        Covariance::new_block_diagonal(SAMPLES, 2)
            .with_tile(0, 0, &variance)?
            .with_tile(1, 1, &variance)
    };
    let full = || -> Result<_, OefpilError> {
        Covariance::new_full(SAMPLES, 2)
            .with_tile(0, 0, &variance)?
            .with_tile(1, 1, &variance)
    };
    let cross = [0.2e-3; SAMPLES];
    // Replace block tile on the diagonal, and block tile off the diagonal of full matrix.
    let replaced = block_diagonal()?
        .with_tile(0, 0, &block)?
        .with_tile(0, 0, &variance)?;
    assert_eq!(block_diagonal()?.to_full(), replaced.to_full());
    assert_eq!(
        fit(&block_diagonal()?, SAMPLES)?.1,
        fit(&replaced, SAMPLES)?.1
    );
    let replaced = full()?
        .with_tile(0, 0, &block)?
        .with_tile(0, 0, &variance)?
        .with_tile(0, 1, &[0.1e-3; SAMPLES * SAMPLES])?
        .with_tile(0, 1, &cross)?;
    let clean = full()?.with_tile(0, 1, &cross)?;
    assert_eq!(clean.to_full(), replaced.to_full());
    assert_eq!(fit(&clean, SAMPLES)?.1, fit(&replaced, SAMPLES)?.1);
    Ok(())
}

#[test]
fn p_value_is_accurate_in_upper_tail() -> Result<(), OefpilError> {
    // Closed forms for even degrees of freedom.
    let closed = |chi_squared: f64, degrees_of_freedom: usize| match degrees_of_freedom {
        2 => (-0.5 * chi_squared).exp(),
        4 => (-0.5 * chi_squared).exp() * (1.0 + 0.5 * chi_squared),
        6 => {
            let x = 0.5 * chi_squared;
            (-x).exp() * (1.0 + x + 0.5 * x * x)
        }
        _ => unreachable!(),
    };
    for degrees_of_freedom in [2, 4, 6] {
        for chi_squared in [0.5, 1.0, 3.0, 8.0, 40.0, 200.0, 1000.0] {
            let report = Report {
                chi_squared,
                degrees_of_freedom,
                ..Report::default()
            };
            let (actual, expected) = (
                report.chi_squared_p_value()?,
                closed(chi_squared, degrees_of_freedom),
            );
            assert!(
                (actual - expected).abs() <= 1e-12 * expected,
                "chi2 {chi_squared} dof {degrees_of_freedom}: {actual} != {expected}"
            );
        }
    }
    let report = Report {
        chi_squared: 0.0,
        degrees_of_freedom: 3,
        ..Report::default()
    };
    assert_eq!(1.0, report.chi_squared_p_value()?);
    let report = Report {
        chi_squared: f64::INFINITY,
        degrees_of_freedom: 3,
        ..Report::default()
    };
    assert_eq!(0.0, report.chi_squared_p_value()?);
    Ok(())
}

#[test]
fn logfile_path_with_nul_byte_is_error() -> Result<(), OefpilError> {
    let covariance = diagonal(4)?;
    let result = Algorithm::default().fit(
        LINE,
        &mut Variable {
            sample: &[1.0, 2.0, 3.0, 4.0, 3.1, 4.9, 7.1, 8.9],
            mean: &mut [],
            covariance: &covariance,
        },
        &mut Parameter {
            mean: &mut [0.9, 2.1],
            covariance: &mut [0f64; 4],
            ..Parameter::default()
        },
        Logfile::Custom(Path::new("log\0file")),
    );
    let Err(OefpilError::Io(error)) = result else {
        panic!("expected I/O error");
    };
    assert_eq!(ErrorKind::InvalidInput, error.kind());
    Ok(())
}

#[test]
fn converging_in_last_allowed_iteration_succeeds() -> Result<(), OefpilError> {
    let covariance = diagonal(4)?;
    let fit = |iteration_limit| {
        Algorithm {
            iteration_limit,
            ..Algorithm::default()
        }
        .fit(
            LINE,
            &mut Variable {
                sample: &[1.0, 2.0, 3.0, 4.0, 3.1, 4.9, 7.1, 8.9],
                mean: &mut [],
                covariance: &covariance,
            },
            &mut Parameter {
                mean: &mut [0.9, 2.1],
                covariance: &mut [0f64; 4],
                ..Parameter::default()
            },
            Logfile::default(),
        )
    };
    let iterations = fit(100)?.iterations;
    assert!(iterations >= 2);
    assert_eq!(fit(iterations)?.iterations, iterations);
    assert!(matches!(
        fit(iterations - 1),
        Err(OefpilError::Algorithm(AlgorithmError::IterationLimit))
    ));
    Ok(())
}

/// Reference values are exact p-values of closed forms evaluated with 80 decimal digits, which exist
/// as `nu / 2` is an integer or half-integer.
#[test]
#[allow(clippy::unreadable_literal)]
fn dcdchi_agrees_with_exact_values() {
    #[rustfmt::skip]
    let table = [
        (0.0, 1.0, 1.0),
        (1e-300, 1.0, 1.0),
        (1e-10, 1.0, 0.9999920211543921),
        (0.5, 1.0, 0.4795001221869535),
        (1.0, 1.0, 0.3173105078629141),
        (3.84, 1.0, 0.050043521248705106),
        (100.0, 1.0, 1.523970604832105e-23),
        (1.0, 2.0, 0.6065306597126334),
        (10.0, 2.0, 0.006737946999085467),
        (1400.0, 2.0, 9.85967654375977e-305),
        (2.5, 3.0, 0.4752910833430206),
        (0.1, 7.0, 0.9999976885812014),
        (7.0, 7.0, 0.42887985755305474),
        (30.0, 7.0, 9.495972508134183e-05),
        (1e-05, 10.0, 1.0),
        (10.0, 10.0, 0.4404932850652124),
        (1000.0, 10.0, 1.8702907209159498e-208),
        (20.0, 30.0, 0.9165415270653372),
        (30.0, 30.0, 0.46565370894400965),
        (45.0, 30.0, 0.038601758266317336),
        (95.0, 100.0, 0.6225792918781735),
        (100.0, 100.0, 0.48119168452795674),
        (150.0, 100.0, 0.0009039320423540091),
        (285.0, 300.0, 0.7242624085345621),
        (300.0, 300.0, 0.4891417702506403),
        (788.7298334620741, 750.0, 0.1585515148958366),
        (1000.0, 1000.0, 0.49405285382923964),
        (1050.0, 1000.0, 0.1324740568214771),
        (910.5572809000084, 1000.0, 0.9797163290773768),
        (1100.0, 1000.0, 0.014614408126295194),
        (8500.0, 5000.0, 1.4560849131328786e-186),
        (4900.0, 5000.0, 0.8413611103322287),
        (9985.857864376268, 10000.0, 0.5379746576481316),
        (10000.0, 10000.0, 0.49811936596618267),
        (10200.0, 10000.0, 0.07932888107761897),
        (100000.0, 100000.0, 0.4994052918952067),
        (101000.0, 100000.0, 0.012868840377233669),
        (1000000.0, 1000000.0, 0.4998119368033945),
        (1008485.2813742385, 1000000.0, 1.0914018031672503e-09),
        (990000.0, 1000000.0, 0.99999999999935),
        (1042426.4068711929, 1000000.0, 1.1176698881326725e-192),
    ];
    let dcdchi = |chisq, nu| {
        let (mut p, mut q, mut error) = (f64::NAN, f64::NAN, 0);
        unsafe {
            oefpil_sys::dcdchi(chisq, nu, &raw mut p, &raw mut q, &raw mut error);
        }
        assert_eq!(0, error);
        q
    };
    for (chisq, nu, q) in table {
        let actual = dcdchi(chisq, nu);
        assert!(
            (actual - q).abs() <= 1e-12 * q,
            "chi2 {chisq} nu {nu}: {actual} != {q}"
        );
    }
    // Close to its median for a billion degrees of freedom instead of zero.
    assert!((dcdchi(1e9, 1e9) - 0.5).abs() < 1e-4);
}

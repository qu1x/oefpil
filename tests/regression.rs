//! Regressions of input validation, input parsing, tile replacement, and statistics.

#![allow(clippy::float_cmp, clippy::suboptimal_flops)]

use oefpil::{
    Algorithm, Covariance, Logfile, Model, OefpilError, Parameter, Report, Variable, VariableError,
    read_matrix,
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

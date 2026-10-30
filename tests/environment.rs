//! Environment variables of the C library, in a test binary of its own as they are process-wide.

#![allow(clippy::suboptimal_flops)]

use oefpil::{Algorithm, Covariance, Logfile, Model, OefpilError, Parameter, Variable, Verbosity};
use std::{env, fs, process};

const LINE: Model = Model {
    fx: |x, p| p[0] + p[1] * x[0],
    dfdx: &[|_, p| p[1]],
    dfdp: &[|_, _| 1.0, |x, _| x[0]],
    implicit: false,
};

fn log_of_fit() -> Result<String, OefpilError> {
    let path = env::temp_dir().join(format!("oefpil-{}-log", process::id()));
    let covariance = Covariance::new_diagonal(4, 2)
        .with_tile(0, 0, &[1e-3; 4])?
        .with_tile(1, 1, &[1e-3; 4])?;
    Algorithm {
        verbosity: Verbosity::ParameterMean,
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
        Logfile::Custom(&path),
    )?;
    let log = fs::read_to_string(&path)?;
    fs::remove_file(&path)?;
    Ok(log)
}

#[test]
fn print_format_is_applied_if_valid() -> Result<(), OefpilError> {
    let default = log_of_fit()?;
    assert!(default.contains("e+00"));
    // SAFETY: This is the only test of this binary and hence there are no other threads.
    unsafe { env::set_var("OEFPIL_PRINTFORMAT", "3f") };
    let log = log_of_fit()?;
    assert!(!log.contains("e+00"));
    assert!(log.contains("2.0") || log.contains("1.9"));
    // Neither crashes nor passes through as format specifier.
    for invalid in ["s", "n", "%s", "3", "3fe", "", "99999999999999999999f"] {
        // SAFETY: See above.
        unsafe { env::set_var("OEFPIL_PRINTFORMAT", invalid) };
        assert_eq!(default, log_of_fit()?, "{invalid}");
    }
    // SAFETY: See above.
    unsafe { env::remove_var("OEFPIL_PRINTFORMAT") };
    Ok(())
}

# Version 0.5.0 (unreleased)

  * Rewrite in pure Rust without system dependencies. `oefpil-sys` is merely used to compare both
    implementations in the test suite. The results agree up to rounding errors.
  * Add optional feature `faer` to speed up large block-diagonal and full covariance matrices.
  * Support `no_std` environments with an allocator by disabling the new default feature `std`,
    which gates `Logfile` and `read_matrix()`.
  * Allow capturing closures in `Model` as `Expression` is now `&dyn Fn(&[f64], &[f64]) -> f64`
    instead of a function pointer. Write `&|x, p| ...` instead of `|x, p| ...`. A `Model` is no
    longer `Sync` and hence cannot be a `static` but it can still be a `const`.
  * Remove `rand` and `rand_distr` from the public API along with feature `random`.
    `Distribution::sample()` takes any `FnMut() -> f64` yielding standard normal values and
    `Distribution` as well as `DistributionError` are always available.
  * Remove `log` argument and logging of `Algorithm::fit()`. Add `Algorithm::fit_with_log()` taking
    any `core::fmt::Write` as log. Open a `Logfile` via `Logfile::open()` and flush it via
    `LogWriter::finish()`.
  * Return `OefpilError::Log` if writing the log fails. The fit still yields its estimates.
  * Add `Algorithm::format` to replace environment variable `OEFPIL_PRINTFORMAT`. Remove environment
    variables `OEFPIL_DEBUG` and `OEFPIL_MAXIT`, the latter is covered by
    `Algorithm::iteration_limit`.
  * Log non-finite model evaluations into the log if not `Verbosity::Silent` instead of always
    writing them to the standard error stream. Omit the workspace sizes of LAPACK from the log.
  * Remove the limit on the sizes of matrices and hence `OefpilError::TryFromInt`.
  * Unwind panics of models through `Algorithm::fit()` instead of aborting the process.
  * Define `Criterion` and `Verbosity` in this crate. They are no longer re-exports.

# Version 0.4.0 (2026-10-31)

  * Support correlations among any number of variables via `Covariance::new_diagonals()` and
    `Covariance::new_full()`, which no longer panic if `variables != 2`.
  * Fix `Covariance::set_tile()` storing asymmetric block tiles among variables transposed instead
    of in documented row-major order. This changes the results of fits using such tiles.
  * Fix `Algorithm::fit()` for diagonal tiles of `Covariance::new_full()`, which yielded wrong
    results or failed.
  * Return `AlgorithmError::NumericalError` if the covariance matrix is not positive definite for
    any tiling mode. Previously, `Covariance::new_diagonal()`, `Covariance::new_block_diagonal()`,
    and `Covariance::new_diagonals()` could report success with meaningless `Report::chi_squared`.
    This includes zero variances and unset tiles. For a variable without errors, use a tiny positive
    variance instead of zero, which yields the same estimates along with a finite chi-squared.
  * Return `AlgorithmError::NumericalError` as soon as a model evaluates to a non-finite value
    rather than possibly `AlgorithmError::IterationLimit`.
  * Succeed if a fit converges in the last allowed iteration instead of returning
    `AlgorithmError::IterationLimit`.
  * Improve accuracy of `Report::chi_squared_p_value()` for many degrees of freedom and fix it
    yielding `0.0` around the expected value for a billion degrees of freedom or more.
  * Fix `Algorithm::fit()` requiring one degree of freedom for a fully known covariance matrix but
    none for `Covariance::with_unknown_scale()` which yielded a NaN covariance matrix, as opposed to
    what `VariableError::InsufficientDegreesOfFreedom` documents.
  * Fix `read_matrix()` yielding wrong fields or panicking for non-square matrices without
    `transpose` and fix the number of found rows and columns in its error messages.
  * Fix `Covariance::set_tile()` keeping off-diagonal fields when replacing a block tile with a
    diagonal tile.
  * Fix `Report::chi_squared_p_value()` losing accuracy in the upper tail and yielding `0.0` for
    small p-values.
  * Fix segmentation fault if `OEFPIL_PRINTFORMAT` is defined. It is ignored unless it comprises up
    to three digits followed by one of `eEfFgG`.
  * Return `OefpilError::TryFromInt` instead of corrupting memory if the size of a matrix exceeds
    what the C library can address.
  * Close custom `Logfile` on conversion errors and adhere to aliasing rules in `Algorithm::fit()`.
  * Return `OefpilError::Io` instead of panicking if the path of a custom `Logfile` has a NUL byte.
  * Correct documentation (e.g., definition of `Report::chi_squared`, aborting panics in models).
  * Update KaTeX.
  * Bump `oefpil-sys`.

# Version 0.3.1 (2025-10-16)

  * Bump KaTeX and fix its integration.
  * Bump dependencies.
  * Drop removed `doc_auto_cfg` feature.

# Version 0.3.0 (2025-02-22)

  * Rust Edition 2024.
  * Bump MSRV.
  * Bump dependencies.

# Version 0.2.0 (2025-01-14)

  * Refactor constructors for sample covariance matrix introducing `with_unknown_scale()`.
  * Redefine model of `sine.rs` example in degrees rather than radians.
  * Improve documentation (e.g., add troubleshooting section, fix index in formula).
  * Update KaTeX.

# Version 0.1.1 (2024-12-15)

  * Clean up documentation (e.g., links).

# Version 0.1.0 (2024-12-12)

  * Initial release.

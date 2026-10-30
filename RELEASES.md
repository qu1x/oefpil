# Version 0.4.0 (unreleased)

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

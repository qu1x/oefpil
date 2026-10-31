# Version 0.4.0 (unreleased)

  * Import respective sources of OEFPIL [307ee4c] into [src] and patch them as follows:
      * Generalize `dpottrf_diags()`, `dpottrf_full()`, `dtrtsv_diags()`, and `dtrtsv_full()` to any
        number of variables.
      * Fix stride of diagonal tiles in `dgetmv_full()`.
      * Report non-positive pivots in `dpottrf_diag()` and `dpottrf_diags()` and failures of any
        tile in `dpottrf_blockdiag()` rather than only of the last one.
      * Check all values in `check_inf_nan()` rather than only the first one and stop iterating at
        once if any is not finite.
      * Remove unconditional print of number of variables.
      * Fix out-of-bounds read when printing parameters if `OEFPIL_DEBUG` is defined.
      * Fix segmentation fault if `OEFPIL_PRINTFORMAT` is defined. It is ignored unless it comprises
        up to three digits followed by one of `eEfFgG`.
      * Clear block tile when setting diagonal tile.
      * Report convergence in the last allowed iteration as such instead of having reached the
        iteration limit.
  * Compute upper incomplete gamma function ratio directly in `dgami()` to keep its accuracy.
  * Compute lower incomplete gamma function ratio directly in `dgami()` and evaluate the prefix of
    both ratios without cancellation to improve accuracy for many degrees of freedom. This also
    fixes a lower ratio of `1.0` for a billion degrees of freedom or more and removes `special`.
  * Use `c_char` for `dpotrf()` and `dtrmv()` to support targets where it is unsigned.
  * Correct documentation of `oefpil_tcm_full_set_tile_full()` which expects a column-major slice.
  * Correct documentation of `Mode` where `Mode::BlockDiagonal` encodes a block tile.

# Version 0.3.1 (2025-10-16)

  * Bump dependencies.

# Version 0.3.0 (2025-02-22)

  * Rust Edition 2024.
  * Bump MSRV.
  * Bump dependencies.

# Version 0.2.0 (2025-01-14)

  * Improve documentation (i.e., use table for system requirements).

# Version 0.1.1 (2024-12-15)

  * Clean up documentation (e.g., license files, copyright notices, system requirements).
  * Fix CI build.

# Version 0.1.0 (2024-12-12)

  * Import respective sources of OEFPIL [923643b] into [src].

[307ee4c]: https://gitlab.com/cmi6014/oefpil/-/tree/307ee4cfe8245a9db79bb411560392a576709b8d
[923643b]: https://gitlab.com/cmi6014/oefpil/-/tree/923643b76ee7cc1aa83750a73e9e5eae3ed6531b
[src]: src

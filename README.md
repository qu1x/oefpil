# oefpil

[![Build][]](https://github.com/qu1x/oefpil/actions/workflows/build.yml)
[![Documentation][]](https://docs.rs/oefpil)
[![Downloads][]](https://crates.io/crates/oefpil)
[![Version][]](https://crates.io/crates/oefpil)
[![Rust][]](https://www.rust-lang.org)

[![License][]](https://opensource.org/licenses)
[![DOI](https://zenodo.org/badge/902104491.svg)](https://doi.org/10.5281/zenodo.15667346)

[Build]: https://github.com/qu1x/oefpil/actions/workflows/build.yml/badge.svg
[Documentation]: https://docs.rs/oefpil/badge.svg
[Downloads]: https://img.shields.io/crates/d/oefpil.svg
[Version]: https://img.shields.io/crates/v/oefpil.svg
[Rust]: https://img.shields.io/badge/rust-v1.85.0-brightgreen.svg
[License]: https://img.shields.io/badge/License-MIT%2FApache--2.0-blue.svg

Optimum Estimate of Function Parameters by Iterated Linearization (OEFPIL)

Algorithm for nonlinear function fitting to data with errors in variables where correlation,
both within variables and among variables, might be present. In principle, OEFPIL can be
employed for fitting both explicit and implicit functions of any number of variables.
Importantly, apart from the parameter estimates, OEFPIL also yields their covariance matrix,
required for further analyses.[^1] Common methods such as ordinary nonlinear least squares are
not capable of treating general uncertainties and correlations in both dependent and independent
variables. A new computation method for nonlinear curve fitting to data with a general
covariance structure is introduced. Numerical simulations show that the new method yields
parameter estimates in agreement with other methods for simple covariance structures. The
obtained uncertainty estimates are in agreement with Monte Carlo studies.[^2]

The most notable features of OEFPIL are as follows:[^2]

  * The OEFPIL algorithm features an improved estimation of uncertainties and covariance
    matrices of the fitted parameters which can be further propagated to physical output
    variables.
  * More general covariance matrices for input data can be applied, although so far they do not
    seem to lead to significant shifts in the physical output variables. However, estimates of
    uncertainties of these are obviously affected.
  * Unlike ordinary least squares OEFPIL is not sensitive to the choice of dependent and
    independent variables.
  * Unlike non-linear least squares and orthogonal distance regression, a reasonable uncertainty
    estimate can be deduced directly without resorting to Monte Carlo studies, evading thus long
    computation times.
  * A possible disadvantage of OEFPIL is a higher sensitivity to initial estimates of parameters
    and the covariance matrix of the input data. However, so far problems have been encountered
    only for exaggerated estimates of the input covariance matrix.

[^1]: R. Šlesinger, A. C. Campbell, Z. Geršlová, V. Šindlář, and G. Wimmer, “OEFPIL: New Method
and Software Tool for Fitting Nonlinear Functions to Correlated Data With Errors in Variables”,
[2023 14th International Conference on Measurement, 126-129
(2023)](https://doi.org/10.23919/MEASUREMENT59122.2023.10164444).

[^2]: A. C. Campbell, Z. Geršlová, V. Šindlář, R. Šlesinger, and G. Wimmer, “New framework for
nanoindentation curve fitting and measurement uncertainty estimation”, [Precision Engineering
85, 166-173 (2024)](https://doi.org/10.1016/j.precisioneng.2023.10.001).

This crate provides a safe API to the [`oefpil-sys`] crate (see its system requirements) which
statically links to the [C/Fortran library].

See the [release history](RELEASES.md) to keep track of the development.

[`oefpil-sys`]: https://crates.io/crates/oefpil-sys
[C/Fortran library]: https://gitlab.com/cmi6014/oefpil

## Licenses

Except as noted below, this work is dual-licensed under either [`MIT`] or [`Apache-2.0`] at your
option. This means you can select the license you prefer. This dual-licensing approach is the
de-facto standard in the Rust ecosystem. Copyrights in this work are retained by their contributors
and no copyright assignment is required to contribute to this work. For full authorship information,
see the individual files and the version control history.

The works in [oefpil-sys](oefpil-sys) carry their own copyright notices and license terms.

[`MIT`]: LICENSE-MIT
[`Apache-2.0`]: LICENSE-APACHE

## Contributions

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in the
work by you, as defined in the [`Apache-2.0`] license, shall be dual-licensed as above, without any
additional terms or conditions.

# oefpil-sys

[![Build][]](https://github.com/qu1x/oefpil/actions/workflows/build.yml)
[![Documentation][]](https://docs.rs/oefpil-sys)
[![Downloads][]](https://crates.io/crates/oefpil-sys)
[![Version][]](https://crates.io/crates/oefpil-sys)
[![Rust][]](https://www.rust-lang.org)
[![License][]](https://opensource.org/licenses)

[Build]: https://github.com/qu1x/oefpil/actions/workflows/build.yml/badge.svg
[Documentation]: https://docs.rs/oefpil-sys/badge.svg
[Downloads]: https://img.shields.io/crates/d/oefpil-sys.svg
[Version]: https://img.shields.io/crates/v/oefpil-sys.svg
[Rust]: https://img.shields.io/badge/rust-v1.82.0-brightgreen.svg
[License]: https://img.shields.io/badge/License-MIT%2FApache--2.0-blue.svg

Rust FFI bindings to statically linked [C/Fortran library](https://gitlab.com/cmi6014/oefpil) OEFPIL

For a safe API, see the [`oefpil`](https://crates.io/crates/oefpil) crate.

See the [release history](RELEASES.md) to keep track of the development.

## System Requirements

By default, this crate dynamically links to the runtime dependency `liblapack` (e.g., package
`liblapack3` on Debian, package `lapack` on Fedora Linux or Arch Linux) and requires a C
compiler as build dependency (e.g., package `clang` or `gcc` on Debian, Fedora Linux, or Arch
Linux). With the `built-in` feature, a subset of `liblapack` and its dependency `libblas`
shipped with this crate is compiled and statically linked. This eliminates the runtime
dependency `liblapack` but requires the GCC Fortran compiler (e.g, `gfortran` on Debian,
`gcc-fortran` on Fedora Linux or Arch Linux) as build dependency which itself depends on and
complements the GCC C compiler such that GCC can compile both C and Fortran sources. It is
attempted to statically link the dependencies of the subset (i.e, `libgfortran` and
`libquadmath`) whereas dynamic linking serves as fallback if no static libraries are found.

## Licenses

This combined work is free, open source, open collaboration, and permissively licensed.

Except where noted (below and/or in individual files), this combined work is dual-licensed under
either [`MIT`] or [`Apache-2.0`] at your option. This means you can select the license you prefer.
This dual-licensing approach is the de-facto standard in the Rust ecosystem. For full authorship
information, see the individual files and/or the commit history.

The dependencies shipped with this combined work are licensed as follows:

  * [oefpil] is licensed under [`MIT`].
  * [blas], [lapack], and [math77_chi2] are licensed under [`BSD-3-Clause`].

[oefpil]: src/oefpil
[blas]: src/blas
[lapack]: src/lapack
[math77_chi2]: src/math77_chi2

[`MIT`]: LICENSE-MIT
[`Apache-2.0`]: LICENSE-APACHE
[`BSD-3-Clause`]: LICENSE-BSD

## Contributions

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in the
work by you, as defined in the [`Apache-2.0`] license, shall be dual-licensed as above, without any
additional terms or conditions.

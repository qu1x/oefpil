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
[Rust]: https://img.shields.io/badge/rust-v1.85.0-brightgreen.svg
[License]: https://img.shields.io/badge/License-MIT%2FApache--2.0-blue.svg

Rust FFI bindings to statically linked [C/Fortran library] OEFPIL

[C/Fortran library]: https://gitlab.com/cmi6014/oefpil

For a safe API, see the [`oefpil`](https://crates.io/crates/oefpil) crate.

See the [release history](RELEASES.md) to keep track of the development.

## System Requirements

By default, this crate dynamically links to the runtime dependency LAPACK and requires a C
compiler as build dependency. With the `built-in` feature enabled (marked with ☑ in the table
below), a subset of LAPACK and its dependency BLAS shipped with this crate is compiled and
statically linked. This eliminates the runtime dependency LAPACK but requires the GCC Fortran
compiler as build dependency which itself depends on and complements the GCC C compiler such
that GCC can compile both C and Fortran sources. It is attempted to statically link the
dependencies of the subset (i.e, the GNU Fortran runtime library and the GCC quad-precision math
library) whereas dynamic linking serves as fallback if no static libraries are found. The
required runtime and build dependencies are satisfied by installing following system packages
where "or" as in `|` has higher precedence than "and" as in `,`:

| Operating System | `built-in` | Runtime Dependencies | Build Dependencies            |
|------------------|:----------:|----------------------|-------------------------------|
| Debian Bookworm  | ☐          | `liblapack3`         | `gcc \| clang, liblapack-dev` |
| Debian Bookworm  | ☑          | &nbsp;               | `gfortran`                    |
| Fedora Linux     | ☐          | `lapack`             | `gcc \| clang, lapack-devel`  |
| Fedora Linux     | ☑          | &nbsp;               | `gcc-gfortran`                |
| Arch Linux       | ☐          | `lapack`             | `gcc \| clang, lapack`        |
| Arch Linux       | ☑          | &nbsp;               | `gcc-fortran`                 |

## Licenses

Except as noted below, this work is dual-licensed under either [`MIT`] or [`Apache-2.0`] at your
option. This means you can select the license you prefer. This dual-licensing approach is the
de-facto standard in the Rust ecosystem. Copyrights in this work are retained by their contributors
and no copyright assignment is required to contribute to this work. For full authorship information,
see the individual files and the version control history.

The works imported from the [C/Fortran library] are licensed as follows:

  * [oefpil](src/oefpil) is licensed under [`MIT`] with copyright notice:

    Copyright © 2020 Czech Metrology Institute

  * [lapack](src/lapack) requests to credit its authors, as suggested, by citing its [Users' Guide].
    It is licensed under **modified** [`BSD-3-Clause`] with copyright notice:

    Copyright © 1992-2013 The University of Tennessee

    Copyright © 1992-2013 The University of Tennessee Research Foundation

    Copyright © 2000-2013 The University of California Berkeley

    Copyright © 2006-2013 The University of Colorado Denver

[`MIT`]: LICENSE-MIT
[`Apache-2.0`]: LICENSE-APACHE
[`BSD-3-Clause`]: LICENSE-BSD
[Users' Guide]: https://www.netlib.org/lapack/lug/lapack_lug.html

## Contributions

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in the
work by you, as defined in the [`Apache-2.0`] license, shall be dual-licensed as above, without any
additional terms or conditions.

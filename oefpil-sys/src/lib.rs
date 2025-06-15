//! Rust FFI bindings to statically linked [C/Fortran library] OEFPIL
//!
//! [C/Fortran library]: https://gitlab.com/cmi6014/oefpil
//!
//! For a safe API, see the [`oefpil`](https://docs.rs/oefpil) crate.
//!
//! # System Requirements
//!
//! By default, this crate dynamically links to the runtime dependency LAPACK and requires a C
//! compiler as build dependency. With the `built-in` feature enabled (marked with ☑ in the table
//! below), a subset of LAPACK and its dependency BLAS shipped with this crate is compiled and
//! statically linked. This eliminates the runtime dependency LAPACK but requires the GCC Fortran
//! compiler as build dependency which itself depends on and complements the GCC C compiler such
//! that GCC can compile both C and Fortran sources. It is attempted to statically link the
//! dependencies of the subset (i.e, the GNU Fortran runtime library and the GCC quad-precision math
//! library) whereas dynamic linking serves as fallback if no static libraries are found. The
//! required runtime and build dependencies are satisfied by installing following system packages
//! where "or" as in `|` has higher precedence than "and" as in `,`:
//!
//! | Operating System | `built-in` | Runtime Dependencies | Build Dependencies            |
//! |------------------|:----------:|----------------------|-------------------------------|
//! | Debian Bookworm  | ☐          | `liblapack3`         | `gcc \| clang, liblapack-dev` |
//! | Debian Bookworm  | ☑          | &nbsp;               | `gfortran`                    |
//! | Fedora Linux     | ☐          | `lapack`             | `gcc \| clang, lapack-devel`  |
//! | Fedora Linux     | ☑          | &nbsp;               | `gcc-gfortran`                |
//! | Arch Linux       | ☐          | `lapack`             | `gcc \| clang, lapack`        |
//! | Arch Linux       | ☑          | &nbsp;               | `gcc-fortran`                 |
//!
//! # Overview
//!
//! The main function of interest is [`oefpil`]. Among other arguments, it expects the convergence
//! [`Criterion`], log [`Verbosity`], log [`FILE`] (e.g., [`stdout_file`], [`stderr_file`]), and a
//! covariance matrix tiled by variables. Among its data fields, a tiled covariance matrix (TCM)
//! comprises metadata about its tilemap and tiling [`Mode`]. The tilemap encodes via
//! [`Mode::Diagonal`] or [`Mode::Full`] which tiles are diagonal or block tiles. The number of
//! `samples` and `variables` define the number of fields per tile and the number of tiles per
//! covariance matrix. A diagonal tile stores `samples` fields whereas a block tile stores
//! `samples.pow(2)` fields. The tiling mode encodes where the tiles are and whether their mode is
//! restricted to be diagonal. For each tiling mode, there are different sets of methods for
//! allocating the tilemap and the data fields and for setting the data fields per tile.
//!
//!   * [`Mode::Diagonal`]: Diagonal tiles on the diagonal.
//!       * [`oefpil_tilemap_diagtiles_new`]: Allocates tilemap.
//!       * [`oefpil_tcm_diag_new`]: Allocates fields.
//!       * [`oefpil_tcm_diag_set_tile_diag`]: Sets fields per diagonal tile.
//!   * [`Mode::BlockDiagonal`]: Diagonal or block tiles on the diagonal.
//!       * [`oefpil_tilemap_diagtiles_new`]: Allocates tilemap.
//!       * [`oefpil_tcm_blockdiag_new`]: Allocates fields.
//!       * [`oefpil_tcm_blockdiag_set_tile_diag`]: Sets fields per diagonal tile.
//!       * [`oefpil_tcm_blockdiag_set_tile_half`]: Sets fields per block tile (row-major
//!         [triangular slice] of lower triangle).
//!       * [`oefpil_tcm_blockdiag_set_tile_full`]: Sets fields per block tile (row-major slice).
//!   * [`Mode::Diagonals`]: Diagonal tiles all over.
//!       * [`oefpil_tilemap_alltiles_new`]: Allocates tilemap.
//!       * [`oefpil_tcm_diags_new`]: Allocates fields.
//!       * [`oefpil_tcm_diags_set_tile_diag`]: Sets fields per diagonal tile.
//!   * [`Mode::Full`]: Diagonal or block tiles all over.
//!       * [`oefpil_tilemap_alltiles_new`]: Allocates tilemap.
//!       * [`oefpil_tcm_full_new`]: Allocates fields.
//!       * [`oefpil_tcm_full_set_tile_diag`]: Sets fields per diagonal tile.
//!       * [`oefpil_tcm_full_set_tile_half`]: Sets fields per block tile (row-major [triangular
//!         slice] of lower triangle).
//!       * [`oefpil_tcm_full_set_tile_full`]: Sets fields per block tile (row-major slice).
//!
//! [triangular slice]: https://en.wikipedia.org/wiki/Triangular_array

pub use libc;

use core::ffi::{c_int, c_long, c_void};
use libc::FILE;

/// Function pointer type passed to [`oefpil`] as 1st argument.
///
/// Arguments:
///
///   * `data`: User-defined structure defining model inclusive number of variables and parameters.
///   * `samples`: Number of samples per variable.
///   * `x`: Sample from independent variables (sample-major).
///   * `p`: Parameters.
///   * `fx`: Evaluated dependent variables.
///   * `dfdx`: Evaluated derivatives in independent variables (sample-major).
///   * `dfdp`: Evaluated derivatives in parameters (sample-major).
pub type Evaluate = Option<
    unsafe extern "C" fn(
        data: *mut c_void,
        samples: c_int,
        x: *const f64,
        p: *const f64,
        fx: *mut f64,
        dfdx: *mut f64,
        dfdp: *mut f64,
    ),
>;

/// Mode of tiled covariance matrix or mode of tile as part of tilemap.
///
/// Valid modes of tile as part of tilemap are:
///
///   * [`Self::None`] for an unset tile,
///   * [`Self::Diagonal`] for a diagonal tile, and
///   * [`Self::Full`] for a block tile.
///
/// Default is [`Self::None`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Mode {
    /// Mode of unset tile.
    #[default]
    None = 0,
    /// Mode of covariance matrix with diagonal tiles on its diagonal (or mode of diagonal tile).
    Diagonal = 1,
    /// Mode of covariance matrix with diagonal or block tiles on its diagonal.
    BlockDiagonal = 2,
    /// Mode of covariance matrix with diagonal tiles all over.
    Diagonals = 3,
    /// Mode of covariance matrix with diagonal or block tiles all over (or mode of block tile).
    Full = 4,
}

/// Convergence criterion.
///
/// Default is [`Self::RelPOrAbsPAndRelXOrAbsX`].
#[derive(Debug, Clone, Copy, Default)]
#[non_exhaustive]
pub enum Criterion {
    /// Convergence in relative change of parameter mean.
    RelP = 0,
    /// Convergence in [`Self::RelP`] or absolute value of parameter mean.
    RelPOrAbsP = 1,
    /// Convergence in relative change of independent variable mean.
    RelX = 2,
    /// Convergence in [`Self::RelX`] or absolute value of independent variable mean.
    RelXOrAbsX = 3,
    /// Convergence in [`Self::RelP`] and [`Self::RelX`].
    RelPAndRelX = 4,
    /// Convergence in [`Self::RelPOrAbsP`] and [`Self::RelXOrAbsX`].
    #[default]
    RelPOrAbsPAndRelXOrAbsX = 5,
    /// Convergence in chi-squared.
    ChiSquared = 6,
}

/// Log verbosity.
///
/// Default is [`Self::Silent`].
#[derive(Debug, Clone, Copy, Default)]
pub enum Verbosity {
    /// No logging.
    #[default]
    Silent = 0,
    /// Logs parameter mean.
    ParameterMean = 1,
    /// Logs covariance matrix among [`Self::ParameterMean`].
    ParameterCovariance = 2,
    /// Logs individual steps among [`Self::ParameterMean`] and [`Self::ParameterCovariance`].
    IndividualSteps = 3,
}

/// Computes the p-value `q` for `chisq` and `nu` degrees of freedom.
///
/// # Safety
///
/// This function is safe as long as the pointers are valid.
#[doc(hidden)]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn dcdchi(chisq: f64, nu: f64, p: *mut f64, q: *mut f64, ierr: *mut c_long) {
    unsafe {
        dgami(0.5 * nu, 0.5 * chisq, p, q, ierr);
    }
}

/// Computes the incomplete gamma function ratios `pans` and `qans` at `a` and `x`.
///
/// # Safety
///
/// This function is safe as long as the pointers are valid.
#[doc(hidden)]
#[unsafe(no_mangle)]
pub unsafe extern "C" fn dgami(a: f64, x: f64, pans: *mut f64, qans: *mut f64, ierr: *mut c_long) {
    unsafe {
        (*pans, *ierr) = if x >= 0.0 && a > 0.0 {
            (special::Gamma::inc_gamma(x, a), false.into())
        } else {
            (f64::NAN, true.into())
        };
        *qans = 1.0 - *pans;
    }
}

unsafe extern "C" {
    /// Returns the standard input file.
    pub safe fn stdin_file() -> *mut FILE;
    /// Returns the standard output file.
    pub safe fn stdout_file() -> *mut FILE;
    /// Returns the standard error file.
    pub safe fn stderr_file() -> *mut FILE;

    /// Computes the Cholesky factorization of a real symmetric positive definite matrix.
    #[doc(hidden)]
    pub unsafe fn dpotrf(uplo: *const i8, n: c_int, a: *mut f64, lda: c_int, info: *mut c_int);

    /// Multiplies lower or upper triangular matrix with vector.
    #[doc(hidden)]
    pub unsafe fn dtrmv(
        uplo: *const i8,
        transa: *const i8,
        diag: *const i8,
        n: c_int,
        a: *const f64,
        lda: c_int,
        x: *mut f64,
        incx: c_int,
    );

    /// Creates initialized tilemap for [`oefpil_tcm_diag_new`] or [`oefpil_tcm_blockdiag_new`].
    ///
    /// Arguments:
    ///
    ///   * Number of `variables` or tiles.
    pub safe fn oefpil_tilemap_diagtiles_new(variables: c_int) -> *mut c_int;

    /// Creates tiled covariance matrix of diagonal tiles on its diagonal.
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tile `map` created with <code>[oefpil_tilemap_diagtiles_new]\(variables\)</code>.
    pub unsafe fn oefpil_tcm_diag_new(
        samples: c_int,
        variables: c_int,
        map: *mut c_int,
    ) -> *mut f64;
    /// Sets `fields` of diagonal tile of `tcm` created with [`oefpil_tcm_diag_new`].
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tiled covariance matrix `tcm` created with
    ///     <code>[oefpil_tcm_diag_new](samples, variables, map)</code>.
    ///   * Tile `map` created with <code>[oefpil_tilemap_diagtiles_new]\(variables\)</code>.
    ///   * Tile `row_column` in `0..variables`.
    ///   * Tile `fields` to copy into diagonal tile at `row_column` of `tcm`. Number of `fields` is
    ///     `samples`.
    pub unsafe fn oefpil_tcm_diag_set_tile_diag(
        samples: c_int,
        variables: c_int,
        tcm: *mut f64,
        map: *mut c_int,
        row_column: c_int,
        fields: *const f64,
    );

    /// Creates tiled covariance matrix of diagonal or block tiles on its diagonal.
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tile `map` created with <code>[oefpil_tilemap_diagtiles_new]\(variables\)</code>.
    pub unsafe fn oefpil_tcm_blockdiag_new(
        samples: c_int,
        variables: c_int,
        map: *mut c_int,
    ) -> *mut f64;
    /// Sets `fields` of diagonal tile of `tcm` created with [`oefpil_tcm_blockdiag_new`].
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tiled covariance matrix `tcm` created with
    ///     <code>[oefpil_tcm_blockdiag_new](samples, variables, map)</code>.
    ///   * Tile `map` created with <code>[oefpil_tilemap_diagtiles_new]\(variables\)</code>.
    ///   * Tile `row_column` in `0..variables`.
    ///   * Tile `fields` to copy into diagonal tile at `row_column` of `tcm`. Number of `fields` is
    ///     `samples`.
    pub unsafe fn oefpil_tcm_blockdiag_set_tile_diag(
        samples: c_int,
        variables: c_int,
        tcm: *mut f64,
        map: *mut c_int,
        row_column: c_int,
        fields: *const f64,
    );
    /// Sets `fields` of symmetric block tile of `tcm` created with [`oefpil_tcm_blockdiag_new`].
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tiled covariance matrix `tcm` created with
    ///     <code>[oefpil_tcm_blockdiag_new](samples, variables, map)</code>.
    ///   * Tile `map` created with <code>[oefpil_tilemap_diagtiles_new]\(variables\)</code>.
    ///   * Tile `row_column` in `0..variables`.
    ///   * Tile `fields` to copy into symmetric block tile at `row_column` of `tcm`. Number of
    ///     `fields` is `samples * (samples + 1) / 2`. `fields` is a row-major [triangular slice] of
    ///     the lower triangle where the upper triangle is automatically constructed.
    ///
    /// [triangular slice]: https://en.wikipedia.org/wiki/Triangular_array
    pub unsafe fn oefpil_tcm_blockdiag_set_tile_half(
        samples: c_int,
        variables: c_int,
        tcm: *mut f64,
        map: *mut c_int,
        row_column: c_int,
        fields: *const f64,
    );
    /// Sets `fields` of symmetric block tile of `tcm` created with [`oefpil_tcm_blockdiag_new`].
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tiled covariance matrix `tcm` created with
    ///     <code>[oefpil_tcm_blockdiag_new](samples, variables, map)</code>.
    ///   * Tile `map` created with <code>[oefpil_tilemap_diagtiles_new]\(variables\)</code>.
    ///   * Tile `row_column` in `0..variables`.
    ///   * Tile `fields` to copy into symmetric block tile at `row_column` of `tcm`. Number of
    ///     `fields` is `samples.pow(2)`. As the block tile is on the diagonal, `fields` must be
    ///     symmetric (not validated) where row-major and column-major order coincide.
    pub unsafe fn oefpil_tcm_blockdiag_set_tile_full(
        samples: c_int,
        variables: c_int,
        tcm: *mut f64,
        map: *mut c_int,
        row_column: c_int,
        fields: *const f64,
    );

    /// Creates initialized tilemap for [`oefpil_tcm_diags_new`] or [`oefpil_tcm_full_new`].
    ///
    /// Arguments:
    ///
    ///   * Number of `variables` or tiles.
    pub safe fn oefpil_tilemap_alltiles_new(bn: c_int) -> *mut c_int;

    /// Creates tiled covariance matrix of diagonal tiles.
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tile `map` created with <code>[oefpil_tilemap_alltiles_new]\(variables\)</code>.
    pub unsafe fn oefpil_tcm_diags_new(
        samples: c_int,
        variables: c_int,
        map: *mut c_int,
    ) -> *mut f64;
    /// Sets `fields` of diagonal tile of `tcm` created with [`oefpil_tcm_diags_new`].
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tiled covariance matrix `tcm` created with
    ///     <code>[oefpil_tcm_diags_new](samples, variables, map)</code>.
    ///   * Tile `map` created with <code>[oefpil_tilemap_diagtiles_new]\(variables\)</code>.
    ///   * Tile `row` in `0..variables`.
    ///   * Tile `column` in `0..variables`.
    ///   * Tile `fields` to copy into diagonal tile at `row` and `column` of `tcm`. Number of
    ///     `fields` is `samples`.
    pub unsafe fn oefpil_tcm_diags_set_tile_diag(
        samples: c_int,
        variables: c_int,
        tcm: *mut f64,
        map: *mut c_int,
        row: c_int,
        column: c_int,
        fields: *const f64,
    );

    /// Creates tiled covariance matrix of diagonal or block tiles.
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tile `map` created with <code>[oefpil_tilemap_alltiles_new]\(variables\)</code>.
    pub unsafe fn oefpil_tcm_full_new(
        samples: c_int,
        variables: c_int,
        map: *mut c_int,
    ) -> *mut f64;
    /// Sets `fields` of diagonal tile of `tcm` created with [`oefpil_tcm_blockdiag_new`].
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tiled covariance matrix `tcm` created with
    ///     <code>[oefpil_tcm_blockdiag_new](samples, variables, map)</code>.
    ///   * Tile `map` created with <code>[oefpil_tilemap_diagtiles_new]\(variables\)</code>.
    ///   * Tile `row` in `0..variables`.
    ///   * Tile `column` in `0..variables`.
    ///   * Tile `fields` to copy into diagonal tile at `row` and `column` of `tcm`. Number of
    ///     `fields` is `samples`.
    pub unsafe fn oefpil_tcm_full_set_tile_diag(
        samples: c_int,
        variables: c_int,
        tcm: *mut f64,
        map: *mut c_int,
        row: c_int,
        column: c_int,
        fields: *const f64,
    );
    /// Sets `fields` of symmetric block tile of `tcm` created with [`oefpil_tcm_blockdiag_new`].
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tiled covariance matrix `tcm` created with
    ///     <code>[oefpil_tcm_blockdiag_new](samples, variables, map)</code>.
    ///   * Tile `map` created with <code>[oefpil_tilemap_diagtiles_new]\(variables\)</code>.
    ///   * Tile `row` in `0..variables`.
    ///   * Tile `column` in `0..variables`.
    ///   * Tile `fields` to copy into symmetric block tile at `row_column` of `tcm`. Number of
    ///     `fields` is `samples * (samples + 1) / 2`. `fields` is a row-major [triangular slice] of
    ///     the lower triangle where the upper triangle is automatically constructed.
    ///
    /// [triangular slice]: https://en.wikipedia.org/wiki/Triangular_array
    pub unsafe fn oefpil_tcm_full_set_tile_half(
        samples: c_int,
        variables: c_int,
        tcm: *mut f64,
        map: *mut c_int,
        row: c_int,
        column: c_int,
        fields: *const f64,
    );
    /// Sets `fields` of row-major block tile of `tcm` created with [`oefpil_tcm_full_new`].
    ///
    /// Arguments:
    ///
    ///   * Number of `samples` per variable.
    ///   * Number of `variables` or tiles.
    ///   * Tiled covariance matrix `tcm` created with
    ///     <code>[oefpil_tcm_full_new](samples, variables, map)</code>.
    ///   * Tile `map` created with <code>[oefpil_tilemap_diagtiles_new]\(variables\)</code>.
    ///   * Tile `row` in `0..variables`.
    ///   * Tile `column` in `0..variables`.
    ///   * Tile `fields` to copy into row-major block tile at `row` and `column` of `tcm`. Number
    ///     of `fields` is `samples.pow(2)`. For block tiles on the diagonal, `fields` must be
    ///     symmetric (not validated) where row-major and column-major order coincide.
    pub unsafe fn oefpil_tcm_full_set_tile_full(
        samples: c_int,
        variables: c_int,
        tcm: *mut f64,
        map: *mut c_int,
        row: c_int,
        column: c_int,
        fields: *const f64,
    );

    /// Fits the initial estimate of the model's parameter to the data sample of its variable.
    ///
    /// Arguments:
    ///
    ///   * Function pointer `evaluate` of signature [`Evaluate`] expecting `data` as 1st argument.
    ///   * `data` passed to `evaluate` as 1st argument, see [`Evaluate`].
    ///   * Whether the model `is_implicit` (`1`) or explicit (`0`).
    ///   * Number of `parameters`.
    ///   * `parameter_mean` vector.
    ///   * `parameter_covariance` matrix of `parameter_mean` vector.
    ///   * Number of `samples` per variable.
    ///   * Number of independent `x_variables`.
    ///   * `x_sample` of `x_variables` (sample-major).
    ///   * `y_sample` if not `is_implicit` else [`core::ptr::null`].
    ///   * `x_mean` of `x_variables` (sample-major).
    ///   * `y_mean` if not `implicit` else [`core::ptr::null`].
    ///   * Tiled `covariance` matrix of `x_sample` and `y_sample`.
    ///   * `covariance_mode` of [`Mode`].
    ///   * `covariance_map` of [`Mode`] per tile.
    ///   * Convergence `iteration_limit`.
    ///   * Convergence `tolerance`.
    ///   * Log `verbosity`, see [`Verbosity`].
    ///   * `logfile`, see [`Verbosity`].
    ///   * `chi_squared` (statistics).
    ///   * Convergence `criterion`, see [`Criterion`].
    ///   * Result `info` with `1` for success, `2` for `iteration_limit`, `3` for numerical error.
    ///   * Number of `iterations` until convergence [`Criterion`] has been reached.
    ///   * `chi_squared_reduced` (statistics).
    ///   * Whether `covariance` is `relative` and rescaled by `chi_squared_reduced` or absolute.
    ///   * `chi_squared_p_value` (statistics).
    pub unsafe fn oefpil(
        evaluate: Evaluate,
        data: *mut c_void,
        is_implicit: c_int,
        parameters: c_int,
        parameter_mean: *mut f64,
        parameter_covariance: *mut f64,
        samples: c_int,
        x_variables: c_int,
        x_sample: *const f64,
        y_sample: *const f64,
        x_mean: *mut f64,
        y_mean: *mut f64,
        covariance: *const f64,
        covariance_mode: c_int,
        covariance_map: *const c_int,
        iteration_limit: c_int,
        tolerance: f64,
        verbosity: c_int,
        logfile: *mut FILE,
        chi_squared: *mut f64,
        criterion: c_int,
        info: *mut c_int,
        iterations: *mut c_int,
        chi_squared_reduced: *mut f64,
        relative: bool,
        chi_squared_p_value: *mut f64,
    );
}

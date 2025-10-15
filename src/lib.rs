//! Optimum Estimate of Function Parameters by Iterated Linearization (OEFPIL)
//!
//! Algorithm for nonlinear function fitting to data with errors in variables where correlation,
//! both within variables and among variables, might be present. In principle, OEFPIL can be
//! employed for fitting both explicit and implicit functions of any number of variables.
//! Importantly, apart from the parameter estimates, OEFPIL also yields their covariance matrix,
//! required for further analyses.[^1] Common methods such as ordinary nonlinear least squares are
//! not capable of treating general uncertainties and correlations in both dependent and independent
//! variables. A new computation method for nonlinear curve fitting to data with a general
//! covariance structure is introduced. Numerical simulations show that the new method yields
//! parameter estimates in agreement with other methods for simple covariance structures. The
//! obtained uncertainty estimates are in agreement with Monte Carlo studies.[^2]
//!
//! The most notable features of OEFPIL are as follows:[^2]
//!
//!   * The OEFPIL algorithm features an improved estimation of uncertainties and covariance
//!     matrices of the fitted parameters which can be further propagated to physical output
//!     variables.
//!   * More general covariance matrices for input data can be applied, although so far they do not
//!     seem to lead to significant shifts in the physical output variables. However, estimates of
//!     uncertainties of these are obviously affected.
//!   * Unlike ordinary least squares OEFPIL is not sensitive to the choice of dependent and
//!     independent variables.
//!   * Unlike non-linear least squares and orthogonal distance regression, a reasonable uncertainty
//!     estimate can be deduced directly without resorting to Monte Carlo studies, evading thus long
//!     computation times.
//!   * A possible disadvantage of OEFPIL is a higher sensitivity to initial estimates of parameters
//!     and the covariance matrix of the input data. However, so far problems have been encountered
//!     only for exaggerated estimates of the input covariance matrix.
//!
//! [^1]: R. Šlesinger, A. C. Campbell, Z. Geršlová, V. Šindlář, and G. Wimmer, “OEFPIL: New Method
//! and Software Tool for Fitting Nonlinear Functions to Correlated Data With Errors in Variables”,
//! [2023 14th International Conference on Measurement, 126-129
//! (2023)](https://doi.org/10.23919/MEASUREMENT59122.2023.10164444).
//!
//! [^2]: A. C. Campbell, Z. Geršlová, V. Šindlář, R. Šlesinger, and G. Wimmer, “New framework for
//! nanoindentation curve fitting and measurement uncertainty estimation”, [Precision Engineering
//! 85, 166-173 (2024)](https://doi.org/10.1016/j.precisioneng.2023.10.001).
//!
//! [`oefpil-sys`]: https://docs.rs/oefpil-sys
//! [C/Fortran library]: https://gitlab.com/cmi6014/oefpil
//!
//! # Introduction
//!
//! This section attempts to introduce some mathematical terminology of probability and statistics
//! (e.g., random variables and their realizations as observation samples) which are used allover in
//! this documentation. Conventionally, to avoid confusion, upper case letters denote random
//! variables whereas the corresponding lower case letters denote their realizations. An omitted
//! subscript or superscript index denotes the respective vector comprising such indexed elements.
//! If all this is unfamiliar to you, jump to the API [overview](#overview) or peek into the
//! [examples](#examples) and come back to look up a particular term or symbol in question.
//!
//! Let $`V\equiv (X,Y)=(X^1,\dots,X^{|X|},Y^1)`$ be a vector of $`|X|`$ independent [random
//! variables] $`X \equiv (X^1,\dots,X^{|X|})`$ and a single dependent random variable $`Y \equiv
//! (Y^1)`$ where $`|V|=|X|+|Y|`$ with $`|\cdot|`$ denoting the respective number of variables. Let
//! $`v`$ be an $`|V|\times|v|`$ sample matrix with $`|v|`$ [multivariate observations]
//! $`v_i=(v^1_i,...,v^I_i)`$ as [realizations] $`v^I_i\sim V^I_i\equiv V^I(\omega^I_i)`$ at
//! $`\omega^I_i\sim\Omega^I`$ from their sample spaces $`\Omega^I`$ where $`I\in[1,|V|]`$ and
//! $`i\in[1,|v|]`$. Equations (1, 2) show a variable per row and an observation per column whereas
//! equation (3) shows the sample matrix $`v`$ split into [dependent and independent] sample vector
//! $`y`$ and sample matrix $`x`$.
//!
//! [dependent and independent]: https://en.wikipedia.org/wiki/Dependent_and_independent_variables
//! [sequence]: https://en.wikipedia.org/wiki/Sequence
//! [random variables]: https://en.wikipedia.org/wiki/Random_variable
//! [multivariate observations]: https://en.wikipedia.org/wiki/Multivariate_statistics
//! [realizations]: https://en.wikipedia.org/wiki/Realization_(probability)
//!
//! ```math
//! \gdef\c#1#2{\textup{Cov}(#1,#2)}\gdef\C#1#2{\textup{Cor}(#1,#2)}
//! \gdef\s{\sigma}\gdef\d{\delta}\gdef\e{\epsilon}
//! \begin{align}
//! v & \equiv \begin{pmatrix}
//! & v^1_1     & \dots  & v^1_{|n|}     & \\[1em]
//! & \vdots    & \ddots & \vdots        & \\[1em]
//! & v^{|V|}_1 & \dots  & v^{|V|}_{|n|} & \end{pmatrix}
//! \\ \nonumber \\
//!   & = \begin{pmatrix}
//! & x^1_1     & \dots  & x^1_{|v|}     & \\[1em]
//! & \vdots    & \ddots & \vdots        & \\[1em]
//! & x^{|X|}_1 & \dots  & x^{|X|}_{|v|} & \\[1em]
//! & y^1_1     & \dots  & y^1_{|v|}     & \end{pmatrix}
//! \\ \nonumber \\
//!   & = \begin{pmatrix}
//! & x & \\[1em]
//! & y & \end{pmatrix}
//! \end{align}
//! ```
//!
//! The observation $`v^I_i=\mu^I_i+\e^I_i`$ deviates by $`\e^I_i`$ from the usually unknown true
//! value $`\mu^I_i`$. The [observational error] $`\e^I_i=(L\d)^I_i`$ is a realization
//! $`\d^I_i\sim\Delta^I_i\equiv\mathcal{N}(0,1)`$ of a standard [normal distribution] possibly
//! correlated by the [Cholesky decomposition] $`L`$ of the [covariance matrix] $`\c{V}{V}=LL^T`$.
//! Equations (4, 5) show the covariance matrix tiled by variables with its tiles $`\c{V^R}{V^C}`$
//! at row $`R`$ and column $`C`$ given by equation (6). The standard deviations $`\s^R_r`$ and
//! $`\s^C_c`$ are possibly correlated within and among variables by the [Pearson correlation
//! coefficients] $`\C{V^R_r}{V^C_c}`$ where $`r,c\in[1,|v|]`$.
//!
//! [observational error]: https://en.wikipedia.org/wiki/Observational_error
//! [normal distribution]: https://en.wikipedia.org/wiki/Normal_distribution
//! [Cholesky decomposition]: https://en.wikipedia.org/wiki/Cholesky_decomposition
//! [standard deviation]: https://en.wikipedia.org/wiki/Standard_deviation
//! [mean]: https://en.wikipedia.org/wiki/Mean#Mean_of_a_probability_distribution
//! [covariance matrix]: https://en.wikipedia.org/wiki/Covariance_matrix
//! [Pearson correlation coefficients]:
//! https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
//!
//! ```math
//! \begin{align}
//! \c{V}{V}     & \equiv \begin{pmatrix}
//! & \c{V^1}{V^1}      & \dots  & \c{V^1}{V^{|V|}}     & \\[1em]
//! & \vdots            & \ddots & \vdots               & \\[1em]
//! & \c{V^{|V|}}{V^1}  & \dots  & \c{V^{|V|}}{V^{|V|}} &
//! \end{pmatrix} \\ \nonumber \\
//!              & = \begin{pmatrix}
//! & \c{X^1}{X^1}      & \dots  & \c{X^1}{X^{|X|}}     & \c{X^1}{Y^1}     & \\[1em]
//! & \vdots            & \ddots & \vdots               & \vdots           & \\[1em]
//! & \c{X^{|X|}}{X^1}  & \dots  & \c{X^{|X|}}{X^{|X|}} & \c{X^{|X|}}{Y^1} & \\[1em]
//! & \c{Y^1}{X^1}      & \dots  & \c{Y^1}{X^{|X|}}     & \c{Y^1}{Y^1}     &
//! \end{pmatrix} \\ \nonumber \\
//! \c{V^R}{V^C} & \equiv \begin{pmatrix}
//! & \C{V^R_1}{V^C_1}\s^R_1\s^C_1         & \dots  & \C{V^R_1}{V^C_{|v|}}\s^R_1\s^C_{|v|}         &
//! \\[1em]
//! & \vdots                               & \ddots & \vdots                                       &
//! \\[1em]
//! & \C{V^R_{|v|}}{V^C_1}\s^R_{|v|}\s^C_1 & \dots  & \C{V^R_{|v|}}{V^C_{|v|}}\s^R_{|v|}\s^C_{|v|} &
//! \end{pmatrix}\ \forall\ R,C \in [1,|V|]
//! \end{align}
//! ```
//!
//! Given an [explicit or implicit model] $`Y=f(X;P)`$ or $`0=f(X;P)`$ with parameter vector
//! $`P\equiv(P^1,...,P^{|P|})`$ and derivative vectors in variables $`df/dX`$ and in parameters
//! $`df/dP`$, its initial estimates of the sample mean vector $`\mu`$ and parameter mean vector
//! $`\lambda`$ are fitted to the sample matrix $`v`$ with covariance matrix $`\c{V}{V}`$. The
//! optimum estimates of the sample mean vector $`\hat\mu`$ and parameter mean vector
//! $`\hat\lambda`$ with covariance matrix $`\c{P}{P}`$ are yielded when the convergence criterion
//! is satisfied. Optionally, the parameter standard deviation vector $`\varsigma`$ and [correlation
//! matrix] $`\C{P}{P}`$ are yielded as well.
//!
//! [explicit or implicit model]: https://en.wikipedia.org/wiki/Explicit_and_implicit_methods
//! [correlation matrix]: https://en.wikipedia.org/wiki/Correlation#Correlation_matrices
//!
//! # Overview
//!
//! This crate provides a safe API to the [`oefpil-sys`] crate (see its system requirements) which
//! statically links to the [C/Fortran library]. The safe Rust API comprises five structures *with
//! public members* of which many can be left at their [`Default`] values:
//!
//!   * [`Algorithm`] for configuring the fitting algorithm (e.g., convergence [`Criterion`]).
//!   * [`Model`] for defining the model comprising the function $`Y=f(X;P)`$ or $`0=f(X;P)`$ and
//!     its derivative vectors in variables $`df/dX`$ and in parameters $`df/dP`$.
//!   * [`Variable`] for assigning the multivariate sample matrix $`v`$, the sample mean matrix
//!     $`\mu`$, and the [`Covariance`] matrix $`\c{V}{V}`$. The initial estimates $`\mu`$ will be
//!     overwritten by the optimum estimates $`\hat\mu`$.
//!   * [`Parameter`] for assigning the initial estimates of parameter mean vector $`\lambda`$ and
//!     output slices for the covariance matrix $`\c{P}{P}`$ and optionally for the deviation vector
//!     $`\varsigma`$ and for the correlation matrix $`\C{P}{P}`$. The initial estimates $`\lambda`$
//!     will be overwritten by the optimum estimates $`\hat\lambda`$.
//!   * [`Report`] for reporting the number of iterations and statistical information.
//!
//! Following feature gates are disabled by default:
//!
//!   * `random` -- Random [`Distribution::sample()`] possibly correlated by
//!     [`Covariance::with_decomposition()`].
//!   * <code>[oefpil-sys]/built-in</code> -- Statically linked Fortran dependencies.
#![cfg_attr(
    not(feature = "random"),
    doc = r"

[`Distribution::sample()`]:
https://docs.rs/oefpil/latest/oefpil/struct.Distribution.html#method.sample
"
)]
//!
//! # Examples
//!
//!   * [`sine.rs`] -- Explicit univariate biparametric sinusoidal model --
//!     Let $`V=(\chi,A)`$ comprise an independent variable $`\chi`$ and an dependent variable $`A`$
//!     and let $`P=(\chi_0,A_0)`$ be the phase and amplitude parameter of the explicit model
//!     $`A=A_0\sin(\chi+\chi_0)`$. Eight observations $`v_i\sim V_i`$ with correlation coefficients
//!     $`\C{V^I_i}{V^I_i}=0.5`$ among variables are sampled where $`i\in[1,8]`$ and $`I\in[1,2]`$.
//!
//! # Troubleshooting
//!
//! As outlined above, this algorithm is more sensitive to the initial estimate of the paramter mean
//! vector $`\lambda`$ and to the sample covariance matrix $`\c{V}{V}`$ than alternative algorithms
//! which also enables it to yield a more accurate estimate of the parameter covariance matrix
//! $`\c{P}{P}`$. Therefore, if the sample covariance matrix $`\c{V}{V}`$ is not well estimated, the
//! algorithm might abort with an [`AlgorithmError`]. The most promising troubleshooting steps are
//! as follows in the order mentioned:
//!
//!   1. Ensure all inputs satisfy [`f64::is_finite()`].
//!   2. Ensure all diagonal fields of the sample covariance matrix $`\c{V}{V}`$ are neither
//!      negative nor zero.
//!   3. Ensure to pass (co)variances and not standard deviations to [`Covariance::set_tile()`].
//!   4. Improve the initial estimate of the parameter mean vector $`\lambda`$.
//!   5. Improve the estimate of the sample covariance matrix $`\c{V}{V}`$. For a diagonal
//!      [`Covariance::with_unknown_scale()`], use a small fraction of the sample means $`v^I_i`$ as
//!      estimates for their standard deviations $`\sigma^I_i`$.
//!   6. Adjust the [`Algorithm`] settings.
//!
//! [`sine.rs`]: ../src/sine/sine.rs.html
//!
//! # Limitation
//!
//! Currently, constructing the [`Covariance`] matrix with correlations among variables (i.e., via
//! [`Covariance::new_diagonals()`] or [`Covariance::new_full()`] is limited to $`|V|=2`$. This
//! covers explicit univariate models $`V=(X^1,Y^1)`$ and implicit bivariate models $`V=(X^1,X^2)`$.
//! In contrast, constructing the [`Covariance`] matrix with correlations within but not among
//! variables (i.e., via [`Covariance::new_diagonal()`] or [`Covariance::new_block_diagonal()`]) is
//! already supported for any number of variables $`|V|`$.
//!
//! # Roadmap
//!
//! Get rid of [limitation](#limitation). On the long run, this crate will provide a pure Rust
//! implementation. The current Rust API won't necessarily need to change much. The [`oefpil-sys`]
//! dependency will become optional. This will allow to compare both implementations. Currently, the
//! [C/Fortran library] allocates working memory and its C API depends on the system type [`FILE`]
//! defined in the C standard library. This might change though. The Rust implementation will work
//! in [`no_std`] environments, probably by leveraging [`faer`] and [`dyn-stack`] internally. Not
//! exposing them would allow a stable release sooner rather than later.
//!
//! [`no_std`]: https://crates.io/categories/no-std
//! [`faer`]: https://docs.rs/faer
//! [`dyn-stack`]: https://docs.rs/dyn-stack

use derive_more::{Display, Error, From};
pub use oefpil_sys::{Criterion, Verbosity};
use oefpil_sys::{
    Evaluate, Mode, dpotrf,
    libc::{FILE, fclose, fflush, fopen},
    oefpil, oefpil_tcm_blockdiag_set_tile_diag, oefpil_tcm_blockdiag_set_tile_full,
    oefpil_tcm_blockdiag_set_tile_half, oefpil_tcm_diag_set_tile_diag,
    oefpil_tcm_diags_set_tile_diag, oefpil_tcm_full_set_tile_diag, oefpil_tcm_full_set_tile_full,
    oefpil_tcm_full_set_tile_half, stderr_file, stdout_file,
};
use std::{
    ffi::{CString, c_int, c_void},
    fs::OpenOptions,
    io::{self, BufRead, BufReader, ErrorKind},
    num::TryFromIntError,
    path::Path,
    ptr,
    slice::{from_raw_parts, from_raw_parts_mut},
};

/// Settings of the fitting algorithm.
///
/// Following environment variables are respected:
///
///   * `OEFPIL_DEBUG` writes intermediate matrices to text files in working directory if defined.
///   * `OEFPIL_PRINTFORMAT` modifies the floating-point formatter (e.g., `"8e"`, `"5f"`).
///   * `OEFPIL_MAXIT` overrides [`Algorithm::iteration_limit`].
#[derive(Debug, Clone, Copy)]
pub struct Algorithm {
    /// Convergence criterion.
    ///
    /// Default is [`Criterion::RelPOrAbsPAndRelXOrAbsX`].
    pub criterion: Criterion,
    /// Convergence tolerance.
    ///
    /// Default is <code>[f64::EPSILON].[powf](f64::powf)(2.0 / 3.0)</code>.
    pub tolerance: f64,
    /// Iteration limit.
    ///
    /// Default is `100`.
    pub iteration_limit: usize,
    /// Log verbosity.
    ///
    /// Default is [`Verbosity::Silent`].
    pub verbosity: Verbosity,
}

impl Default for Algorithm {
    fn default() -> Self {
        Self {
            criterion: Criterion::default(),
            tolerance: f64::EPSILON.powf(2.0 / 3.0),
            iteration_limit: 100,
            verbosity: Verbosity::default(),
        }
    }
}

impl Algorithm {
    /// Fits the initial estimates of the [`Model`]'s [`Parameter`] to its [`Variable`] sample.
    ///
    /// The initial estimates of the [`Model`]'s [`Variable::mean`] matrix and [`Parameter::mean`]
    /// vector are fitted to the [`Variable::sample`] matrix yielding the optimum estimates of
    /// [`Variable::mean`] and [`Parameter::mean`] along with the [`Parameter::covariance`] matrix
    /// and optionally the [`Parameter::deviation`] vector and [`Parameter::correlation`] matrix.
    ///
    /// Returns the fitting [`Report`] on success.
    ///
    /// # Errors
    ///
    /// Returns [`OefpilError`] on failure.
    pub fn fit(
        &self,
        model: Model,
        variable: &mut Variable,
        parameter: &mut Parameter,
        logfile: Logfile,
    ) -> Result<Report, OefpilError> {
        model.validate()?;
        variable.validate(model)?;
        parameter.validate(model)?;
        let x_variables = model.dfdx.len();
        let samples = variable.covariance.samples;
        let parameters = parameter.mean.len();
        let x_fields = samples * x_variables;
        let mut x_mean;
        let mut y_mean;
        let x_mean = if variable.mean.is_empty() {
            x_mean = variable.sample[..x_fields].to_vec();
            x_mean.as_mut_ptr()
        } else {
            variable.mean.as_mut_ptr()
        };
        let y_mean = if variable.mean.is_empty() {
            y_mean = variable.sample[x_fields..].to_vec();
            y_mean.as_mut_ptr()
        } else {
            variable.mean[x_fields..].as_mut_ptr()
        };
        let evaluate: Evaluate = Some(evaluate);
        let mut data = Data {
            model,
            x: &mut vec![0.0; x_variables],
        };
        let mut info = 0;
        let mut iterations = 0;
        let mut report = Report {
            degrees_of_freedom: samples
                .checked_sub(parameters)
                .filter(|&dof| dof >= usize::from(!variable.covariance.unknown_scale))
                .ok_or(VariableError::InsufficientDegreesOfFreedom)?,
            ..Report::default()
        };
        unsafe {
            let file = logfile.open()?;
            oefpil(
                evaluate,
                ptr::from_mut(&mut data).cast(),
                model.implicit.into(),
                parameters.try_into()?,
                parameter.mean.as_mut_ptr(),
                parameter.covariance.as_mut_ptr(),
                samples.try_into()?,
                x_variables.try_into()?,
                variable.sample[..x_fields].as_ptr(),
                variable.sample[x_fields..].as_ptr(),
                x_mean,
                y_mean,
                variable.covariance.fields.as_ptr(),
                variable.covariance.mode as c_int,
                variable.covariance.tiles.as_ptr(),
                self.iteration_limit.try_into()?,
                self.tolerance,
                self.verbosity as i32,
                file,
                ptr::from_mut(&mut report.chi_squared),
                self.criterion as i32,
                &raw mut info,
                &raw mut iterations,
                &raw mut report.chi_squared_reduced,
                variable.covariance.unknown_scale,
                ptr::null_mut(),
            );
            let log = if matches!(logfile, Logfile::Custom(_path)) {
                fclose(file)
            } else {
                fflush(file)
            };
            if log != 0 {
                Err(io::Error::last_os_error())?;
            }
        };
        report.iterations = iterations.try_into()?;
        if !parameter.deviation.is_empty() {
            for row_column in 0..parameters {
                let diagonal = row_column * parameters + row_column;
                parameter.deviation[row_column] = parameter.covariance[diagonal].sqrt();
            }
        }
        if !parameter.correlation.is_empty() {
            for row in 0..parameters {
                for column in 0..parameters {
                    let field = row * parameters + column;
                    parameter.correlation[field] = parameter.covariance[field]
                        / (parameter.deviation[row] * parameter.deviation[column]);
                }
            }
        }
        match info {
            1 => Ok(report),
            2 => Err(AlgorithmError::IterationLimit.into()),
            3 => Err(AlgorithmError::NumericalError.into()),
            _ => unreachable!("unknown error"),
        }
    }
}

struct Data<'a> {
    model: Model<'a>,
    x: &'a mut [f64],
}

#[allow(clippy::similar_names)]
unsafe extern "C" fn evaluate(
    data: *mut c_void,
    samples: c_int,
    x_sample: *const f64,
    parameter: *const f64,
    fx: *mut f64,
    dfdx: *mut f64,
    dfdp: *mut f64,
) {
    let Data { model, x } = unsafe { &mut *data.cast::<Data>() };
    let samples: usize = samples.try_into().unwrap();
    let x_variables = model.dfdx.len();
    let parameters = model.dfdp.len();
    let x_sample = unsafe { from_raw_parts(x_sample, samples * x_variables) };
    let p = unsafe { from_raw_parts(parameter, parameters) };
    let fx = unsafe { from_raw_parts_mut(fx, samples) };
    let dfdx = unsafe { from_raw_parts_mut(dfdx, samples * x_variables) };
    let dfdp = unsafe { from_raw_parts_mut(dfdp, samples * parameters) };
    for sample in 0..samples {
        for x_variable in 0..x_variables {
            x[x_variable] = x_sample[x_variable * samples + sample];
        }
        fx[sample] = (model.fx)(x, p);
        for x_variable in 0..x_variables {
            dfdx[x_variable * samples + sample] = (model.dfdx[x_variable])(x, p);
        }
        for parameter in 0..parameters {
            dfdp[parameter * samples + sample] = (model.dfdp[parameter])(x, p);
        }
    }
}

/// Logfile of either standard output/error or custom file path.
#[derive(Debug, Clone, Copy, Default)]
pub enum Logfile<'a> {
    /// Standard output file descriptor/handle.
    #[default]
    StdOut,
    /// Standard error file descriptor/handle.
    StdErr,
    /// Custom file path to be opened.
    Custom(&'a Path),
}

impl<'a> Logfile<'a> {
    fn open(&'a self) -> Result<*mut FILE, io::Error> {
        unsafe {
            match self {
                Self::StdOut => Ok(stdout_file()),
                Self::StdErr => Ok(stderr_file()),
                Self::Custom(path) => {
                    let path = CString::new(path.as_os_str().as_encoded_bytes()).unwrap();
                    let file = fopen(path.as_c_str().as_ptr(), c"w".as_ptr());
                    if file.is_null() {
                        Err(io::Error::last_os_error())
                    } else {
                        Ok(file)
                    }
                }
            }
        }
    }
}

/// Model comprising the function $`Y=f(X;P)`$ or $`0=f(X;P)`$ and its derivative vectors in
/// variables $`df/dX`$ and parameters $`df/dP`$.
#[derive(Debug, Clone, Copy)]
#[allow(clippy::similar_names)]
pub struct Model<'a> {
    /// Explicit or implicit function $`Y=f(X;P)`$ or $`0=f(X;P)`$.
    pub fx: Expression,
    /// Derivative vector $`df/dX`$ in independent variables $`X`$.
    pub dfdx: &'a [Expression],
    /// Derivative vector $`df/dP`$ in parameters $`P`$.
    pub dfdp: &'a [Expression],
    /// Whether the model is explicit or `implicit`.
    pub implicit: bool,
}

impl Model<'_> {
    /// Number of independent as well as dependent variables $`|V|\equiv|X|+|Y|`$.
    #[must_use]
    #[inline]
    pub fn variables(&self) -> usize {
        self.dfdx.len() + usize::from(!self.implicit)
    }
    /// Ensures non-empty derivative slices.
    ///
    /// Called in [`Algorithm::fit()`].
    ///
    /// # Errors
    ///
    /// Returns [`ModelError`] on failure.
    pub fn validate(&self) -> Result<(), ModelError> {
        if self.dfdx.is_empty() {
            return Err(ModelError::EmptyVariableDerivative);
        }
        if self.dfdp.is_empty() {
            return Err(ModelError::EmptyParameterDerivative);
        }
        Ok(())
    }
}

/// Mathematical expression in independent variables `x` and parameters `p`.
pub type Expression = fn(x: &[f64], p: &[f64]) -> f64;

/// Distribution comprising sample matrix $`v`$, sample mean matrix $`\mu`$, and sample
/// [`Covariance`] matrix $`\gdef\c#1#2{\textup{Cov}(#1,#2)}\c{V}{V}`$.
#[derive(Debug)]
#[cfg(feature = "random")]
pub struct Distribution<'a> {
    /// Sample matrix $`v`$ as output.
    ///
    /// Expects a row-major slice, that is, for multivariate models the samples of the first
    /// independent variable are appended by the samples of the second independent variable, etc.,
    /// followed by the samples of the dependent variable if the model is explicit.
    pub sample: &'a mut [f64],
    /// Sample mean (i.e., true value) matrix $`\hat\mu`$ as input.
    ///
    /// Expects a row-major slice, that is, for multivariate models the means of the first
    /// independent variable are appended by the means of the second independent variable, etc.,
    /// followed by the means of the dependent variable if the model is explicit.
    pub mean: &'a [f64],
    /// Sample covariance matrix $`\c{V}{V}`$.
    ///
    /// The covariance matrix is tiled by variables, that is, for multivariate models the
    /// covariances of the first independent variable are appended by the covariances of the second
    /// independent variable, etc., followed by the dependent variable if the model is explicit.
    pub covariance: &'a Covariance,
}

#[cfg(feature = "random")]
impl Distribution<'_> {
    /// Samples `random` observations possibly correlated by [`Covariance::with_decomposition()`].
    ///
    /// This facilitates a [Monte Carlo simulation] to confirm the [`Parameter::covariance`] matrix.
    ///
    /// [Monte Carlo simulation]: https://en.wikipedia.org/wiki/Monte_Carlo_method
    ///
    /// # Errors
    ///
    /// Returns [`OefpilError`] on failure.
    pub fn sample(&mut self, random: &mut dyn rand::RngCore) -> Result<(), OefpilError> {
        self.validate()?;
        let length = self.sample.len().try_into()?;
        // Uncorrelated observational unit errors (mean of `0` and deviation of `1`).
        self.sample
            .iter_mut()
            .zip(rand::Rng::sample_iter(random, rand_distr::StandardNormal))
            .for_each(|(sample, random): (&mut f64, f64)| *sample = random);
        // Correlated observational errors.
        unsafe {
            oefpil_sys::dtrmv(
                c"L".as_ptr(),
                c"N".as_ptr(),
                c"N".as_ptr(),
                length,
                self.covariance.decomposition.as_ptr(),
                length,
                self.sample.as_mut_ptr(),
                1,
            );
        };
        // Observational sample.
        self.sample
            .iter_mut()
            .zip(self.mean)
            .for_each(|(sample, mean)| {
                *sample += mean;
            });
        Ok(())
    }
    /// Ensures slice lengths agree and [`Covariance::with_decomposition()`].
    ///
    /// Called in [`Self::sample()`].
    ///
    /// # Errors
    ///
    /// Returns [`DistributionError`] on failure.
    pub fn validate(&self) -> Result<(), DistributionError> {
        if self.covariance.decomposition.is_empty() {
            return Err(DistributionError::MissingDecomposition);
        }
        let fields = self.covariance.order();
        if self.sample.len() != fields {
            return Err(DistributionError::MismatchingSample);
        }
        if self.mean.len() != fields {
            return Err(DistributionError::MismatchingMean);
        }
        Ok(())
    }
}

/// Variable comprising sample matrix $`v`$, sample mean matrix $`\mu`$, and sample [`Covariance`]
/// matrix $`\gdef\c#1#2{\textup{Cov}(#1,#2)}\c{V}{V}`$.
#[derive(Debug)]
pub struct Variable<'a> {
    /// Sample matrix $`v`$.
    ///
    /// Expects a row-major slice, that is, for multivariate models the samples of the first
    /// independent variable are appended by the samples of the second independent variable, etc.,
    /// followed by the samples of the dependent variable if the model is explicit.
    pub sample: &'a [f64],
    /// Sample mean (i.e., true value) matrix $`\mu`$ as input and $`\hat\mu`$ as output.
    ///
    /// An empty slice implies [`Self::sample`] whereas a non-empty slice serves as input for the
    /// initial estimates and as output for the optimum estimates.
    ///
    /// Expects a row-major slice, that is, for multivariate models the means of the first
    /// independent variable are appended by the means of the second independent variable, etc.,
    /// followed by the means of the dependent variable if the model is explicit.
    pub mean: &'a mut [f64],
    /// Sample covariance matrix $`\c{V}{V}`$.
    ///
    /// The covariance matrix is tiled by variables, that is, for multivariate models the
    /// covariances of the first independent variable are appended by the covariances of the second
    /// independent variable, etc., followed by the dependent variable if the model is explicit.
    pub covariance: &'a Covariance,
}

impl Variable<'_> {
    /// Ensures slice lengths agree with `model`.
    ///
    /// Called in [`Algorithm::fit()`].
    ///
    /// # Errors
    ///
    /// Returns [`VariableError`] on failure.
    pub fn validate(&self, model: Model) -> Result<(), VariableError> {
        if self.covariance.variables != model.variables() {
            return Err(VariableError::MismatchingCovariance);
        }
        let fields = self.covariance.samples * self.covariance.variables;
        if self.sample.len() != fields {
            return Err(VariableError::MismatchingSample);
        }
        if !self.mean.is_empty() && self.mean.len() != fields {
            return Err(VariableError::MismatchingMean);
        }
        Ok(())
    }
}

/// $`|V||v|\times|V||v|`$ covariance matrix $`\gdef\c#1#2{\textup{Cov}(#1,#2)}\c{V}{V}`$ with
/// $`|v|\times|v|`$ tiles $`\c{V^R}{V^C}`$ at row $`R`$ column $`C`$ where $`R,C\in[1,|V|]`$.
///
/// For certain problems, many fields are zero. For instance, when a problem has only variances, the
/// only non-zero fields are on the diagonal. When a problem has only variances and covariances
/// within variables, the non-zero fields make up a block-diagonal matrix. When a problem has
/// covariances among variables but not within variables, the non-zero fields are on the diagonal of
/// each tile. The different tiling modes allow to reduce memory allocations for many fields which
/// are known to be zero for certain problems. Only non-zero tiles are allocated and a tilemap
/// encodes whether a tile is a diagonal or block tile. The algorithm leverages the tiling modes as
/// well, allowing it to make mathematical shortcuts.
///
/// There is a constructor for each tiling mode:
///
///   * [`Self::new_diagonal()`]
///   * [`Self::new_block_diagonal()`]
///   * [`Self::new_diagonals()`]
///   * [`Self::new_full()`]
///
/// All constructors require the number of samples $`|v|`$ and the number of variables $`|V|`$. The
/// fields of the matrix are set per tile via [`Self::set_tile()`] or [`Self::with_tile()`] which
/// deduce whether a tile is a diagonal or block tile by the number of provided fields. When the
/// (co)variances are partially or fully unknown, communicate it via [`Self::with_unknown_scale()`].
#[derive(Debug, Clone)]
pub struct Covariance {
    mode: Mode,
    tiles: Vec<c_int>,
    fields: Vec<f64>,
    decomposition: Vec<f64>,
    samples: usize,
    variables: usize,
    unknown_scale: bool,
}

impl Covariance {
    /// Constructs a diagonal matrix (i.e, diagonal tiles on its diagonal).
    ///
    /// Allocates $`|V|`$ diagonal tiles on its diagonal with $`|v|`$ fields per tile.
    ///
    /// ```math
    /// \gdef\c#1#2{\textup{Cov}(#1,#2)}
    /// \c{V}{V} = \begin{pmatrix}
    /// & \diagdown &        &           & \\[1em]
    /// &           & \ddots &           & \\[1em]
    /// &           &        & \diagdown & \end{pmatrix}
    /// ```
    ///
    /// # Panics
    ///
    ///   * Asserts `samples` and `variables` to be non-zero.
    ///
    /// # Example
    ///
    /// ```
    /// use oefpil::Covariance;
    ///
    /// # fn main() -> Result<(), oefpil::OefpilError> {
    /// let covariance = Covariance::new_diagonal(3, 2)
    ///     .with_tile(0, 0, &[1.0, 2.0, 3.0])?
    ///     .with_tile(1, 1, &[4.0, 5.0, 6.0])?;
    /// let full = [
    ///     1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ///     0.0, 2.0, 0.0, 0.0, 0.0, 0.0,
    ///     0.0, 0.0, 3.0, 0.0, 0.0, 0.0,
    ///     0.0, 0.0, 0.0, 4.0, 0.0, 0.0,
    ///     0.0, 0.0, 0.0, 0.0, 5.0, 0.0,
    ///     0.0, 0.0, 0.0, 0.0, 0.0, 6.0,
    /// ];
    /// assert_eq!(full, covariance.to_full().as_slice());
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn new_diagonal(samples: usize, variables: usize) -> Self {
        assert!(samples > 0);
        assert!(variables > 0);
        let tiles = vec![Mode::None as c_int; variables];
        let fields = vec![0.0; samples * variables];
        Self {
            mode: Mode::Diagonal,
            tiles,
            fields,
            decomposition: Vec::new(),
            samples,
            variables,
            unknown_scale: false,
        }
    }
    /// Constructs a block-diagonal matrix (i.e, diagonal and block tiles on its diagonal).
    ///
    /// Allocates $`|V|`$ block tiles on its diagonal with $`|v|^2`$ fields per tile. Although the
    /// fields of block tiles are allocated, a tile can still be encoded as diagonal tile leveraging
    /// mathematical shortcuts.
    ///
    /// ```math
    /// \c{V}{V} = \begin{pmatrix}
    /// & \boxed\diagdown &        &                 & \\
    /// &                 & \ddots &                 & \\
    /// &                 &        & \boxed\diagdown & \end{pmatrix}
    /// ```
    ///
    /// # Panics
    ///
    ///   * Asserts `samples` and `variables` to be non-zero.
    ///
    /// # Example
    ///
    /// ```
    /// use oefpil::Covariance;
    ///
    /// # fn main() -> Result<(), oefpil::OefpilError> {
    /// let covariance = Covariance::new_block_diagonal(3, 3)
    ///     .with_tile(0, 0, &[1.0, 2.0, 3.0])?
    ///     .with_tile(1, 1, &[4.0,
    ///                        7.0, 5.0,
    ///                        8.0, 9.0, 6.0])?
    ///     .with_tile(2, 2, &[4.0, 7.0, 8.0,
    ///                        7.0, 5.0, 9.0,
    ///                        8.0, 9.0, 6.0])?;
    /// let full = [
    ///     1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ///     0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ///     0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ///     0.0, 0.0, 0.0, 4.0, 7.0, 8.0, 0.0, 0.0, 0.0,
    ///     0.0, 0.0, 0.0, 7.0, 5.0, 9.0, 0.0, 0.0, 0.0,
    ///     0.0, 0.0, 0.0, 8.0, 9.0, 6.0, 0.0, 0.0, 0.0,
    ///     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 7.0, 8.0,
    ///     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 5.0, 9.0,
    ///     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 9.0, 6.0,
    /// ];
    /// assert_eq!(full, covariance.to_full().as_slice());
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn new_block_diagonal(samples: usize, variables: usize) -> Self {
        assert!(samples > 0);
        assert!(variables > 0);
        let tiles = vec![Mode::None as c_int; variables];
        let fields = vec![0.0; samples.pow(2) * variables];
        Self {
            mode: Mode::BlockDiagonal,
            tiles,
            fields,
            decomposition: Vec::new(),
            samples,
            variables,
            unknown_scale: false,
        }
    }
    /// Constructs a matrix of solely diagonal tiles (i.e, diagonal tiles allover).
    ///
    /// Allocates $`|V|^2`$ diagonal tiles with $`|v|`$ fields per tile.
    ///
    /// ```math
    /// \c{V}{V} = \begin{pmatrix}
    /// & \diagdown & \dots  & \diagdown & \\
    /// & \vdots    & \ddots & \vdots    & \\
    /// & \diagdown & \dots  & \diagdown & \end{pmatrix}
    /// ```
    ///
    /// # Panics
    ///
    ///   * Asserts `samples` and `variables` to be non-zero.
    ///   * Currently, the implementation is limited to `variables == 2` and will panic otherwise.
    ///
    /// # Example
    ///
    /// ```
    /// use oefpil::Covariance;
    ///
    /// # fn main() -> Result<(), oefpil::OefpilError> {
    /// let covariance = Covariance::new_diagonals(3, 2)
    ///     .with_tile(0, 0, &[1.0, 2.0, 3.0])?
    ///     .with_tile(0, 1, &[7.0, 8.0, 9.0])?
    ///     .with_tile(1, 1, &[4.0, 5.0, 6.0])?;
    /// let full = [
    ///     1.0, 0.0, 0.0, 7.0, 0.0, 0.0,
    ///     0.0, 2.0, 0.0, 0.0, 8.0, 0.0,
    ///     0.0, 0.0, 3.0, 0.0, 0.0, 9.0,
    ///     7.0, 0.0, 0.0, 4.0, 0.0, 0.0,
    ///     0.0, 8.0, 0.0, 0.0, 5.0, 0.0,
    ///     0.0, 0.0, 9.0, 0.0, 0.0, 6.0,
    /// ];
    /// assert_eq!(full, covariance.to_full().as_slice());
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn new_diagonals(samples: usize, variables: usize) -> Self {
        assert!(samples > 0);
        assert!(variables > 0);
        assert_eq!(
            variables, 2,
            "currently only implemented for `variables == 2`"
        );
        let tiles = vec![Mode::None as c_int; variables.pow(2)];
        let fields = vec![0.0; samples * variables.pow(2)];
        Self {
            mode: Mode::Diagonals,
            tiles,
            fields,
            decomposition: Vec::new(),
            samples,
            variables,
            unknown_scale: false,
        }
    }
    /// Constructs a full matrix (i.e, diagonal and block tiles allover).
    ///
    /// Allocates $`|V|^2`$ block tiles with $`|v|^2`$ fields per tile. Although the fields of block
    /// tiles are allocated, a tile can still be encoded as diagonal tile leveraging mathematical
    /// shortcuts.
    ///
    /// ```math
    /// \c{V}{V} = \begin{pmatrix}
    /// & \boxed\diagdown & \dots  & \boxed\diagdown & \\
    /// & \vdots          & \ddots & \vdots          & \\
    /// & \boxed\diagdown & \dots  & \boxed\diagdown & \end{pmatrix}
    /// ```
    ///
    /// # Panics
    ///
    ///   * Asserts `samples` and `variables` to be non-zero.
    ///   * Currently, the implementation is limited to `variables == 2` and will panic otherwise.
    ///
    /// # Example
    ///
    /// ```
    /// use oefpil::Covariance;
    ///
    /// # fn main() -> Result<(), oefpil::OefpilError> {
    /// let covariance_1 = Covariance::new_full(3, 2)
    ///     .with_tile(0, 0, &[1.0, 2.0, 3.0])?
    ///     .with_tile(0, 1, &[4.0,
    ///                        7.0, 5.0,
    ///                        8.0, 9.0, 6.0])?
    ///     .with_tile(1, 1, &[4.0, 5.0, 6.0])?;
    /// let covariance_2 = Covariance::new_full(3, 2)
    ///     .with_tile(0, 0, &[1.0, 2.0, 3.0])?
    ///     .with_tile(0, 1, &[4.0, 7.0, 8.0,
    ///                        7.0, 5.0, 9.0,
    ///                        8.0, 9.0, 6.0])?
    ///     .with_tile(1, 1, &[4.0, 5.0, 6.0])?;
    /// let full = [
    ///     1.0, 0.0, 0.0, 4.0, 7.0, 8.0,
    ///     0.0, 2.0, 0.0, 7.0, 5.0, 9.0,
    ///     0.0, 0.0, 3.0, 8.0, 9.0, 6.0,
    ///     4.0, 7.0, 8.0, 4.0, 0.0, 0.0,
    ///     7.0, 5.0, 9.0, 0.0, 5.0, 0.0,
    ///     8.0, 9.0, 6.0, 0.0, 0.0, 6.0,
    /// ];
    /// assert_eq!(full, covariance_1.to_full().as_slice());
    /// assert_eq!(full, covariance_2.to_full().as_slice());
    /// # Ok(())
    /// # }
    /// ```
    #[must_use]
    pub fn new_full(samples: usize, variables: usize) -> Self {
        assert!(samples > 0);
        assert!(variables > 0);
        assert_eq!(
            variables, 2,
            "currently only implemented for `variables == 2`"
        );
        let tiles = vec![Mode::None as c_int; variables.pow(2)];
        let fields = vec![0.0; (samples * variables).pow(2)];
        Self {
            mode: Mode::Full,
            tiles,
            fields,
            decomposition: Vec::new(),
            samples,
            variables,
            unknown_scale: false,
        }
    }
    /// Sets `fields` of diagonal or block tile at `row` and `column`.
    ///
    /// See [`Self::with_tile()`] for its method chaining variant.
    ///
    /// Setting a tile besides the matrix diagonal will automatically set the transposed tile such
    /// that the matrix is always symmetric. Whether the tile is encoded as diagonal or block tile
    /// is deduced by the length of `fields`. A diagonal tile has `samples` fields whereas a block
    /// tile has `samples.pow(2)` row-major fields. A block tile on the matrix diagonal is validated
    /// to be symmetric unless only the lower triangle inclusive the diagonal is provided as
    /// row-major [triangular slice] of length `samples * (samples + 1) / 2` where the upper
    /// triangle is automatically constructed.
    ///
    /// Clears decomposition computed via [`Self::with_decomposition()`].
    ///
    /// [triangular slice]: https://en.wikipedia.org/wiki/Triangular_array
    ///
    /// # Errors
    ///
    /// Returns [`OefpilError`] on failure.
    pub fn set_tile(
        &mut self,
        row: usize,
        column: usize,
        fields: &[f64],
    ) -> Result<(), OefpilError> {
        self.decomposition.clear();
        unsafe {
            if !(row < self.variables && column < self.variables) {
                return Err(CovarianceError::OutOfBoundTileIndex.into());
            }
            let samples = i32::try_from(self.samples)?;
            let variables = i32::try_from(self.variables)?;
            let row = i32::try_from(row)?;
            let column = i32::try_from(column)?;
            match self.mode {
                Mode::None => unreachable!("invalid mode"),
                Mode::Diagonal => {
                    if row != column {
                        return Err(CovarianceError::NonDiagonalTileIndex.into());
                    }
                    if fields.len() != self.samples {
                        return Err(CovarianceError::MismatchingTile.into());
                    }
                    oefpil_tcm_diag_set_tile_diag(
                        samples,
                        variables,
                        self.fields.as_mut_ptr(),
                        self.tiles.as_mut_ptr(),
                        row,
                        fields.as_ptr(),
                    );
                }
                Mode::BlockDiagonal => {
                    if row != column {
                        return Err(CovarianceError::NonDiagonalTileIndex.into());
                    }
                    let oefpil_tcm_blockdiag_set_tile = if fields.len() == self.samples {
                        oefpil_tcm_blockdiag_set_tile_diag
                    } else if fields.len() == self.samples * (self.samples + 1) / 2 {
                        oefpil_tcm_blockdiag_set_tile_half
                    } else if fields.len() == self.samples.pow(2) {
                        if !self.is_symmetric(fields) {
                            return Err(CovarianceError::AsymmetricTile.into());
                        }
                        oefpil_tcm_blockdiag_set_tile_full
                    } else {
                        return Err(CovarianceError::MismatchingTile.into());
                    };
                    oefpil_tcm_blockdiag_set_tile(
                        samples,
                        variables,
                        self.fields.as_mut_ptr(),
                        self.tiles.as_mut_ptr(),
                        row,
                        fields.as_ptr(),
                    );
                }
                Mode::Diagonals => {
                    if fields.len() != self.samples {
                        return Err(CovarianceError::MismatchingTile.into());
                    }
                    oefpil_tcm_diags_set_tile_diag(
                        samples,
                        variables,
                        self.fields.as_mut_ptr(),
                        self.tiles.as_mut_ptr(),
                        row,
                        column,
                        fields.as_ptr(),
                    );
                }
                Mode::Full => {
                    let oefpil_tcm_full_set_tile = if fields.len() == self.samples {
                        oefpil_tcm_full_set_tile_diag
                    } else if fields.len() == self.samples * (self.samples + 1) / 2 {
                        oefpil_tcm_full_set_tile_half
                    } else if fields.len() == self.samples.pow(2) {
                        if row == column && !self.is_symmetric(fields) {
                            return Err(CovarianceError::AsymmetricTile.into());
                        }
                        oefpil_tcm_full_set_tile_full
                    } else {
                        return Err(CovarianceError::MismatchingTile.into());
                    };
                    oefpil_tcm_full_set_tile(
                        samples,
                        variables,
                        self.fields.as_mut_ptr(),
                        self.tiles.as_mut_ptr(),
                        row,
                        column,
                        fields.as_ptr(),
                    );
                }
            }
        }
        Ok(())
    }
    /// Method chaining variant of [`Self::set_tile()`].
    ///
    /// # Errors
    ///
    /// Returns [`OefpilError`] on failure.
    pub fn with_tile(
        mut self,
        row: usize,
        column: usize,
        fields: &[f64],
    ) -> Result<Self, OefpilError> {
        self.set_tile(row, column, fields).map(|()| self)
    }
    /// Assumes this sample covariance matrix is (already) multiplied by an unknown scalar.
    ///
    /// Causes [`Algorithm::fit()`] to multiply the [`Parameter::covariance`] matrix by
    /// [`Report::chi_squared_reduced`] which is equivalent to rescaling this sample covariance
    /// matrix such that $`\chi_{\nu}^2\approx 1`$. This effectively pretends [Pearson's $`\chi^2`$
    /// test] of [goodness of fit] to be successful. This is common practice to correct for over- or
    /// under-dispersion due to partially or fully unknown (co)variances. Note however that in this
    /// way, the parameter uncertainties do no longer scale with the sample uncertainties and hence
    /// this correction is in conflict with uncertainty propagation.
    ///
    /// [Pearson's $`\chi^2`$ test]: https://en.wikipedia.org/wiki/Pearson%27s_chi-squared_test
    /// [goodness of fit]: https://en.wikipedia.org/wiki/Goodness_of_fit
    ///
    /// # Validity
    ///
    /// This correction is as valid as the definition of [`Report::degrees_of_freedom`].
    #[must_use]
    #[inline]
    pub const fn with_unknown_scale(mut self) -> Self {
        self.unknown_scale = true;
        self
    }
    /// Computes and stores [Cholesky decomposition] $`L`$ where $`\c{V}{V}=LL^T`$.
    ///
    /// Required by [`Distribution::sample()`].
    #[cfg_attr(
        not(feature = "random"),
        doc = r"

[`Distribution::sample()`]:
https://docs.rs/oefpil/latest/oefpil/struct.Distribution.html#method.sample
"
    )]
    ///
    /// [Cholesky decomposition]: https://en.wikipedia.org/wiki/Cholesky_decomposition
    ///
    /// # Errors
    ///
    /// Returns [`OefpilError`] on failure.
    pub fn with_decomposition(mut self) -> Result<Self, OefpilError> {
        self.decomposition = self.to_full();
        let order = self.order().try_into()?;
        let mut info = 0;
        unsafe {
            dpotrf(
                c"L".as_ptr(),
                order,
                self.decomposition.as_mut_ptr(),
                order,
                &raw mut info,
            );
        };
        match info {
            0 => Ok(self),
            _ => Err(AlgorithmError::NumericalError.into()),
        }
    }
    /// Creates full matrix from tiled covariance matrix.
    ///
    /// As the covariance matrix is symmetric, row-major and column-major order coincide.
    #[must_use]
    pub fn to_full(&self) -> Vec<f64> {
        let order = self.order();
        let mut fields = vec![0f64; order * order];
        match self.mode {
            Mode::None => unreachable!("invalid mode"),
            Mode::Diagonal => {
                for row_column in 0..self.fields.len() {
                    fields[row_column * order + row_column] = self.fields[row_column];
                }
            }
            Mode::BlockDiagonal => {
                for variable in 0..self.variables {
                    let target = variable * (order + 1) * self.samples;
                    let source = variable * self.samples.pow(2);
                    for row in 0..self.samples {
                        for column in 0..self.samples {
                            let target = target + order * row + column;
                            let source = source + row * self.samples + column;
                            fields[target] = self.fields[source];
                        }
                    }
                }
            }
            Mode::Diagonals => {
                for row in 0..self.variables {
                    let target = row * order * self.samples;
                    for column in 0..self.variables {
                        let target = target + column * self.samples;
                        let source = (row * self.variables + column) * self.samples;
                        for sample in 0..self.samples {
                            let target = target + sample * order + sample;
                            let source = source + sample;
                            fields[target] = self.fields[source];
                        }
                    }
                }
            }
            Mode::Full => fields.clone_from(&self.fields),
        }
        fields
    }
    /// Number of observations $`|v|`$ (i.e., the number of samples per variable $`V^I`$).
    #[must_use]
    #[inline]
    pub const fn samples(&self) -> usize {
        self.samples
    }
    /// Number of independent as well as dependent variables $`|V|\equiv|X|+|Y|`$.
    #[must_use]
    #[inline]
    pub const fn variables(&self) -> usize {
        self.variables
    }
    /// Order $`|V||v|`$.
    #[must_use]
    #[inline]
    pub const fn order(&self) -> usize {
        self.samples * self.variables
    }
    #[must_use]
    fn is_symmetric(&self, fields: &[f64]) -> bool {
        for row in 0..self.samples {
            for column in 0..self.samples {
                #[allow(clippy::float_cmp)]
                if fields[row * self.samples + column] != fields[column * self.samples + row] {
                    return false;
                }
            }
        }
        true
    }
}

/// Parameter comprising mean vector $`\lambda`$, deviation vector $`\varsigma`$, covariance matrix
/// $`\gdef\c#1#2{\textup{Cov}(#1,#2)}\c{P}{P}`$, and correlation matrix
/// $`\gdef\C#1#2{\textup{Cor}(#1,#2)}\C{P}{P}`$.
#[allow(clippy::too_long_first_doc_paragraph)]
#[derive(Debug, Default)]
pub struct Parameter<'a> {
    /// Mean vector $`\lambda`$ as input and $`\hat\lambda`$ as output.
    ///
    /// Serves as input for the initial estimates and as output for the  optimum estimates.
    pub mean: &'a mut [f64],
    /// Standard deviation vector $`\varsigma`$ as output.
    ///
    /// ```math
    /// \varsigma_i\equiv\sqrt{\c{P_i}{P_i}}\quad \forall\ i \in [1,|P|]
    /// ````
    ///
    /// This vector is optional, indicated by an empty slice.
    pub deviation: &'a mut [f64],
    /// Covariance matrix $`\c{P}{P}`$ as output.
    ///
    /// ```math
    /// \c{P}{P} \equiv \begin{pmatrix}
    /// & \c{P_1}{P_1}     & \dots  & \c{P_1}{P_{|P|}}     & \\[1em]
    /// & \vdots           & \ddots & \vdots               & \\[1em]
    /// & \c{P_{|P|}}{P_1} & \dots  & \c{P_{|P|}}{P_{|P|}} &
    /// \end{pmatrix}
    /// ```
    ///
    /// As a covariance matrix is symmetric, row-major and column-major order coincide.
    pub covariance: &'a mut [f64],
    /// Correlation matrix $`\C{P}{P}`$ as output.
    ///
    /// ```math
    /// \C{P_r}{P_c}\equiv\frac{\c{P_r}{P_c}}{\varsigma_r\varsigma_c}\quad \forall\ r,c \in [1,|P|]
    /// ```
    ///
    /// As a correlation matrix is symmetric, row-major and column-major order coincide.
    ///
    /// This matrix is optional, indicated by an empty slice. If non-empty, [`Self::deviation`] must
    /// be non-empty as well.
    pub correlation: &'a mut [f64],
}

impl Parameter<'_> {
    /// Ensures slice lengths agree with `model`.
    ///
    /// Called in [`Algorithm::fit()`].
    ///
    /// # Errors
    ///
    /// Returns [`ParameterError`] on failure.
    pub fn validate(&self, model: Model) -> Result<(), ParameterError> {
        if self.mean.len() != model.dfdp.len() {
            return Err(ParameterError::MismatchingMean);
        }
        if self.covariance.len() != self.mean.len().pow(2) {
            return Err(ParameterError::MismatchingCovariance);
        }
        if !self.deviation.is_empty() && self.deviation.len() != self.mean.len() {
            return Err(ParameterError::MismatchingDeviation);
        }
        if !self.correlation.is_empty() && self.correlation.len() != self.covariance.len() {
            return Err(ParameterError::MismatchingCorrelation);
        }
        if !self.correlation.is_empty() && self.deviation.is_empty() {
            return Err(ParameterError::CorrelationRequiresDeviation);
        }
        Ok(())
    }
}

/// Fitting report on success.
#[derive(Debug, Clone, Copy)]
pub struct Report {
    /// $`\chi^2\equiv\sum_{I,i} (v^I_i-\hat\mu^I_i)^2`$
    pub chi_squared: f64,
    /// $`\chi^2_\nu\equiv\chi^2/\nu`$
    ///
    /// # Validity
    ///
    /// See [`Report::degrees_of_freedom`] for its validity.
    pub chi_squared_reduced: f64,
    /// $`\nu\equiv|v|-|P|`$
    ///
    /// # Validity
    ///
    /// Valid whenever [`Covariance::samples()`] is sufficiently large and [`Model`] is sufficiently
    /// linear.[^1]
    ///
    /// [^1]: R. Andrae, T. Schulze-Hartung, and P. Melchior, “Dos and don'ts of reduced
    /// chi-squared”, [Instrumentation and Methods for Astrophysics
    /// (2010)](https://doi.org/10.48550/arXiv.1012.3754).
    pub degrees_of_freedom: usize,
    /// Iterations until convergence.
    pub iterations: usize,
}

impl Report {
    /// Computes [p-value] from [incomplete gamma function].
    ///
    /// ```math
    /// p = Q(s=\frac{1}{2}\nu,x=\frac{1}{2}\chi^2) = \frac{\Gamma(s,x)}{\Gamma(s)}
    ///                                             = 1-P(s,x)
    ///                                             = 1-\frac{\gamma(s,x)}{\Gamma(s)}
    /// ```
    ///
    /// [p-value]: https://en.wikipedia.org/wiki/P-value
    /// [incomplete gamma function]: https://en.wikipedia.org/wiki/Incomplete_gamma_function
    ///
    /// # Validity
    ///
    /// See [`Report::degrees_of_freedom`] for its validity.
    ///
    /// # Errors
    ///
    /// Returns [`OefpilError`] on failure.
    pub fn chi_squared_p_value(&self) -> Result<f64, OefpilError> {
        let mut p = f64::NAN;
        let mut q = f64::NAN;
        let mut error = 0;
        unsafe {
            oefpil_sys::dcdchi(
                self.chi_squared,
                f64::from(i32::try_from(self.degrees_of_freedom)?),
                &raw mut p,
                &raw mut q,
                &raw mut error,
            );
        }
        match error {
            0 => Ok(q),
            _ => Err(AlgorithmError::NumericalError.into()),
        }
    }
}

impl Default for Report {
    fn default() -> Self {
        Self {
            chi_squared: f64::NAN,
            chi_squared_reduced: f64::NAN,
            degrees_of_freedom: 0,
            iterations: 0,
        }
    }
}

/// OEFPIL errors.
#[derive(Debug, Display, Error, From)]
#[non_exhaustive]
pub enum OefpilError {
    /// I/O error.
    #[display("I/O error")]
    Io(#[error(source)] io::Error),
    /// Conversation error between [`usize`] and [`i32`].
    #[display("conversation error")]
    TryFromInt(#[error(source)] TryFromIntError),
    /// Invalid model.
    #[display("invalid model")]
    Model(#[error(source)] ModelError),
    /// Invalid variable.
    #[display("invalid variable")]
    Variable(#[error(source)] VariableError),
    /// Invalid distribution.
    #[display("invalid distribution")]
    #[cfg(feature = "random")]
    Distribution(#[error(source)] DistributionError),
    /// Invalid covariance.
    #[display("invalid covariance")]
    Covariance(#[error(source)] CovarianceError),
    /// Invalid parameter.
    #[display("invalid parameter")]
    Parameter(#[error(source)] ParameterError),
    /// Algorithm error.
    #[display("algorithm error")]
    Algorithm(#[error(source)] AlgorithmError),
}

/// Model errors.
#[derive(Debug, Display, Error, PartialEq, Eq, Clone, Copy)]
#[non_exhaustive]
pub enum ModelError {
    /// Empty variable derivative vector [`Model::dfdx`].
    #[display("empty variable derivative vector")]
    EmptyVariableDerivative,
    /// Empty parameter derivative vector [`Model::dfdp`].
    #[display("empty parameter derivative vector")]
    EmptyParameterDerivative,
}

/// Variable errors.
#[derive(Debug, Display, Error, PartialEq, Eq, Clone, Copy)]
#[non_exhaustive]
pub enum VariableError {
    /// Mismatching length of multivariate [`Variable::sample`] matrix.
    #[display("mismatching length of multivariate sample matrix")]
    MismatchingSample,
    /// Mismatching length of multivariate [`Variable::mean`] matrix.
    #[display("mismatching length of multivariate mean matrix")]
    MismatchingMean,
    /// Mismatching [`Covariance::order()`].
    #[display("mismatching covariance order")]
    MismatchingCovariance,
    /// Insufficient degrees of freedom $`\nu\equiv|v|-|P|`$.
    ///
    /// Requirements are:
    ///
    ///   * $`\nu\ge 0`$ for a fully known [`Covariance`],
    ///   * $`\nu\ge 1`$ for a [`Covariance::with_unknown_scale()`].
    ///
    /// There must be at least one parameter to fit against a sample of at least one observation,
    /// leading to the requirement $`\nu\ge 0`$. For a [`Covariance::with_unknown_scale()`], the
    /// [`Parameter::covariance`] matrix is being multiplied by $`\chi_\nu^2\equiv\chi^2/\nu`$ which
    /// would be infinite without another observation due to division by zero, leading to the
    /// requirement $`\nu\ge 1`$.
    #[display("insufficient degrees of freedom")]
    InsufficientDegreesOfFreedom,
}

/// Distribution errors.
#[derive(Debug, Display, Error, PartialEq, Eq, Clone, Copy)]
#[non_exhaustive]
#[cfg(feature = "random")]
pub enum DistributionError {
    /// Mismatching length of multivariate [`Distribution::sample`] matrix.
    #[display("mismatching length of multivariate sample matrix")]
    MismatchingSample,
    /// Mismatching length of multivariate [`Distribution::mean`] matrix.
    #[display("mismatching length of multivariate mean matrix")]
    MismatchingMean,
    /// [`Distribution::sample()`] requires [`Covariance::with_decomposition()`].
    #[display("missing decomposition")]
    MissingDecomposition,
}

/// Model errors.
#[derive(Debug, Display, Error, PartialEq, Eq, Clone, Copy)]
#[non_exhaustive]
pub enum CovarianceError {
    /// Out of bound tile index.
    #[display("out of bound tile index")]
    OutOfBoundTileIndex,
    /// Non-diagonal tile index.
    #[display("non-diagonal tile index")]
    NonDiagonalTileIndex,
    /// Mismatching tile length.
    #[display("mismatching tile length")]
    MismatchingTile,
    /// Asymmetric tile on diagonal of block-diagonal or full covariance matrix.
    #[display("asymmetric tile on diagonal of block-diagonal or full covariance matrix")]
    AsymmetricTile,
}

/// Model errors.
#[derive(Debug, Display, Error, PartialEq, Eq, Clone, Copy)]
#[non_exhaustive]
pub enum ParameterError {
    /// Mismatching length of [`Parameter::mean`] vector.
    #[display("mismatching length of mean vector")]
    MismatchingMean,
    /// Mismatching length [`Parameter::deviation`] vector.
    #[display("mismatching length of deviation vector")]
    MismatchingDeviation,
    /// Mismatching length [`Parameter::covariance`] matrix.
    #[display("mismatching length of covariance matrix")]
    MismatchingCovariance,
    /// Mismatching length [`Parameter::correlation`] matrix.
    #[display("mismatching length of correlation matrix")]
    MismatchingCorrelation,
    /// Non-empty [`Parameter::correlation`] requires non-empty  [`Parameter::deviation`].
    #[display("non-empty correlation matrix slice requires non-empty deviation vector slice")]
    CorrelationRequiresDeviation,
}

/// Algorithm errors.
#[derive(Debug, Display, Error, PartialEq, Eq, Clone, Copy)]
#[non_exhaustive]
pub enum AlgorithmError {
    /// Iteration limit reached before satisfying convergence [`Criterion`].
    #[display("iteration limit reached before satisfying convergence criterion")]
    IterationLimit,
    /// Numerical error probably due to infinite values.
    #[display("numerical error probably due to infinite values")]
    NumericalError,
}

/// Reads matrix from text file at `path` into row-major [`Vec`].
///
/// The `rows` are [`newline-separated`] and the `columns` are [`whitespace-separated`].
///
/// If `transpose`, reads the text file of `rows` lines and `columns` words into column-major
/// [`Vec`] effectively transposing the matrix as [`Vec`] is supposed to be in row-major order.
///
/// [`newline-separated`]: `str::lines()`
/// [`whitespace-separated`]: `str::split_whitespace()`
///
/// # Errors
///
/// In addition to standard I/O errors, returns errors of [`ErrorKind::InvalidData`] if the
/// requested number of `rows` or number of `columns` differ from the number of lines or the number
/// of words, respectively or a word cannot be parsed as [`f64`].
pub fn read_matrix(
    path: &Path,
    rows: usize,
    columns: usize,
    transpose: bool,
) -> io::Result<Vec<f64>> {
    let file = BufReader::new(OpenOptions::new().read(true).open(path)?);
    let mut matrix = vec![0f64; rows * columns];
    let mut lines = file.lines();
    for row in 0..rows {
        let line = lines.next().ok_or_else(|| {
            io::Error::new(
                ErrorKind::InvalidData,
                format!("expected {rows} row(s), found {} row(s)", row + 1),
            )
        })??;
        let mut words = line.split_whitespace();
        for column in 0..columns {
            let word = words.next().ok_or_else(|| {
                io::Error::new(
                    ErrorKind::InvalidData,
                    format!(
                        "expected {columns} column(s), found {} column(s)",
                        column + 1
                    ),
                )
            })?;
            let field = if transpose {
                column * rows + row
            } else {
                row * rows + column
            };
            matrix[field] = word
                .parse()
                .map_err(|err| io::Error::new(ErrorKind::InvalidData, err))?;
        }
        if words.next().is_some() {
            Err(io::Error::new(
                ErrorKind::InvalidData,
                format!("expected {columns} column(s), found more"),
            ))?;
        }
    }
    if lines.next().is_some() {
        Err(io::Error::new(
            ErrorKind::InvalidData,
            format!("expected {rows} rows(s), found more"),
        ))?;
    }
    Ok(matrix)
}

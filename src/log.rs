//! Logging into a [`core::fmt::Write`] sink.

use crate::Verbosity;
use alloc::string::String;
use core::fmt::{self, Write as _};

/// Notation of floating-point numbers in the log.
///
/// The variants print like the respective conversion of C's `printf()`.
///
/// Default is [`Self::Exponential`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Notation {
    /// Exponential notation like `1.50e+00`, conversion `e`.
    #[default]
    Exponential,
    /// Exponential notation like `1.50E+00`, conversion `E`.
    UpperExponential,
    /// Fixed-point notation like `1.50`, conversion `f`.
    Fixed,
    /// Fixed-point or exponential notation by magnitude without trailing zeros, conversion `g`.
    General,
    /// Fixed-point or exponential notation by magnitude without trailing zeros, conversion `G`.
    UpperGeneral,
}

/// Format of floating-point numbers in the log.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Format {
    /// Notation.
    ///
    /// Default is [`Notation::Exponential`].
    pub notation: Notation,
    /// Number of digits after the decimal point or number of significant digits if
    /// [`Notation::General`] or [`Notation::UpperGeneral`].
    ///
    /// Default is `8`.
    pub precision: usize,
}

impl Default for Format {
    fn default() -> Self {
        Self {
            notation: Notation::default(),
            precision: 8,
        }
    }
}

/// Log sink discarding everything for use with [`Verbosity::Silent`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Discard;

impl fmt::Write for Discard {
    #[inline]
    fn write_str(&mut self, _text: &str) -> fmt::Result {
        Ok(())
    }
}

/// Appends the magnitude of `x` in exponential notation with an exponent of at least two digits.
fn push_exponential(out: &mut String, x: f64, precision: usize) {
    let text = alloc::format!("{x:.precision$e}");
    let (mantissa, exponent) = text.split_once('e').unwrap();
    let (sign, digits) = exponent
        .strip_prefix('-')
        .map_or(('+', exponent), |digits| ('-', digits));
    out.push_str(mantissa);
    out.push('e');
    out.push(sign);
    if digits.len() < 2 {
        out.push('0');
    }
    out.push_str(digits);
}

/// Appends the magnitude of `x` in general notation.
fn push_general(out: &mut String, x: f64, precision: usize) {
    let precision = precision.max(1);
    let text = alloc::format!("{:.*e}", precision - 1, x);
    let exponent = text
        .split_once('e')
        .and_then(|(_mantissa, exponent)| exponent.parse::<isize>().ok())
        .unwrap();
    let mut body = String::new();
    #[allow(clippy::cast_possible_wrap, clippy::cast_sign_loss)]
    let (number, suffix) = if exponent >= -4 && exponent < precision as isize {
        let decimals = (precision as isize - 1 - exponent) as usize;
        let _ = write!(body, "{x:.decimals$}");
        (body.as_str(), "")
    } else {
        push_exponential(&mut body, x, precision - 1);
        body.split_at(body.find('e').unwrap())
    };
    let number = if number.contains('.') {
        number.trim_end_matches('0').trim_end_matches('.')
    } else {
        number
    };
    out.push_str(number);
    out.push_str(suffix);
}

/// Formats `x` as C's `printf()` with `flags` of `space_for_sign` and `width`.
pub fn c_float(
    x: f64,
    notation: Notation,
    precision: usize,
    space_for_sign: bool,
    width: usize,
) -> String {
    let mut out = String::new();
    if x.is_sign_negative() {
        out.push('-');
    } else if space_for_sign {
        out.push(' ');
    }
    let magnitude = x.abs();
    if x.is_nan() {
        out.push_str("nan");
    } else if x.is_infinite() {
        out.push_str("inf");
    } else {
        match notation {
            Notation::Exponential | Notation::UpperExponential => {
                push_exponential(&mut out, magnitude, precision);
            }
            Notation::Fixed => {
                let _ = write!(out, "{magnitude:.precision$}");
            }
            Notation::General | Notation::UpperGeneral => {
                push_general(&mut out, magnitude, precision);
            }
        }
    }
    if matches!(
        notation,
        Notation::UpperExponential | Notation::UpperGeneral
    ) {
        out.make_ascii_uppercase();
    }
    let padding = width.saturating_sub(out.len());
    if padding > 0 {
        out.insert_str(0, &" ".repeat(padding));
    }
    out
}

/// Log of a fit, which continues if writing fails but remembers the failure.
pub struct Log<'a> {
    sink: &'a mut dyn fmt::Write,
    verbosity: Verbosity,
    format: Format,
    failed: bool,
}

impl<'a> Log<'a> {
    pub fn new(sink: &'a mut dyn fmt::Write, verbosity: Verbosity, format: Format) -> Self {
        Self {
            sink,
            verbosity,
            format,
            failed: false,
        }
    }
    /// Whether writing failed.
    pub const fn failed(&self) -> bool {
        self.failed
    }
    /// Whether the verbosity is at least `level`.
    pub fn at(&self, level: Verbosity) -> bool {
        self.verbosity >= level
    }
    pub const fn format(&self) -> Format {
        self.format
    }
    pub fn print(&mut self, text: fmt::Arguments<'_>) {
        if !self.failed && self.sink.write_fmt(text).is_err() {
            self.failed = true;
        }
    }
    /// Prints column-major `m` by `n` matrix `a` below optional `label`.
    pub fn matrix(&mut self, m: usize, n: usize, a: &[f64], label: Option<&str>) {
        self.matrix_with_format(m, n, a, label, self.format);
    }
    pub fn matrix_with_format(
        &mut self,
        m: usize,
        n: usize,
        a: &[f64],
        label: Option<&str>,
        format: Format,
    ) {
        if let Some(label) = label {
            self.print(format_args!("{label}\n"));
        }
        for i in 0..m {
            for j in 0..n {
                let field = c_float(a[j * m + i], format.notation, format.precision, true, 10);
                self.print(format_args!(" {field}"));
            }
            self.print(format_args!("\n"));
        }
    }
}

#[cfg(test)]
#[allow(clippy::cast_sign_loss, clippy::unnecessary_cast)]
mod tests {
    use super::{Notation, c_float};
    use alloc::{ffi::CString, format, string::String};
    use core::ffi::c_char;

    unsafe extern "C" {
        fn snprintf(buffer: *mut c_char, size: usize, format: *const c_char, ...) -> i32;
    }

    fn c_printf(format: &str, x: f64) -> String {
        let format = CString::new(format).unwrap();
        let mut buffer = [0 as c_char; 512];
        let length = unsafe { snprintf(buffer.as_mut_ptr(), buffer.len(), format.as_ptr(), x) };
        let length = usize::try_from(length).unwrap();
        String::from_utf8(buffer[..length].iter().map(|&c| c as u8).collect()).unwrap()
    }

    #[test]
    fn agrees_with_printf() {
        let values = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            1.5,
            0.1,
            0.25,
            0.125,
            2.5,
            3.5,
            1e-4,
            9.9999e-5,
            1e-5,
            999_999.5,
            1e5,
            1e6,
            1e22,
            1e100,
            1e-100,
            1.25e-300,
            123.456,
            -123.456,
            0.000_123_456,
            99.999_999_999,
            1.797_693_134_862_315_7e308,
            5e-324,
            2.225_073_858_507_201_4e-308,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
            -f64::NAN,
        ];
        for (conversion, notation) in [
            ('e', Notation::Exponential),
            ('E', Notation::UpperExponential),
            ('f', Notation::Fixed),
            ('g', Notation::General),
            ('G', Notation::UpperGeneral),
        ] {
            for precision in [0, 1, 3, 8, 15] {
                for x in values {
                    assert_eq!(
                        c_printf(&format!("% 10.{precision}{conversion}"), x),
                        c_float(x, notation, precision, true, 10),
                        "% 10.{precision}{conversion} of {x}",
                    );
                    assert_eq!(
                        c_printf(&format!("%.{precision}{conversion}"), x),
                        c_float(x, notation, precision, false, 0),
                        "%.{precision}{conversion} of {x}",
                    );
                }
            }
        }
    }
}

//! Ratios of the incomplete gamma function.

use special::Gamma;

/// Computes the p-value `q` of `chisq` and `nu` degrees of freedom along with `p = 1 - q`.
///
/// Both are NaN if `chisq` is negative or `nu` is not positive.
pub fn dcdchi(chisq: f64, nu: f64) -> (f64, f64) {
    dgami(0.5 * nu, 0.5 * chisq)
}

/// Computes the ratios `p` and `q` of the lower and upper incomplete gamma function at `a` and `x`.
fn dgami(a: f64, x: f64) -> (f64, f64) {
    if x >= 0.0 && a > 0.0 {
        if x < a + 1.0 {
            let p = x.inc_gamma(a);
            (p, 1.0 - p)
        } else {
            let q = upper_inc_gamma(a, x);
            (1.0 - q, q)
        }
    } else {
        (f64::NAN, f64::NAN)
    }
}

/// Computes the upper incomplete gamma function ratio at `a` and `x` where `x >= a + 1`.
///
/// Evaluates its continued fraction directly rather than the complement of the lower ratio as the
/// latter loses all accuracy once the upper ratio approaches the machine epsilon. Returns NaN if
/// the continued fraction does not converge.
fn upper_inc_gamma(a: f64, x: f64) -> f64 {
    const TINY: f64 = f64::MIN_POSITIVE / f64::EPSILON;
    if x.is_infinite() {
        return 0.0;
    }
    // Modified Lentz's method.
    let mut denominator = x + 1.0 - a;
    let mut upper = 1.0 / TINY;
    let mut lower = 1.0 / denominator;
    let mut fraction = lower;
    for term in 1..=1_000_000u32 {
        let term = f64::from(term);
        let numerator = term * (a - term);
        denominator += 2.0;
        lower = libm::fma(numerator, lower, denominator);
        lower = 1.0 / if lower.abs() < TINY { TINY } else { lower };
        upper = denominator + numerator / upper;
        upper = if upper.abs() < TINY { TINY } else { upper };
        let delta = upper * lower;
        fraction *= delta;
        if (delta - 1.0).abs() < f64::EPSILON {
            let (ln_gamma, _sign) = Gamma::ln_gamma(a);
            return libm::exp(libm::fma(a, libm::log(x), -x) - ln_gamma) * fraction;
        }
    }
    f64::NAN
}

#[cfg(test)]
mod tests {
    use super::dcdchi;

    #[test]
    fn agrees_with_reference() {
        for nu in [1.0, 2.0, 3.0, 7.0, 20.0, 150.0, 1000.0] {
            for chisq in [0.0, 1e-3, 0.5, 1.0, 3.0, 8.0, 40.0, 200.0, 1000.0, 5000.0] {
                let (mut p, mut q, mut error) = (f64::NAN, f64::NAN, 0);
                unsafe {
                    oefpil_sys::dcdchi(chisq, nu, &raw mut p, &raw mut q, &raw mut error);
                }
                assert_eq!(0, error);
                let (actual_p, actual_q) = dcdchi(chisq, nu);
                assert!((actual_p - p).abs() <= 1e-11, "{actual_p} != {p}");
                assert!((actual_q - q).abs() <= 1e-11 * q, "{actual_q} != {q}");
            }
        }
    }

    #[test]
    fn invalid_arguments_are_nan() {
        assert!(dcdchi(1.0, 0.0).1.is_nan());
        assert!(dcdchi(-1.0, 1.0).1.is_nan());
        assert!(dcdchi(f64::NAN, 1.0).1.is_nan());
    }
}

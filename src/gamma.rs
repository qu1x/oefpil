//! Ratios of the incomplete gamma function.

#![allow(clippy::suboptimal_flops)]

/// Half the logarithm of two times pi.
const LN_SQRT_2PI: f64 = 0.918_938_533_204_672_7;

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
            let p = lower_inc_gamma(a, x);
            (p, 1.0 - p)
        } else {
            let q = upper_inc_gamma(a, x);
            (1.0 - q, q)
        }
    } else {
        (f64::NAN, f64::NAN)
    }
}

/// Computes the difference between `ln(Gamma(a))` and its Stirling approximation.
fn stirling_error(a: f64) -> f64 {
    if a < 16.0 {
        libm::lgamma(a) - ((a - 0.5) * libm::log(a) - a + LN_SQRT_2PI)
    } else {
        let square = a * a;
        (1.0 / 12.0
            - (1.0 / 360.0
                - (1.0 / 1260.0 - (1.0 / 1680.0 - 1.0 / 1188.0 / square) / square) / square)
                / square)
            / a
    }
}

/// Computes the logarithm of `x^a * exp(-x) / Gamma(a)`.
///
/// Evaluates it in terms of the relative deviation of `x` from `a` rather than directly as the
/// latter subtracts terms growing with `a * ln(a)` and hence loses accuracy for large `a`.
fn ln_prefix(a: f64, x: f64) -> f64 {
    let deviation = (x - a) / a;
    // Logarithm of `1 + deviation` minus `deviation`.
    let difference = if deviation.abs() < 0.1 {
        let mut series = 0.05;
        for order in (2..20).rev() {
            series = 1.0 / f64::from(order) - deviation * series;
        }
        -deviation * deviation * series
    } else {
        libm::log(x / a) - deviation
    };
    a * difference + 0.5 * libm::log(a) - LN_SQRT_2PI - stirling_error(a)
}

/// Computes the lower incomplete gamma function ratio at `a` and `x` where `x < a + 1`.
///
/// Evaluates its series. Returns NaN if the series does not converge.
fn lower_inc_gamma(a: f64, x: f64) -> f64 {
    let (mut term, mut sum, mut denominator) = (1.0 / a, 1.0 / a, a);
    for _ in 0..1_000_000u32 {
        denominator += 1.0;
        term *= x / denominator;
        sum += term;
        if term <= sum * f64::EPSILON / 16.0 {
            return libm::exp(ln_prefix(a, x)) * sum;
        }
    }
    f64::NAN
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
            return libm::exp(ln_prefix(a, x)) * fraction;
        }
    }
    f64::NAN
}

#[cfg(test)]
mod tests {
    use super::dcdchi;

    /// Reference values are exact p-values of closed forms evaluated with 80 decimal digits, which
    /// exist as `nu / 2` is an integer or half-integer.
    #[test]
    #[allow(clippy::excessive_precision, clippy::unreadable_literal)]
    fn agrees_with_exact_values() {
        #[rustfmt::skip]
        let table = [
            (0.0, 1.0, 1.0),
            (1e-300, 1.0, 1.0),
            (1e-10, 1.0, 0.9999920211543921),
            (0.5, 1.0, 0.4795001221869535),
            (1.0, 1.0, 0.3173105078629141),
            (3.84, 1.0, 0.050043521248705106),
            (100.0, 1.0, 1.523970604832105e-23),
            (1.0, 2.0, 0.6065306597126334),
            (10.0, 2.0, 0.006737946999085467),
            (1400.0, 2.0, 9.85967654375977e-305),
            (2.5, 3.0, 0.4752910833430206),
            (0.1, 7.0, 0.9999976885812014),
            (7.0, 7.0, 0.42887985755305474),
            (30.0, 7.0, 9.495972508134183e-05),
            (1e-05, 10.0, 1.0),
            (10.0, 10.0, 0.4404932850652124),
            (1000.0, 10.0, 1.8702907209159498e-208),
            (20.0, 30.0, 0.9165415270653372),
            (30.0, 30.0, 0.46565370894400965),
            (45.0, 30.0, 0.038601758266317336),
            (95.0, 100.0, 0.6225792918781735),
            (100.0, 100.0, 0.48119168452795674),
            (150.0, 100.0, 0.0009039320423540091),
            (285.0, 300.0, 0.7242624085345621),
            (300.0, 300.0, 0.4891417702506403),
            (788.7298334620741, 750.0, 0.1585515148958366),
            (1000.0, 1000.0, 0.49405285382923964),
            (1050.0, 1000.0, 0.1324740568214771),
            (910.5572809000084, 1000.0, 0.9797163290773768),
            (1100.0, 1000.0, 0.014614408126295194),
            (8500.0, 5000.0, 1.4560849131328786e-186),
            (4900.0, 5000.0, 0.8413611103322287),
            (9985.857864376268, 10000.0, 0.5379746576481316),
            (10000.0, 10000.0, 0.49811936596618267),
            (10200.0, 10000.0, 0.07932888107761897),
            (100000.0, 100000.0, 0.4994052918952067),
            (101000.0, 100000.0, 0.012868840377233669),
            (1000000.0, 1000000.0, 0.4998119368033945),
            (1008485.2813742385, 1000000.0, 1.0914018031672503e-09),
            (990000.0, 1000000.0, 0.99999999999935),
            (1042426.4068711929, 1000000.0, 1.1176698881326725e-192),
        ];
        for (chisq, nu, q) in table {
            let (actual_p, actual_q) = dcdchi(chisq, nu);
            assert!(
                (actual_q - q).abs() <= 1e-12 * q,
                "chi2 {chisq} nu {nu}: {actual_q} != {q}"
            );
            assert!(
                (actual_p - (1.0 - q)).abs() <= 1e-12,
                "chi2 {chisq} nu {nu}: {actual_p} != {}",
                1.0 - q
            );
        }
    }

    /// Neither crate can depend on the other, so `oefpil-sys` has its own copy of these formulas.
    #[test]
    fn agrees_with_oefpil_sys() {
        for nu in [1.0, 2.0, 3.0, 7.0, 20.0, 150.0, 1000.0, 1e4, 1e6] {
            for factor in [
                0.0, 1e-3, 0.1, 0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 1.5, 2.0, 5.0,
            ] {
                let chisq = factor * nu;
                let (mut p, mut q, mut error) = (f64::NAN, f64::NAN, 0);
                unsafe {
                    oefpil_sys::dcdchi(chisq, nu, &raw mut p, &raw mut q, &raw mut error);
                }
                assert_eq!(0, error);
                let (actual_p, actual_q) = dcdchi(chisq, nu);
                assert!((actual_p - p).abs() <= 1e-12, "{actual_p} != {p}");
                assert!((actual_q - q).abs() <= 1e-12 * q, "{actual_q} != {q}");
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

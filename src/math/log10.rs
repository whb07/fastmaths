//! log10(x) implementation.
//!
//! Computes log10 via ln(x) scaled by log10(e) with high/low constants. This
//! mirrors glibc-style error-compensated scaling for ≤1 ULP accuracy.

use super::{fma_internal, ln, scalbn_internal};

const NEAR1_MIN: f64 = 0.75;
const NEAR1_MAX: f64 = 1.5;

const LOG10_E_HI: f64 = 0.434_294_462_203_979_5;
const LOG10_E_LO: f64 = 1.969_927_232_448_043_2e-08;
const LOG10_E_FULL: f64 = 0.434_294_481_903_251_8;
const LOG10_2_HI: f64 = 0.301_029_920_578_002_93;
const LOG10_2_LO: f64 = 7.508_597_826_832_997e-08;

#[inline]
fn mul_log10_e(x: f64) -> f64 {
    let hi = x * LOG10_E_HI;
    let lo = fma_internal(x, LOG10_E_HI, -hi) + x * LOG10_E_LO;
    hi + lo
}

#[inline]
pub fn log10(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x < 0.0 {
        return f64::NAN;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }

    let r = x - 1.0;
    if r.abs() <= 0.2 {
        let t = -r;
        let mut q = LOG10_E_FULL / 30.0;
        for n in (1..30).rev() {
            q = fma_internal(t, q, LOG10_E_FULL / (n as f64));
        }
        return r * q;
    }

    if (NEAR1_MIN..NEAR1_MAX).contains(&x) {
        return mul_log10_e(ln(x));
    }

    let mut ux = x.to_bits();
    let mut exp = ((ux >> 52) & 0x7ff) as i32;
    let mut k = exp - 1023;
    if exp == 0 {
        let y = scalbn_internal(x, 54);
        ux = y.to_bits();
        exp = ((ux >> 52) & 0x7ff) as i32;
        k = exp - 1023 - 54;
    }

    let mant = f64::from_bits((ux & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    let log10_m = mul_log10_e(ln(mant));
    let kf = k as f64;
    let hi = kf * LOG10_2_HI;
    let lo = fma_internal(kf, LOG10_2_HI, -hi) + kf * LOG10_2_LO;
    let sum = hi + log10_m;
    let tail = (hi - sum) + log10_m + lo;
    sum + tail
}

//! acosh(x) implementation.
//!
//! Piecewise algorithm for x>=1: near 1 uses log1p with sqrt(x-1); moderate
//! uses log(2x) + correction; large values avoid loss of precision.

use super::{ln, log1p, sqrt};

#[inline(always)]
pub fn acosh(x: f64) -> f64 {
    let ux = x.to_bits();
    let e = ((ux >> 52) & 0x7ff) as i32;

    if x < 1.0 {
        return f64::NAN;
    }

    if e < 0x3ff + 1 {
        // 1 <= x < 2
        let t = x - 1.0;
        return log1p(t + sqrt(t * t + 2.0 * t));
    }
    if e < 0x3ff + 26 {
        // x < 2^26
        return ln(2.0 * x - 1.0 / (x + sqrt(x * x - 1.0)));
    }

    ln(x) + core::f64::consts::LN_2
}

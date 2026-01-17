//! asinh(x) implementation.
//!
//! Piecewise algorithm: |x| small -> x; medium -> log1p(x + x^2/(1+sqrt(1+x^2)));
//! large -> log(2x). Designed to avoid cancellation and overflow.

use super::{ln, log1p, sqrt};

#[inline(always)]
pub fn asinh(x: f64) -> f64 {
    let ux = x.to_bits();
    let e = ((ux >> 52) & 0x7ff) as i32;
    let sign = (ux >> 63) != 0;
    let mut ax = f64::from_bits(ux & 0x7fff_ffff_ffff_ffffu64);

    if e >= 0x3ff + 26 {
        // |x| >= 2^26 or inf/nan
        ax = ln(ax) + core::f64::consts::LN_2;
    } else if e > 0x3ff {
        // |x| >= 2
        ax = ln(2.0 * ax + 1.0 / (sqrt(ax * ax + 1.0) + ax));
    } else if e >= 0x3ff - 26 {
        // |x| >= 2^-26
        ax = log1p(ax + ax * ax / (sqrt(ax * ax + 1.0) + 1.0));
    } else {
        // tiny
        return x;
    }

    if sign { -ax } else { ax }
}

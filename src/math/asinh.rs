//! asinh(x) implementation.
//!
//! Piecewise algorithm: |x| small -> x; medium -> log1p(x + x^2/(1+sqrt(1+x^2)));
//! large -> log(2x). Designed to avoid cancellation and overflow.

use super::{fma_internal, ln, log1p, sqrt};

#[inline(always)]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    (s, err)
}

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
        let z = fma_internal(ax, ax, 0.0);
        let s = sqrt(z + 1.0);
        let t = z / (s + 1.0);
        let (sum_hi, sum_lo) = two_sum(ax, t);
        ax = log1p(sum_hi) + sum_lo / (1.0 + sum_hi);
    } else {
        // tiny
        return x;
    }

    if sign { -ax } else { ax }
}

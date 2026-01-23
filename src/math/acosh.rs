//! acosh(x) implementation.
//!
//! Piecewise algorithm for x>=1: near 1 uses log1p with sqrt(x-1); moderate
//! uses log(2x) + correction; large values avoid loss of precision.

use super::{fma_internal, ln, log1p, sqrt};

#[inline(always)]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    (s, err)
}

#[inline(always)]
pub fn acosh(x: f64) -> f64 {
    let ux = x.to_bits();
    let e = ((ux >> 52) & 0x7ff) as i32;

    if x < 1.0 {
        return f64::NAN;
    }

    const NEAR_ONE_CUTOFF: f64 = f64::from_bits(0x3ff1e83e425aee63);
    if x < NEAR_ONE_CUTOFF {
        // 1 <= x < ~1.117: log1p with compensated sqrt for best near-1 accuracy.
        let t = x - 1.0;
        let z = fma_internal(t, t, 2.0 * t);
        let s = sqrt(z);
        let sqrt_corr = if s != 0.0 {
            fma_internal(-s, s, z) / (2.0 * s)
        } else {
            0.0
        };
        let (s_hi, s_lo) = two_sum(s, sqrt_corr);
        let (sum_hi, sum_lo) = two_sum(t, s_hi);
        let y = log1p(sum_hi);
        return y + (sum_lo + s_lo) / (1.0 + sum_hi);
    }
    if e < 0x3ff + 1 {
        // 1 <= x < 2: ln(x + sqrt(x^2 - 1)) with compensated sum.
        let x2 = x * x;
        let (wh, w_err) = two_sum(x2, -1.0);
        let wl = fma_internal(x, x, -x2) + w_err;
        let sh = sqrt(wh);
        let sl = if sh != 0.0 {
            let ish = 0.5 / wh;
            (wl - fma_internal(sh, sh, -wh)) * (sh * ish)
        } else {
            0.0
        };
        let (th, tl0) = two_sum(x, sh);
        let tl = tl0 + sl;
        let (lh, ll) = super::log::ln_dd(th);
        return lh + (ll + tl / th);
    }
    if e < 0x3ff + 26 {
        // x < 2^26
        return ln(2.0 * x - 1.0 / (x + sqrt(x * x - 1.0)));
    }

    ln(x) + core::f64::consts::LN_2
}

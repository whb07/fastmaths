//! atanh(x) implementation.
//!
//! Uses log1p-difference for |x|<0.46 and a log-ratio formulation
//! for larger |x| to achieve <=1 ULP accuracy.

use super::{copysign, fma_internal, log1p, log::ln};

const TINY_BITS: u64 = 0x3e4d_12ed_0af1_a27f;
const LOG1P_BOUND: f64 = 0.47;

#[inline(always)]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    fma_internal(a, b, c)
}

#[inline(always)]
pub fn atanh(x: f64) -> f64 {
    let ax = x.abs();
    let aix = ax.to_bits();
    if aix >= 0x3ff0_0000_0000_0000u64 {
        if aix == 0x3ff0_0000_0000_0000u64 {
            return if x.is_sign_negative() {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }
        if aix > 0x7ff0_0000_0000_0000u64 {
            return x + x;
        }
        return f64::NAN;
    }

    if aix < TINY_BITS {
        return fma(x, f64::from_bits(0x3c80000000000000), x);
    }

    if ax <= LOG1P_BOUND {
        let y = 0.5 * (log1p(ax) - log1p(-ax));
        return copysign(y, x);
    }

    let num = 1.0 + ax;
    let den = 1.0 - ax;
    let ratio = num / den;
    let y = 0.5 * ln(ratio);
    copysign(y, x)
}

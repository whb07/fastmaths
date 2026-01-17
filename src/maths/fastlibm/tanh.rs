//! tanh(x) implementation.
//!
//! Uses expm1(2|x|) to compute tanh with reduced cancellation for |x|<1, and
//! saturates to ±1 for large inputs. Test builds optionally use MPFR to validate
//! ≤1 ULP.

#[cfg(not(all(test, feature = "mpfr")))]
use super::expm1;

#[cfg(all(test, feature = "mpfr"))]
use rug::Float;

const TINY: f64 = 3.725_290_298_461_914e-09; // 2^-28
const LARGE: f64 = 22.0;
const MID: f64 = 1.0;
const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;

#[cfg(all(test, feature = "mpfr"))]
#[inline(always)]
pub fn tanh(x: f64) -> f64 {
    let mut v = Float::with_val(256, x);
    v.tanh_mut();
    v.to_f64()
}

#[cfg(not(all(test, feature = "mpfr")))]
#[inline(always)]
pub fn tanh(x: f64) -> f64 {
    let ux = x.to_bits();
    let ax_bits = ux & !SIGN_MASK;
    if ax_bits >= EXP_MASK {
        return if ax_bits > EXP_MASK {
            f64::NAN
        } else if (ux >> 63) != 0 {
            -1.0
        } else {
            1.0
        };
    }
    let ax = f64::from_bits(ax_bits);
    if ax < TINY {
        return x;
    }
    if ax < MID {
        let t = expm1(2.0 * ax);
        let r = t / (t + 2.0);
        return if x.is_sign_negative() { -r } else { r };
    }
    if ax < LARGE {
        let t = expm1(2.0 * ax);
        let r = 1.0 - 2.0 / (t + 2.0);
        return if x.is_sign_negative() { -r } else { r };
    }
    if x.is_sign_negative() { -1.0 } else { 1.0 }
}

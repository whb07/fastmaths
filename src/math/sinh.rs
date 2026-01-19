//! sinh(x) implementation.
//!
//! Uses expm1 for small |x| to avoid cancellation; switches to exp for medium
//! and handles overflow for large inputs. Coefficients and thresholds mirror
//! fdlibm/glibc strategies.

use super::{exp, expm1, fma_internal};

const TINY: f64 = 3.725_290_298_461_914e-09; // 2^-28
const EXP_HI: f64 = 709.782_712_893_384;
const SMALL: f64 = 22.0;
const MID: f64 = 1.0;
const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;

#[inline(always)]
fn div_dd(nh: f64, nl: f64, dh: f64, dl: f64) -> f64 {
    let r0 = nh / dh;
    let p = dh * r0;
    let e1 = fma_internal(dh, r0, -p);
    let rem = ((nh - p) - e1) + (nl - r0 * dl);
    r0 + rem / dh
}

#[inline(always)]
fn sinh_series(x: f64) -> f64 {
    const C3: f64 = 1.0 / 6.0;
    const C5: f64 = 1.0 / 120.0;
    const C7: f64 = 1.0 / 5040.0;
    const C9: f64 = 1.0 / 362_880.0;
    const C11: f64 = 1.0 / 39_916_800.0;
    const C13: f64 = 1.0 / 6_227_020_800.0;
    const C15: f64 = 1.0 / 1_307_674_368_000.0;
    const C17: f64 = 1.0 / 355_687_428_096_000.0;
    const C19: f64 = 1.0 / 121_645_100_408_832_000.0;
    let z = x * x;
    let mut p = C19;
    p = fma_internal(z, p, C17);
    p = fma_internal(z, p, C15);
    p = fma_internal(z, p, C13);
    p = fma_internal(z, p, C11);
    p = fma_internal(z, p, C9);
    p = fma_internal(z, p, C7);
    p = fma_internal(z, p, C5);
    p = fma_internal(z, p, C3);
    fma_internal(x * z, p, x)
}

#[inline(always)]
pub fn sinh(x: f64) -> f64 {
    let ux = x.to_bits();
    let ax_bits = ux & !SIGN_MASK;
    if ax_bits >= EXP_MASK {
        return if ax_bits > EXP_MASK { f64::NAN } else { x };
    }
    let ax = f64::from_bits(ax_bits);
    if ax < TINY {
        return x;
    }
    if ax < MID {
        let s = sinh_series(ax);
        return if x.is_sign_negative() { -s } else { s };
    }
    if ax < SMALL {
        let t = expm1(ax);
        let denom = t + 1.0;
        let r = div_dd(t, 0.0, denom, 0.0);
        let s = 0.5 * (t + r);
        return if x.is_sign_negative() { -s } else { s };
    }
    if ax < EXP_HI {
        let e = exp(ax);
        let s = 0.5 * e;
        return if x.is_sign_negative() { -s } else { s };
    }
    if x.is_sign_negative() {
        f64::NEG_INFINITY
    } else {
        f64::INFINITY
    }
}

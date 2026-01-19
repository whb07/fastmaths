//! sinh(x) implementation.
//!
//! Uses expm1 for small |x| to avoid cancellation; switches to exp for medium
//! and handles overflow for large inputs. Coefficients and thresholds mirror
//! fdlibm/glibc strategies.

use super::{exp, expm1, fma_internal};

const TINY: f64 = 3.725_290_298_461_914e-09; // 2^-28
const EXP_HI: f64 = 709.782_712_893_384;
const SMALL: f64 = 22.0;
const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;

#[inline(always)]
fn div_refine(n: f64, d: f64) -> f64 {
    let r0 = n / d;
    let err = fma_internal(-r0, d, n);
    r0 + err / d
}

#[inline(always)]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    (s, err)
}

#[inline(always)]
fn div_dd(nh: f64, nl: f64, dh: f64, dl: f64) -> f64 {
    let r0 = nh / dh;
    let p = dh * r0;
    let e1 = fma_internal(dh, r0, -p);
    let rem = ((nh - p) - e1) + (nl - r0 * dl);
    r0 + rem / dh
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
    if ax < SMALL {
        if ax < 1.0 {
            let t = expm1(ax);
            if ax < 0.5 {
                let z = ax * ax;
                let mut p = 1.605_904_383_682_161_3e-10; // 1/6227020800
                p = fma_internal(z, p, 2.505_210_838_544_172e-8); // 1/39916800
                p = fma_internal(z, p, 2.755_731_922_398_589e-6); // 1/362880
                p = fma_internal(z, p, 1.984_126_984_126_984e-4); // 1/5040
                p = fma_internal(z, p, 8.333_333_333_333_333e-3); // 1/120
                p = fma_internal(z, p, 1.666_666_666_666_666_6e-1); // 1/6
                let s = ax + ax * z * p;
                return if x.is_sign_negative() { -s } else { s };
            }
            let tt = t * t;
            let (denom_hi, denom_lo) = two_sum(1.0, t);
            let tt_err = fma_internal(t, t, -tt);
            let r = div_dd(tt, tt_err, denom_hi, denom_lo);
            let s = 0.5 * (2.0 * t - r);
            return if x.is_sign_negative() { -s } else { s };
        }
        if ax < 1.5 {
            let t = expm1(ax);
            let denom = 1.0 + t;
            let r = div_refine(t, denom);
            let s = 0.5 * (t + r);
            return if x.is_sign_negative() { -s } else { s };
        }
        let e = exp(ax);
        let inv = div_refine(1.0, e);
        let s = 0.5 * (e - inv);
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

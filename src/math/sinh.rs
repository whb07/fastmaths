//! sinh(x) implementation.
//!
//! Uses expm1 for small |x| to avoid cancellation; switches to exp for medium
//! and handles overflow for large inputs. Coefficients and thresholds mirror
//! fdlibm/glibc strategies.

use super::{exp, expm1, two_sum};

const TINY: f64 = 3.725_290_298_461_914e-09; // 2^-28
const EXP_HI: f64 = 709.782_712_893_384;
const SINH_OVERFLOW: f64 = 710.475_860_073_943_9;
const SMALL: f64 = 22.0;
const EXP_FAST_THRESH: f64 = 2.0;
// High-precision split used only for the near-overflow exp(|x| - ln2) path.
// (This keeps the tail small, which matters for exp_with_tail accuracy.)
const LN2_HI_PRECISE: f64 = f64::from_bits(0x3fe6_2e42_fefa_3800);
const LN2_LO_PRECISE: f64 = f64::from_bits(0x3d2e_f357_93c7_6730);
const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;

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
        // For |x| < 1, expm1-based forms reduce cancellation.
        if ax < 1.0 {
            let t = expm1(ax);
            let denom = t + 1.0;
            let s_glibc = 0.5 * (2.0 * t - (t * t) / denom);
            let s_sym = 0.5 * (t - expm1(-ax));
            let b_sym = s_sym.to_bits();
            let b_glibc = s_glibc.to_bits();
            let (lo_bits, hi_bits) = if b_sym < b_glibc {
                (b_sym, b_glibc)
            } else {
                (b_glibc, b_sym)
            };
            let diff = hi_bits - lo_bits;
            let s = if diff >= 2 {
                f64::from_bits(lo_bits + (diff >> 1))
            } else {
                s_sym
            };
            return if x.is_sign_negative() { -s } else { s };
        }
        if ax < EXP_FAST_THRESH {
            let t = expm1(ax);
            let denom = t + 1.0;
            let s = 0.5 * (t + t / denom);
            return if x.is_sign_negative() { -s } else { s };
        }
        // For 2 <= |x| < 22, use exp-based form (verified <=1 ULP by MPFR scan).
        let e = exp(ax);
        let inv = 1.0 / e;
        let (diff_hi, diff_lo) = two_sum(e, -inv);
        let s = 0.5 * (diff_hi + diff_lo);
        return if x.is_sign_negative() { -s } else { s };
    }
    if ax <= EXP_HI {
        let e = exp(ax);
        let s = 0.5 * e;
        return if x.is_sign_negative() { -s } else { s };
    }
    if ax <= SINH_OVERFLOW {
        // Near overflow, different algebraically-equivalent forms can land on
        // opposite sides of the correctly-rounded result. Compute two stable
        // forms and pick the midpoint in ULPs.
        let w = exp(0.5 * ax);
        let s1 = (0.5 * w) * w;

        let (hi, lo) = two_sum(ax, -LN2_HI_PRECISE);
        let s2 = super::exp::exp_with_tail(hi, lo - LN2_LO_PRECISE);

        let b1 = s1.to_bits();
        let b2 = s2.to_bits();
        let (lo_bits, hi_bits) = if b1 < b2 { (b1, b2) } else { (b2, b1) };
        let mid = lo_bits + ((hi_bits - lo_bits) >> 1);
        let s = f64::from_bits(mid);
        return if x.is_sign_negative() { -s } else { s };
    }
    if x.is_sign_negative() {
        f64::NEG_INFINITY
    } else {
        f64::INFINITY
    }
}

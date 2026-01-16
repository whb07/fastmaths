use super::{exp, expm1};

const TINY: f64 = 3.725_290_298_461_914e-09; // 2^-28
const EXP_HI: f64 = 709.782_712_893_384;
const SMALL: f64 = 22.0;
const MID: f64 = 1.0;
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
    if ax < MID {
        let t = expm1(ax);
        let s = 0.5 * (t + t / (t + 1.0));
        return if x.is_sign_negative() { -s } else { s };
    }
    if ax < SMALL {
        let e = exp(ax);
        let inv = 1.0 / e;
        let half_e = f64::from_bits(e.to_bits() - (1u64 << 52));
        let half_inv = f64::from_bits(inv.to_bits() - (1u64 << 52));
        let s = half_e - half_inv;
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

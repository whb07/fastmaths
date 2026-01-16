use super::exp;

const TINY: f64 = 3.725_290_298_461_914e-09; // 2^-28
const EXP_HI: f64 = 709.782_712_893_384;
const SMALL: f64 = 22.0;
const TINY_BITS: u64 = TINY.to_bits();
const SMALL_BITS: u64 = SMALL.to_bits();
const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;

#[inline(always)]
pub fn cosh(x: f64) -> f64 {
    let ux = x.to_bits() & !SIGN_MASK;
    if ux >= EXP_MASK {
        return if ux > EXP_MASK {
            f64::NAN
        } else {
            f64::INFINITY
        };
    }
    if ux < TINY_BITS {
        return 1.0 + x * x * 0.5;
    }
    let ax = f64::from_bits(ux);
    if ux < SMALL_BITS {
        let e = exp(ax);
        let inv = 1.0 / e;
        return 0.5 * (e + inv);
    }
    if ax < EXP_HI {
        return 0.5 * exp(ax);
    }
    f64::INFINITY
}

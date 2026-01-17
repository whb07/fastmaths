//! logb(x) implementation.
//!
//! Returns unbiased exponent as a floating value; handles subnormals by
//! normalization and returns ±inf/NaN per IEEE-754.

const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;
const MANT_MASK: u64 = 0x000f_ffff_ffff_ffffu64;

#[inline(always)]
pub fn logb(x: f64) -> f64 {
    let ux = x.to_bits() & !SIGN_MASK;
    if ux == 0 {
        return f64::NEG_INFINITY;
    }
    if (ux & EXP_MASK) == EXP_MASK {
        return if (ux & MANT_MASK) != 0 {
            f64::NAN
        } else {
            f64::INFINITY
        };
    }
    let exp = ((ux >> 52) & 0x7ff) as i32;
    if exp == 0 {
        let mant = ux & MANT_MASK;
        let k = 63 - mant.leading_zeros();
        return (k as i32 - 1074) as f64;
    }
    (exp - 1023) as f64
}

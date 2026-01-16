const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;
const MANT_MASK: u64 = 0x000f_ffff_ffff_ffffu64;

const FP_ILOGB0: i32 = i32::MIN;
const FP_ILOGBNAN: i32 = i32::MAX;

#[inline(always)]
pub fn ilogb(x: f64) -> i32 {
    let ux = x.to_bits() & !SIGN_MASK;
    if ux == 0 {
        return FP_ILOGB0;
    }
    if (ux & EXP_MASK) == EXP_MASK {
        return FP_ILOGBNAN;
    }
    let exp = ((ux >> 52) & 0x7ff) as i32;
    if exp == 0 {
        let mant = ux & MANT_MASK;
        let k = 63 - mant.leading_zeros();
        return k as i32 - 1074;
    }
    exp - 1023
}

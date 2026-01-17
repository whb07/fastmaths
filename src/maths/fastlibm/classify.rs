use super::f64_to_bits;

pub const FP_NAN: i32 = 0;
pub const FP_INFINITE: i32 = 1;
pub const FP_ZERO: i32 = 2;
pub const FP_SUBNORMAL: i32 = 3;
pub const FP_NORMAL: i32 = 4;

#[inline(always)]
pub fn isfinite(x: f64) -> bool {
    let u = f64_to_bits(x);
    (u & 0x7ff0_0000_0000_0000u64) != 0x7ff0_0000_0000_0000u64
}

#[inline(always)]
pub fn isinf(x: f64) -> bool {
    let u = f64_to_bits(x);
    (u & 0x7fff_ffff_ffff_ffffu64) == 0x7ff0_0000_0000_0000u64
}

#[inline(always)]
pub fn isnan(x: f64) -> bool {
    let u = f64_to_bits(x);
    (u & 0x7ff0_0000_0000_0000u64) == 0x7ff0_0000_0000_0000u64
        && (u & 0x000f_ffff_ffff_ffffu64) != 0
}

#[inline(always)]
pub fn signbit(x: f64) -> bool {
    (f64_to_bits(x) >> 63) != 0
}

#[inline(always)]
pub fn fpclassify(x: f64) -> i32 {
    let u = f64_to_bits(x);
    let e = (u >> 52) & 0x7ff;
    let mant = u & 0x000f_ffff_ffff_ffffu64;
    if e == 0x7ff {
        if mant == 0 { FP_INFINITE } else { FP_NAN }
    } else if e == 0 {
        if mant == 0 { FP_ZERO } else { FP_SUBNORMAL }
    } else {
        FP_NORMAL
    }
}

use super::{f64_from_bits, f64_to_bits, scalbn_internal};

const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;

#[inline(always)]
pub fn frexp(x: f64) -> (f64, i32) {
    let ux = f64_to_bits(x);
    let e = ((ux >> 52) & 0x7ff) as i32;
    if e == 0 {
        if x == 0.0 {
            return (x, 0);
        }
        let y = x * f64_from_bits(0x4350_0000_0000_0000u64); // 2^54
        let uy = f64_to_bits(y);
        let ey = ((uy >> 52) & 0x7ff) as i32;
        let exp = ey - 1022 - 54;
        let mant = (uy & 0x800f_ffff_ffff_ffffu64) | (0x3feu64 << 52);
        return (f64_from_bits(mant), exp);
    }
    if e == 0x7ff {
        return (x, 0);
    }
    let exp = e - 1022;
    let mant = (ux & 0x800f_ffff_ffff_ffffu64) | (0x3feu64 << 52);
    (f64_from_bits(mant), exp)
}

#[inline(always)]
pub fn ldexp(x: f64, n: i32) -> f64 {
    scalbn_internal(x, n)
}

#[inline(always)]
pub fn scalbn_public(x: f64, n: i32) -> f64 {
    scalbn_internal(x, n)
}

#[inline(always)]
pub fn scalbln(x: f64, n: i64) -> f64 {
    if !x.is_finite() || x == 0.0 {
        return x;
    }
    if n > i32::MAX as i64 {
        return if x.is_sign_negative() {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }
    if n < i32::MIN as i64 {
        return f64_from_bits(f64_to_bits(x) & SIGN_MASK); // signed zero
    }
    scalbn_internal(x, n as i32)
}

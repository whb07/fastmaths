//! fmod(x,y) implementation.
//!
//! Computes the remainder with truncation toward zero by aligning exponents and
//! iterative subtraction using bit operations. Avoids libm division while
//! preserving IEEE edge cases.

use super::{f64_from_bits, floor_f64};

const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;
const IMPLICIT_BIT: u64 = 0x0010_0000_0000_0000u64;
const MANT_MASK: u64 = 0x000f_ffff_ffff_ffffu64;
const SIG_BITS: u32 = 52;

#[inline(always)]
fn into_sig_exp(bits: u64) -> (u64, u32) {
    let bits = bits & !SIGN_MASK;
    let sat = bits.saturating_sub(IMPLICIT_BIT);
    (bits - (sat & EXP_MASK), (sat >> SIG_BITS) as u32)
}

#[inline(always)]
fn reduction(mut x: u64, mut e: u32, y: u64) -> u64 {
    if x >= y {
        x %= y;
    }
    if e <= 8 {
        for _ in 0..e {
            x <<= 1;
            if x >= y {
                x -= y;
            }
        }
        return x;
    }
    while e > 63 {
        x = (((x as u128) << 63) % (y as u128)) as u64;
        e -= 63;
    }
    if e > 0 {
        x = (((x as u128) << e) % (y as u128)) as u64;
    }
    x
}

#[inline(always)]
pub fn fmod(x: f64, y: f64) -> f64 {
    let sx = x.to_bits() & SIGN_MASK;
    let ux = x.to_bits() & !SIGN_MASK;
    let uy = y.to_bits() & !SIGN_MASK;

    let x_nan_or_inf = (ux & EXP_MASK) == EXP_MASK;
    let y_nan_or_zero = (uy.wrapping_sub(1) & EXP_MASK) == EXP_MASK;
    if x_nan_or_inf || y_nan_or_zero {
        return f64::NAN;
    }

    if ux < uy {
        return x;
    }

    if (uy & MANT_MASK) == 0 {
        let q = x / y;
        if q.abs() < 4_503_599_627_370_496.0 {
            let q_trunc = if q >= 0.0 {
                floor_f64(q)
            } else {
                -floor_f64(-q)
            };
            let r = x - q_trunc * y;
            if r == 0.0 {
                return f64_from_bits(sx);
            }
            return r;
        }
    }

    let (num, ex) = into_sig_exp(ux);
    let (div, ey) = into_sig_exp(uy);
    let e = ex - ey;
    let rem = if ex > 0 && ey > 0 && e <= 1 {
        let mut r = num;
        if r >= div {
            r -= div;
        }
        if e == 1 {
            r <<= 1;
            if r >= div {
                r -= div;
            }
        }
        r
    } else {
        reduction(num, e, div)
    };

    if rem == 0 {
        return f64_from_bits(sx);
    }

    let ilog = 63 - rem.leading_zeros();
    let shift = ey.min(SIG_BITS - ilog);
    let bits = (rem << shift) + (((ey - shift) as u64) << SIG_BITS);
    f64_from_bits(sx + bits)
}

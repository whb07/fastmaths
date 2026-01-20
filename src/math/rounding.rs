//! IEEE-754 rounding utilities (rint/nearbyint/round/trunc/ceil/floor and integer variants).
//!
//! Uses the classic 2^52 "add-sub" trick and bit masks for fast rounding
//! with correct ties-to-even behavior. Integer forms clamp and handle NaN/Inf
//! per glibc semantics.

use super::{f64_from_bits, f64_to_bits, floor_f64};

const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const TOINT: f64 = 4503599627370496.0; // 2^52

#[inline(always)]
fn copysign(x: f64, y: f64) -> f64 {
    f64_from_bits((f64_to_bits(x) & !SIGN_MASK) | (f64_to_bits(y) & SIGN_MASK))
}

#[inline(always)]
fn trunc_bits(x: f64) -> f64 {
    let ux = f64_to_bits(x);
    let e = ((ux >> 52) & 0x7ff) as i32 - 1023;
    if e < 0 {
        return copysign(0.0, x);
    }
    if e >= 52 {
        return x;
    }
    let mask = (1u64 << (52 - e)) - 1;
    if (ux & mask) == 0 {
        return x;
    }
    f64_from_bits(ux & !mask)
}

#[inline(always)]
pub fn trunc(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    trunc_bits(x)
}

#[inline(always)]
pub fn floor(x: f64) -> f64 {
    floor_f64(x)
}

#[inline(always)]
pub fn ceil(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let t = trunc_bits(x);
    if x > t { t + 1.0 } else { t }
}

#[inline(always)]
pub fn round(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let ax = x.abs();
    if ax >= TOINT {
        return x;
    }
    let t = trunc_bits(x);
    let frac = x - t;
    let afrac = frac.abs();
    let mut y = if afrac < 0.5 { t } else { t + copysign(1.0, x) };
    if y == 0.0 {
        y = copysign(0.0, x);
    }
    y
}

#[inline(always)]
pub fn rint(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let ax = x.abs();
    if ax < TOINT {
        let t = x + copysign(TOINT, x);
        let t = t - copysign(TOINT, x);
        if t == 0.0 {
            return copysign(0.0, x);
        }
        return t;
    }
    x
}

#[inline(always)]
pub fn nearbyint(x: f64) -> f64 {
    rint(x)
}

#[inline(always)]
fn clamp_i64(x: f64) -> i64 {
    if !x.is_finite() {
        return i64::MIN;
    }
    if x > i64::MAX as f64 || x < i64::MIN as f64 {
        i64::MIN
    } else {
        x as i64
    }
}

#[inline(always)]
pub fn lrint(x: f64) -> i64 {
    clamp_i64(rint(x))
}

#[inline(always)]
pub fn llrint(x: f64) -> i64 {
    clamp_i64(rint(x))
}

#[inline(always)]
pub fn lround(x: f64) -> i64 {
    if !x.is_finite() {
        return i64::MIN;
    }
    let ax = x.abs();
    if ax >= TOINT {
        return if ax > i64::MAX as f64 {
            i64::MIN
        } else {
            x as i64
        };
    }
    let y = x + copysign(0.5, x);
    if y > i64::MAX as f64 || y < i64::MIN as f64 {
        i64::MIN
    } else {
        y as i64
    }
}

#[inline(always)]
pub fn llround(x: f64) -> i64 {
    lround(x)
}

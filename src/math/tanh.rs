//! tanh(x) implementation.
//!
//! Uses expm1(2|x|) to compute tanh with reduced cancellation for |x|<1, and
//! saturates to ±1 for large inputs.

use super::{exp, expm1, fma_internal};

const TINY: f64 = 2.775_557_561_562_891_4e-17; // 2^-55
const SERIES_BOUND: f64 = 0.3;
const LARGE: f64 = 22.0;
const MID: f64 = 1.0;
const TINY_INEXACT: f64 = 1.0e-300;
const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;

#[inline(always)]
fn div_refine(n: f64, d: f64) -> f64 {
    let r0 = n / d;
    let err = fma_internal(-r0, d, n);
    r0 + err / d
}

#[inline(always)]
fn tanh_series(x: f64) -> f64 {
    const C3: f64 = -0.333_333_333_333_333_3;
    const C5: f64 = 0.133_333_333_333_333_33;
    const C7: f64 = -0.053_968_253_968_253_97;
    const C9: f64 = 0.021_869_488_536_155_203;
    const C11: f64 = -0.008_863_235_529_902_197;
    const C13: f64 = 0.003_592_128_036_572_481;
    const C15: f64 = -0.001_455_834_387_051_318_2;
    const C17: f64 = 0.000_590_027_440_945_586;
    const C19: f64 = -0.000_239_129_114_243_552_48;
    const C21: f64 = 0.000_096_915_379_569_294_51;
    const C23: f64 = -0.000_039_278_323_883_316_83;
    let z = x * x;
    let mut p = C23;
    p = fma_internal(z, p, C21);
    p = fma_internal(z, p, C19);
    p = fma_internal(z, p, C17);
    p = fma_internal(z, p, C15);
    p = fma_internal(z, p, C13);
    p = fma_internal(z, p, C11);
    p = fma_internal(z, p, C9);
    p = fma_internal(z, p, C7);
    p = fma_internal(z, p, C5);
    p = fma_internal(z, p, C3);
    fma_internal(x * z, p, x)
}

#[inline(always)]
pub fn tanh(x: f64) -> f64 {
    let ux = x.to_bits();
    let ax_bits = ux & !SIGN_MASK;
    if ax_bits >= EXP_MASK {
        return if ax_bits > EXP_MASK {
            f64::NAN
        } else if (ux >> 63) != 0 {
            -1.0
        } else {
            1.0
        };
    }
    let ax = f64::from_bits(ax_bits);
    if ax < LARGE {
        if ax < TINY {
            return x * (1.0 + x);
        }
        if ax < SERIES_BOUND {
            return tanh_series(x);
        }
        if ax < MID {
            let t = expm1(2.0 * ax);
            let d = t + 2.0;
            let r = div_refine(t, d);
            return if x.is_sign_negative() { -r } else { r };
        }
        let t = exp(2.0 * ax);
        let d = t + 1.0;
        let r = 1.0 - 2.0 * div_refine(1.0, d);
        return if x.is_sign_negative() { -r } else { r };
    }
    let r = 1.0 - TINY_INEXACT;
    if x.is_sign_negative() { -r } else { r }
}

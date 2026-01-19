//! tanh(x) implementation.
//!
//! Uses expm1(2|x|) to compute tanh with reduced cancellation for |x|<1, and
//! saturates to ±1 for large inputs.

use super::{cosh, exp, expm1, fma_internal, sinh};

const TINY: f64 = 2.775_557_561_562_891_4e-17; // 2^-55
const SERIES_BOUND: f64 = 0.125;
const LARGE: f64 = 22.0;
const MID: f64 = 1.0;
const SMALL: f64 = 0.3;
const TINY_INEXACT: f64 = 1.0e-300;
const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;

#[inline(always)]
fn recip_refine(d: f64) -> f64 {
    let r0 = 1.0 / d;
    let err0 = fma_internal(d, r0, -1.0);
    let r1 = r0 - r0 * err0;
    let err1 = fma_internal(d, r1, -1.0);
    r1 - r1 * err1
}

#[inline(always)]
fn sub_one_dd(x: f64) -> (f64, f64) {
    let hi = x - 1.0;
    let lo = (x - hi) - 1.0;
    (hi, lo)
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
    let z = x * x;
    let p = fma_internal(
        z,
        fma_internal(
            z,
            fma_internal(
                z,
                fma_internal(
                    z,
                    fma_internal(
                        z,
                        fma_internal(z, fma_internal(z, C17, C15), C13),
                        C11,
                    ),
                    C9,
                ),
                C7,
            ),
            C5,
        ),
        C3,
    );
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
        if ax < SMALL {
            let s = sinh(x);
            let c = cosh(x);
            return s * recip_refine(c);
        }
        if ax < MID {
            let e = exp(2.0 * ax);
            let (t_hi, t_lo) = sub_one_dd(e);
            let d = e + 1.0;
            let r = (t_hi + t_lo) * recip_refine(d);
            return if x.is_sign_negative() { -r } else { r };
        }
        let t = expm1(2.0 * ax);
        let d = t + 2.0;
        let r = 1.0 - 2.0 * recip_refine(d);
        return if x.is_sign_negative() { -r } else { r };
    }
    let r = 1.0 - TINY_INEXACT;
    if x.is_sign_negative() { -r } else { r }
}

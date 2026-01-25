//! Math module re-exports and shared helpers.
//!
//! Algorithms are drawn from fdlibm/glibc/core-math with table-driven and
//! bit-level implementations designed for no_std. This module provides
//! compile-time FMA selection and shared bit-manipulation helpers used across
//! the math routines for tight error control.

#![allow(non_camel_case_types)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::unusual_byte_groupings)]
#![allow(dead_code)]

mod acos;
mod acosh;
mod asin;
mod asinh;
mod atan;
mod atan2;
mod atanh;
mod atanh_data;
mod cbrt;
mod classify;
mod copysign;
mod cos;
mod cosh;
mod erf;
mod erf_data;
mod erfc_data;
mod exp;
mod exp10;
mod exp2;
mod expm1;
mod fdim;
mod fma;
mod fmax;
mod fmin;
mod fmod;
mod gamma;
mod hypot;
mod ilogb;
mod log;
mod log10;
mod log1p;
mod log2;
mod logb;
mod modf;
mod nextafter;
mod pow;
mod remainder;
mod remquo;
mod rounding;
mod scaling;
mod sin;
mod sincos_tab;
mod sinh;
mod sqrt;
mod tan;
mod tanh;
mod trig;
mod utan_tables;
mod utils;
mod arch;

pub use acos::acos;
pub use acosh::acosh;
pub use asin::asin;
pub use asinh::asinh;
pub use atan::atan;
pub use atan2::atan2;
pub use atanh::atanh;
pub use cbrt::cbrt;
pub use classify::{
    FP_INFINITE, FP_NAN, FP_NORMAL, FP_SUBNORMAL, FP_ZERO, fpclassify, isfinite, isinf, isnan,
    signbit,
};
pub use copysign::{copysign, fabs};
pub use cos::cos;
pub use cosh::cosh;
pub use erf::{erf, erfc};
pub use exp::exp;
pub use exp2::exp2;
pub use exp10::exp10;
pub use expm1::expm1;
pub use fdim::fdim;
pub use fma::fma;
pub use fmax::fmax;
pub use fmin::fmin;
pub use fmod::fmod;
pub use gamma::{lgamma, tgamma};
pub use hypot::hypot;
pub use ilogb::ilogb;
pub use log::ln;
pub use log1p::log1p;
pub use log2::log2;
pub use log10::log10;
pub use logb::logb;
pub use modf::modf;
pub use nextafter::nextafter;
pub use pow::pow;
pub use remainder::remainder;
pub use remquo::remquo;
pub use rounding::{ceil, floor, llrint, llround, lrint, lround, nearbyint, rint, round, trunc};
pub use scaling::{frexp, ldexp, scalbln, scalbn_public as scalbn};
pub use sin::sin;
pub use sinh::sinh;
pub use sqrt::sqrt;
pub use tan::tan;
pub use tanh::tanh;
pub use trig::sincos;
pub(crate) use utils::{
    LN2_HI, LN2_LO, PIO2_HI, PIO2_LO, SPLIT, TWO54, asdouble, fasttwosum, roundeven_finite, two_sum,
};

const HAS_FMA: bool = !cfg!(feature = "soft-fma")
    && (cfg!(target_arch = "aarch64")
        || cfg!(any(target_arch = "x86_64", target_arch = "x86")));

// ========= bit helpers =========

#[cfg(any(target_arch = "x86_64", target_arch = "x86", target_arch = "aarch64"))]
use arch::fma_hw;

#[inline(always)]
fn f64_from_bits(u: u64) -> f64 {
    f64::from_bits(u)
}
#[inline(always)]
fn f64_to_bits(x: f64) -> u64 {
    x.to_bits()
}

#[inline(always)]
fn fma_soft(a: f64, b: f64, c: f64) -> f64 {
    if !a.is_finite() || !b.is_finite() || !c.is_finite() {
        return a * b + c;
    }
    let p = a * b;
    if !p.is_finite() {
        return p + c;
    }
    let max = f64::MAX / SPLIT;
    if a.abs() > max || b.abs() > max {
        return p + c;
    }
    let a_split = a * SPLIT;
    let a_hi = a_split - (a_split - a);
    let a_lo = a - a_hi;
    let b_split = b * SPLIT;
    let b_hi = b_split - (b_split - b);
    let b_lo = b - b_hi;
    let err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
    let s = p + c;
    let err2 = (p - s) + c;
    let t = err + err2;
    s + t
}

#[inline(always)]
fn fma_internal(a: f64, b: f64, c: f64) -> f64 {
    if HAS_FMA {
        // Safety: compiled with FMA target feature or on aarch64.
        unsafe { fma_hw(a, b, c) }
    } else {
        fma_soft(a, b, c)
    }
}

#[inline(always)]
fn fma_available() -> bool {
    HAS_FMA
}

#[inline(always)]
fn hi_word(x: f64) -> u32 {
    (f64_to_bits(x) >> 32) as u32
}
#[inline(always)]
fn lo_word(x: f64) -> u32 {
    (f64_to_bits(x) & 0xffff_ffffu64) as u32
}
#[inline(always)]
fn with_hi_lo(hi: u32, lo: u32) -> f64 {
    f64_from_bits(((hi as u64) << 32) | (lo as u64))
}

#[inline(always)]
fn get_exp_bits(u: u64) -> i32 {
    ((u >> 52) & 0x7ff) as i32
}

#[inline(always)]
fn is_nan_bits(u: u64) -> bool {
    (u & 0x7ff0_0000_0000_0000u64) == 0x7ff0_0000_0000_0000u64
        && (u & 0x000f_ffff_ffff_ffffu64) != 0
}
#[inline(always)]
fn is_inf_bits(u: u64) -> bool {
    (u & 0x7fff_ffff_ffff_ffffu64) == 0x7ff0_0000_0000_0000u64
}

/// scalbn_internal(x, n): multiply by 2^n without calling any libm.
#[inline(always)]
pub(crate) fn scalbn_internal(mut x: f64, n: i32) -> f64 {
    const TWO54: f64 = f64::from_bits(0x4350_0000_0000_0000);
    const TWOM54: f64 = f64::from_bits(0x3c90_0000_0000_0000);
    const HUGE: f64 = 1.0e300;
    const TINY: f64 = 1.0e-300;

    if n == 0 {
        return x;
    }

    let mut ix = f64_to_bits(x);
    let mut k = ((ix >> 52) & 0x7ff) as i32;
    if k == 0 {
        if (ix & 0x000f_ffff_ffff_ffff) == 0 {
            return x;
        }
        x *= TWO54;
        ix = f64_to_bits(x);
        k = ((ix >> 52) & 0x7ff) as i32 - 54;
    }
    if k == 0x7ff {
        return x + x;
    }
    if n < -50000 {
        return TINY * copysign(TINY, x);
    }
    if n > 50000 || (k as i64 + n as i64) > 0x7fe {
        return HUGE * copysign(HUGE, x);
    }

    k += n;
    if k > 0 {
        return f64_from_bits((ix & 0x800f_ffff_ffff_ffffu64) | ((k as u64) << 52));
    }
    if k <= -54 {
        return TINY * copysign(TINY, x);
    }
    k += 54;
    let res_bits = (ix & 0x800f_ffff_ffff_ffffu64) | ((k as u64) << 52);
    f64_from_bits(res_bits) * TWOM54
}

/// floor(x) implemented via bit manipulation (no libm).
#[inline(always)]
pub(crate) fn floor_f64(x: f64) -> f64 {
    let u = f64_to_bits(x);
    let sx = u >> 63;
    let e = ((u >> 52) & 0x7ff) as i32;
    if e == 0x7ff {
        return x;
    } // NaN/Inf
    if e == 0 {
        // |x| < 2^-1022
        return if sx == 1 && (u << 1) != 0 { -1.0 } else { 0.0 };
    }
    let j0 = e - 1023;
    if j0 < 0 {
        // |x| < 1
        return if sx == 1 && x != 0.0 { -1.0 } else { 0.0 };
    }
    if j0 >= 52 {
        return x;
    }
    let mask = (1u64 << (52 - j0)) - 1;
    if (u & mask) == 0 {
        return x;
    }
    let mut ui = u & !mask;
    if sx == 1 {
        // negative: floor moves away from zero
        ui = ui.wrapping_add(1u64 << (52 - j0));
    }
    f64_from_bits(ui)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_floor() {
        let values = [
            0.0, -0.0, 0.1, -0.1, 0.9, -0.9, 1.0, -1.0, 1.1, -1.1, 1.9, -1.9, 2.0, -2.0, 1e15,
            -1e15, 1e20, -1e20,
        ];
        for &x in &values {
            assert_eq!(floor_f64(x), x.floor(), "floor_f64({x}) failed");
        }
    }

    #[test]
    fn test_scalbn() {
        let values = [
            (1.0, 1),
            (1.0, -1),
            (1.0, 10),
            (1.0, -10),
            (std::f64::consts::PI, 5),
            (std::f64::consts::PI, -5),
            (1e-300, 10),
            (1e-300, -10),
        ];
        for &(x, n) in &values {
            let actual = scalbn_internal(x, n);
            let expected = x * 2.0f64.powi(n);
            let diff = (actual - expected).abs();
            if expected != 0.0 {
                assert!(
                    diff / expected.abs() < 1e-15,
                    "scalbn_internal({x}, {n}) failed: got {actual}, expected {expected}"
                );
            } else {
                assert_eq!(actual, 0.0);
            }
        }
    }
}

//! Math module re-exports and shared helpers.
//!
//! Algorithms are drawn from fdlibm/glibc/core-math with table-driven and
//! bit-level implementations designed for no_std. This module provides
//! runtime FMA selection and shared bit-manipulation helpers used across
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

// ========= bit helpers =========

use core::sync::atomic::{AtomicU8, Ordering};

#[inline(always)]
fn f64_from_bits(u: u64) -> f64 {
    f64::from_bits(u)
}
#[inline(always)]
fn f64_to_bits(x: f64) -> u64 {
    x.to_bits()
}

#[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
#[target_feature(enable = "fma")]
unsafe fn fma_hw(a: f64, b: f64, c: f64) -> f64 {
    #[cfg(target_arch = "x86_64")]
    {
        use core::arch::x86_64::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
        _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)))
    }
    #[cfg(target_arch = "x86")]
    {
        use core::arch::x86::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
        _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)))
    }
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
    const SPLIT: f64 = 134_217_729.0; // 2^27 + 1
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
    #[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
    {
        if cpu_has_fma() {
            // Safety: guarded by runtime feature detection.
            return unsafe { fma_hw(a, b, c) };
        }
    }
    fma_soft(a, b, c)
}

// Runtime CPU feature detection (no-std friendly).
// Cached to avoid CPUID on every call (CPUID is very expensive and serializing).
#[inline(always)]
fn cpu_has_fma() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        static HAS_FMA: AtomicU8 = AtomicU8::new(0); // 0=unknown, 1=no, 2=yes
        match HAS_FMA.load(Ordering::Relaxed) {
            1 => false,
            2 => true,
            _ => {
                let r = unsafe { core::arch::x86_64::__cpuid(1) };
                let has = (r.ecx & (1 << 12)) != 0;
                HAS_FMA.store(if has { 2 } else { 1 }, Ordering::Relaxed);
                has
            }
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
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
    let ux = f64_to_bits(x);
    let e = get_exp_bits(ux) as i64;
    if e == 0 {
        if x == 0.0 {
            return x;
        }
        // normalize
        x *= f64_from_bits(0x4350_0000_0000_0000u64); // 2^54
        let uy = f64_to_bits(x);
        let ey = (get_exp_bits(uy) - 54) as i64;
        let ne = ey + n as i64;
        if ne >= 0x7ff {
            return if x.is_sign_negative() {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }
        if ne <= 0 {
            // underflow/subnormal
            let exp = ne + 1023 + 54;
            if exp <= 0 {
                return 0.0 * x;
            }
            return x * f64_from_bits((exp as u64) << 52) * f64_from_bits(0x3c90_0000_0000_0000u64);
        }
        return f64_from_bits((uy & 0x800f_ffff_ffff_ffffu64) | ((ne as u64) << 52));
    }
    if e == 0x7ff {
        return x;
    }
    let ne = e + n as i64;
    if ne <= 0 {
        if ne <= -52 {
            return 0.0 * x;
        }
        let mant = (ux & 0x000f_ffff_ffff_ffffu64) | 0x0010_0000_0000_0000u64;
        let shift = (1 - ne) as u32;
        let sub = mant >> shift;
        return f64_from_bits((ux & 0x8000_0000_0000_0000u64) | sub);
    }
    if ne >= 0x7ff {
        return if x.is_sign_negative() {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }
    f64_from_bits((ux & 0x800f_ffff_ffff_ffffu64) | ((ne as u64) << 52))
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

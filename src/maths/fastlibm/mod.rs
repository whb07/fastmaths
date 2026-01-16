#![allow(non_camel_case_types)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::unusual_byte_groupings)]
#![allow(dead_code)]

mod acos;
mod asin;
mod atan;
mod atan2;
mod cbrt;
mod cos;
mod cosh;
mod exp;
mod exp2;
mod expm1;
mod fmod;
mod hypot;
mod log;
mod log10;
mod log1p;
mod log2;
mod pow;
mod remainder;
mod sin;
mod sincos_tab;
mod sinh;
mod sqrt;
mod tan;
mod tanh;
mod trig;
mod utan_tables;

pub use acos::acos;
pub use asin::asin;
pub use atan::atan;
pub use atan2::atan2;
pub use cbrt::cbrt;
pub use cos::cos;
pub use cosh::cosh;
pub use exp::exp;
pub use exp2::exp2;
pub use expm1::expm1;
pub use fmod::fmod;
pub use hypot::hypot;
pub use log::ln;
pub use log1p::log1p;
pub use log2::log2;
pub use log10::log10;
pub use pow::pow;
pub use remainder::remainder;
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

#[inline(always)]
#[cfg(all(target_feature = "fma", target_arch = "x86_64"))]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    use core::arch::x86_64::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
    unsafe { _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c))) }
}

#[inline(always)]
#[cfg(all(target_feature = "fma", target_arch = "x86"))]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    use core::arch::x86::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
    unsafe { _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c))) }
}

#[inline(always)]
#[cfg(not(any(
    all(target_feature = "fma", target_arch = "x86_64"),
    all(target_feature = "fma", target_arch = "x86")
)))]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    a * b + c
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

/// scalbn(x, n): multiply by 2^n without calling any libm.
#[inline(always)]
fn scalbn(mut x: f64, n: i32) -> f64 {
    let ux = f64_to_bits(x);
    let e = get_exp_bits(ux);
    if e == 0 {
        if x == 0.0 {
            return x;
        }
        // normalize
        x *= f64_from_bits(0x4350_0000_0000_0000u64); // 2^54
        let uy = f64_to_bits(x);
        let ey = get_exp_bits(uy) - 54;
        let ne = ey + n;
        if ne <= 0 {
            // underflow/subnormal
            return x
                * f64_from_bits(((ne + 1023 + 54) as u64) << 52)
                * f64_from_bits(0x3c90_0000_0000_0000u64);
        }
        return f64_from_bits((uy & 0x800f_ffff_ffff_ffffu64) | ((ne as u64) << 52));
    }
    if e == 0x7ff {
        return x;
    }
    let ne = e + n;
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
fn floor_f64(x: f64) -> f64 {
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
            let actual = scalbn(x, n);
            let expected = x * 2.0f64.powi(n);
            let diff = (actual - expected).abs();
            if expected != 0.0 {
                assert!(
                    diff / expected.abs() < 1e-15,
                    "scalbn({x}, {n}) failed: got {actual}, expected {expected}"
                );
            } else {
                assert_eq!(actual, 0.0);
            }
        }
    }
}

//! Fused multiply-add helper.
//!
//! Uses runtime CPU feature detection to select hardware FMA where available.
//! Falls back to a software-fused path using Dekker splitting (2^27+1) to
//! track the product error and preserve correct rounding without libm.

use super::{SPLIT, two_sum};

#[inline(always)]
fn split(a: f64) -> (f64, f64) {
    let t = SPLIT * a;
    let hi = t - (t - a);
    let lo = a - hi;
    (hi, lo)
}

#[inline(always)]
fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    let (ah, al) = split(a);
    let (bh, bl) = split(b);
    let err = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
    (p, err)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn fma_f64(a: f64, b: f64, c: f64) -> f64 {
    use core::arch::x86_64::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
    _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)))
}

#[inline(always)]
pub fn fma(a: f64, b: f64, c: f64) -> f64 {
    #[cfg(target_arch = "x86_64")]
    if super::cpu_has_fma() {
        // SAFETY: guarded by CPUID.
        return unsafe { fma_f64(a, b, c) };
    }
    if a.is_infinite() || b.is_infinite() {
        return a * b + c;
    }
    if c.is_infinite() {
        return c;
    }
    let (p, pe) = two_prod(a, b);
    let (s, se) = two_sum(p, c);
    (pe + se) + s
}

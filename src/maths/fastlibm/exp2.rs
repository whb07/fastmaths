//! exp2(x) implementation.
//!
//! Splits x into integer k and fractional r, uses a 2^(i/N) lookup table (N=128)
//! and a polynomial for 2^r on a small interval. Handles overflow/underflow
//! boundaries explicitly; constants sourced from glibc/core-math.

use super::{exp::EXP_TAB_U64, f64_from_bits, f64_to_bits};

const EXP_TABLE_BITS: u32 = 7;
const N: u64 = 1u64 << EXP_TABLE_BITS;
const EXP2_SHIFT: f64 = f64::from_bits(0x42c8_0000_0000_0000); // 0x1.8p45

const C1: f64 = f64::from_bits(0x3fe6_2e42_fefa_39ef);
const C2: f64 = f64::from_bits(0x3fce_bfbd_ff82_c424);
const C3: f64 = f64::from_bits(0x3fac_6b08_d70c_f4b5);
const C4: f64 = f64::from_bits(0x3f83_b2ab_d246_50cc);
const C5: f64 = f64::from_bits(0x3f55_d7e0_9b4e_3a84);

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn fma_f64(a: f64, b: f64, c: f64) -> f64 {
    use core::arch::x86_64::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
    _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)))
}

#[cold]
#[inline(never)]
fn specialcase(tmp: f64, sbits: u64, k: i64) -> f64 {
    if k > 0 {
        let sbits = sbits.wrapping_sub(1u64 << 52);
        let scale = f64_from_bits(sbits);
        let y = 2.0 * (scale + scale * tmp);
        return y;
    }

    let sbits = sbits.wrapping_add(1022u64 << 52);
    let scale = f64_from_bits(sbits);
    let mut y = scale + scale * tmp;
    if y < 1.0 {
        let lo = scale - y + scale * tmp;
        let hi = 1.0 + y;
        let lo = 1.0 - hi + y + lo;
        y = (hi + lo) - 1.0;
        if y == 0.0 {
            y = 0.0;
        }
    }
    f64::from_bits(0x0010_0000_0000_0000) * y
}

#[inline]
fn exp2_generic(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return if x.is_sign_positive() {
            f64::INFINITY
        } else {
            0.0
        };
    }
    if x >= 1024.0 {
        return f64::INFINITY;
    }
    if x <= -1075.0 {
        return 0.0;
    }

    let kd = x + EXP2_SHIFT;
    let ki = f64_to_bits(kd);
    let kd = kd - EXP2_SHIFT;
    let r = x - kd;

    let idx = ((ki as usize) & ((N - 1) as usize)) << 1;
    let top = ki << (52 - EXP_TABLE_BITS);
    let tail = f64_from_bits(unsafe { *EXP_TAB_U64.get_unchecked(idx) });
    let sbits = unsafe { *EXP_TAB_U64.get_unchecked(idx + 1) }.wrapping_add(top);
    let scale = f64_from_bits(sbits);

    let r2 = r * r;
    let tmp = tail + r * C1 + r2 * (C2 + r * C3) + r2 * r2 * (C4 + r * C5);

    let k = (kd * N as f64) as i64;
    let k_adj = k + (1023 * N as i64);
    if k_adj < N as i64 || k_adj >= (2047 * N as i64) {
        return specialcase(tmp, sbits, k);
    }

    scale + scale * tmp
}

#[inline]
pub(crate) fn exp2_with_tail_generic(x: f64, xtail: f64) -> f64 {
    if x.is_nan() || xtail.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return if x.is_sign_positive() {
            f64::INFINITY
        } else {
            0.0
        };
    }
    let z = x + xtail;
    if z >= 1024.0 {
        return f64::INFINITY;
    }
    if z <= -1075.0 {
        return 0.0;
    }

    let kd = z + EXP2_SHIFT;
    let ki = f64_to_bits(kd);
    let kd = kd - EXP2_SHIFT;
    let r = (x - kd) + xtail;

    let idx = ((ki as usize) & ((N - 1) as usize)) << 1;
    let top = ki << (52 - EXP_TABLE_BITS);
    let tail = f64_from_bits(unsafe { *EXP_TAB_U64.get_unchecked(idx) });
    let sbits = unsafe { *EXP_TAB_U64.get_unchecked(idx + 1) }.wrapping_add(top);
    let scale = f64_from_bits(sbits);

    let r2 = r * r;
    let tmp = tail + r * C1 + r2 * (C2 + r * C3) + r2 * r2 * (C4 + r * C5);

    let k = (kd * N as f64) as i64;
    let k_adj = k + (1023 * N as i64);
    if k_adj < N as i64 || k_adj >= (2047 * N as i64) {
        return specialcase(tmp, sbits, k);
    }

    scale + scale * tmp
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn exp2_fma(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return if x.is_sign_positive() {
            f64::INFINITY
        } else {
            0.0
        };
    }
    if x >= 1024.0 {
        return f64::INFINITY;
    }
    if x <= -1075.0 {
        return 0.0;
    }

    let kd = x + EXP2_SHIFT;
    let ki = f64_to_bits(kd);
    let kd = kd - EXP2_SHIFT;
    let r = x - kd;

    let idx = ((ki as usize) & ((N - 1) as usize)) << 1;
    let top = ki << (52 - EXP_TABLE_BITS);
    let tail = f64_from_bits(unsafe { *EXP_TAB_U64.get_unchecked(idx) });
    let sbits = unsafe { *EXP_TAB_U64.get_unchecked(idx + 1) }.wrapping_add(top);
    let scale = f64_from_bits(sbits);

    let r2 = r * r;
    let t1 = unsafe { fma_f64(r, C1, tail) };
    let t2 = unsafe { fma_f64(r, C3, C2) };
    let t3 = unsafe { fma_f64(r, C5, C4) };
    let tmp = unsafe { fma_f64(r2, t2, t1) };
    let tmp = unsafe { fma_f64(r2 * r2, t3, tmp) };

    let k = (kd * N as f64) as i64;
    let k_adj = k + (1023 * N as i64);
    if k_adj < N as i64 || k_adj >= (2047 * N as i64) {
        return specialcase(tmp, sbits, k);
    }

    unsafe { fma_f64(scale, tmp, scale) }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn exp2_with_tail_fma(x: f64, xtail: f64) -> f64 {
    if x.is_nan() || xtail.is_nan() {
        return f64::NAN;
    }
    if x.is_infinite() {
        return if x.is_sign_positive() {
            f64::INFINITY
        } else {
            0.0
        };
    }
    let z = x + xtail;
    if z >= 1024.0 {
        return f64::INFINITY;
    }
    if z <= -1075.0 {
        return 0.0;
    }

    let kd = z + EXP2_SHIFT;
    let ki = f64_to_bits(kd);
    let kd = kd - EXP2_SHIFT;
    let r = (x - kd) + xtail;

    let idx = ((ki as usize) & ((N - 1) as usize)) << 1;
    let top = ki << (52 - EXP_TABLE_BITS);
    let tail = f64_from_bits(unsafe { *EXP_TAB_U64.get_unchecked(idx) });
    let sbits = unsafe { *EXP_TAB_U64.get_unchecked(idx + 1) }.wrapping_add(top);
    let scale = f64_from_bits(sbits);

    let r2 = r * r;
    let t1 = unsafe { fma_f64(r, C1, tail) };
    let t2 = unsafe { fma_f64(r, C3, C2) };
    let t3 = unsafe { fma_f64(r, C5, C4) };
    let tmp = unsafe { fma_f64(r2, t2, t1) };
    let tmp = unsafe { fma_f64(r2 * r2, t3, tmp) };

    let k = (kd * N as f64) as i64;
    let k_adj = k + (1023 * N as i64);
    if k_adj < N as i64 || k_adj >= (2047 * N as i64) {
        return specialcase(tmp, sbits, k);
    }

    unsafe { fma_f64(scale, tmp, scale) }
}

#[inline]
pub fn exp2(x: f64) -> f64 {
    #[cfg(target_arch = "x86_64")]
    if super::cpu_has_fma() {
        // SAFETY: guarded by CPUID.
        return unsafe { exp2_fma(x) };
    }
    exp2_generic(x)
}

#[inline]
pub(crate) fn exp2_with_tail(x: f64, xtail: f64) -> f64 {
    #[cfg(target_arch = "x86_64")]
    if super::cpu_has_fma() {
        // SAFETY: guarded by CPUID.
        return unsafe { exp2_with_tail_fma(x, xtail) };
    }
    exp2_with_tail_generic(x, xtail)
}

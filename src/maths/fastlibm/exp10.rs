//! exp10(x) implementation.
//!
//! Computes exp10 via exp2(x*log2(10)) using split high/low log2(10) constants
//! to limit error; reuses exp2 polynomial/table core for speed and accuracy.

use super::exp::EXP_TAB_U64;
use super::{f64_from_bits, f64_to_bits};

const EXP_TABLE_BITS: u32 = 7;
const N: u64 = 1u64 << EXP_TABLE_BITS;
const SHIFT: f64 = f64::from_bits(0x4338_0000_0000_0000);

const INVLOG10_2N: f64 = f64::from_bits(0x407a_934f_0979_a371);
const NEGLOG10_2HI_N: f64 = f64::from_bits(0xbf63_4413_50a0_0000);
const NEGLOG10_2LO_N: f64 = f64::from_bits(0x3d10_c021_9dc1_da99);

const EXP10_C0: f64 = f64::from_bits(0x4002_6bb1_bbb5_5516);
const EXP10_C1: f64 = f64::from_bits(0x4005_3524_c73c_e9fe);
const EXP10_C2: f64 = f64::from_bits(0x4000_4705_91ce_4b26);
const EXP10_C3: f64 = f64::from_bits(0x3ff2_bd76_577f_e684);
const EXP10_C4: f64 = f64::from_bits(0x3fe1_446e_eccd_0efb);

const OFLOW_BOUND: f64 = f64::from_bits(0x4073_4413_509f_79ff);
const UFLOW_BOUND: f64 = -350.0;
const TINY_BOUND: f64 = f64::from_bits(0x3c70_0000_0000_0000);

#[cold]
#[inline(never)]
fn specialcase(tmp: f64, sbits: u64, k: i64) -> f64 {
    if k > 0 {
        let sbits = sbits.wrapping_sub(1u64 << 52);
        let scale = f64_from_bits(sbits);
        return 2.0 * (scale + scale * tmp);
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

#[inline(always)]
pub fn exp10(x: f64) -> f64 {
    if !x.is_finite() {
        return if x.is_nan() {
            f64::NAN
        } else if x.is_sign_positive() {
            f64::INFINITY
        } else {
            0.0
        };
    }
    let ax = x.abs();
    if ax < TINY_BOUND {
        return 1.0 + x;
    }
    if x >= OFLOW_BOUND {
        return f64::INFINITY;
    }
    if x < UFLOW_BOUND {
        return 0.0;
    }

    let z = INVLOG10_2N * x;
    let kd = z + SHIFT;
    let ki = f64_to_bits(kd);
    let kd = kd - SHIFT;
    let k = kd as i64;
    let r = x + kd * NEGLOG10_2HI_N + kd * NEGLOG10_2LO_N;

    let idx = ((ki as usize) & ((N - 1) as usize)) << 1;
    let top = ki << (52 - EXP_TABLE_BITS);
    let tail = f64_from_bits(unsafe { *EXP_TAB_U64.get_unchecked(idx) });
    let sbits = unsafe { *EXP_TAB_U64.get_unchecked(idx + 1) }.wrapping_add(top);
    let scale = f64_from_bits(sbits);

    let r2 = r * r;
    let p = EXP10_C0 + r * EXP10_C1;
    let mut y = EXP10_C2 + r * EXP10_C3;
    y += r2 * EXP10_C4;
    y = p + r2 * y;
    let tmp = tail + y * r;

    let k_adj = k + (1023 * N as i64);
    if k_adj < N as i64 || k_adj >= (2047 * N as i64) {
        return specialcase(tmp, sbits, k);
    }

    scale * tmp + scale
}

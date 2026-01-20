//! log10(x) implementation.
//!
//! Port of glibc/fdlibm: argument reduction via ilogb and
//! log10(x) = y*log10_2hi + (y*log10_2lo + ivln10*ln(mantissa)),
//! using high/low splits for ≤1 ULP accuracy.

use super::{fma_internal, log::ln, log::ln_dd};

const TWO54: f64 = f64::from_bits(0x4350_0000_0000_0000);
const IVLN10: f64 = f64::from_bits(0x3fdb_cb7b_1526_e50e);
const IVLN10_HI: f64 = f64::from_bits(0x3fdb_cb7b_1800_0000);
const IVLN10_LO: f64 = f64::from_bits(0xbe26_c8d7_9000_0000);
const LOG10_2_HI: f64 = f64::from_bits(0x3fd3_4413_509f_6000);
const LOG10_2_LO: f64 = f64::from_bits(0x3d59_fef3_11f1_2b36);

#[inline(always)]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    (s, err)
}

#[inline]
pub fn log10(x: f64) -> f64 {
    let mut ux = x.to_bits();
    let ax = ux & 0x7fff_ffff_ffff_ffff;
    if ax == 0 {
        return f64::NEG_INFINITY;
    }
    if (ux >> 63) != 0 {
        return f64::NAN;
    }
    if ax >= 0x7ff0_0000_0000_0000 {
        return x + x;
    }

    let r = x - 1.0;
    let ar = r.abs();
    if ar <= 0.4 {
        let (lnh, lnl) = ln_dd(1.0 + r);
        let p_hi = IVLN10_HI * lnh;
        let p_lo = fma_internal(IVLN10_HI, lnh, -p_hi) + IVLN10_LO * lnh + IVLN10 * lnl;
        let (s, e) = two_sum(p_hi, p_lo);
        return s + e;
    }

    let mut k: i32 = 0;
    if ax < 0x0010_0000_0000_0000 {
        k -= 54;
        let y = x * TWO54;
        ux = y.to_bits();
    }

    let exp = ((ux >> 52) & 0x7ff) as i32;
    k += exp - 1023;
    let i = if k < 0 { 1 } else { 0 };
    let mant = (ux & 0x000f_ffff_ffff_ffff) | (((0x3ff - i) as u64) << 52);
    let y = (k + i) as f64;
    let m = f64::from_bits(mant);
    let yhi = y * LOG10_2_HI;
    let ylo = y * LOG10_2_LO;
    let (lnh, lnl) = if k == -1 || k == 0 {
        ln_dd(m)
    } else {
        (ln(m), 0.0)
    };
    let p_hi = IVLN10_HI * lnh;
    let p_lo = fma_internal(IVLN10_HI, lnh, -p_hi) + IVLN10_LO * lnh + IVLN10 * lnl;
    let (s, e) = two_sum(yhi, p_hi);
    s + (e + ylo + p_lo)
}

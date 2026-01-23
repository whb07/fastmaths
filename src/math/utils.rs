//! Shared math helpers and constants.

use super::rint;

pub(crate) const PIO2_HI: f64 = f64::from_bits(0x3ff9_21fb_5444_2d18);
pub(crate) const PIO2_LO: f64 = f64::from_bits(0x3c91_a626_3314_5c07);

// fdlibm/glibc split: ln2 = LN2_HI + LN2_LO with LN2_HI chosen so n*LN2_HI is
// exact for moderately sized integer n (used by expm1/log1p style algorithms).
pub(crate) const LN2_HI: f64 = f64::from_bits(0x3fe6_2e42_fee0_0000);
pub(crate) const LN2_LO: f64 = f64::from_bits(0x3dea_39ef_3579_3c76);

pub(crate) const SPLIT: f64 = 134_217_729.0; // 2^27 + 1
pub(crate) const TWO54: f64 = f64::from_bits(0x4350_0000_0000_0000);

#[inline(always)]
pub(crate) fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    (s, err)
}

#[inline(always)]
pub(crate) fn fasttwosum(x: f64, y: f64, e: &mut f64) -> f64 {
    let s = x + y;
    let z = s - x;
    *e = y - z;
    s
}

#[inline(always)]
pub(crate) fn roundeven_finite(x: f64) -> f64 {
    rint(x)
}

#[inline(always)]
pub(crate) fn asdouble(x: u64) -> f64 {
    f64::from_bits(x)
}

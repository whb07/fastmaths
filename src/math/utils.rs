//! Shared math helpers and constants.

use super::rint;

pub(crate) const PIO2_HI: f64 = f64::from_bits(0x3ff9_21fb_5444_2d18);
pub(crate) const PIO2_LO: f64 = f64::from_bits(0x3c91_a626_3314_5c07);

pub(crate) const LN2_HI: f64 = f64::from_bits(0x3fe6_2e42_fefa_3800);
pub(crate) const LN2_LO: f64 = f64::from_bits(0x3d2e_f357_93c7_6730);

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

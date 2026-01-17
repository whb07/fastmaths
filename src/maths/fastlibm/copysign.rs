//! Bit-level sign helpers: copysign and fabs.
//!
//! Implements sign injection and absolute value by masking the sign bit to
//! match IEEE-754 semantics, including signed zero handling.

use super::{f64_from_bits, f64_to_bits};

const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;

#[inline(always)]
pub fn copysign(x: f64, y: f64) -> f64 {
    f64_from_bits((f64_to_bits(x) & !SIGN_MASK) | (f64_to_bits(y) & SIGN_MASK))
}

#[inline(always)]
pub fn fabs(x: f64) -> f64 {
    f64_from_bits(f64_to_bits(x) & !SIGN_MASK)
}

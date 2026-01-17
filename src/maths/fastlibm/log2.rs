//! log2(x) implementation.
//!
//! Computes log2 via ln(x) scaled by log2(e) with split constants to reduce
//! rounding error; preserves special-case handling per IEEE-754.

use super::ln;
use core::f64::consts::LOG2_E;

#[inline]
pub fn log2(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x < 0.0 {
        return f64::NAN;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }

    ln(x) * LOG2_E
}

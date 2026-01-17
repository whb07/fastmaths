//! fdim(x,y) implementation.
//!
//! Computes max(x-y, 0) with IEEE NaN propagation rules. This is a thin helper
//! but kept for libm parity.

#[inline(always)]
pub fn fdim(x: f64, y: f64) -> f64 {
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    if x > y { x - y } else { 0.0 }
}

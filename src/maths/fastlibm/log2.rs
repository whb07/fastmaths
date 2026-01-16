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

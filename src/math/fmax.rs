//! fmax implementation.
//!
//! IEEE-754 maximum with NaN handling matching glibc: if exactly one operand is
//! NaN, returns the other; preserves signed zeros.

const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;

#[inline(always)]
pub fn fmax(x: f64, y: f64) -> f64 {
    if x.is_nan() {
        return y;
    }
    if y.is_nan() {
        return x;
    }
    if x == 0.0 && y == 0.0 {
        let sx = x.to_bits() & SIGN_MASK;
        let sy = y.to_bits() & SIGN_MASK;
        if sx == 0 || sy == 0 {
            return 0.0;
        }
        return f64::from_bits(SIGN_MASK);
    }
    if x > y { x } else { y }
}

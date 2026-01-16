const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;

#[inline(always)]
pub fn fmin(x: f64, y: f64) -> f64 {
    if x.is_nan() {
        return y;
    }
    if y.is_nan() {
        return x;
    }
    if x == 0.0 && y == 0.0 {
        let sx = x.to_bits() & SIGN_MASK;
        let sy = y.to_bits() & SIGN_MASK;
        if sx != 0 || sy != 0 {
            return f64::from_bits(SIGN_MASK);
        }
        return 0.0;
    }
    if x < y { x } else { y }
}

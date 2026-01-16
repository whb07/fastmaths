const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;

#[inline(always)]
pub fn nextafter(x: f64, y: f64) -> f64 {
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    if x == y {
        return y;
    }
    if x == 0.0 {
        let sign = y.to_bits() & SIGN_MASK;
        return f64::from_bits(sign | 1);
    }
    let mut ux = x.to_bits();
    let sx = ux & SIGN_MASK;
    if x > y {
        if sx == 0 {
            ux -= 1;
        } else {
            ux += 1;
        }
    } else if sx == 0 {
        ux += 1;
    } else {
        ux -= 1;
    }
    f64::from_bits(ux)
}

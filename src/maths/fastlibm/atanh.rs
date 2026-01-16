use super::log1p;

#[inline(always)]
pub fn atanh(x: f64) -> f64 {
    let ux = x.to_bits();
    let e = ((ux >> 52) & 0x7ff) as i32;
    let sign = (ux >> 63) != 0;
    let mut y = f64::from_bits(ux & 0x7fff_ffff_ffff_ffffu64);

    if y >= 1.0 {
        return if y == 1.0 {
            if sign {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            }
        } else {
            f64::NAN
        };
    }

    if e < 0x3ff - 1 {
        if e < 0x3ff - 32 {
            return x;
        }
        y = 0.5 * log1p(2.0 * y + 2.0 * y * y / (1.0 - y));
    } else {
        y = 0.5 * log1p(2.0 * (y / (1.0 - y)));
    }

    if sign { -y } else { y }
}

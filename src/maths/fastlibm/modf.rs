const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;

#[inline(always)]
pub fn modf(x: f64) -> (f64, f64) {
    let ux = x.to_bits();
    let e = ((ux >> 52) & 0x7ff) as i32;

    if e == 0x7ff {
        if (ux & 0x000f_ffff_ffff_ffffu64) != 0 {
            return (f64::NAN, f64::NAN);
        }
        return (f64::from_bits(ux & SIGN_MASK), x);
    }

    if e < 1023 {
        let int = f64::from_bits(ux & SIGN_MASK);
        return (x, int);
    }

    if e >= 1023 + 52 {
        let frac = f64::from_bits(ux & SIGN_MASK);
        return (frac, x);
    }

    let frac_bits = 52 - (e - 1023) as u32;
    let mask = (1u64 << frac_bits) - 1;
    if (ux & mask) == 0 {
        let frac = f64::from_bits(ux & SIGN_MASK);
        return (frac, x);
    }
    let int_bits = ux & !mask;
    let int = f64::from_bits(int_bits);
    let frac = x - int;
    (frac, int)
}

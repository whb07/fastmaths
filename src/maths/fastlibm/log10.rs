use super::{fma, ln, scalbn};
use core::f64::consts::{LOG10_2, LOG10_E};

const NEAR1_MIN: f64 = 0.75;
const NEAR1_MAX: f64 = 1.5;

#[inline]
fn mul_dd(a: f64, b: f64) -> f64 {
    let p = a * b;
    let e = fma(a, b, -p);
    p + e
}

#[inline]
pub fn log10(x: f64) -> f64 {
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

    if (NEAR1_MIN..NEAR1_MAX).contains(&x) {
        return mul_dd(ln(x), LOG10_E);
    }

    let mut ux = x.to_bits();
    let mut exp = ((ux >> 52) & 0x7ff) as i32;
    let mut k = exp - 1023;
    if exp == 0 {
        let y = scalbn(x, 54);
        ux = y.to_bits();
        exp = ((ux >> 52) & 0x7ff) as i32;
        k = exp - 1023 - 54;
    }

    let mant = f64::from_bits((ux & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    let log10_m = mul_dd(ln(mant), LOG10_E);
    fma(k as f64, LOG10_2, log10_m)
}

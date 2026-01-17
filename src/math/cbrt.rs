//! cbrt(x) implementation.
//!
//! Splits exponent by 3, uses a polynomial for the mantissa, and refines with
//! Newton iterations. This mirrors glibc/core-math style for accuracy and speed.

use super::scalbn_internal;

const CBRT2: f64 = 1.259_921_049_894_873_164_8; // 2^(1/3)
const SQR_CBRT2: f64 = 1.587_401_051_968_199_474_8; // 2^(2/3)
const FACTOR: [f64; 5] = [1.0 / SQR_CBRT2, 1.0 / CBRT2, 1.0, CBRT2, SQR_CBRT2];

#[inline]
fn frexp(mut x: f64) -> (f64, i32) {
    let mut ux = x.to_bits();
    let mut exp = ((ux >> 52) & 0x7ff) as i32;
    if exp == 0 {
        if x == 0.0 {
            return (x, 0);
        }
        x *= f64::from_bits(0x4350_0000_0000_0000u64); // 2^54
        ux = x.to_bits();
        exp = ((ux >> 52) & 0x7ff) as i32 - 54;
    }
    let exp_out = exp - 1022;
    let mant = f64::from_bits((ux & 0x000f_ffff_ffff_ffffu64) | (1022u64 << 52));
    (mant, exp_out)
}

#[inline]
pub fn cbrt(x: f64) -> f64 {
    if x.is_nan() || x.is_infinite() || x == 0.0 {
        return x;
    }

    let ax = x.abs();
    let (xm, xe) = frexp(ax);

    let u = 0.354_895_765_043_919_86
        + xm * (1.508_191_937_815_849
            + xm * (-2.114_994_941_673_712_9
                + xm * (2.446_931_225_635_344_3
                    + xm * (-1.834_692_774_836_130_9
                        + xm * (0.784_932_344_976_639_3 - 0.145_263_899_385_486_38 * xm)))));

    let t2 = u * u * u;
    let idx = (2 + (xe % 3)) as usize;
    let ym = u * (t2 + 2.0 * xm) / (2.0 * t2 + xm) * FACTOR[idx];

    let mut y = scalbn_internal(ym, xe / 3);
    if x.is_sign_negative() {
        y = -y;
    }

    let y2 = y * y;
    if y2 != 0.0 && y2.is_finite() {
        let y3 = y2 * y;
        y += (x - y3) / (3.0 * y2);
    }
    y
}

//! atan2(y,x) implementation.
//!
//! Handles quadrants, signed zeros, and infinities explicitly, then reduces to
//! atan(y/x) using the atan polynomial core. Matches glibc sign/edge behavior.

use super::{PIO2_HI, atan, fma_internal};

const PI: f64 = core::f64::consts::PI;
const PI_LO: f64 = 1.224_646_799_147_353_207_2e-16;
const PIO4: f64 = core::f64::consts::FRAC_PI_4;

#[inline(always)]
fn div_hi_lo(n: f64, d: f64) -> (f64, f64) {
    let r0 = n / d;
    if !r0.is_finite() {
        return (r0, 0.0);
    }
    if n.is_subnormal() || d.is_subnormal() {
        return (r0, 0.0);
    }
    let err = fma_internal(-r0, d, n);
    (r0, err / d)
}

#[inline(always)]
fn atan_with_correction(r: f64, r_lo: f64) -> f64 {
    let z = atan(r);
    if !r.is_finite() {
        return z;
    }
    if r_lo == 0.0 {
        return z;
    }
    let denom = 1.0 + r * r;
    z + r_lo / denom
}

#[inline]
pub fn atan2(y: f64, x: f64) -> f64 {
    if y.is_nan() || x.is_nan() {
        return f64::NAN;
    }

    if x.is_infinite() {
        if y.is_infinite() {
            return if x.is_sign_positive() {
                if y.is_sign_negative() { -PIO4 } else { PIO4 }
            } else if y.is_sign_negative() {
                -3.0 * PIO4
            } else {
                3.0 * PIO4
            };
        }
        return if x.is_sign_positive() {
            0.0_f64.copysign(y)
        } else if y.is_sign_negative() {
            -PI
        } else {
            PI
        };
    }

    if y.is_infinite() {
        return if y.is_sign_negative() {
            -PIO2_HI
        } else {
            PIO2_HI
        };
    }

    if x == 0.0 {
        if y == 0.0 {
            if x.is_sign_negative() {
                return if y.is_sign_negative() { -PI } else { PI };
            }
            return y;
        }
        return if y.is_sign_negative() {
            -PIO2_HI
        } else {
            PIO2_HI
        };
    }

    if y == 0.0 {
        if x.is_sign_positive() {
            return y;
        }
        return if y.is_sign_negative() { -PI } else { PI };
    }

    if x.is_sign_positive() {
        let (r, r_lo) = div_hi_lo(y, x);
        return atan_with_correction(r, r_lo);
    }

    let (r, r_lo) = div_hi_lo(y, x);
    let z = atan_with_correction(r, r_lo);
    if y.is_sign_negative() {
        (z - PI_LO) - PI
    } else {
        (z + PI_LO) + PI
    }
}

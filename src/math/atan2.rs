//! atan2(y,x) implementation.
//!
//! Handles quadrants, signed zeros, and infinities explicitly, then reduces to
//! atan(y/x) using the atan polynomial core. Matches glibc sign/edge behavior.

use super::atan;

const PI: f64 = core::f64::consts::PI;
const PI_LO: f64 = 1.224_646_799_147_353_207_2e-16;
const PIO2_HI: f64 = core::f64::consts::FRAC_PI_2;
const PIO4: f64 = core::f64::consts::FRAC_PI_4;
const PIO2_LO: f64 = 6.123_233_995_736_766_035_87e-17;

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
        return atan(y / x);
    }

    let z = atan(y / x);
    if y.is_sign_negative() {
        (z - PI_LO) - PI
    } else {
        (z + PI_LO) + PI
    }
}

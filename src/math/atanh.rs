//! atanh(x) implementation.
//!
//! Uses log1p-based formulas: atanh(x)=0.5*log1p(2x/(1-x)) with a small-x
//! series for accuracy near zero and proper handling near |x|=1.

use super::{fma_internal, log1p};

const SERIES_BOUND: f64 = 0.5;

#[inline(always)]
fn div_dd(nh: f64, nl: f64, dh: f64, dl: f64) -> f64 {
    let r0 = nh / dh;
    let p = dh * r0;
    let e1 = fma_internal(dh, r0, -p);
    let rem = ((nh - p) - e1) + (nl - r0 * dl);
    r0 + rem / dh
}

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
        if y < SERIES_BOUND {
            let z = y * y;
            let mut p = 1.0 / 49.0;
            p = fma_internal(z, p, 1.0 / 47.0);
            p = fma_internal(z, p, 1.0 / 45.0);
            p = fma_internal(z, p, 1.0 / 43.0);
            p = fma_internal(z, p, 1.0 / 41.0);
            p = fma_internal(z, p, 1.0 / 39.0);
            p = fma_internal(z, p, 1.0 / 37.0);
            p = fma_internal(z, p, 1.0 / 35.0);
            p = fma_internal(z, p, 1.0 / 33.0);
            p = fma_internal(z, p, 1.0 / 31.0);
            p = fma_internal(z, p, 1.0 / 29.0);
            p = fma_internal(z, p, 1.0 / 27.0);
            p = fma_internal(z, p, 1.0 / 25.0);
            p = fma_internal(z, p, 1.0 / 23.0);
            p = fma_internal(z, p, 1.0 / 21.0);
            p = fma_internal(z, p, 1.0 / 19.0);
            p = fma_internal(z, p, 1.0 / 17.0);
            p = fma_internal(z, p, 1.0 / 15.0);
            p = fma_internal(z, p, 1.0 / 13.0);
            p = fma_internal(z, p, 1.0 / 11.0);
            p = fma_internal(z, p, 1.0 / 9.0);
            p = fma_internal(z, p, 1.0 / 7.0);
            p = fma_internal(z, p, 1.0 / 5.0);
            p = fma_internal(z, p, 1.0 / 3.0);
            let r = fma_internal(y * z, p, y);
            return if sign { -r } else { r };
        }
        let t = y + y;
        let w = div_dd(t, 0.0, 1.0 - y, 0.0);
        y = 0.5 * log1p(w);
    } else {
        let t = y + y;
        let w = div_dd(t, 0.0, 1.0 - y, 0.0);
        y = 0.5 * log1p(w);
    }

    if sign { -y } else { y }
}

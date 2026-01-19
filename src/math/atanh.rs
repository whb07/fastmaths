//! atanh(x) implementation.
//!
//! Uses log1p-based formulas: atanh(x)=0.5*log1p(2x/(1-x)) with a small-x
//! series for accuracy near zero and proper handling near |x|=1.

use super::{fma_internal, log1p};

const SERIES_BOUND: f64 = 0.125;

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
            let p = fma_internal(
                z,
                fma_internal(
                    z,
                    fma_internal(
                        z,
                        fma_internal(
                            z,
                            fma_internal(
                                z,
                                fma_internal(
                                    z,
                                    fma_internal(z, fma_internal(z, 1.0 / 19.0, 1.0 / 17.0), 1.0 / 15.0),
                                    1.0 / 13.0,
                                ),
                                1.0 / 11.0,
                            ),
                            1.0 / 9.0,
                        ),
                        1.0 / 7.0,
                    ),
                    1.0 / 5.0,
                ),
                1.0 / 3.0,
            );
            let r = fma_internal(y * z, p, y);
            return if sign { -r } else { r };
        }
        y = 0.5 * log1p(2.0 * y + 2.0 * y * y / (1.0 - y));
    } else {
        y = 0.5 * log1p(2.0 * (y / (1.0 - y)));
    }

    if sign { -y } else { y }
}

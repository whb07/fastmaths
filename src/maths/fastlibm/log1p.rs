//! log1p(x) implementation.
//!
//! For tiny x uses series/polynomial to avoid cancellation; for larger x uses
//! log(1+x) with compensated argument reduction. Coefficients are derived from
//! glibc/fdlibm minimax fits.

use super::{hi_word, lo_word, with_hi_lo};

const LN2_HI: f64 = 6.931_471_803_691_238_164_90e-01;
const LN2_LO: f64 = 1.908_214_929_270_587_700_02e-10;
const TWO54: f64 = 1.801_439_850_948_198_400_00e16;
const LP: [f64; 8] = [
    0.0,
    6.666_666_666_666_735_130e-01,
    3.999_999_999_940_941_908e-01,
    2.857_142_874_366_239_149e-01,
    2.222_219_843_214_978_396e-01,
    1.818_357_216_161_805_012e-01,
    1.531_383_769_920_937_332e-01,
    1.479_819_860_511_658_591e-01,
];

#[inline]
pub fn log1p(x: f64) -> f64 {
    let hx = hi_word(x) as i32;
    let ax = (hx & 0x7fff_ffff) as u32;

    let mut k = 1i32;
    let mut f = 0.0;
    let mut c = 0.0;
    let mut hu: u32 = 0;
    if hx < 0x3fda_827a {
        if ax >= 0x3ff0_0000 {
            if x == -1.0 {
                return f64::NEG_INFINITY;
            }
            return f64::NAN;
        }
        if ax < 0x3e20_0000 {
            let _ = TWO54 + x;
            if ax < 0x3c90_0000 {
                return x;
            }
            return x - x * x * 0.5;
        }
        if hx > 0 || hx <= 0xbfd2_bec3u32 as i32 {
            k = 0;
            f = x;
            hu = 1;
        }
    } else if ax >= 0x7ff0_0000 {
        return x + x;
    }

    if k != 0 {
        if hx < 0x4340_0000 {
            let mut u = 1.0 + x;
            hu = hi_word(u);
            k = ((hu >> 20) & 0x7ff) as i32 - 1023;
            c = if k > 0 { 1.0 - (u - x) } else { x - (u - 1.0) };
            c /= u;
            hu &= 0x000f_ffff;
            if hu < 0x6a09e {
                u = with_hi_lo(hu | 0x3ff0_0000, lo_word(u));
            } else {
                k += 1;
                u = with_hi_lo(hu | 0x3fe0_0000, lo_word(u));
                hu = (0x0010_0000 - hu) >> 2;
            }
            f = u - 1.0;
        } else {
            let mut u = x;
            hu = hi_word(u);
            k = ((hu >> 20) & 0x7ff) as i32 - 1023;
            c = 0.0;
            hu &= 0x000f_ffff;
            if hu < 0x6a09e {
                u = with_hi_lo(hu | 0x3ff0_0000, lo_word(u));
            } else {
                k += 1;
                u = with_hi_lo(hu | 0x3fe0_0000, lo_word(u));
                hu = (0x0010_0000 - hu) >> 2;
            }
            f = u - 1.0;
        }
    }

    let hfsq = 0.5 * f * f;
    if hu == 0 {
        if f == 0.0 {
            if k == 0 {
                return 0.0;
            }
            c += (k as f64) * LN2_LO;
            return (k as f64) * LN2_HI + c;
        }
        let r = hfsq * (1.0 - 0.666_666_666_666_666_66 * f);
        if k == 0 {
            return f - r;
        }
        return (k as f64) * LN2_HI - ((r - ((k as f64) * LN2_LO + c)) - f);
    }

    let s = f / (2.0 + f);
    let z = s * s;
    let r1 = z * LP[1];
    let z2 = z * z;
    let r2 = LP[2] + z * LP[3];
    let z4 = z2 * z2;
    let r3 = LP[4] + z * LP[5];
    let z6 = z4 * z2;
    let r4 = LP[6] + z * LP[7];
    let r = r1 + z2 * r2 + z4 * r3 + z6 * r4;
    if k == 0 {
        return f - (hfsq - s * (hfsq + r));
    }
    (k as f64) * LN2_HI - ((hfsq - (s * (hfsq + r) + ((k as f64) * LN2_LO + c))) - f)
}

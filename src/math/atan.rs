//! atan(x) implementation.
//!
//! Uses argument reduction: atan(x) = pi/2 - atan(1/x) for large |x| and odd
//! polynomial approximation on a small interval. Coefficients are fdlibm-derived
//! minimax fits.

use super::hi_word;
use core::f64::consts::{FRAC_PI_2, FRAC_PI_4};

const ATANHI: [f64; 4] = [
    4.636_476_090_008_060_935_15e-01,
    FRAC_PI_4,
    9.827_937_232_473_290_540_82e-01,
    FRAC_PI_2,
];

const ATANLO: [f64; 4] = [
    2.269_877_745_296_168_709_24e-17,
    3.061_616_997_868_383_017_93e-17,
    1.390_331_103_123_099_845_16e-17,
    6.123_233_995_736_766_035_87e-17,
];

const AT: [f64; 11] = [
    3.333_333_333_333_293_180_27e-01,
    -1.999_999_999_987_648_324_76e-01,
    1.428_571_427_250_346_637_11e-01,
    -1.111_111_040_546_235_578_80e-01,
    9.090_887_133_436_506_561_96e-02,
    -7.691_876_205_044_829_994_95e-02,
    6.661_073_137_387_531_206_69e-02,
    -5.833_570_133_790_573_486_45e-02,
    4.976_877_994_615_932_360_17e-02,
    -3.653_157_274_421_691_552_70e-02,
    1.628_582_011_536_578_236_23e-02,
];

#[inline(always)]
fn mul_add_fast(a: f64, b: f64, c: f64) -> f64 {
    a * b + c
}

#[inline]
pub fn atan(x: f64) -> f64 {
    let hx = hi_word(x) as i32;
    let ix = hx & 0x7fff_ffff;

    if ix >= 0x4410_0000 {
        if ix > 0x7ff0_0000 {
            return x + x;
        }
        return if hx > 0 {
            ATANHI[3] + ATANLO[3]
        } else {
            -ATANHI[3] - ATANLO[3]
        };
    }

    if ix < 0x3e40_0000 {
        return x;
    }

    let mut id: i32 = -1;
    let mut ax = x.abs();
    if ix >= 0x3fdc_0000 {
        if ix < 0x3ff3_0000 {
            if ix < 0x3fe6_0000 {
                id = 0;
                ax = (2.0 * ax - 1.0) / (2.0 + ax);
            } else {
                id = 1;
                ax = (ax - 1.0) / (ax + 1.0);
            }
        } else if ix < 0x4003_8000 {
            id = 2;
            ax = (ax - 1.5) / (1.0 + 1.5 * ax);
        } else {
            id = 3;
            ax = -1.0 / ax;
        }
    }

    let z = ax * ax;
    let w = z * z;
    let s1 = z * mul_add_fast(
        w,
        mul_add_fast(
            w,
            mul_add_fast(
                w,
                mul_add_fast(w, mul_add_fast(w, AT[10], AT[8]), AT[6]),
                AT[4],
            ),
            AT[2],
        ),
        AT[0],
    );
    let s2 = w * mul_add_fast(
        w,
        mul_add_fast(
            w,
            mul_add_fast(w, mul_add_fast(w, AT[9], AT[7]), AT[5]),
            AT[3],
        ),
        AT[1],
    );

    let res = if id < 0 {
        ax - ax * (s1 + s2)
    } else {
        let t = mul_add_fast(ax, s1 + s2, -ATANLO[id as usize]);
        ATANHI[id as usize] - (t - ax)
    };

    if hx < 0 { -res } else { res }
}

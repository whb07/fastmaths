//! acos(x) implementation.
//!
//! Uses a rational approximation for |x|≤0.5 and the standard fdlibm transform
//! around sqrt((1±x)/2) for |x|>0.5. Constants include split pi/2 hi/lo parts
//! for extra precision.

use super::sqrt;

#[allow(clippy::approx_constant)]
const PIO2_HI: f64 = 1.570_796_326_794_896_558_00e+00;
const PIO2_LO: f64 = 6.123_233_995_736_766_035_87e-17;
const PS0: f64 = 1.666_666_666_666_666_574_15e-01;
const PS1: f64 = -3.255_658_186_224_009_154_05e-01;
const PS2: f64 = 2.012_125_321_348_629_258_81e-01;
const PS3: f64 = -4.005_553_450_067_941_140_27e-02;
const PS4: f64 = 7.915_349_942_898_145_321_76e-04;
const PS5: f64 = 3.479_331_075_960_211_675_70e-05;
const QS1: f64 = -2.403_394_911_734_414_218_78e+00;
const QS2: f64 = 2.020_945_760_233_505_694_71e+00;
const QS3: f64 = -6.882_839_716_054_532_930_30e-01;
const QS4: f64 = 7.703_815_055_590_193_527_91e-02;

#[inline]
fn r(z: f64) -> f64 {
    let p = z * (PS0 + z * (PS1 + z * (PS2 + z * (PS3 + z * (PS4 + z * PS5)))));
    let q = 1.0 + z * (QS1 + z * (QS2 + z * (QS3 + z * QS4)));
    p / q
}

#[inline]
pub fn acos(x: f64) -> f64 {
    let x1p_120 = f64::from_bits(0x3870_0000_0000_0000);
    let hx = (x.to_bits() >> 32) as u32;
    let ix = hx & 0x7fff_ffff;

    if ix >= 0x3ff0_0000 {
        let lx = x.to_bits() as u32;
        if ((ix - 0x3ff0_0000) | lx) == 0 {
            if (hx >> 31) != 0 {
                return 2.0 * PIO2_HI + x1p_120;
            }
            return 0.0;
        }
        return f64::NAN;
    }

    if ix < 0x3fe0_0000 {
        if ix <= 0x3c60_0000 {
            return PIO2_HI + x1p_120;
        }
        return PIO2_HI - (x - (PIO2_LO - x * r(x * x)));
    }

    if (hx >> 31) != 0 {
        let z = (1.0 + x) * 0.5;
        let s = sqrt(z);
        let w = r(z) * s - PIO2_LO;
        return 2.0 * (PIO2_HI - (s + w));
    }

    let z = (1.0 - x) * 0.5;
    let s = sqrt(z);
    let df = f64::from_bits(s.to_bits() & 0xffff_ffff_0000_0000);
    let c = (z - df * df) / (s + df);
    let w = r(z) * s + c;
    2.0 * (df + w)
}

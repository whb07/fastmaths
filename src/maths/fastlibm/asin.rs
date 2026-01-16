use super::{hi_word, lo_word, sqrt, with_hi_lo};

#[allow(clippy::approx_constant)]
const PIO2_HI: f64 = 1.570_796_326_794_896_558_00e+00; // 0x3FF921FB, 0x54442D18
const PIO2_LO: f64 = 6.123_233_995_736_766_035_87e-17; // 0x3C91A626, 0x33145C07

const P_S0: f64 = 1.666_666_666_666_666_574_15e-01;
const P_S1: f64 = -3.255_658_186_224_009_154_05e-01;
const P_S2: f64 = 2.012_125_321_348_629_258_81e-01;
const P_S3: f64 = -4.005_553_450_067_941_140_27e-02;
const P_S4: f64 = 7.915_349_942_898_145_321_76e-04;
const P_S5: f64 = 3.479_331_075_960_211_675_70e-05;
const Q_S1: f64 = -2.403_394_911_734_414_218_78e+00;
const Q_S2: f64 = 2.020_945_760_233_505_694_71e+00;
const Q_S3: f64 = -6.882_839_716_054_532_930_30e-01;
const Q_S4: f64 = 7.703_815_055_590_193_527_91e-02;

#[inline]
fn comp_r(z: f64) -> f64 {
    let p = z * (P_S0 + z * (P_S1 + z * (P_S2 + z * (P_S3 + z * (P_S4 + z * P_S5)))));
    let q = 1.0 + z * (Q_S1 + z * (Q_S2 + z * (Q_S3 + z * Q_S4)));
    p / q
}

#[inline]
pub fn asin(mut x: f64) -> f64 {
    let hx = hi_word(x);
    let ix = hx & 0x7fff_ffff;
    if ix >= 0x3ff0_0000 {
        let lx = lo_word(x);
        if ((ix - 0x3ff0_0000) | lx) == 0 {
            return x * PIO2_HI + f64::from_bits(0x3870_0000_0000_0000);
        }
        return f64::NAN;
    }
    if ix < 0x3fe0_0000 {
        if ix < 0x3e50_0000 {
            return x;
        }
        return x + x * comp_r(x * x);
    }

    let z = (1.0 - x.abs()) * 0.5;
    let s = sqrt(z);
    let r = comp_r(z);
    if ix >= 0x3fef_3333 {
        x = PIO2_HI - (2.0 * (s + s * r) - PIO2_LO);
    } else {
        let f = with_hi_lo(hi_word(s), 0);
        let c = (z - f * f) / (s + f);
        x = 0.5 * PIO2_HI - (2.0 * s * r - (PIO2_LO - 2.0 * c) - (0.5 * PIO2_HI - 2.0 * f));
    }
    if (hx >> 31) != 0 { -x } else { x }
}

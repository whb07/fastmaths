#![no_std]

#[cfg(test)]
extern crate std;

pub mod maths;

pub use maths::fastlibm;

#[cfg(test)]
mod tests {
    use super::fastlibm;
    use libloading::Library;
    #[cfg(feature = "mpfr")]
    use rug::{Float, ops::Pow};
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6, PI, TAU};
    use std::path::Path;
    use std::string::String;
    use std::vec;
    use std::vec::Vec;
    use std::{eprintln, format};

    const MAX_ULP_TOL: f64 = 1.0;
    const DERIVED_ULP_TOL: f64 = 1.0;
    const PROPTEST_ULP_TOL: f64 = 1.0;
    #[cfg(feature = "mpfr")]
    const MPFR_PREC: u32 = 256;
    #[cfg(feature = "mpfr")]
    const MPFR_TRIG_LIMIT: f64 = 1.0e6;

    fn ulp_size(x: f64) -> f64 {
        if x == 0.0 {
            return f64::from_bits(1);
        }
        if x.is_nan() || x.is_infinite() {
            return f64::NAN;
        }
        let next = if x.is_sign_negative() {
            x.next_down()
        } else {
            x.next_up()
        };
        (next - x).abs()
    }

    fn ulp_error(actual: f64, expected: f64) -> f64 {
        let diff = (actual - expected).abs();
        if diff == 0.0 {
            return 0.0;
        }
        let ulp = ulp_size(expected);
        if !ulp.is_finite() || ulp == 0.0 {
            return f64::INFINITY;
        }
        diff / ulp
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_exp_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.exp_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_exp2_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.exp2_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_expm1_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.exp_m1_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_ln_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.ln_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_log2_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.log2_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_log10_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.log10_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_log1p_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.ln_1p_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_sin_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.sin_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_cos_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.cos_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_tan_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.tan_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_asin_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.asin_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_acos_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.acos_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_atan_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.atan_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_atan2_f64(y: f64, x: f64) -> f64 {
        let mut vy = Float::with_val(MPFR_PREC, y);
        let vx = Float::with_val(MPFR_PREC, x);
        vy.atan2_mut(&vx);
        vy.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_sqrt_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.sqrt_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_cbrt_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.cbrt_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_hypot_f64(x: f64, y: f64) -> f64 {
        let mut vx = Float::with_val(MPFR_PREC, x);
        let vy = Float::with_val(MPFR_PREC, y);
        vx.hypot_mut(&vy);
        vx.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_pow_f64(x: f64, y: f64) -> f64 {
        let base = Float::with_val(MPFR_PREC, x);
        let exp = Float::with_val(MPFR_PREC, y);
        base.pow(exp).to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_sinh_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.sinh_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_cosh_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.cosh_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_tanh_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.tanh_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_fmod_f64(x: f64, y: f64) -> f64 {
        if x.is_nan() || y.is_nan() || y == 0.0 || x.is_infinite() {
            return f64::NAN;
        }
        let mut vx = Float::with_val(MPFR_PREC, x);
        let vy = Float::with_val(MPFR_PREC, y);
        vx %= vy;
        vx.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_remainder_f64(x: f64, y: f64) -> f64 {
        if x.is_nan() || y.is_nan() || y == 0.0 || x.is_infinite() {
            return f64::NAN;
        }
        let mut vx = Float::with_val(MPFR_PREC, x);
        let vy = Float::with_val(MPFR_PREC, y);
        vx.remainder_mut(&vy);
        vx.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn ln_reference(x: f64) -> f64 {
        mpfr_ln_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn ln_reference(x: f64) -> f64 {
        x.ln()
    }

    #[cfg(feature = "mpfr")]
    fn exp2_reference(x: f64) -> f64 {
        mpfr_exp2_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn exp2_reference(x: f64) -> f64 {
        x.exp2()
    }

    #[cfg(feature = "mpfr")]
    fn expm1_reference(x: f64) -> f64 {
        mpfr_expm1_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn expm1_reference(x: f64) -> f64 {
        x.exp_m1()
    }

    #[cfg(feature = "mpfr")]
    fn log2_reference(x: f64) -> f64 {
        mpfr_log2_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn log2_reference(x: f64) -> f64 {
        x.log2()
    }

    #[cfg(feature = "mpfr")]
    fn log10_reference(x: f64) -> f64 {
        mpfr_log10_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn log10_reference(x: f64) -> f64 {
        x.log10()
    }

    #[cfg(feature = "mpfr")]
    fn log1p_reference(x: f64) -> f64 {
        mpfr_log1p_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn log1p_reference(x: f64) -> f64 {
        x.ln_1p()
    }

    #[cfg(feature = "mpfr")]
    fn sin_reference(x: f64) -> f64 {
        if x.abs() <= MPFR_TRIG_LIMIT {
            mpfr_sin_f64(x)
        } else {
            x.sin()
        }
    }

    #[cfg(not(feature = "mpfr"))]
    fn sin_reference(x: f64) -> f64 {
        x.sin()
    }

    #[cfg(feature = "mpfr")]
    fn cos_reference(x: f64) -> f64 {
        if x.abs() <= MPFR_TRIG_LIMIT {
            mpfr_cos_f64(x)
        } else {
            x.cos()
        }
    }

    #[cfg(not(feature = "mpfr"))]
    fn cos_reference(x: f64) -> f64 {
        x.cos()
    }

    #[cfg(feature = "mpfr")]
    fn tan_reference(x: f64) -> f64 {
        if x.abs() <= MPFR_TRIG_LIMIT {
            mpfr_tan_f64(x)
        } else {
            x.tan()
        }
    }

    #[cfg(not(feature = "mpfr"))]
    fn tan_reference(x: f64) -> f64 {
        x.tan()
    }

    #[cfg(feature = "mpfr")]
    fn asin_reference(x: f64) -> f64 {
        mpfr_asin_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn asin_reference(x: f64) -> f64 {
        x.asin()
    }

    #[cfg(feature = "mpfr")]
    fn acos_reference(x: f64) -> f64 {
        mpfr_acos_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn acos_reference(x: f64) -> f64 {
        x.acos()
    }

    #[cfg(feature = "mpfr")]
    fn atan_reference(x: f64) -> f64 {
        mpfr_atan_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn atan_reference(x: f64) -> f64 {
        x.atan()
    }

    #[cfg(feature = "mpfr")]
    fn atan2_reference(y: f64, x: f64) -> f64 {
        mpfr_atan2_f64(y, x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn atan2_reference(y: f64, x: f64) -> f64 {
        y.atan2(x)
    }

    #[cfg(feature = "mpfr")]
    fn sinh_reference(x: f64) -> f64 {
        mpfr_sinh_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn sinh_reference(x: f64) -> f64 {
        x.sinh()
    }

    #[cfg(feature = "mpfr")]
    fn cosh_reference(x: f64) -> f64 {
        mpfr_cosh_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn cosh_reference(x: f64) -> f64 {
        x.cosh()
    }

    #[cfg(feature = "mpfr")]
    fn tanh_reference(x: f64) -> f64 {
        mpfr_tanh_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn tanh_reference(x: f64) -> f64 {
        x.tanh()
    }

    #[cfg(feature = "mpfr")]
    fn fmod_reference(x: f64, y: f64) -> f64 {
        mpfr_fmod_f64(x, y)
    }

    #[cfg(not(feature = "mpfr"))]
    fn fmod_reference(x: f64, y: f64) -> f64 {
        x % y
    }

    #[cfg(feature = "mpfr")]
    fn remainder_reference(x: f64, y: f64) -> f64 {
        mpfr_remainder_f64(x, y)
    }

    #[cfg(not(feature = "mpfr"))]
    fn remainder_reference(x: f64, y: f64) -> f64 {
        let q = (x / y).round();
        x - q * y
    }

    #[cfg(feature = "mpfr")]
    fn sqrt_reference(x: f64) -> f64 {
        mpfr_sqrt_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn sqrt_reference(x: f64) -> f64 {
        x.sqrt()
    }

    #[cfg(feature = "mpfr")]
    fn cbrt_reference(x: f64) -> f64 {
        mpfr_cbrt_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn cbrt_reference(x: f64) -> f64 {
        x.cbrt()
    }

    #[cfg(feature = "mpfr")]
    fn hypot_reference(x: f64, y: f64) -> f64 {
        mpfr_hypot_f64(x, y)
    }

    #[cfg(not(feature = "mpfr"))]
    fn hypot_reference(x: f64, y: f64) -> f64 {
        x.hypot(y)
    }

    #[cfg(feature = "mpfr")]
    fn pow_reference(x: f64, y: f64) -> f64 {
        mpfr_pow_f64(x, y)
    }

    #[cfg(not(feature = "mpfr"))]
    fn pow_reference(x: f64, y: f64) -> f64 {
        x.powf(y)
    }

    fn assert_ulp_eq(actual: f64, expected: f64, max_ulps: f64, context: &str) {
        if actual.is_nan() && expected.is_nan() {
            return;
        }
        if actual == expected {
            return;
        }
        if actual.is_infinite() || expected.is_infinite() {
            assert_eq!(
                actual, expected,
                "{context}: expected {expected}, got {actual}"
            );
            return;
        }
        let ulps = ulp_error(actual, expected);
        assert!(
            ulps <= max_ulps,
            "{context}: expected {expected}, got {actual} (ulps={ulps})"
        );
    }

    fn assert_ulp_eq_exp(actual: f64, x: f64, max_ulps: f64, context: &str) {
        let expected_std = x.exp();
        if expected_std.is_infinite() || expected_std.is_nan() {
            assert_ulp_eq(actual, expected_std, max_ulps, context);
            return;
        }
        let ulps_std = ulp_error(actual, expected_std);
        if ulps_std <= max_ulps {
            return;
        }

        #[cfg(feature = "mpfr")]
        {
            let expected_mpfr = mpfr_exp_f64(x);
            let ulps_mpfr = ulp_error(actual, expected_mpfr);
            if ulps_mpfr <= max_ulps {
                return;
            }
            panic!(
                "{context}: expected {expected_std} (std) / {expected_mpfr} (mpfr), got {actual} (ulps_std={ulps_std}, ulps_mpfr={ulps_mpfr})"
            );
        }

        #[cfg(not(feature = "mpfr"))]
        {
            panic!("{context}: expected {expected_std}, got {actual} (ulps={ulps_std})");
        }
    }

    fn push_unique(values: &mut Vec<f64>, x: f64) {
        if !values.iter().any(|v| v.to_bits() == x.to_bits()) {
            values.push(x);
        }
    }

    fn exp_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            f64::from_bits(1),
            -f64::from_bits(1),
            f64::MIN_POSITIVE,
            -f64::MIN_POSITIVE,
            1e-300,
            -1e-300,
            1e-200,
            -1e-200,
            1e-100,
            -1e-100,
            std::f64::consts::LN_2,
            -std::f64::consts::LN_2,
            std::f64::consts::LN_2 / 128.0,
            -std::f64::consts::LN_2 / 128.0,
            std::f64::consts::LN_2 / 256.0,
            -std::f64::consts::LN_2 / 256.0,
            -745.133_219_101_941_1,
            -744.0,
            -720.0,
            -709.78,
            -100.0,
            -20.0,
            -10.0,
            -1.0,
            -0.5,
            -1e-12,
            -1e-6,
            -1e-16,
            0.0,
            1e-16,
            1e-12,
            1e-6,
            0.5,
            1.0,
            2.0,
            10.0,
            20.0,
            100.0,
            700.0,
            709.0,
            709.5,
            709.7,
            709.782_712_893_384,
            709.8,
            710.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for &x in &[
            0.0f64.next_up(),
            0.0f64.next_down(),
            1.0f64.next_up(),
            1.0f64.next_down(),
            (-1.0f64).next_up(),
            (-1.0f64).next_down(),
        ] {
            push_unique(&mut inputs, x);
        }
        for i in -200..=200 {
            push_unique(&mut inputs, (i as f64) * 0.25);
        }
        for i in -100..=100 {
            push_unique(&mut inputs, (i as f64) * 7.5);
        }
        for i in -70..=70 {
            push_unique(&mut inputs, (i as f64) * 10.0);
        }
        inputs
    }

    fn exp_special_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            0.0,
            -0.0,
            1.0,
            -1.0,
            0.5,
            -0.5,
            std::f64::consts::LN_2,
            -std::f64::consts::LN_2,
            std::f64::consts::LN_2 / 128.0,
            -std::f64::consts::LN_2 / 128.0,
            -100.0,
            100.0,
            -700.0,
            700.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        inputs
    }

    fn ln_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let min_sub = f64::from_bits(1);
        let max_sub = f64::from_bits(0x000f_ffff_ffff_ffff);
        let specials = [
            min_sub,
            max_sub,
            f64::MIN_POSITIVE,
            f64::from_bits(0x03fe_ffff_ffff_ffff),
            f64::from_bits(0x3ff0000000000001),
            1e-308,
            1e-300,
            1e-200,
            1e-100,
            1e-50,
            1e-20,
            1e-10,
            1e-5,
            0.1,
            0.5,
            0.9,
            0.999_999_999_999,
            1.0,
            1.000_000_000_001,
            1.5,
            2.0,
            10.0,
            1e5,
            1e10,
            1e100,
            f64::MAX,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for &x in &[
            1.0f64.next_up(),
            1.0f64.next_down(),
            f64::MIN_POSITIVE.next_up(),
            f64::MIN_POSITIVE.next_down(),
        ] {
            if x > 0.0 {
                push_unique(&mut inputs, x);
            }
        }
        for i in -60..=60 {
            let x = 2f64.powi(i);
            push_unique(&mut inputs, x);
        }
        for i in -100..=100 {
            let x = 1.0 + (i as f64) * 1e-6;
            if x > 0.0 {
                push_unique(&mut inputs, x);
            }
        }
        inputs
    }

    fn trig_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            0.0,
            1e-12,
            -1e-12,
            1e-6,
            -1e-6,
            0.5,
            -0.5,
            1.0,
            -1.0,
            PI / 7.0,
            -PI / 7.0,
            FRAC_PI_2,
            FRAC_PI_2 + 1e-15,
            FRAC_PI_2 - 1e-15,
            PI,
            PI + 1e-15,
            PI - 1e-15,
            3.0 * FRAC_PI_2,
            TAU,
            10.0,
            -10.0,
            1e6,
            -1e6,
            1e12,
            -1e12,
            1e20,
            -1e20,
            1e100,
            -1e100,
            1e300,
            -1e300,
            (1u64 << 53) as f64,
            (1u64 << 62) as f64,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for &base in &[0.0, FRAC_PI_2, PI, 3.0 * PI / 2.0, TAU] {
            push_unique(&mut inputs, base);
            push_unique(&mut inputs, base.next_up());
            push_unique(&mut inputs, base.next_down());
            push_unique(&mut inputs, -base);
        }
        for i in -200..=200 {
            push_unique(&mut inputs, (i as f64) * 0.25);
        }
        for i in -100..=100 {
            push_unique(&mut inputs, (i as f64) * 2.5);
        }
        for i in 1..=64 {
            let x = (i as f64) * PI / 32.0;
            push_unique(&mut inputs, x);
            push_unique(&mut inputs, -x);
        }
        inputs
    }

    fn exp2_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            -1074.0,
            -1022.0,
            -100.0,
            -10.0,
            -1.0,
            -0.5,
            -1e-6,
            0.0,
            1e-6,
            0.5,
            1.0,
            2.0,
            10.0,
            100.0,
            1023.0,
            1024.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -16..=16 {
            push_unique(&mut inputs, (i as f64) * 0.125);
        }
        inputs
    }

    fn expm1_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            -1e-12,
            -1e-6,
            -1e-3,
            -0.1,
            -1.0,
            -10.0,
            0.0,
            1e-12,
            1e-6,
            1e-3,
            0.1,
            1.0,
            10.0,
            50.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        inputs
    }

    fn tan_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            0.0,
            -0.0,
            1e-12,
            -1e-12,
            FRAC_PI_6,
            FRAC_PI_4,
            PI / 3.0,
            FRAC_PI_2 - 1e-12,
            FRAC_PI_2 + 1e-12,
            -FRAC_PI_2 + 1e-12,
            -FRAC_PI_2 - 1e-12,
            PI,
            2.0 * PI,
            10.0,
            -10.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -64..=64 {
            push_unique(&mut inputs, (i as f64) * PI / 32.0);
        }
        inputs
    }

    fn atan_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            0.0,
            -0.0,
            1e-12,
            -1e-12,
            0.5,
            -0.5,
            1.0,
            -1.0,
            10.0,
            -10.0,
            1e6,
            -1e6,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        inputs
    }

    fn atan2_inputs() -> Vec<(f64, f64)> {
        vec![
            (0.0, 0.0),
            (-0.0, 0.0),
            (0.0, -0.0),
            (-0.0, -0.0),
            (1.0, 0.0),
            (-1.0, 0.0),
            (0.0, 1.0),
            (0.0, -1.0),
            (1.0, 1.0),
            (1.0, -1.0),
            (-1.0, 1.0),
            (-1.0, -1.0),
            (1e-12, 1.0),
            (-1e-12, 1.0),
            (1.0, 1e-12),
            (1.0, -1e-12),
            (1e6, 1.0),
            (-1e6, 1.0),
            (1.0, 1e6),
            (1.0, -1e6),
            (f64::INFINITY, 1.0),
            (-f64::INFINITY, 1.0),
            (1.0, f64::INFINITY),
            (1.0, f64::NEG_INFINITY),
        ]
    }

    fn hypot_inputs() -> Vec<(f64, f64)> {
        vec![
            (0.0, 0.0),
            (3.0, 4.0),
            (1e-300, 1e-300),
            (1e-200, -1e-200),
            (1e-50, 1e-60),
            (1e100, -1e100),
            (1e300, 1e300),
            (f64::INFINITY, 1.0),
            (f64::NAN, 1.0),
        ]
    }

    fn pow_inputs() -> Vec<(f64, f64)> {
        vec![
            (2.0, 3.0),
            (2.0, -3.0),
            (10.0, 0.5),
            (0.5, 2.0),
            (-2.0, 3.0),
            (-2.0, 4.0),
            (-2.0, 0.5),
            (0.0, 2.0),
            (0.0, -2.0),
            (1e-300, 2.0),
            (1e300, 2.0),
            (-1.0, 1e6),
        ]
    }

    fn sqrt_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            f64::NAN,
            f64::INFINITY,
            0.0,
            -0.0,
            1.0,
            2.0,
            4.0,
            1e-300,
            1e300,
            -1.0,
            -0.5,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        inputs
    }

    fn cbrt_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            0.0,
            -0.0,
            1.0,
            -1.0,
            8.0,
            -8.0,
            1e-300,
            -1e-300,
            1e300,
            -1e300,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        inputs
    }

    fn log1p_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            -1.0,
            -0.999_999_999_999,
            -0.9,
            -0.5,
            -0.25,
            -1e-12,
            -1e-6,
            -0.0,
            0.0,
            1e-12,
            1e-6,
            0.25,
            0.5,
            1.0,
            2.0,
            10.0,
            1e6,
            1e12,
            1e20,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for &x in &[
            (-1.0f64).next_up(),
            (-1.0f64).next_down(),
            0.0f64.next_up(),
            0.0f64.next_down(),
        ] {
            push_unique(&mut inputs, x);
        }
        for i in -100..=100 {
            push_unique(&mut inputs, (i as f64) * 1e-3);
        }
        inputs
    }

    fn asin_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            -1.0,
            -0.999_999_999_999,
            -0.75,
            -0.5,
            -0.25,
            -1e-12,
            -1e-6,
            -0.0,
            0.0,
            1e-12,
            1e-6,
            0.25,
            0.5,
            0.75,
            0.999_999_999_999,
            1.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -100..=100 {
            let x = (i as f64) / 100.0;
            if (-1.0..=1.0).contains(&x) {
                push_unique(&mut inputs, x);
            }
        }
        inputs
    }

    fn sinh_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            -700.0, -100.0, -50.0, -20.0, -10.0, -1.0, -0.5, -1e-6, -1e-12, -0.0, 0.0, 1e-12, 1e-6,
            0.5, 1.0, 10.0, 20.0, 50.0, 100.0, 700.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -100..=100 {
            push_unique(&mut inputs, (i as f64) * 0.25);
        }
        inputs
    }

    fn cosh_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            0.0, 1e-12, 1e-6, 0.5, 1.0, 2.0, 10.0, 20.0, 50.0, 100.0, 700.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
            push_unique(&mut inputs, -x);
        }
        for i in -100..=100 {
            push_unique(&mut inputs, (i as f64) * 0.25);
        }
        inputs
    }

    fn tanh_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            -50.0, -20.0, -10.0, -3.0, -1.0, -0.5, -1e-6, -1e-12, -0.0, 0.0, 1e-12, 1e-6, 0.5, 1.0,
            3.0, 10.0, 20.0, 50.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -100..=100 {
            push_unique(&mut inputs, (i as f64) * 0.1);
        }
        inputs
    }

    fn fmod_inputs() -> Vec<(f64, f64)> {
        vec![
            (0.0, 1.0),
            (-0.0, 1.0),
            (1.0, 0.5),
            (1.5, 0.5),
            (-1.5, 0.5),
            (5.3, 2.1),
            (-5.3, 2.1),
            (1e20, 3.0),
            (1e-10, 1e-12),
            (1e-10, 1e-6),
            (1.0, 1.0),
            (2.0, 3.0),
            (-2.0, 3.0),
        ]
    }

    fn remainder_inputs() -> Vec<(f64, f64)> {
        vec![
            (0.0, 1.0),
            (-0.0, 1.0),
            (1.0, 0.5),
            (1.5, 0.5),
            (-1.5, 0.5),
            (5.3, 2.1),
            (-5.3, 2.1),
            (1e20, 3.0),
            (1e-10, 1e-12),
            (1.0, 1.0),
            (2.0, 3.0),
            (-2.0, 3.0),
        ]
    }

    fn glibc_libm_path() -> Option<String> {
        if std::env::var("FASTLIBM_GLIBC_TEST").is_err() {
            return None;
        }
        let path = std::env::var("FASTLIBM_GLIBC_LIBM")
            .unwrap_or_else(|_| String::from("/tmp/maths/glibc-build/math/libm.so"));
        if !Path::new(&path).exists() {
            eprintln!("glibc libm not found at {path}");
            return None;
        }
        Some(path)
    }

    fn glibc_libm_path_dist() -> Option<String> {
        if std::env::var("FASTLIBM_GLIBC_DIST").is_err() {
            return None;
        }
        let path = std::env::var("FASTLIBM_GLIBC_LIBM")
            .unwrap_or_else(|_| String::from("/tmp/maths/glibc-build/math/libm.so"));
        if !Path::new(&path).exists() {
            eprintln!("glibc libm not found at {path}");
            return None;
        }
        Some(path)
    }

    fn assert_ulp_eq_glibc(actual: f64, expected: f64, max_ulps: f64, context: &str) {
        if actual == 0.0 && expected == 0.0 {
            assert_eq!(
                actual.to_bits(),
                expected.to_bits(),
                "{context}: zero sign mismatch"
            );
            return;
        }
        assert_ulp_eq(actual, expected, max_ulps, context);
    }

    fn rand_u64(state: &mut u64) -> u64 {
        const A: u64 = 6364136223846793005;
        const C: u64 = 1442695040888963407;
        *state = state.wrapping_mul(A).wrapping_add(C);
        *state
    }

    fn rand_f64_unit(state: &mut u64) -> f64 {
        let bits = rand_u64(state) >> 11;
        (bits as f64) / ((1u64 << 53) as f64)
    }

    fn rand_range(state: &mut u64, min: f64, max: f64) -> f64 {
        min + (max - min) * rand_f64_unit(state)
    }

    fn rand_f64_pos(state: &mut u64) -> f64 {
        let exp = (rand_u64(state) % 0x7fe) + 1;
        let mant = rand_u64(state) & 0x000f_ffff_ffff_ffff;
        f64::from_bits((exp << 52) | mant)
    }

    #[test]
    fn exp_special_cases() {
        let nan = f64::NAN;
        let pos_inf = f64::INFINITY;
        let neg_inf = f64::NEG_INFINITY;

        assert!(fastlibm::exp(nan).is_nan());
        assert_eq!(fastlibm::exp(pos_inf), f64::INFINITY);
        assert_eq!(fastlibm::exp(neg_inf), 0.0);
        assert_eq!(fastlibm::exp(0.0).to_bits(), 1.0f64.to_bits());
        assert_eq!(fastlibm::exp(-0.0).to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn exp_matches_std_ulps() {
        let inputs = exp_inputs();

        for &x in &inputs {
            let actual = fastlibm::exp(x);
            let context = format!("exp({x})");
            assert_ulp_eq_exp(actual, x, MAX_ULP_TOL, &context);
        }
    }

    #[test]
    fn ln_special_cases() {
        let nan = f64::NAN;
        let pos_inf = f64::INFINITY;

        assert!(fastlibm::ln(nan).is_nan());
        assert_eq!(fastlibm::ln(pos_inf), f64::INFINITY);
        assert_eq!(fastlibm::ln(0.0), f64::NEG_INFINITY);
        assert_eq!(fastlibm::ln(-0.0), (-0.0f64).ln());
        assert!(fastlibm::ln(-1.0).is_nan());
    }

    #[test]
    fn ln_matches_std_ulps() {
        let inputs = ln_inputs();

        for &x in &inputs {
            let expected = x.ln();
            let actual = fastlibm::ln(x);
            let context = format!("ln({x})");
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &context);
        }
    }

    #[test]
    fn log1p_special_cases() {
        assert!(fastlibm::log1p(f64::NAN).is_nan());
        assert_eq!(fastlibm::log1p(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastlibm::log1p(0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(fastlibm::log1p(-0.0).to_bits(), (-0.0f64).to_bits());
        assert_eq!(fastlibm::log1p(-1.0), f64::NEG_INFINITY);
        assert!(fastlibm::log1p(-1.5).is_nan());
    }

    #[test]
    fn log1p_matches_reference_ulps() {
        for &x in &log1p_inputs() {
            let expected = log1p_reference(x);
            let actual = fastlibm::log1p(x);
            if expected.is_nan() {
                assert!(actual.is_nan(), "log1p({x}) expected NaN, got {actual}");
            } else if expected.is_infinite() {
                assert_eq!(actual, expected, "log1p({x}) expected {expected}");
            } else {
                assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("log1p({x})"));
            }
        }
    }

    #[test]
    fn exp2_special_cases() {
        assert!(fastlibm::exp2(f64::NAN).is_nan());
        assert_eq!(fastlibm::exp2(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastlibm::exp2(f64::NEG_INFINITY), 0.0);
        assert_eq!(fastlibm::exp2(0.0).to_bits(), 1.0f64.to_bits());
        assert_eq!(fastlibm::exp2(-0.0).to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn exp2_matches_reference_ulps() {
        for &x in &exp2_inputs() {
            let actual = fastlibm::exp2(x);
            let expected = exp2_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("exp2({x})"));
        }
    }

    #[test]
    fn expm1_special_cases() {
        assert!(fastlibm::expm1(f64::NAN).is_nan());
        assert_eq!(fastlibm::expm1(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastlibm::expm1(f64::NEG_INFINITY), -1.0);
        assert_eq!(fastlibm::expm1(0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(fastlibm::expm1(-0.0).to_bits(), (-0.0f64).to_bits());
    }

    #[test]
    fn expm1_matches_reference_ulps() {
        for &x in &expm1_inputs() {
            let actual = fastlibm::expm1(x);
            let expected = expm1_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("expm1({x})"));
        }
    }

    #[test]
    fn log2_log10_special_cases() {
        assert!(fastlibm::log2(f64::NAN).is_nan());
        assert!(fastlibm::log10(f64::NAN).is_nan());
        assert_eq!(fastlibm::log2(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastlibm::log10(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastlibm::log2(0.0), f64::NEG_INFINITY);
        assert_eq!(fastlibm::log10(0.0), f64::NEG_INFINITY);
        assert!(fastlibm::log2(-1.0).is_nan());
        assert!(fastlibm::log10(-1.0).is_nan());
    }

    #[test]
    fn log2_log10_matches_reference_ulps() {
        for &x in &ln_inputs() {
            let actual = fastlibm::log2(x);
            let expected = log2_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("log2({x})"));

            let actual = fastlibm::log10(x);
            let expected = log10_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("log10({x})"));
        }
    }

    #[test]
    fn tan_special_cases() {
        assert!(fastlibm::tan(f64::NAN).is_nan());
        assert!(fastlibm::tan(f64::INFINITY).is_nan());
        assert!(fastlibm::tan(f64::NEG_INFINITY).is_nan());
    }

    #[test]
    fn tan_matches_reference_ulps() {
        for &x in &tan_inputs() {
            let actual = fastlibm::tan(x);
            let expected = tan_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("tan({x})"));
        }
    }

    #[test]
    fn asin_acos_special_cases() {
        assert!(fastlibm::asin(f64::NAN).is_nan());
        assert!(fastlibm::acos(f64::NAN).is_nan());
        assert_eq!(fastlibm::asin(1.0), FRAC_PI_2);
        assert_eq!(fastlibm::asin(-1.0), -FRAC_PI_2);
        assert_eq!(fastlibm::acos(1.0), 0.0);
        assert_eq!(fastlibm::acos(-1.0), PI);
        assert!(fastlibm::asin(1.1).is_nan());
        assert!(fastlibm::acos(-1.1).is_nan());
    }

    #[test]
    fn asin_acos_matches_reference_ulps() {
        for &x in &asin_inputs() {
            let actual = fastlibm::asin(x);
            let expected = asin_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("asin({x})"));

            let actual = fastlibm::acos(x);
            let expected = acos_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("acos({x})"));
        }
    }

    #[test]
    fn atan_matches_reference_ulps() {
        for &x in &atan_inputs() {
            let actual = fastlibm::atan(x);
            let expected = atan_reference(x);
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &format!("atan({x})"));
        }
    }

    #[test]
    fn atan2_matches_reference_ulps() {
        for &(y, x) in &atan2_inputs() {
            let actual = fastlibm::atan2(y, x);
            let expected = atan2_reference(y, x);
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &format!("atan2({y},{x})"));
        }
    }

    #[test]
    fn sinh_cosh_tanh_special_cases() {
        assert!(fastlibm::sinh(f64::NAN).is_nan());
        assert!(fastlibm::cosh(f64::NAN).is_nan());
        assert!(fastlibm::tanh(f64::NAN).is_nan());

        assert_eq!(fastlibm::sinh(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastlibm::sinh(f64::NEG_INFINITY), f64::NEG_INFINITY);
        assert_eq!(fastlibm::cosh(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastlibm::cosh(f64::NEG_INFINITY), f64::INFINITY);
        assert_eq!(fastlibm::tanh(f64::INFINITY), 1.0);
        assert_eq!(fastlibm::tanh(f64::NEG_INFINITY), -1.0);
    }

    #[test]
    fn sinh_cosh_tanh_matches_reference_ulps() {
        for &x in &sinh_inputs() {
            let actual = fastlibm::sinh(x);
            let expected = sinh_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("sinh({x})"));
        }
        for &x in &cosh_inputs() {
            let actual = fastlibm::cosh(x);
            let expected = cosh_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("cosh({x})"));
        }
        for &x in &tanh_inputs() {
            let actual = fastlibm::tanh(x);
            let expected = tanh_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("tanh({x})"));
        }
    }

    #[test]
    fn hypot_matches_reference_ulps() {
        for &(x, y) in &hypot_inputs() {
            let actual = fastlibm::hypot(x, y);
            let expected = hypot_reference(x, y);
            assert_ulp_eq(
                actual,
                expected,
                DERIVED_ULP_TOL,
                &format!("hypot({x},{y})"),
            );
        }
    }

    #[test]
    fn fmod_special_cases() {
        assert!(fastlibm::fmod(f64::NAN, 1.0).is_nan());
        assert!(fastlibm::fmod(1.0, f64::NAN).is_nan());
        assert!(fastlibm::fmod(f64::INFINITY, 1.0).is_nan());
        assert!(fastlibm::fmod(1.0, 0.0).is_nan());
        assert_eq!(fastlibm::fmod(0.0, 1.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(fastlibm::fmod(-0.0, 1.0).to_bits(), (-0.0f64).to_bits());
    }

    #[test]
    fn fmod_matches_reference_ulps() {
        for &(x, y) in &fmod_inputs() {
            let actual = fastlibm::fmod(x, y);
            let expected = fmod_reference(x, y);
            if expected.is_nan() {
                assert!(actual.is_nan(), "fmod({x},{y}) expected NaN");
            } else {
                assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("fmod({x},{y})"));
            }
        }
    }

    #[test]
    fn remainder_special_cases() {
        assert!(fastlibm::remainder(f64::NAN, 1.0).is_nan());
        assert!(fastlibm::remainder(1.0, f64::NAN).is_nan());
        assert!(fastlibm::remainder(f64::INFINITY, 1.0).is_nan());
        assert!(fastlibm::remainder(1.0, 0.0).is_nan());
        assert_eq!(fastlibm::remainder(0.0, 1.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(
            fastlibm::remainder(-0.0, 1.0).to_bits(),
            (-0.0f64).to_bits()
        );
        assert_eq!(fastlibm::remainder(1.0, f64::INFINITY), 1.0);
    }

    #[test]
    fn remainder_matches_reference_ulps() {
        for &(x, y) in &remainder_inputs() {
            let actual = fastlibm::remainder(x, y);
            let expected = remainder_reference(x, y);
            if expected.is_nan() {
                assert!(actual.is_nan(), "remainder({x},{y}) expected NaN");
            } else {
                assert_ulp_eq(
                    actual,
                    expected,
                    DERIVED_ULP_TOL,
                    &format!("remainder({x},{y})"),
                );
            }
        }
    }

    #[test]
    fn pow_matches_reference_ulps() {
        for &(x, y) in &pow_inputs() {
            let actual = fastlibm::pow(x, y);
            let expected = pow_reference(x, y);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("pow({x},{y})"));
        }
    }

    #[test]
    fn sqrt_matches_reference_ulps() {
        for &x in &sqrt_inputs() {
            let actual = fastlibm::sqrt(x);
            let expected = sqrt_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("sqrt({x})"));
        }
    }

    #[test]
    fn cbrt_matches_reference_ulps() {
        for &x in &cbrt_inputs() {
            let actual = fastlibm::cbrt(x);
            let expected = cbrt_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("cbrt({x})"));
        }
    }

    #[test]
    fn sin_cos_special_cases() {
        let nan = f64::NAN;
        let pos_inf = f64::INFINITY;
        let neg_inf = f64::NEG_INFINITY;

        assert!(fastlibm::sin(nan).is_nan());
        assert!(fastlibm::cos(nan).is_nan());
        assert!(fastlibm::sin(pos_inf).is_nan());
        assert!(fastlibm::cos(pos_inf).is_nan());
        assert!(fastlibm::sin(neg_inf).is_nan());
        assert!(fastlibm::cos(neg_inf).is_nan());

        let neg_zero = -0.0f64;
        assert_eq!(fastlibm::sin(neg_zero).to_bits(), neg_zero.to_bits());
        assert_eq!(fastlibm::cos(neg_zero).to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn sin_cos_known_angles() {
        let inputs = [
            0.0,
            FRAC_PI_6,
            FRAC_PI_4,
            PI / 3.0,
            FRAC_PI_2,
            PI,
            2.0 * PI,
            TAU,
            -FRAC_PI_2,
            -PI,
        ];

        for &x in &inputs {
            let sin_expected = x.sin();
            let cos_expected = x.cos();
            let sin_actual = fastlibm::sin(x);
            let cos_actual = fastlibm::cos(x);
            assert_ulp_eq(sin_actual, sin_expected, MAX_ULP_TOL, &format!("sin({x})"));
            assert_ulp_eq(cos_actual, cos_expected, MAX_ULP_TOL, &format!("cos({x})"));
        }
    }

    #[test]
    fn sin_cos_matches_std_ulps() {
        let inputs = trig_inputs();

        for &x in &inputs {
            let sin_expected = x.sin();
            let cos_expected = x.cos();
            let sin_actual = fastlibm::sin(x);
            let cos_actual = fastlibm::cos(x);
            assert_ulp_eq(sin_actual, sin_expected, MAX_ULP_TOL, &format!("sin({x})"));
            assert_ulp_eq(cos_actual, cos_expected, MAX_ULP_TOL, &format!("cos({x})"));
        }
    }

    #[test]
    fn sincos_matches_std_ulps() {
        let inputs = trig_inputs();

        for &x in &inputs {
            let (sin_actual, cos_actual) = fastlibm::sincos(x);
            let sin_expected = x.sin();
            let cos_expected = x.cos();
            assert_ulp_eq(
                sin_actual,
                sin_expected,
                MAX_ULP_TOL,
                &format!("sincos sin({x})"),
            );
            assert_ulp_eq(
                cos_actual,
                cos_expected,
                MAX_ULP_TOL,
                &format!("sincos cos({x})"),
            );
        }
    }

    #[test]
    fn sin_cos_symmetry() {
        let inputs = [
            -10.0, -3.0, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 3.0, 10.0, 1e6, 1e12, 1e20,
        ];

        for &x in &inputs {
            let sin_pos = fastlibm::sin(x);
            let sin_neg = fastlibm::sin(-x);
            let cos_pos = fastlibm::cos(x);
            let cos_neg = fastlibm::cos(-x);

            assert_ulp_eq(
                sin_neg,
                -sin_pos,
                MAX_ULP_TOL,
                &format!("sin symmetry at {x}"),
            );
            assert_ulp_eq(
                cos_neg,
                cos_pos,
                MAX_ULP_TOL,
                &format!("cos symmetry at {x}"),
            );
        }
    }

    #[test]
    fn exp_matches_glibc_ulps() {
        let Some(path) = glibc_libm_path() else {
            return;
        };
        let lib = unsafe { Library::new(&path).expect("load glibc libm") };
        let exp: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            unsafe { lib.get(b"exp").expect("load exp") };

        for &x in &exp_inputs() {
            let expected = unsafe { exp(x) };
            let actual = fastlibm::exp(x);
            let context = format!("glibc exp({x})");
            assert_ulp_eq_glibc(actual, expected, MAX_ULP_TOL, &context);
        }
    }

    #[test]
    fn ln_matches_glibc_ulps() {
        let Some(path) = glibc_libm_path() else {
            return;
        };
        let lib = unsafe { Library::new(&path).expect("load glibc libm") };
        let log: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            unsafe { lib.get(b"log").expect("load log") };

        for &x in &ln_inputs() {
            let expected = unsafe { log(x) };
            let actual = fastlibm::ln(x);
            let context = format!("glibc log({x})");
            assert_ulp_eq_glibc(actual, expected, MAX_ULP_TOL, &context);
        }
    }

    #[test]
    fn sin_cos_match_glibc_ulps() {
        let Some(path) = glibc_libm_path() else {
            return;
        };
        let lib = unsafe { Library::new(&path).expect("load glibc libm") };
        let sin: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            unsafe { lib.get(b"sin").expect("load sin") };
        let cos: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            unsafe { lib.get(b"cos").expect("load cos") };

        for &x in &trig_inputs() {
            let sin_expected = unsafe { sin(x) };
            let cos_expected = unsafe { cos(x) };
            let sin_actual = fastlibm::sin(x);
            let cos_actual = fastlibm::cos(x);
            assert_ulp_eq_glibc(
                sin_actual,
                sin_expected,
                MAX_ULP_TOL,
                &format!("glibc sin({x})"),
            );
            assert_ulp_eq_glibc(
                cos_actual,
                cos_expected,
                MAX_ULP_TOL,
                &format!("glibc cos({x})"),
            );
        }
    }

    #[test]
    fn compare_glibc_fastlibm() {
        let Some(path) = glibc_libm_path() else {
            return;
        };
        let lib = unsafe { Library::new(&path).expect("load glibc libm") };
        unsafe {
            let g_exp: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
                lib.get(b"exp").unwrap();
            let g_log: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
                lib.get(b"log").unwrap();
            let g_sin: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
                lib.get(b"sin").unwrap();
            let g_cos: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
                lib.get(b"cos").unwrap();

            let test_inputs = [1.0, 2.0, PI, 1e10, -6.5684415251369026e19, 0.0, -0.0];

            std::println!("| Input | Func | glibc (bits) | fastlibm (bits) | ULP Delta |");
            std::println!("| :--- | :--- | :--- | :--- | :--- |");
            for &x in &test_inputs {
                type CFn = unsafe extern "C" fn(f64) -> f64;
                type RustFn = fn(f64) -> f64;
                type FnSpec = (&'static str, CFn, RustFn);

                let fns: [FnSpec; 4] = [
                    ("exp", *g_exp, fastlibm::exp),
                    ("log", *g_log, fastlibm::ln),
                    ("sin", *g_sin, fastlibm::sin),
                    ("cos", *g_cos, fastlibm::cos),
                ];
                for (name, gf, ff) in fns {
                    let gv = gf(x);
                    let fv = ff(x);
                    let delta = ulp_error(fv, gv);
                    std::println!(
                        "| {:e} | {} | {:016x} | {:016x} | {:.4} |",
                        x,
                        name,
                        gv.to_bits(),
                        fv.to_bits(),
                        delta
                    );
                }
            }
        }
    }

    #[test]
    fn glibc_distribution_ulps() {
        let Some(path) = glibc_libm_path_dist() else {
            return;
        };
        let lib = unsafe { Library::new(&path).expect("load glibc libm") };

        type CFn = unsafe extern "C" fn(f64) -> f64;
        type CFn2 = unsafe extern "C" fn(f64, f64) -> f64;

        let exp: libloading::Symbol<CFn> = unsafe { lib.get(b"exp").unwrap() };
        let exp2: libloading::Symbol<CFn> = unsafe { lib.get(b"exp2").unwrap() };
        let expm1: libloading::Symbol<CFn> = unsafe { lib.get(b"expm1").unwrap() };
        let log: libloading::Symbol<CFn> = unsafe { lib.get(b"log").unwrap() };
        let log2: libloading::Symbol<CFn> = unsafe { lib.get(b"log2").unwrap() };
        let log10: libloading::Symbol<CFn> = unsafe { lib.get(b"log10").unwrap() };
        let log1p: libloading::Symbol<CFn> = unsafe { lib.get(b"log1p").unwrap() };
        let sin: libloading::Symbol<CFn> = unsafe { lib.get(b"sin").unwrap() };
        let cos: libloading::Symbol<CFn> = unsafe { lib.get(b"cos").unwrap() };
        let tan: libloading::Symbol<CFn> = unsafe { lib.get(b"tan").unwrap() };
        let asin: libloading::Symbol<CFn> = unsafe { lib.get(b"asin").unwrap() };
        let acos: libloading::Symbol<CFn> = unsafe { lib.get(b"acos").unwrap() };
        let atan: libloading::Symbol<CFn> = unsafe { lib.get(b"atan").unwrap() };
        let atan2: libloading::Symbol<CFn2> = unsafe { lib.get(b"atan2").unwrap() };
        let sinh: libloading::Symbol<CFn> = unsafe { lib.get(b"sinh").unwrap() };
        let cosh: libloading::Symbol<CFn> = unsafe { lib.get(b"cosh").unwrap() };
        let tanh: libloading::Symbol<CFn> = unsafe { lib.get(b"tanh").unwrap() };
        let hypot: libloading::Symbol<CFn2> = unsafe { lib.get(b"hypot").unwrap() };
        let fmod: libloading::Symbol<CFn2> = unsafe { lib.get(b"fmod").unwrap() };
        let remainder: libloading::Symbol<CFn2> = unsafe { lib.get(b"remainder").unwrap() };
        let pow: libloading::Symbol<CFn2> = unsafe { lib.get(b"pow").unwrap() };
        let sqrt: libloading::Symbol<CFn> = unsafe { lib.get(b"sqrt").unwrap() };
        let cbrt: libloading::Symbol<CFn> = unsafe { lib.get(b"cbrt").unwrap() };

        let mut state = 0x1234_5678_9abc_def0u64;
        let samples = 256usize;

        for _ in 0..samples {
            let x = rand_range(&mut state, -100.0, 100.0);
            let g = unsafe { exp(x) };
            let f = fastlibm::exp(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist exp({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -100.0, 100.0);
            let g = unsafe { exp2(x) };
            let f = fastlibm::exp2(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist exp2({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1.0, 1.0);
            let g = unsafe { expm1(x) };
            let f = fastlibm::expm1(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist expm1({x})"));
        }

        for _ in 0..samples {
            let x = rand_f64_pos(&mut state);
            let g = unsafe { log(x) };
            let f = fastlibm::ln(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist ln({x})"));
        }

        for _ in 0..samples {
            let x = rand_f64_pos(&mut state);
            let g = unsafe { log2(x) };
            let f = fastlibm::log2(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist log2({x})"));
        }

        for _ in 0..samples {
            let x = rand_f64_pos(&mut state);
            let g = unsafe { log10(x) };
            let f = fastlibm::log10(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist log10({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -0.9, 1e6);
            let g = unsafe { log1p(x) };
            let f = fastlibm::log1p(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist log1p({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let g = unsafe { sin(x) };
            let f = fastlibm::sin(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist sin({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let g = unsafe { cos(x) };
            let f = fastlibm::cos(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist cos({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let g = unsafe { tan(x) };
            let f = fastlibm::tan(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist tan({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1.0, 1.0);
            let g = unsafe { asin(x) };
            let f = fastlibm::asin(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist asin({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1.0, 1.0);
            let g = unsafe { acos(x) };
            let f = fastlibm::acos(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist acos({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let g = unsafe { atan(x) };
            let f = fastlibm::atan(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist atan({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            let g = unsafe { sinh(x) };
            let f = fastlibm::sinh(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist sinh({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            let g = unsafe { cosh(x) };
            let f = fastlibm::cosh(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist cosh({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            let g = unsafe { tanh(x) };
            let f = fastlibm::tanh(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist tanh({x})"));
        }

        for _ in 0..samples {
            let y = rand_range(&mut state, -1e6, 1e6);
            let x = rand_range(&mut state, -1e6, 1e6);
            if x == 0.0 && y == 0.0 {
                continue;
            }
            let g = unsafe { atan2(y, x) };
            let f = fastlibm::atan2(y, x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist atan2({y},{x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e200, 1e200);
            let y = rand_range(&mut state, -1e200, 1e200);
            let g = unsafe { hypot(x, y) };
            let f = fastlibm::hypot(x, y);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist hypot({x},{y})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let mut y = rand_range(&mut state, 1e-6, 1e6);
            if rand_u64(&mut state) & 1 == 0 {
                y = -y;
            }
            let g = unsafe { fmod(x, y) };
            let f = fastlibm::fmod(x, y);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist fmod({x},{y})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let mut y = rand_range(&mut state, 1e-6, 1e6);
            if rand_u64(&mut state) & 1 == 0 {
                y = -y;
            }
            let g = unsafe { remainder(x, y) };
            let f = fastlibm::remainder(x, y);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist remainder({x},{y})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, 0.1, 10.0);
            let y = rand_range(&mut state, -10.0, 10.0);
            let g = unsafe { pow(x, y) };
            let f = fastlibm::pow(x, y);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist pow({x},{y})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, 0.0, 1e300);
            let g = unsafe { sqrt(x) };
            let f = fastlibm::sqrt(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist sqrt({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e300, 1e300);
            let g = unsafe { cbrt(x) };
            let f = fastlibm::cbrt(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist cbrt({x})"));
        }
    }

    #[cfg(feature = "mpfr")]
    #[test]
    fn glibc_distribution_accuracy() {
        let Some(path) = glibc_libm_path_dist() else {
            return;
        };
        let lib = unsafe { Library::new(&path).expect("load glibc libm") };

        type CFn = unsafe extern "C" fn(f64) -> f64;
        type CFn2 = unsafe extern "C" fn(f64, f64) -> f64;

        let exp: libloading::Symbol<CFn> = unsafe { lib.get(b"exp").unwrap() };
        let log: libloading::Symbol<CFn> = unsafe { lib.get(b"log").unwrap() };
        let sin: libloading::Symbol<CFn> = unsafe { lib.get(b"sin").unwrap() };
        let cos: libloading::Symbol<CFn> = unsafe { lib.get(b"cos").unwrap() };
        let tan: libloading::Symbol<CFn> = unsafe { lib.get(b"tan").unwrap() };
        let exp2: libloading::Symbol<CFn> = unsafe { lib.get(b"exp2").unwrap() };
        let expm1: libloading::Symbol<CFn> = unsafe { lib.get(b"expm1").unwrap() };
        let log2: libloading::Symbol<CFn> = unsafe { lib.get(b"log2").unwrap() };
        let log10: libloading::Symbol<CFn> = unsafe { lib.get(b"log10").unwrap() };
        let log1p: libloading::Symbol<CFn> = unsafe { lib.get(b"log1p").unwrap() };
        let atan: libloading::Symbol<CFn> = unsafe { lib.get(b"atan").unwrap() };
        let atan2: libloading::Symbol<CFn2> = unsafe { lib.get(b"atan2").unwrap() };
        let asin: libloading::Symbol<CFn> = unsafe { lib.get(b"asin").unwrap() };
        let acos: libloading::Symbol<CFn> = unsafe { lib.get(b"acos").unwrap() };
        let sinh: libloading::Symbol<CFn> = unsafe { lib.get(b"sinh").unwrap() };
        let cosh: libloading::Symbol<CFn> = unsafe { lib.get(b"cosh").unwrap() };
        let tanh: libloading::Symbol<CFn> = unsafe { lib.get(b"tanh").unwrap() };
        let hypot: libloading::Symbol<CFn2> = unsafe { lib.get(b"hypot").unwrap() };
        let fmod: libloading::Symbol<CFn2> = unsafe { lib.get(b"fmod").unwrap() };
        let remainder: libloading::Symbol<CFn2> = unsafe { lib.get(b"remainder").unwrap() };
        let pow: libloading::Symbol<CFn2> = unsafe { lib.get(b"pow").unwrap() };
        let sqrt: libloading::Symbol<CFn> = unsafe { lib.get(b"sqrt").unwrap() };
        let cbrt: libloading::Symbol<CFn> = unsafe { lib.get(b"cbrt").unwrap() };

        let mut state = 0xdead_beef_cafe_f00du64;
        let samples = 128usize;

        for _ in 0..samples {
            let x = rand_range(&mut state, -100.0, 100.0);
            let mp = mpfr_exp_f64(x);
            let g = unsafe { exp(x) };
            let f = fastlibm::exp(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr exp ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast exp ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_f64_pos(&mut state);
            let mp = mpfr_ln_f64(x);
            let g = unsafe { log(x) };
            let f = fastlibm::ln(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr ln ulp {ulp_f} > 1 at {x}");
            assert!(ulp_f <= ulp_g, "fast ln ulp {ulp_f} > glibc {ulp_g} at {x}");
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let mp = mpfr_sin_f64(x);
            let g = unsafe { sin(x) };
            let f = fastlibm::sin(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr sin ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast sin ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let mp = mpfr_cos_f64(x);
            let g = unsafe { cos(x) };
            let f = fastlibm::cos(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr cos ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast cos ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let mp = mpfr_tan_f64(x);
            let g = unsafe { tan(x) };
            let f = fastlibm::tan(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr tan ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast tan ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -100.0, 100.0);
            let mp = mpfr_exp2_f64(x);
            let g = unsafe { exp2(x) };
            let f = fastlibm::exp2(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr exp2 ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast exp2 ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1.0, 1.0);
            let mp = mpfr_expm1_f64(x);
            let g = unsafe { expm1(x) };
            let f = fastlibm::expm1(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr expm1 ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast expm1 ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_f64_pos(&mut state);
            let mp = mpfr_log2_f64(x);
            let g = unsafe { log2(x) };
            let f = fastlibm::log2(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr log2 ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast log2 ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_f64_pos(&mut state);
            let mp = mpfr_log10_f64(x);
            let g = unsafe { log10(x) };
            let f = fastlibm::log10(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr log10 ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast log10 ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -0.9, 1e6);
            let mp = mpfr_log1p_f64(x);
            let g = unsafe { log1p(x) };
            let f = fastlibm::log1p(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr log1p ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast log1p ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, 0.1, 10.0);
            let y = rand_range(&mut state, -10.0, 10.0);
            let mp = mpfr_pow_f64(x, y);
            let g = unsafe { pow(x, y) };
            let f = fastlibm::pow(x, y);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr pow ulp {ulp_f} > 1 at {x},{y}");
            assert!(
                ulp_f <= ulp_g,
                "fast pow ulp {ulp_f} > glibc {ulp_g} at {x},{y}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let mp = mpfr_atan_f64(x);
            let g = unsafe { atan(x) };
            let f = fastlibm::atan(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr atan ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast atan ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1.0, 1.0);
            let mp = mpfr_asin_f64(x);
            let g = unsafe { asin(x) };
            let f = fastlibm::asin(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr asin ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast asin ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1.0, 1.0);
            let mp = mpfr_acos_f64(x);
            let g = unsafe { acos(x) };
            let f = fastlibm::acos(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr acos ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast acos ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let y = rand_range(&mut state, -1e6, 1e6);
            let x = rand_range(&mut state, -1e6, 1e6);
            if x == 0.0 && y == 0.0 {
                continue;
            }
            let mp = mpfr_atan2_f64(y, x);
            let g = unsafe { atan2(y, x) };
            let f = fastlibm::atan2(y, x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr atan2 ulp {ulp_f} > 1 at {y},{x}");
            assert!(
                ulp_f <= ulp_g,
                "fast atan2 ulp {ulp_f} > glibc {ulp_g} at {y},{x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            let mp = mpfr_sinh_f64(x);
            let g = unsafe { sinh(x) };
            let f = fastlibm::sinh(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr sinh ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast sinh ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            let mp = mpfr_cosh_f64(x);
            let g = unsafe { cosh(x) };
            let f = fastlibm::cosh(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr cosh ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast cosh ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            let mp = mpfr_tanh_f64(x);
            let g = unsafe { tanh(x) };
            let f = fastlibm::tanh(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr tanh ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast tanh ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e200, 1e200);
            let y = rand_range(&mut state, -1e200, 1e200);
            let mp = mpfr_hypot_f64(x, y);
            let g = unsafe { hypot(x, y) };
            let f = fastlibm::hypot(x, y);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr hypot ulp {ulp_f} > 1 at {x},{y}");
            assert!(
                ulp_f <= ulp_g,
                "fast hypot ulp {ulp_f} > glibc {ulp_g} at {x},{y}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let mut y = rand_range(&mut state, 1e-6, 1e6);
            if rand_u64(&mut state) & 1 == 0 {
                y = -y;
            }
            let mp = mpfr_fmod_f64(x, y);
            let g = unsafe { fmod(x, y) };
            let f = fastlibm::fmod(x, y);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr fmod ulp {ulp_f} > 1 at {x},{y}");
            assert!(
                ulp_f <= ulp_g,
                "fast fmod ulp {ulp_f} > glibc {ulp_g} at {x},{y}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let mut y = rand_range(&mut state, 1e-6, 1e6);
            if rand_u64(&mut state) & 1 == 0 {
                y = -y;
            }
            let mp = mpfr_remainder_f64(x, y);
            let g = unsafe { remainder(x, y) };
            let f = fastlibm::remainder(x, y);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr remainder ulp {ulp_f} > 1 at {x},{y}");
            assert!(
                ulp_f <= ulp_g,
                "fast remainder ulp {ulp_f} > glibc {ulp_g} at {x},{y}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, 0.0, 1e300);
            let mp = mpfr_sqrt_f64(x);
            let g = unsafe { sqrt(x) };
            let f = fastlibm::sqrt(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr sqrt ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast sqrt ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e300, 1e300);
            let mp = mpfr_cbrt_f64(x);
            let g = unsafe { cbrt(x) };
            let f = fastlibm::cbrt(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr cbrt ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast cbrt ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }
    }

    use proptest::prelude::*;
    proptest! {
        #[test]
        fn ptest_exp_special(x in proptest::sample::select(exp_special_inputs())) {
            let actual = fastlibm::exp(x);
            assert_ulp_eq_exp(actual, x, PROPTEST_ULP_TOL, &format!("exp special({x})"));
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_exp(x in -745.0..709.78_f64) {
            let actual = fastlibm::exp(x);
            assert_ulp_eq_exp(actual, x, PROPTEST_ULP_TOL, &format!("exp({x})"));
        }

        #[test]
        fn ptest_ln(x in proptest::num::f64::POSITIVE) {
            if x.is_finite() && x > 0.0 {
                let actual = fastlibm::ln(x);
                let expected = ln_reference(x);
                assert_ulp_eq(
                    actual,
                    expected,
                    PROPTEST_ULP_TOL,
                    &format!("ln({x})"),
                );
            }
        }

        #[test]
        fn ptest_sin(x in -1e20..1e20_f64) {
            let actual = fastlibm::sin(x);
            let expected = sin_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format!("sin({x})"),
            );
        }

        #[test]
        fn ptest_cos(x in -1e20..1e20_f64) {
            let actual = fastlibm::cos(x);
            let expected = cos_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format!("cos({x})"),
            );
        }

        #[test]
        fn ptest_sincos(x in -1e20..1e20_f64) {
            let (s_actual, c_actual) = fastlibm::sincos(x);
            assert_ulp_eq(
                s_actual,
                sin_reference(x),
                PROPTEST_ULP_TOL,
                &format!("sincos sin({x})"),
            );
            assert_ulp_eq(
                c_actual,
                cos_reference(x),
                PROPTEST_ULP_TOL,
                &format!("sincos cos({x})"),
            );
        }

        #[test]
        fn ptest_exp2(x in -1074.0..1024.0_f64) {
            let actual = fastlibm::exp2(x);
            let expected = exp2_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("exp2({x})"));
        }

        #[test]
        fn ptest_expm1(x in -50.0..50.0_f64) {
            let actual = fastlibm::expm1(x);
            let expected = expm1_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format!("expm1({x})"),
            );
        }

        #[test]
        fn ptest_log2(x in proptest::num::f64::POSITIVE) {
            if x.is_finite() && x > 0.0 {
                let actual = fastlibm::log2(x);
                let expected = log2_reference(x);
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("log2({x})"));
            }
        }

        #[test]
        fn ptest_log10(x in proptest::num::f64::POSITIVE) {
            if x.is_finite() && x > 0.0 {
                let actual = fastlibm::log10(x);
                let expected = log10_reference(x);
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("log10({x})"));
            }
        }

        #[test]
        fn ptest_log1p(x in -0.999_999_999_999_f64..1e6_f64) {
            let actual = fastlibm::log1p(x);
            let expected = log1p_reference(x);
            if expected.is_nan() {
                assert!(actual.is_nan(), "log1p({x}) expected NaN");
            } else {
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("log1p({x})"));
            }
        }

        #[test]
        fn ptest_tan(x in -1e6..1e6_f64) {
            let actual = fastlibm::tan(x);
            let expected = tan_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format!("tan({x})"),
            );
        }

        #[test]
        fn ptest_atan(x in -1e6..1e6_f64) {
            let actual = fastlibm::atan(x);
            let expected = atan_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("atan({x})"));
        }

        #[test]
        fn ptest_asin(x in -1.0..1.0_f64) {
            let actual = fastlibm::asin(x);
            let expected = asin_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("asin({x})"));
        }

        #[test]
        fn ptest_acos(x in -1.0..1.0_f64) {
            let actual = fastlibm::acos(x);
            let expected = acos_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("acos({x})"));
        }

        #[test]
        fn ptest_atan2(y in -1e6..1e6_f64, x in -1e6..1e6_f64) {
            let actual = fastlibm::atan2(y, x);
            let expected = atan2_reference(y, x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("atan2({y},{x})"));
        }

        #[test]
        fn ptest_hypot(x in -1e200..1e200_f64, y in -1e200..1e200_f64) {
            let actual = fastlibm::hypot(x, y);
            let expected = hypot_reference(x, y);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("hypot({x},{y})"));
        }

        #[test]
        fn ptest_sinh(x in -100.0..100.0_f64) {
            let actual = fastlibm::sinh(x);
            let expected = sinh_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("sinh({x})"));
        }

        #[test]
        fn ptest_cosh(x in -100.0..100.0_f64) {
            let actual = fastlibm::cosh(x);
            let expected = cosh_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("cosh({x})"));
        }

        #[test]
        fn ptest_tanh(x in -20.0..20.0_f64) {
            let actual = fastlibm::tanh(x);
            let expected = tanh_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("tanh({x})"));
        }

        #[test]
        fn ptest_pow(x in -10.0..10.0_f64, y in -10.0..10.0_f64) {
            let actual = fastlibm::pow(x, y);
            let expected = pow_reference(x, y);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format!("pow({x},{y})"),
            );
        }

        #[test]
        fn ptest_sqrt(x in -1e300..1e300_f64) {
            let actual = fastlibm::sqrt(x);
            let expected = sqrt_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("sqrt({x})"));
        }

        #[test]
        fn ptest_cbrt(x in -1e300..1e300_f64) {
            let actual = fastlibm::cbrt(x);
            let expected = cbrt_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format!("cbrt({x})"),
            );
        }

        #[test]
        fn ptest_fmod(
            x in -1e6..1e6_f64,
            y in prop_oneof![1e-6..1e6_f64, -1e6..-1e-6_f64],
        ) {
            let actual = fastlibm::fmod(x, y);
            let expected = fmod_reference(x, y);
            if expected.is_nan() {
                assert!(actual.is_nan(), "fmod({x},{y}) expected NaN");
            } else {
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("fmod({x},{y})"));
            }
        }

        #[test]
        fn ptest_remainder(
            x in -1e6..1e6_f64,
            y in prop_oneof![1e-6..1e6_f64, -1e6..-1e-6_f64],
        ) {
            let actual = fastlibm::remainder(x, y);
            let expected = remainder_reference(x, y);
            if expected.is_nan() {
                assert!(actual.is_nan(), "remainder({x},{y}) expected NaN");
            } else {
                assert_ulp_eq(
                    actual,
                    expected,
                    PROPTEST_ULP_TOL,
                    &format!("remainder({x},{y})"),
                );
            }
        }
    }
}

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
    use rug::Float;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6, PI, TAU};
    use std::path::Path;
    use std::string::String;
    use std::vec::Vec;
    use std::{eprintln, format};

    const MAX_ULP_TOL: f64 = 0.6;
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
    fn mpfr_ln_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.ln_mut();
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
    fn ln_reference(x: f64) -> f64 {
        mpfr_ln_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn ln_reference(x: f64) -> f64 {
        x.ln()
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

    fn assert_ulp_eq_exp(actual: f64, x: f64, context: &str) {
        let expected_std = x.exp();
        if expected_std.is_infinite() || expected_std.is_nan() {
            assert_ulp_eq(actual, expected_std, MAX_ULP_TOL, context);
            return;
        }
        let ulps_std = ulp_error(actual, expected_std);
        if ulps_std <= MAX_ULP_TOL {
            return;
        }

        #[cfg(feature = "mpfr")]
        {
            let expected_mpfr = mpfr_exp_f64(x);
            let ulps_mpfr = ulp_error(actual, expected_mpfr);
            if ulps_mpfr <= MAX_ULP_TOL {
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
            assert_ulp_eq_exp(actual, x, &context);
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

    use proptest::prelude::*;
    proptest! {
        #[test]
        fn ptest_exp_special(x in proptest::sample::select(exp_special_inputs())) {
            let actual = fastlibm::exp(x);
            let expected = x.exp();
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &format!("exp special({x})"));
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_exp(x in -745.0..709.78_f64) {
            let actual = fastlibm::exp(x);
            assert_ulp_eq_exp(actual, x, &format!("exp({x})"));
        }

        #[test]
        fn ptest_ln(x in proptest::num::f64::POSITIVE) {
            if x.is_finite() && x > 0.0 {
                let actual = fastlibm::ln(x);
                let expected = ln_reference(x);
                assert_ulp_eq(actual, expected, MAX_ULP_TOL, &format!("ln({x})"));
            }
        }

        #[test]
        fn ptest_sin(x in -1e20..1e20_f64) {
            let actual = fastlibm::sin(x);
            let expected = sin_reference(x);
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &format!("sin({x})"));
        }

        #[test]
        fn ptest_cos(x in -1e20..1e20_f64) {
            let actual = fastlibm::cos(x);
            let expected = cos_reference(x);
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &format!("cos({x})"));
        }

        #[test]
        fn ptest_sincos(x in -1e20..1e20_f64) {
            let (s_actual, c_actual) = fastlibm::sincos(x);
            assert_ulp_eq(
                s_actual,
                sin_reference(x),
                MAX_ULP_TOL,
                &format!("sincos sin({x})"),
            );
            assert_ulp_eq(
                c_actual,
                cos_reference(x),
                MAX_ULP_TOL,
                &format!("sincos cos({x})"),
            );
        }
    }
}

pub mod maths;

pub use maths::fastlibm;

#[cfg(test)]
mod tests {
    use super::fastlibm;
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6, PI, TAU};

    fn to_ordered_i64(x: f64) -> i64 {
        let bits = x.to_bits() as i64;
        if bits < 0 { i64::MIN - bits } else { bits }
    }

    fn ulp_diff(a: f64, b: f64) -> u64 {
        let a = to_ordered_i64(a);
        let b = to_ordered_i64(b);
        a.abs_diff(b)
    }

    fn assert_ulp_eq(actual: f64, expected: f64, max_ulps: u64, context: &str) {
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
        let ulps = ulp_diff(actual, expected);
        assert!(
            ulps <= max_ulps,
            "{context}: expected {expected}, got {actual} (ulps={ulps})"
        );
    }

    fn push_unique(values: &mut Vec<f64>, x: f64) {
        if !values.iter().any(|v| v.to_bits() == x.to_bits()) {
            values.push(x);
        }
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
            -709.78,
            -100.0,
            -10.0,
            -1.0,
            -0.5,
            -1e-12,
            -1e-6,
            0.0,
            1e-12,
            1e-6,
            0.5,
            1.0,
            2.0,
            10.0,
            100.0,
            700.0,
            709.0,
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

        for &x in &inputs {
            let expected = x.exp();
            let actual = fastlibm::exp(x);
            let context = format!("exp({x})");
            assert_ulp_eq(actual, expected, 1, &context);
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
        let mut inputs = Vec::new();
        let min_sub = f64::from_bits(1);
        let max_sub = f64::from_bits(0x000f_ffff_ffff_ffff);
        let specials = [
            min_sub,
            max_sub,
            f64::MIN_POSITIVE,
            f64::from_bits(0x3fefffffff_fffff),
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

        for &x in &inputs {
            let expected = x.ln();
            let actual = fastlibm::ln(x);
            let context = format!("ln({x})");
            assert_ulp_eq(actual, expected, 1, &context);
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
            assert_ulp_eq(sin_actual, sin_expected, 1, &format!("sin({x})"));
            assert_ulp_eq(cos_actual, cos_expected, 1, &format!("cos({x})"));
        }
    }

    #[test]
    fn sin_cos_matches_std_ulps() {
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

        for &x in &inputs {
            let sin_expected = x.sin();
            let cos_expected = x.cos();
            let sin_actual = fastlibm::sin(x);
            let cos_actual = fastlibm::cos(x);
            assert_ulp_eq(sin_actual, sin_expected, 1, &format!("sin({x})"));
            assert_ulp_eq(cos_actual, cos_expected, 1, &format!("cos({x})"));
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

            assert_ulp_eq(sin_neg, -sin_pos, 1, &format!("sin symmetry at {x}"));
            assert_ulp_eq(cos_neg, cos_pos, 1, &format!("cos symmetry at {x}"));
        }
    }
}

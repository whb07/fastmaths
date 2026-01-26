#![no_std]

#[cfg(test)]
extern crate std;

mod math;

pub use self::math::*;

#[cfg(test)]
mod tests {
    use crate as fastmaths;
    use libloading::Library;
    #[cfg(feature = "mpfr")]
    use rug::{Float, ops::Pow};
    use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, FRAC_PI_6, PI, TAU};
    use std::path::Path;
    use std::string::String;
    #[cfg(not(feature = "mpfr"))]
    use std::sync::OnceLock;
    use std::vec;
    use std::vec::Vec;
    use std::{eprintln, format};

    const MAX_ULP_TOL: f64 = 1.0;
    const DERIVED_ULP_TOL: f64 = 1.0;
    const PROPTEST_ULP_TOL: f64 = 1.0;
    const COMPOSED_ULP_TOL: f64 = 4.0;
    #[cfg(feature = "mpfr")]
    const TANH_ULP_TOL: f64 = 1.0;
    #[cfg(not(feature = "mpfr"))]
    const TANH_ULP_TOL: f64 = 3.0;
    #[cfg(feature = "mpfr")]
    const ATANH_ULP_TOL: f64 = 1.0;
    #[cfg(not(feature = "mpfr"))]
    const ATANH_ULP_TOL: f64 = 2.0;
    #[cfg(feature = "mpfr")]
    const MPFR_PREC: u32 = 256;
    #[cfg(feature = "mpfr")]
    const MPFR_TRIG_LIMIT: f64 = 1.0e6;
    const ROUNDTRIP_EXP_MIN_ABS: f64 = 0.5;

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

    fn classify_f64(x: f64) -> &'static str {
        if x.is_nan() {
            return "nan";
        }
        if x.is_infinite() {
            return "inf";
        }
        if x == 0.0 {
            return if x.is_sign_negative() { "-0" } else { "+0" };
        }
        if x.is_subnormal() {
            return "subnormal";
        }
        "normal"
    }

    fn f64_details(x: f64) -> String {
        let bits = x.to_bits();
        let sign = bits >> 63;
        let exp = (bits >> 52) & 0x7ff;
        let mant = bits & ((1u64 << 52) - 1);
        format!(
            "x={x} bits=0x{bits:016x} sign={sign} exp=0x{exp:03x} mant=0x{mant:013x} class={}",
            classify_f64(x)
        )
    }

    fn format_case(func: &str, x: f64, bucket: &str) -> String {
        format!("{func}({x}) [bucket={bucket}] {}", f64_details(x))
    }

    fn format_case2(func: &str, x: f64, y: f64, bucket: &str) -> String {
        format!(
            "{func}({x},{y}) [bucket={bucket}] x:{} y:{}",
            f64_details(x),
            f64_details(y)
        )
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
    fn mpfr_floor_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.floor_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_ceil_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.ceil_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_trunc_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.trunc_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_round_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.round_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_rint_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.round_even_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_fma_f64(x: f64, y: f64, z: f64) -> f64 {
        let mut a = Float::with_val(MPFR_PREC, x);
        let b = Float::with_val(MPFR_PREC, y);
        let c = Float::with_val(MPFR_PREC, z);
        a.mul_add_mut(&b, &c);
        a.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_frexp_f64(x: f64) -> (f64, i32) {
        let mut v = Float::with_val(MPFR_PREC, x);
        let exp = v.frexp_mut();
        (v.to_f64(), exp)
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_ldexp_f64(x: f64, n: i32) -> f64 {
        let v = Float::with_val(MPFR_PREC, x);
        v.as_shl(n).to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_scalbln_f64(x: f64, n: i64) -> f64 {
        if n > i32::MAX as i64 {
            return if x.is_sign_negative() {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }
        if n < i32::MIN as i64 {
            return 0.0_f64.copysign(x);
        }
        mpfr_ldexp_f64(x, n as i32)
    }

    fn clamp_f64_to_i64(x: f64) -> i64 {
        #[cfg(any(target_arch = "aarch64", target_arch = "arm"))]
        {
            if x.is_nan() {
                return 0;
            }
            if x.is_infinite() {
                return if x.is_sign_negative() {
                    i64::MIN
                } else {
                    i64::MAX
                };
            }
            if x > i64::MAX as f64 {
                return i64::MAX;
            }
            if x < i64::MIN as f64 {
                return i64::MIN;
            }
            return x as i64;
        }
        #[cfg(not(any(target_arch = "aarch64", target_arch = "arm")))]
        {
            if !x.is_finite() {
                return i64::MIN;
            }
            if x > i64::MAX as f64 || x < i64::MIN as f64 {
                i64::MIN
            } else {
                x as i64
            }
        }
    }

    fn remquo_sig_exp(bits: u64) -> (u64, i32) {
        let exp = ((bits >> 52) & 0x7ff) as i32;
        let mant = bits & 0x000f_ffff_ffff_ffffu64;
        if exp == 0 {
            if mant == 0 {
                return (0, 0);
            }
            return (mant, -1074);
        }
        (mant | (1u64 << 52), exp - 1023 - 52)
    }

    fn remquo_pow2_mod(mut exp: u32, modulus: u128) -> u128 {
        let mut result = 1u128 % modulus;
        let mut base = 2u128 % modulus;
        while exp != 0 {
            if (exp & 1) != 0 {
                result = (result * base) % modulus;
            }
            base = (base * base) % modulus;
            exp >>= 1;
        }
        result
    }

    fn remquo_quotient_mod8(x: f64, y: f64) -> i32 {
        let ax_bits = x.to_bits() & 0x7fff_ffff_ffff_ffffu64;
        let ay_bits = y.to_bits() & 0x7fff_ffff_ffff_ffffu64;
        let (sigx, ex) = remquo_sig_exp(ax_bits);
        let (sigy, ey) = remquo_sig_exp(ay_bits);
        if sigx == 0 || sigy == 0 {
            return 0;
        }
        let shift = ex - ey;
        if shift >= 0 {
            let d = sigy as u128;
            let d16 = d << 4;
            let pow2_mod_d = remquo_pow2_mod(shift as u32, d);
            let pow2_mod_d16 = remquo_pow2_mod(shift as u32, d16);
            let n_mod = ((sigx as u128) % d) * pow2_mod_d % d;
            let n_mod16 = ((sigx as u128) % d16) * pow2_mod_d16 % d16;
            let r = n_mod;
            let q_mod16 = ((n_mod16 + d16 - r) / d) & 0x0f;
            let mut q_mod8 = (q_mod16 & 0x7) as u8;
            let twice_r = r << 1;
            if twice_r > d || (twice_r == d && (q_mod16 & 1) != 0) {
                q_mod8 = q_mod8.wrapping_add(1) & 0x7;
            }
            return q_mod8 as i32;
        }
        let s = (-shift) as u32;
        if s > 60 {
            return 0;
        }
        let d = (sigy as u128) << s;
        let n = sigx as u128;
        let q0 = n / d;
        let r = n % d;
        let mut q = q0;
        let twice_r = r << 1;
        if twice_r > d || (twice_r == d && (q0 & 1) != 0) {
            q += 1;
        }
        (q as u8 & 0x7) as i32
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
    fn mpfr_asinh_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.asinh_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_acosh_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.acosh_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_atanh_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.atanh_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_erf_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.erf_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_erfc_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.erfc_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_exp10_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.exp10_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_lgamma_f64(x: f64) -> f64 {
        let v = Float::with_val(MPFR_PREC, x);
        let (lg, _) = v.ln_abs_gamma();
        lg.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_tgamma_f64(x: f64) -> f64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.gamma_mut();
        v.to_f64()
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_logb_f64(x: f64) -> f64 {
        if x == 0.0 {
            return f64::NEG_INFINITY;
        }
        if x.is_infinite() {
            return f64::INFINITY;
        }
        if x.is_nan() {
            return f64::NAN;
        }
        let v = Float::with_val(MPFR_PREC, x);
        let Some((int, exp)) = v.to_integer_exp() else {
            return f64::NAN;
        };
        if int == 0 {
            return f64::NEG_INFINITY;
        }
        let bits = int.significant_bits() as i32;
        (exp + bits - 1) as f64
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_ilogb_i32(x: f64) -> i32 {
        if x == 0.0 {
            return i32::MIN;
        }
        if x.is_infinite() || x.is_nan() {
            return i32::MAX;
        }
        let v = Float::with_val(MPFR_PREC, x);
        let Some((int, exp)) = v.to_integer_exp() else {
            return i32::MAX;
        };
        if int == 0 {
            return i32::MIN;
        }
        let bits = int.significant_bits() as i32;
        exp + bits - 1
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_nextafter_f64(x: f64, y: f64) -> f64 {
        const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
        if x.is_nan() || y.is_nan() {
            return f64::NAN;
        }
        if x == y {
            return y;
        }
        if x == 0.0 {
            let sign = y.to_bits() & SIGN_MASK;
            return f64::from_bits(sign | 1);
        }
        let mut ux = x.to_bits();
        let sx = ux & SIGN_MASK;
        if x > y {
            if sx == 0 {
                ux -= 1;
            } else {
                ux += 1;
            }
        } else if sx == 0 {
            ux += 1;
        } else {
            ux -= 1;
        }
        f64::from_bits(ux)
    }

    #[cfg(feature = "mpfr")]
    fn mpfr_modf_f64(x: f64) -> (f64, f64) {
        let mut v = Float::with_val(MPFR_PREC, x);
        let mut fract = Float::new(MPFR_PREC);
        v.trunc_fract_mut(&mut fract);
        (fract.to_f64(), v.to_f64())
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
    fn floor_reference(x: f64) -> f64 {
        mpfr_floor_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn floor_reference(x: f64) -> f64 {
        x.floor()
    }

    #[cfg(feature = "mpfr")]
    fn ceil_reference(x: f64) -> f64 {
        mpfr_ceil_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn ceil_reference(x: f64) -> f64 {
        x.ceil()
    }

    #[cfg(feature = "mpfr")]
    fn trunc_reference(x: f64) -> f64 {
        mpfr_trunc_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn trunc_reference(x: f64) -> f64 {
        x.trunc()
    }

    #[cfg(feature = "mpfr")]
    fn round_reference(x: f64) -> f64 {
        mpfr_round_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn round_reference(x: f64) -> f64 {
        x.round()
    }

    #[cfg(feature = "mpfr")]
    fn rint_reference(x: f64) -> f64 {
        mpfr_rint_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn rint_reference(x: f64) -> f64 {
        if let Some(f) = glibc_sym_f64(b"rint") {
            // Safety: symbol is loaded from libm with correct signature.
            return unsafe { f(x) };
        }
        x.round_ties_even()
    }

    #[cfg(feature = "mpfr")]
    fn nearbyint_reference(x: f64) -> f64 {
        mpfr_rint_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn nearbyint_reference(x: f64) -> f64 {
        if let Some(f) = glibc_sym_f64(b"nearbyint") {
            // Safety: symbol is loaded from libm with correct signature.
            return unsafe { f(x) };
        }
        rint_reference(x)
    }

    #[cfg(feature = "mpfr")]
    fn lrint_reference(x: f64) -> i64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.round_even_mut();
        clamp_f64_to_i64(v.to_f64())
    }

    #[cfg(not(feature = "mpfr"))]
    fn lrint_reference(x: f64) -> i64 {
        if let Some(f) = glibc_sym_i64(b"lrint") {
            // Safety: symbol is loaded from libm with correct signature.
            return unsafe { f(x) };
        }
        clamp_f64_to_i64(x.round_ties_even())
    }

    #[cfg(feature = "mpfr")]
    fn llrint_reference(x: f64) -> i64 {
        lrint_reference(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn llrint_reference(x: f64) -> i64 {
        if let Some(f) = glibc_sym_i64(b"llrint") {
            // Safety: symbol is loaded from libm with correct signature.
            return unsafe { f(x) };
        }
        lrint_reference(x)
    }

    #[cfg(feature = "mpfr")]
    fn lround_reference(x: f64) -> i64 {
        let mut v = Float::with_val(MPFR_PREC, x);
        v.round_mut();
        clamp_f64_to_i64(v.to_f64())
    }

    #[cfg(not(feature = "mpfr"))]
    fn lround_reference(x: f64) -> i64 {
        if let Some(f) = glibc_sym_i64(b"lround") {
            // Safety: symbol is loaded from libm with correct signature.
            return unsafe { f(x) };
        }
        clamp_f64_to_i64(x.round())
    }

    #[cfg(feature = "mpfr")]
    fn llround_reference(x: f64) -> i64 {
        lround_reference(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn llround_reference(x: f64) -> i64 {
        if let Some(f) = glibc_sym_i64(b"llround") {
            // Safety: symbol is loaded from libm with correct signature.
            return unsafe { f(x) };
        }
        clamp_f64_to_i64(x.round())
    }

    #[cfg(feature = "mpfr")]
    fn fma_reference(x: f64, y: f64, z: f64) -> f64 {
        mpfr_fma_f64(x, y, z)
    }

    #[cfg(not(feature = "mpfr"))]
    fn fma_reference(x: f64, y: f64, z: f64) -> f64 {
        x.mul_add(y, z)
    }

    fn copysign_reference(x: f64, y: f64) -> f64 {
        x.copysign(y)
    }

    fn fabs_reference(x: f64) -> f64 {
        x.abs()
    }

    #[cfg(feature = "mpfr")]
    fn frexp_reference(x: f64) -> (f64, i32) {
        mpfr_frexp_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn frexp_reference(x: f64) -> (f64, i32) {
        fastmaths::frexp(x)
    }

    #[cfg(feature = "mpfr")]
    fn ldexp_reference(x: f64, n: i32) -> f64 {
        mpfr_ldexp_f64(x, n)
    }

    #[cfg(not(feature = "mpfr"))]
    fn ldexp_reference(x: f64, n: i32) -> f64 {
        fastmaths::ldexp(x, n)
    }

    #[cfg(feature = "mpfr")]
    fn scalbn_reference(x: f64, n: i32) -> f64 {
        mpfr_ldexp_f64(x, n)
    }

    #[cfg(not(feature = "mpfr"))]
    fn scalbn_reference(x: f64, n: i32) -> f64 {
        fastmaths::scalbn(x, n)
    }

    #[cfg(feature = "mpfr")]
    fn scalbln_reference(x: f64, n: i64) -> f64 {
        mpfr_scalbln_f64(x, n)
    }

    #[cfg(not(feature = "mpfr"))]
    fn scalbln_reference(x: f64, n: i64) -> f64 {
        fastmaths::scalbln(x, n)
    }

    #[cfg(feature = "mpfr")]
    fn remquo_reference(x: f64, y: f64) -> (f64, i32) {
        if x.is_nan() || y.is_nan() || y == 0.0 || x.is_infinite() {
            return (f64::NAN, 0);
        }
        let mut q = Float::with_val(MPFR_PREC, x);
        let vy = Float::with_val(MPFR_PREC, y);
        q /= &vy;
        q.round_even_mut();
        let mut r = Float::with_val(MPFR_PREC, x);
        r -= &vy * &q;
        let mut quo = remquo_quotient_mod8(x, y);
        if (x.is_sign_negative() ^ y.is_sign_negative()) && quo != 0 {
            quo = -quo;
        }
        (r.to_f64(), quo)
    }

    #[cfg(not(feature = "mpfr"))]
    fn remquo_reference(x: f64, y: f64) -> (f64, i32) {
        fastmaths::remquo(x, y)
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
        if let Some(f) = glibc_sym_f64(b"tanh") {
            // Safety: symbol is loaded from libm with correct signature.
            return unsafe { f(x) };
        }
        x.tanh()
    }

    #[cfg(feature = "mpfr")]
    fn asinh_reference(x: f64) -> f64 {
        mpfr_asinh_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn asinh_reference(x: f64) -> f64 {
        if let Some(f) = glibc_sym_f64(b"asinh") {
            // Safety: symbol is loaded from libm with correct signature.
            return unsafe { f(x) };
        }
        x.asinh()
    }

    #[cfg(feature = "mpfr")]
    fn acosh_reference(x: f64) -> f64 {
        mpfr_acosh_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn acosh_reference(x: f64) -> f64 {
        if let Some(f) = glibc_sym_f64(b"acosh") {
            // Safety: symbol is loaded from libm with correct signature.
            return unsafe { f(x) };
        }
        x.acosh()
    }

    #[cfg(feature = "mpfr")]
    fn atanh_reference(x: f64) -> f64 {
        mpfr_atanh_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn atanh_reference(x: f64) -> f64 {
        if let Some(f) = glibc_sym_f64(b"atanh") {
            // Safety: symbol is loaded from libm with correct signature.
            return unsafe { f(x) };
        }
        x.atanh()
    }

    #[cfg(feature = "mpfr")]
    fn erf_reference(x: f64) -> f64 {
        mpfr_erf_f64(x)
    }

    #[cfg(feature = "mpfr")]
    fn erfc_reference(x: f64) -> f64 {
        mpfr_erfc_f64(x)
    }

    #[cfg(feature = "mpfr")]
    fn exp10_reference(x: f64) -> f64 {
        mpfr_exp10_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn exp10_reference(x: f64) -> f64 {
        10.0f64.powf(x)
    }

    #[cfg(feature = "mpfr")]
    fn lgamma_reference(x: f64) -> f64 {
        mpfr_lgamma_f64(x)
    }

    #[cfg(feature = "mpfr")]
    fn tgamma_reference(x: f64) -> f64 {
        mpfr_tgamma_f64(x)
    }

    #[cfg(feature = "mpfr")]
    fn logb_reference(x: f64) -> f64 {
        mpfr_logb_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn logb_reference(x: f64) -> f64 {
        fastmaths::logb(x)
    }

    #[cfg(feature = "mpfr")]
    fn ilogb_reference(x: f64) -> i32 {
        mpfr_ilogb_i32(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn ilogb_reference(x: f64) -> i32 {
        fastmaths::ilogb(x)
    }

    #[cfg(feature = "mpfr")]
    fn nextafter_reference(x: f64, y: f64) -> f64 {
        mpfr_nextafter_f64(x, y)
    }

    #[cfg(not(feature = "mpfr"))]
    fn nextafter_reference(x: f64, y: f64) -> f64 {
        fastmaths::nextafter(x, y)
    }

    #[cfg(feature = "mpfr")]
    fn modf_reference(x: f64) -> (f64, f64) {
        mpfr_modf_f64(x)
    }

    #[cfg(not(feature = "mpfr"))]
    fn modf_reference(x: f64) -> (f64, f64) {
        fastmaths::modf(x)
    }

    fn fdim_reference(x: f64, y: f64) -> f64 {
        if x.is_nan() || y.is_nan() {
            return f64::NAN;
        }
        if x > y { x - y } else { 0.0 }
    }

    fn fmax_reference(x: f64, y: f64) -> f64 {
        if x.is_nan() {
            return y;
        }
        if y.is_nan() {
            return x;
        }
        if x == 0.0 && y == 0.0 {
            let sx = x.to_bits() & 0x8000_0000_0000_0000u64;
            let sy = y.to_bits() & 0x8000_0000_0000_0000u64;
            if sx == 0 || sy == 0 {
                0.0
            } else {
                f64::from_bits(0x8000_0000_0000_0000u64)
            }
        } else if x > y {
            x
        } else {
            y
        }
    }

    fn fmin_reference(x: f64, y: f64) -> f64 {
        if x.is_nan() {
            return y;
        }
        if y.is_nan() {
            return x;
        }
        if x == 0.0 && y == 0.0 {
            let sx = x.to_bits() & 0x8000_0000_0000_0000u64;
            let sy = y.to_bits() & 0x8000_0000_0000_0000u64;
            if sx != 0 || sy != 0 {
                f64::from_bits(0x8000_0000_0000_0000u64)
            } else {
                0.0
            }
        } else if x < y {
            x
        } else {
            y
        }
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
        if let Some(f) = glibc_sym_f64_f64(b"remainder") {
            // Safety: symbol is loaded from libm with correct signature.
            return unsafe { f(x, y) };
        }
        if !x.is_finite() || !y.is_finite() || y == 0.0 {
            return f64::NAN;
        }
        let ay = y.abs();
        let mut r = (x % (y + y)).abs();
        if r + r > ay {
            r -= ay;
            if r + r >= ay {
                r -= ay;
            } else if r == 0.0 {
                r = 0.0;
            }
        }
        if x.is_sign_negative() { -r } else { r }
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
            -700.0,
            -100.0,
            -50.0,
            -20.0,
            -10.0,
            -1.0,
            -0.5,
            -1e-6,
            -1e-12,
            -0.0,
            0.0,
            1e-12,
            1e-6,
            0.5,
            0.839_264_735_768_179_8,
            1.0,
            10.0,
            20.0,
            50.0,
            100.0,
            700.0,
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

    fn asinh_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            -1e20, -1e10, -100.0, -10.0, -2.0, -1.0, -1e-6, -1e-12, -0.0, 0.0, 1e-12, 1e-6, 1.0,
            2.0, 10.0, 100.0, 1e10, 1e20,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -100..=100 {
            push_unique(&mut inputs, (i as f64) * 0.1);
        }
        inputs
    }

    fn acosh_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            1.0,
            1.0 + 1e-12,
            1.0 + 1e-6,
            1.125,
            1.111_059_132_586_820_9,
            1.5,
            2.0,
            10.0,
            100.0,
            1e6,
            1e20,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in 0..=200 {
            push_unique(&mut inputs, 1.0 + (i as f64) * 0.01);
        }
        inputs
    }

    fn atanh_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            -0.999_999_999_999,
            -0.99,
            -0.9,
            -0.5,
            -1e-6,
            -1e-12,
            -0.0,
            0.0,
            1e-12,
            1e-6,
            0.5,
            0.9,
            0.99,
            0.999_999_999_999,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -99..=99 {
            push_unique(&mut inputs, (i as f64) / 100.0);
        }
        inputs
    }

    #[cfg(feature = "mpfr")]
    fn erf_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            -6.0, -3.0, -2.0, -1.0, -0.5, -1e-6, -1e-12, -0.0, 0.0, 1e-12, 1e-6, 0.5, 1.0, 2.0,
            3.0, 6.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -60..=60 {
            push_unique(&mut inputs, (i as f64) * 0.1);
        }
        inputs
    }

    #[cfg(feature = "mpfr")]
    fn erfc_inputs() -> Vec<f64> {
        erf_inputs()
    }

    fn exp10_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            -308.0, -300.0, -100.0, -10.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 2.0, 10.0, 100.0, 300.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -20..=20 {
            push_unique(&mut inputs, (i as f64) * 0.5);
        }
        inputs
    }

    #[cfg(feature = "mpfr")]
    fn lgamma_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            -10.5, -9.5, -4.5, -3.5, -2.5, -1.5, -0.5, -4.0, -3.0, -2.0, -1.0, -0.0, 0.0, 0.1, 0.5,
            1.0, 1.5, 2.0, 3.0, 10.0, 50.0, 100.0, 170.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -20..=20 {
            let base = i as f64;
            push_unique(&mut inputs, base + 0.25);
            push_unique(&mut inputs, base + 0.5);
            push_unique(&mut inputs, base + 0.75);
        }
        inputs
    }

    #[cfg(feature = "mpfr")]
    fn tgamma_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            -10.5, -9.5, -4.5, -3.5, -2.5, -1.5, -0.5, -4.0, -3.0, -2.0, -1.0, -0.0, 0.0, 0.1, 0.5,
            1.0, 1.5, 2.0, 3.0, 10.0, 50.0, 100.0, 170.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -20..=20 {
            let base = i as f64;
            push_unique(&mut inputs, base + 0.25);
            push_unique(&mut inputs, base + 0.5);
            push_unique(&mut inputs, base + 0.75);
        }
        inputs
    }

    fn logb_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            f64::MIN_POSITIVE,
            -f64::MIN_POSITIVE,
            f64::from_bits(1),
            -f64::from_bits(1),
            1e-300,
            -1e-300,
            1e-10,
            -1e-10,
            1.0,
            -1.0,
            2.0,
            -2.0,
            1024.0,
            -1024.0,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -100..=100 {
            let x = (i as f64) * 0.25;
            if x != 0.0 {
                push_unique(&mut inputs, x);
            }
        }
        inputs
    }

    fn ilogb_inputs() -> Vec<f64> {
        logb_inputs()
    }

    fn modf_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            0.0,
            -0.0,
            0.5,
            -0.5,
            1.5,
            -1.5,
            10.25,
            -10.25,
            1e20,
            -1e20,
            f64::MIN_POSITIVE,
            -f64::MIN_POSITIVE,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -100..=100 {
            push_unique(&mut inputs, (i as f64) * 0.1);
        }
        inputs
    }

    fn rounding_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            -3.5,
            -2.5,
            -1.5,
            -1.0,
            -0.9,
            -0.5,
            -0.1,
            -0.0,
            0.0,
            0.1,
            0.5,
            0.9,
            1.0,
            1.5,
            2.5,
            3.5,
            1e6,
            -1e6,
            2f64.powi(52),
            -2f64.powi(52),
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -100..=100 {
            push_unique(&mut inputs, (i as f64) * 0.25);
        }
        inputs
    }

    fn scaling_inputs() -> Vec<f64> {
        let mut inputs = Vec::new();
        let specials = [
            0.0,
            -0.0,
            f64::MIN_POSITIVE,
            -f64::MIN_POSITIVE,
            1e-300,
            -1e-300,
            1e-10,
            -1e-10,
            1.0,
            -1.0,
            2.0,
            -2.0,
            1024.0,
            -1024.0,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ];
        for &x in &specials {
            push_unique(&mut inputs, x);
        }
        for i in -100..=100 {
            push_unique(&mut inputs, (i as f64) * 0.5);
        }
        inputs
    }

    fn scalbn_inputs() -> Vec<(f64, i32)> {
        vec![
            (1.0, 0),
            (1.0, 1),
            (1.0, -1),
            (1e-300, 10),
            (1e-300, -10),
            (1e300, -10),
            (-2.5, 3),
            (-2.5, -3),
        ]
    }

    fn scalbln_inputs() -> Vec<(f64, i64)> {
        vec![
            (1.0, 0),
            (1.0, 1),
            (1.0, -1),
            (1e-300, 10),
            (1e-300, -10),
            (1e300, -10),
            (-2.5, 3),
            (-2.5, -3),
            (1.0, i64::from(i32::MAX)),
            (1.0, i64::from(i32::MIN)),
        ]
    }

    fn remquo_inputs() -> Vec<(f64, f64)> {
        vec![
            (5.3, 2.0),
            (-5.3, 2.0),
            (5.3, -2.0),
            (-5.3, -2.0),
            (1.0, 0.5),
            (10.0, 3.0),
            (1e-10, 1e-3),
            (1e10, 3.0),
        ]
    }

    fn fma_inputs() -> Vec<(f64, f64, f64)> {
        vec![
            (0.0, 0.0, 0.0),
            (1.0, 2.0, 3.0),
            (-1.0, 2.0, 3.0),
            (1e300, 1e-300, 1.0),
            (1e200, 1e200, f64::NEG_INFINITY),
            (1.2345, 6.789, -3.21),
        ]
    }

    fn fdim_inputs() -> Vec<(f64, f64)> {
        vec![
            (0.0, 0.0),
            (-0.0, 0.0),
            (1.0, 2.0),
            (2.0, 1.0),
            (-1.0, -2.0),
            (5.3, 2.1),
            (-5.3, 2.1),
            (1e20, 3.0),
        ]
    }

    fn fmax_inputs() -> Vec<(f64, f64)> {
        vec![
            (0.0, -0.0),
            (-0.0, 0.0),
            (1.0, 2.0),
            (-1.0, -2.0),
            (f64::NAN, 1.0),
            (1.0, f64::NAN),
            (f64::NAN, f64::NAN),
        ]
    }

    fn fmin_inputs() -> Vec<(f64, f64)> {
        vec![
            (0.0, -0.0),
            (-0.0, 0.0),
            (1.0, 2.0),
            (-1.0, -2.0),
            (f64::NAN, 1.0),
            (1.0, f64::NAN),
            (f64::NAN, f64::NAN),
        ]
    }

    fn nextafter_inputs() -> Vec<(f64, f64)> {
        vec![
            (0.0, 1.0),
            (0.0, -1.0),
            (-0.0, 1.0),
            (1.0, 2.0),
            (1.0, 0.0),
            (-1.0, -2.0),
            (-1.0, 0.0),
            (1e-300, 0.0),
            (-1e-300, 0.0),
        ]
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
        if std::env::var("FASTMATHS_GLIBC_TEST").is_err() {
            return None;
        }
        let path = std::env::var("FASTMATHS_GLIBC_LIBM")
            .unwrap_or_else(|_| String::from("/tmp/maths/glibc-build/math/libm.so"));
        if !Path::new(&path).exists() {
            eprintln!("glibc libm not found at {path}");
            return None;
        }
        Some(path)
    }

    fn glibc_libm_path_dist() -> Option<String> {
        if std::env::var("FASTMATHS_GLIBC_DIST").is_err() {
            return None;
        }
        let path = std::env::var("FASTMATHS_GLIBC_LIBM")
            .unwrap_or_else(|_| String::from("/tmp/maths/glibc-build/math/libm.so"));
        if !Path::new(&path).exists() {
            eprintln!("glibc libm not found at {path}");
            return None;
        }
        Some(path)
    }

    #[cfg(not(feature = "mpfr"))]
    fn glibc_lib_any() -> Option<&'static Library> {
        static LIB: OnceLock<Option<Library>> = OnceLock::new();
        LIB.get_or_init(|| {
            if let Ok(path) = std::env::var("FASTMATHS_GLIBC_LIBM") {
                if Path::new(&path).exists() {
                    if let Ok(lib) = unsafe { Library::new(&path) } {
                        return Some(lib);
                    }
                }
            }
            let candidates = [
                "/lib/x86_64-linux-gnu/libm.so.6",
                "/usr/lib/x86_64-linux-gnu/libm.so.6",
                "libm.so.6",
            ];
            for path in candidates {
                if let Ok(lib) = unsafe { Library::new(path) } {
                    return Some(lib);
                }
            }
            None
        })
        .as_ref()
    }

    #[cfg(not(feature = "mpfr"))]
    fn glibc_sym_f64(name: &'static [u8]) -> Option<unsafe extern "C" fn(f64) -> f64> {
        let lib = glibc_lib_any()?;
        unsafe {
            lib.get::<unsafe extern "C" fn(f64) -> f64>(name)
                .ok()
                .map(|s| *s)
        }
    }

    #[cfg(not(feature = "mpfr"))]
    fn glibc_sym_f64_f64(name: &'static [u8]) -> Option<unsafe extern "C" fn(f64, f64) -> f64> {
        let lib = glibc_lib_any()?;
        unsafe {
            lib.get::<unsafe extern "C" fn(f64, f64) -> f64>(name)
                .ok()
                .map(|s| *s)
        }
    }

    #[cfg(not(feature = "mpfr"))]
    fn glibc_sym_i64(name: &'static [u8]) -> Option<unsafe extern "C" fn(f64) -> i64> {
        let lib = glibc_lib_any()?;
        unsafe {
            lib.get::<unsafe extern "C" fn(f64) -> i64>(name)
                .ok()
                .map(|s| *s)
        }
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

        assert!(fastmaths::exp(nan).is_nan());
        assert_eq!(fastmaths::exp(pos_inf), f64::INFINITY);
        assert_eq!(fastmaths::exp(neg_inf), 0.0);
        assert_eq!(fastmaths::exp(0.0).to_bits(), 1.0f64.to_bits());
        assert_eq!(fastmaths::exp(-0.0).to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn exp_matches_std_ulps() {
        let inputs = exp_inputs();

        for &x in &inputs {
            let actual = fastmaths::exp(x);
            let context = format!("exp({x})");
            assert_ulp_eq_exp(actual, x, MAX_ULP_TOL, &context);
        }
    }

    #[test]
    fn ln_special_cases() {
        let nan = f64::NAN;
        let pos_inf = f64::INFINITY;

        assert!(fastmaths::ln(nan).is_nan());
        assert_eq!(fastmaths::ln(pos_inf), f64::INFINITY);
        assert_eq!(fastmaths::ln(0.0), f64::NEG_INFINITY);
        assert_eq!(fastmaths::ln(-0.0), (-0.0f64).ln());
        assert!(fastmaths::ln(-1.0).is_nan());
    }

    #[test]
    fn ln_matches_std_ulps() {
        let inputs = ln_inputs();

        for &x in &inputs {
            let expected = x.ln();
            let actual = fastmaths::ln(x);
            let context = format!("ln({x})");
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &context);
        }
    }

    #[test]
    fn log1p_special_cases() {
        assert!(fastmaths::log1p(f64::NAN).is_nan());
        assert_eq!(fastmaths::log1p(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::log1p(0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(fastmaths::log1p(-0.0).to_bits(), (-0.0f64).to_bits());
        assert_eq!(fastmaths::log1p(-1.0), f64::NEG_INFINITY);
        assert!(fastmaths::log1p(-1.5).is_nan());
    }

    #[test]
    fn log1p_matches_reference_ulps() {
        for &x in &log1p_inputs() {
            let expected = log1p_reference(x);
            let actual = fastmaths::log1p(x);
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
    fn rounding_special_cases() {
        assert!(fastmaths::floor(f64::NAN).is_nan());
        assert!(fastmaths::ceil(f64::NAN).is_nan());
        assert!(fastmaths::trunc(f64::NAN).is_nan());
        assert!(fastmaths::round(f64::NAN).is_nan());
        assert!(fastmaths::rint(f64::NAN).is_nan());
        assert!(fastmaths::nearbyint(f64::NAN).is_nan());

        assert_eq!(fastmaths::floor(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::ceil(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::trunc(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::round(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::rint(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::nearbyint(f64::INFINITY), f64::INFINITY);

        assert_eq!(fastmaths::floor(f64::NEG_INFINITY), f64::NEG_INFINITY);
        assert_eq!(fastmaths::ceil(f64::NEG_INFINITY), f64::NEG_INFINITY);
        assert_eq!(fastmaths::trunc(f64::NEG_INFINITY), f64::NEG_INFINITY);
        assert_eq!(fastmaths::round(f64::NEG_INFINITY), f64::NEG_INFINITY);
        assert_eq!(fastmaths::rint(f64::NEG_INFINITY), f64::NEG_INFINITY);
        assert_eq!(fastmaths::nearbyint(f64::NEG_INFINITY), f64::NEG_INFINITY);

        assert_eq!(fastmaths::trunc(-0.0).to_bits(), (-0.0_f64).to_bits());
        assert_eq!(fastmaths::ceil(-0.3).to_bits(), (-0.0_f64).to_bits());
        assert_eq!(fastmaths::round(-0.3).to_bits(), (-0.0_f64).to_bits());
        assert_eq!(fastmaths::rint(-0.3).to_bits(), (-0.0_f64).to_bits());
    }

    #[test]
    fn rounding_matches_reference_ulps() {
        for &x in &rounding_inputs() {
            assert_ulp_eq(
                fastmaths::floor(x),
                floor_reference(x),
                MAX_ULP_TOL,
                &format!("floor({x})"),
            );
            assert_ulp_eq(
                fastmaths::ceil(x),
                ceil_reference(x),
                MAX_ULP_TOL,
                &format!("ceil({x})"),
            );
            assert_ulp_eq(
                fastmaths::trunc(x),
                trunc_reference(x),
                MAX_ULP_TOL,
                &format!("trunc({x})"),
            );
            assert_ulp_eq(
                fastmaths::round(x),
                round_reference(x),
                MAX_ULP_TOL,
                &format!("round({x})"),
            );
            assert_ulp_eq(
                fastmaths::rint(x),
                rint_reference(x),
                MAX_ULP_TOL,
                &format!("rint({x})"),
            );
            assert_ulp_eq(
                fastmaths::nearbyint(x),
                nearbyint_reference(x),
                MAX_ULP_TOL,
                &format!("nearbyint({x})"),
            );
        }
    }

    #[test]
    fn int_rounding_matches_reference() {
        for &x in &rounding_inputs() {
            assert_eq!(fastmaths::lrint(x), lrint_reference(x), "lrint({x})");
            assert_eq!(fastmaths::llrint(x), llrint_reference(x), "llrint({x})");
            assert_eq!(fastmaths::lround(x), lround_reference(x), "lround({x})");
            assert_eq!(fastmaths::llround(x), llround_reference(x), "llround({x})");
        }
    }

    #[test]
    fn copysign_fabs_special_cases() {
        assert_eq!(fastmaths::fabs(-0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(
            fastmaths::copysign(1.0, -0.0).to_bits(),
            (-1.0f64).to_bits()
        );
        assert!(fastmaths::fabs(f64::NAN).is_nan());
    }

    #[test]
    fn copysign_fabs_matches_reference() {
        let inputs = [
            (1.0, 1.0),
            (1.0, -1.0),
            (-1.0, 1.0),
            (-1.0, -1.0),
            (0.0, -1.0),
            (-0.0, 1.0),
            (1e-300, -1.0),
            (-1e-300, 1.0),
            (1e6, -1e6),
        ];
        for &(x, y) in &inputs {
            let actual = fastmaths::copysign(x, y);
            let expected = copysign_reference(x, y);
            assert_eq!(
                actual.to_bits(),
                expected.to_bits(),
                "copysign({x}, {y}) expected {expected}, got {actual}"
            );
        }
        for &(x, _) in &inputs {
            let actual = fastmaths::fabs(x);
            let expected = fabs_reference(x);
            if expected.is_nan() {
                assert!(actual.is_nan(), "fabs({x}) expected NaN");
            } else {
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "fabs({x}) expected {expected}, got {actual}"
                );
            }
        }
    }

    #[test]
    fn fma_matches_reference_ulps() {
        for &(x, y, z) in &fma_inputs() {
            let actual = fastmaths::fma(x, y, z);
            let expected = fma_reference(x, y, z);
            if expected.is_nan() {
                assert!(actual.is_nan(), "fma({x}, {y}, {z}) expected NaN");
            } else {
                assert_ulp_eq(
                    actual,
                    expected,
                    MAX_ULP_TOL,
                    &format!("fma({x}, {y}, {z})"),
                );
            }
        }
    }

    #[test]
    fn scaling_special_cases() {
        let (m, e) = fastmaths::frexp(0.0);
        assert_eq!(m.to_bits(), 0.0f64.to_bits());
        assert_eq!(e, 0);
        let (m, e) = fastmaths::frexp(f64::INFINITY);
        assert_eq!(m, f64::INFINITY);
        assert_eq!(e, 0);
        assert_eq!(fastmaths::ldexp(f64::INFINITY, 5), f64::INFINITY);
        assert_eq!(fastmaths::scalbn(f64::INFINITY, -5), f64::INFINITY);
    }

    #[test]
    fn scaling_matches_reference_ulps() {
        for &x in &scaling_inputs() {
            let (m_a, e_a) = fastmaths::frexp(x);
            let (m_e, e_e) = frexp_reference(x);
            assert_ulp_eq(m_a, m_e, MAX_ULP_TOL, &format!("frexp({x}) mantissa"));
            assert_eq!(e_a, e_e, "frexp({x}) exponent");
        }
        for &(x, n) in &scalbn_inputs() {
            let actual = fastmaths::scalbn(x, n);
            let expected = scalbn_reference(x, n);
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &format!("scalbn({x},{n})"));
            let actual = fastmaths::ldexp(x, n);
            let expected = ldexp_reference(x, n);
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &format!("ldexp({x},{n})"));
        }
        for &(x, n) in &scalbln_inputs() {
            let actual = fastmaths::scalbln(x, n);
            let expected = scalbln_reference(x, n);
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &format!("scalbln({x},{n})"));
        }
    }

    #[test]
    fn remquo_matches_reference_ulps() {
        for &(x, y) in &remquo_inputs() {
            let (actual_r, actual_q) = fastmaths::remquo(x, y);
            let (expected_r, expected_q) = remquo_reference(x, y);
            if expected_r.is_nan() {
                assert!(actual_r.is_nan(), "remquo({x}, {y}) expected NaN");
            } else {
                assert_ulp_eq(
                    actual_r,
                    expected_r,
                    MAX_ULP_TOL,
                    &format!("remquo({x}, {y})"),
                );
            }
            assert_eq!(actual_q, expected_q, "remquo({x}, {y}) quotient");
        }
    }

    #[test]
    fn exp2_special_cases() {
        assert!(fastmaths::exp2(f64::NAN).is_nan());
        assert_eq!(fastmaths::exp2(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::exp2(f64::NEG_INFINITY), 0.0);
        assert_eq!(fastmaths::exp2(0.0).to_bits(), 1.0f64.to_bits());
        assert_eq!(fastmaths::exp2(-0.0).to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn exp2_matches_reference_ulps() {
        for &x in &exp2_inputs() {
            let actual = fastmaths::exp2(x);
            let expected = exp2_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("exp2({x})"));
        }
    }

    #[test]
    fn expm1_special_cases() {
        assert!(fastmaths::expm1(f64::NAN).is_nan());
        assert_eq!(fastmaths::expm1(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::expm1(f64::NEG_INFINITY), -1.0);
        assert_eq!(fastmaths::expm1(0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(fastmaths::expm1(-0.0).to_bits(), (-0.0f64).to_bits());
    }

    #[test]
    fn expm1_matches_reference_ulps() {
        for &x in &expm1_inputs() {
            let actual = fastmaths::expm1(x);
            let expected = expm1_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("expm1({x})"));
        }
    }

    #[test]
    fn expm1_regression_for_sinh_case_seed() {
        let x = 1.7335082797141748_f64;
        let actual = fastmaths::expm1(x);
        let expected = expm1_reference(x);
        assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("expm1({x})"));
    }

    #[test]
    fn log2_log10_special_cases() {
        assert!(fastmaths::log2(f64::NAN).is_nan());
        assert!(fastmaths::log10(f64::NAN).is_nan());
        assert_eq!(fastmaths::log2(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::log10(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::log2(0.0), f64::NEG_INFINITY);
        assert_eq!(fastmaths::log10(0.0), f64::NEG_INFINITY);
        assert!(fastmaths::log2(-1.0).is_nan());
        assert!(fastmaths::log10(-1.0).is_nan());
    }

    #[test]
    fn log2_log10_matches_reference_ulps() {
        for &x in &ln_inputs() {
            let actual = fastmaths::log2(x);
            let expected = log2_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("log2({x})"));

            let actual = fastmaths::log10(x);
            let expected = log10_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("log10({x})"));
        }
    }

    #[test]
    fn log10_regression_case_seed() {
        let x = 1.017_091_338_276_825_4_f64;
        let actual = fastmaths::log10(x);
        let expected = log10_reference(x);
        assert_ulp_eq(
            actual,
            expected,
            DERIVED_ULP_TOL,
            &format_case("log10", x, "seed"),
        );
    }

    #[test]
    fn tan_special_cases() {
        assert!(fastmaths::tan(f64::NAN).is_nan());
        assert!(fastmaths::tan(f64::INFINITY).is_nan());
        assert!(fastmaths::tan(f64::NEG_INFINITY).is_nan());
    }

    #[test]
    fn tan_matches_reference_ulps() {
        for &x in &tan_inputs() {
            let actual = fastmaths::tan(x);
            let expected = tan_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("tan({x})"));
        }
    }

    #[test]
    fn asin_acos_special_cases() {
        assert!(fastmaths::asin(f64::NAN).is_nan());
        assert!(fastmaths::acos(f64::NAN).is_nan());
        assert_eq!(fastmaths::asin(1.0), FRAC_PI_2);
        assert_eq!(fastmaths::asin(-1.0), -FRAC_PI_2);
        assert_eq!(fastmaths::acos(1.0), 0.0);
        assert_eq!(fastmaths::acos(-1.0), PI);
        assert!(fastmaths::asin(1.1).is_nan());
        assert!(fastmaths::acos(-1.1).is_nan());
    }

    #[test]
    fn asin_acos_matches_reference_ulps() {
        for &x in &asin_inputs() {
            let actual = fastmaths::asin(x);
            let expected = asin_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("asin({x})"));

            let actual = fastmaths::acos(x);
            let expected = acos_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("acos({x})"));
        }
    }

    #[test]
    fn atan_matches_reference_ulps() {
        for &x in &atan_inputs() {
            let actual = fastmaths::atan(x);
            let expected = atan_reference(x);
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &format!("atan({x})"));
        }
    }

    #[test]
    fn atan2_matches_reference_ulps() {
        for &(y, x) in &atan2_inputs() {
            let actual = fastmaths::atan2(y, x);
            let expected = atan2_reference(y, x);
            assert_ulp_eq(actual, expected, MAX_ULP_TOL, &format!("atan2({y},{x})"));
        }
    }

    #[test]
    fn sinh_cosh_tanh_special_cases() {
        assert!(fastmaths::sinh(f64::NAN).is_nan());
        assert!(fastmaths::cosh(f64::NAN).is_nan());
        assert!(fastmaths::tanh(f64::NAN).is_nan());

        assert_eq!(fastmaths::sinh(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::sinh(f64::NEG_INFINITY), f64::NEG_INFINITY);
        assert_eq!(fastmaths::cosh(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::cosh(f64::NEG_INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::tanh(f64::INFINITY), 1.0);
        assert_eq!(fastmaths::tanh(f64::NEG_INFINITY), -1.0);
    }

    #[test]
    fn sinh_cosh_tanh_matches_reference_ulps() {
        for &x in &sinh_inputs() {
            let actual = fastmaths::sinh(x);
            let expected = sinh_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                DERIVED_ULP_TOL,
                &format_case("sinh", x, "unit"),
            );
        }
        for &x in &cosh_inputs() {
            let actual = fastmaths::cosh(x);
            let expected = cosh_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("cosh({x})"));
        }
        for &x in &tanh_inputs() {
            let actual = fastmaths::tanh(x);
            let expected = tanh_reference(x);
            assert_ulp_eq(actual, expected, TANH_ULP_TOL, &format!("tanh({x})"));
        }
    }

    #[test]
    fn sinh_regression_case_seed() {
        let x = 0.8697034726629151_f64;
        let actual = fastmaths::sinh(x);
        let expected = sinh_reference(x);
        assert_ulp_eq(
            actual,
            expected,
            DERIVED_ULP_TOL,
            &format_case("sinh", x, "seed"),
        );
    }

    #[test]
    fn sinh_regression_near_overflow_threshold() {
        // This exercises the high-|x| sinh path where exp(|x|) would overflow
        // but sinh(x) is still finite.
        let x = -710.475_860_073_934_f64;
        let actual = fastmaths::sinh(x);
        let expected = sinh_reference(x);
        assert_ulp_eq(
            actual,
            expected,
            DERIVED_ULP_TOL,
            &format_case("sinh", x, "regression"),
        );
    }

    #[test]
    fn sinh_regression_rounding_corner_cases() {
        // These inputs were found by MPFR proptest as 2-ULP misses at various
        // points during development. Keep them as explicit unit tests to catch
        // future regressions in the |x|<1 rounding logic.
        for &x in &[
            0.06177710537659864_f64,
            0.7699349981625943_f64,
            0.7782559419902234_f64,
            0.815019628039047_f64,
            0.8640551460732968_f64,
        ] {
            let actual = fastmaths::sinh(x);
            let expected = sinh_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                DERIVED_ULP_TOL,
                &format_case("sinh", x, "regression"),
            );

            // Odd symmetry: also validate the negative side explicitly.
            let xn = -x;
            let actual = fastmaths::sinh(xn);
            let expected = sinh_reference(xn);
            assert_ulp_eq(
                actual,
                expected,
                DERIVED_ULP_TOL,
                &format_case("sinh", xn, "regression"),
            );
        }
    }

    #[test]
    fn sinh_regression_medium_range() {
        // MPFR proptest regressions in the 1<=|x|<22 path.
        for &x in &[1.408992436517082_f64, 2.0854509241613868_f64] {
            let actual = fastmaths::sinh(x);
            let expected = sinh_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                DERIVED_ULP_TOL,
                &format_case("sinh", x, "regression"),
            );

            let xn = -x;
            let actual = fastmaths::sinh(xn);
            let expected = sinh_reference(xn);
            assert_ulp_eq(
                actual,
                expected,
                DERIVED_ULP_TOL,
                &format_case("sinh", xn, "regression"),
            );
        }
    }

    #[test]
    fn asinh_acosh_atanh_special_cases() {
        assert!(fastmaths::asinh(f64::NAN).is_nan());
        assert!(fastmaths::acosh(f64::NAN).is_nan());
        assert!(fastmaths::atanh(f64::NAN).is_nan());
        assert_eq!(fastmaths::asinh(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::asinh(f64::NEG_INFINITY), f64::NEG_INFINITY);
        assert_eq!(fastmaths::acosh(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::atanh(0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(fastmaths::atanh(-0.0).to_bits(), (-0.0f64).to_bits());
        assert!(fastmaths::acosh(0.5).is_nan());
        assert!(fastmaths::atanh(1.0).is_infinite());
        assert!(fastmaths::atanh(-1.0).is_infinite());
    }

    #[test]
    fn atanh_regression_case_seed() {
        let x = -0.4789704365236613_f64;
        let actual = fastmaths::atanh(x);
        let expected = atanh_reference(x);
        assert_ulp_eq(actual, expected, ATANH_ULP_TOL, &format!("atanh({x})"));
    }

    #[test]
    fn acosh_regression_case_seed() {
        let x = 1.1226630563945177_f64;
        let actual = fastmaths::acosh(x);
        let expected = acosh_reference(x);
        assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("acosh({x})"));
    }

    #[test]
    fn acosh_sinh_regression_cases() {
        let x = 1.8399999999999999_f64;
        let actual = fastmaths::acosh(x);
        let expected = acosh_reference(x);
        assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("acosh({x})"));

        let x = -0.7398078960390417_f64;
        let actual = fastmaths::sinh(x);
        let expected = sinh_reference(x);
        assert_ulp_eq(
            actual,
            expected,
            DERIVED_ULP_TOL,
            &format_case("sinh", x, "regression"),
        );

        let x = -0.8648745336167547_f64;
        let actual = fastmaths::sinh(x);
        let expected = sinh_reference(x);
        assert_ulp_eq(
            actual,
            expected,
            DERIVED_ULP_TOL,
            &format_case("sinh", x, "regression"),
        );

        let x = 0.8392647357681798_f64;
        let actual = fastmaths::sinh(x);
        let expected = sinh_reference(x);
        assert_ulp_eq(
            actual,
            expected,
            DERIVED_ULP_TOL,
            &format_case("sinh", x, "regression"),
        );

        let x = -0.723199985761056_f64;
        let actual = fastmaths::sinh(x);
        let expected = sinh_reference(x);
        assert_ulp_eq(
            actual,
            expected,
            DERIVED_ULP_TOL,
            &format_case("sinh", x, "regression"),
        );
    }

    #[test]
    fn asinh_acosh_atanh_matches_reference_ulps() {
        for &x in &asinh_inputs() {
            let actual = fastmaths::asinh(x);
            let expected = asinh_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("asinh({x})"));
        }
        for &x in &acosh_inputs() {
            let actual = fastmaths::acosh(x);
            let expected = acosh_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("acosh({x})"));
        }
        for &x in &atanh_inputs() {
            let actual = fastmaths::atanh(x);
            let expected = atanh_reference(x);
            assert_ulp_eq(actual, expected, ATANH_ULP_TOL, &format!("atanh({x})"));
        }
    }

    #[test]
    fn erf_erfc_special_cases() {
        assert!(fastmaths::erf(f64::NAN).is_nan());
        assert!(fastmaths::erfc(f64::NAN).is_nan());
        assert_eq!(fastmaths::erf(f64::INFINITY), 1.0);
        assert_eq!(fastmaths::erf(f64::NEG_INFINITY), -1.0);
        assert_eq!(fastmaths::erfc(f64::INFINITY), 0.0);
        assert_eq!(fastmaths::erfc(f64::NEG_INFINITY), 2.0);
    }

    #[test]
    fn erf_erfc_matches_reference_ulps() {
        #[cfg(feature = "mpfr")]
        {
            for &x in &erf_inputs() {
                let actual = fastmaths::erf(x);
                let expected = erf_reference(x);
                assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("erf({x})"));
            }
            for &x in &erfc_inputs() {
                let actual = fastmaths::erfc(x);
                let expected = erfc_reference(x);
                assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("erfc({x})"));
            }
        }
    }

    #[test]
    fn exp10_special_cases() {
        assert!(fastmaths::exp10(f64::NAN).is_nan());
        assert_eq!(fastmaths::exp10(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::exp10(f64::NEG_INFINITY), 0.0);
        assert_eq!(fastmaths::exp10(0.0).to_bits(), 1.0f64.to_bits());
        assert_eq!(fastmaths::exp10(-0.0).to_bits(), 1.0f64.to_bits());
    }

    #[test]
    fn exp10_matches_reference_ulps() {
        for &x in &exp10_inputs() {
            let actual = fastmaths::exp10(x);
            let expected = exp10_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("exp10({x})"));
        }
    }

    #[test]
    fn lgamma_tgamma_special_cases() {
        assert!(fastmaths::lgamma(f64::NAN).is_nan());
        assert_eq!(fastmaths::lgamma(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::lgamma(0.0), f64::INFINITY);
        assert_eq!(fastmaths::lgamma(-0.0), f64::INFINITY);
        assert_eq!(fastmaths::lgamma(-1.0), f64::INFINITY);
        assert_eq!(fastmaths::lgamma(-2.0), f64::INFINITY);

        assert!(fastmaths::tgamma(f64::NAN).is_nan());
        assert_eq!(fastmaths::tgamma(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::tgamma(1.0), 1.0);
        assert_eq!(fastmaths::tgamma(2.0), 1.0);
        assert_eq!(fastmaths::tgamma(0.5), core::f64::consts::PI.sqrt());
        assert_eq!(fastmaths::tgamma(0.0), f64::INFINITY);
        assert_eq!(fastmaths::tgamma(-0.0), f64::NEG_INFINITY);
        assert!(fastmaths::tgamma(-1.0).is_nan());
        assert!(fastmaths::tgamma(-2.0).is_nan());
    }

    #[test]
    fn lgamma_tgamma_matches_reference_ulps() {
        #[cfg(feature = "mpfr")]
        {
            for &x in &lgamma_inputs() {
                let actual = fastmaths::lgamma(x);
                let expected = lgamma_reference(x);
                assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("lgamma({x})"));
            }
            for &x in &tgamma_inputs() {
                let actual = fastmaths::tgamma(x);
                let expected = tgamma_reference(x);
                assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("tgamma({x})"));
            }
        }
    }

    #[test]
    fn logb_ilogb_special_cases() {
        assert!(fastmaths::logb(f64::NAN).is_nan());
        assert_eq!(fastmaths::logb(0.0), f64::NEG_INFINITY);
        assert_eq!(fastmaths::logb(f64::INFINITY), f64::INFINITY);
        assert_eq!(fastmaths::ilogb(0.0), i32::MIN);
        assert_eq!(fastmaths::ilogb(f64::INFINITY), i32::MAX);
        assert_eq!(fastmaths::ilogb(f64::NAN), i32::MAX);
    }

    #[test]
    fn logb_ilogb_matches_reference_ulps() {
        for &x in &logb_inputs() {
            let actual = fastmaths::logb(x);
            let expected = logb_reference(x);
            if expected.is_nan() {
                assert!(actual.is_nan(), "logb({x}) expected NaN, got {actual}");
            } else {
                assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("logb({x})"));
            }
        }
        for &x in &ilogb_inputs() {
            let actual = fastmaths::ilogb(x);
            let expected = ilogb_reference(x);
            assert_eq!(
                actual, expected,
                "ilogb({x}) expected {expected}, got {actual}"
            );
        }
    }

    #[test]
    fn modf_special_cases() {
        let (frac, int) = fastmaths::modf(f64::INFINITY);
        assert_eq!(int, f64::INFINITY);
        assert_eq!(frac.to_bits(), 0.0f64.to_bits());
        let (frac, int) = fastmaths::modf(f64::NEG_INFINITY);
        assert_eq!(int, f64::NEG_INFINITY);
        assert_eq!(frac.to_bits(), (-0.0f64).to_bits());
        let (frac, int) = fastmaths::modf(f64::NAN);
        assert!(frac.is_nan());
        assert!(int.is_nan());
    }

    #[test]
    fn modf_matches_reference_ulps() {
        for &x in &modf_inputs() {
            let (frac, int) = fastmaths::modf(x);
            let (frac_e, int_e) = modf_reference(x);
            assert_ulp_eq(frac, frac_e, DERIVED_ULP_TOL, &format!("modf frac({x})"));
            assert_ulp_eq(int, int_e, DERIVED_ULP_TOL, &format!("modf int({x})"));
        }
    }

    #[test]
    fn fdim_fmax_fmin_special_cases() {
        assert!(fastmaths::fdim(f64::NAN, 1.0).is_nan());
        assert_eq!(fastmaths::fmax(0.0, -0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(fastmaths::fmin(0.0, -0.0).to_bits(), (-0.0f64).to_bits());
    }

    #[test]
    fn fdim_fmax_fmin_matches_reference_ulps() {
        for &(x, y) in &fdim_inputs() {
            let actual = fastmaths::fdim(x, y);
            let expected = fdim_reference(x, y);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("fdim({x},{y})"));
        }
        for &(x, y) in &fmax_inputs() {
            let actual = fastmaths::fmax(x, y);
            let expected = fmax_reference(x, y);
            if actual == 0.0 && expected == 0.0 {
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "fmax({x},{y}) sign mismatch"
                );
            } else {
                assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("fmax({x},{y})"));
            }
        }
        for &(x, y) in &fmin_inputs() {
            let actual = fastmaths::fmin(x, y);
            let expected = fmin_reference(x, y);
            if actual == 0.0 && expected == 0.0 {
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "fmin({x},{y}) sign mismatch"
                );
            } else {
                assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("fmin({x},{y})"));
            }
        }
    }

    #[test]
    fn nextafter_special_cases() {
        assert!(fastmaths::nextafter(f64::NAN, 1.0).is_nan());
        assert!(fastmaths::nextafter(1.0, f64::NAN).is_nan());
        assert_eq!(fastmaths::nextafter(1.0, 1.0), 1.0);
    }

    #[test]
    fn nextafter_matches_reference_ulps() {
        for &(x, y) in &nextafter_inputs() {
            let actual = fastmaths::nextafter(x, y);
            let expected = nextafter_reference(x, y);
            if actual.is_nan() {
                assert!(expected.is_nan(), "nextafter({x},{y}) expected NaN");
            } else {
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "nextafter({x},{y}) expected {expected:?}, got {actual:?}"
                );
            }
        }
    }

    #[test]
    fn hypot_matches_reference_ulps() {
        for &(x, y) in &hypot_inputs() {
            let actual = fastmaths::hypot(x, y);
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
        assert!(fastmaths::fmod(f64::NAN, 1.0).is_nan());
        assert!(fastmaths::fmod(1.0, f64::NAN).is_nan());
        assert!(fastmaths::fmod(f64::INFINITY, 1.0).is_nan());
        assert!(fastmaths::fmod(1.0, 0.0).is_nan());
        assert_eq!(fastmaths::fmod(0.0, 1.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(fastmaths::fmod(-0.0, 1.0).to_bits(), (-0.0f64).to_bits());
    }

    #[test]
    fn fmod_matches_reference_ulps() {
        for &(x, y) in &fmod_inputs() {
            let actual = fastmaths::fmod(x, y);
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
        assert!(fastmaths::remainder(f64::NAN, 1.0).is_nan());
        assert!(fastmaths::remainder(1.0, f64::NAN).is_nan());
        assert!(fastmaths::remainder(f64::INFINITY, 1.0).is_nan());
        assert!(fastmaths::remainder(1.0, 0.0).is_nan());
        assert_eq!(fastmaths::remainder(0.0, 1.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(
            fastmaths::remainder(-0.0, 1.0).to_bits(),
            (-0.0f64).to_bits()
        );
        assert_eq!(fastmaths::remainder(1.0, f64::INFINITY), 1.0);
    }

    #[test]
    fn remainder_matches_reference_ulps() {
        for &(x, y) in &remainder_inputs() {
            let actual = fastmaths::remainder(x, y);
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
            let actual = fastmaths::pow(x, y);
            let expected = pow_reference(x, y);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("pow({x},{y})"));
        }
    }

    #[test]
    fn pow_regression_negative_tiny_base_negative_odd_int_exp_overflow() {
        let x = -5.048_709_793_414_476e-29_f64;
        let y = -11.0_f64;
        let actual = fastmaths::pow(x, y);
        let expected = pow_reference(x, y);
        assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("pow({x},{y})"));
    }

    #[test]
    fn sqrt_matches_reference_ulps() {
        for &x in &sqrt_inputs() {
            let actual = fastmaths::sqrt(x);
            let expected = sqrt_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("sqrt({x})"));
        }
    }

    #[test]
    fn cbrt_matches_reference_ulps() {
        for &x in &cbrt_inputs() {
            let actual = fastmaths::cbrt(x);
            let expected = cbrt_reference(x);
            assert_ulp_eq(actual, expected, DERIVED_ULP_TOL, &format!("cbrt({x})"));
        }
    }

    #[test]
    fn sin_cos_special_cases() {
        let nan = f64::NAN;
        let pos_inf = f64::INFINITY;
        let neg_inf = f64::NEG_INFINITY;

        assert!(fastmaths::sin(nan).is_nan());
        assert!(fastmaths::cos(nan).is_nan());
        assert!(fastmaths::sin(pos_inf).is_nan());
        assert!(fastmaths::cos(pos_inf).is_nan());
        assert!(fastmaths::sin(neg_inf).is_nan());
        assert!(fastmaths::cos(neg_inf).is_nan());

        let neg_zero = -0.0f64;
        assert_eq!(fastmaths::sin(neg_zero).to_bits(), neg_zero.to_bits());
        assert_eq!(fastmaths::cos(neg_zero).to_bits(), 1.0f64.to_bits());
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
            let sin_actual = fastmaths::sin(x);
            let cos_actual = fastmaths::cos(x);
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
            let sin_actual = fastmaths::sin(x);
            let cos_actual = fastmaths::cos(x);
            assert_ulp_eq(sin_actual, sin_expected, MAX_ULP_TOL, &format!("sin({x})"));
            assert_ulp_eq(cos_actual, cos_expected, MAX_ULP_TOL, &format!("cos({x})"));
        }
    }

    #[test]
    fn sincos_matches_std_ulps() {
        let inputs = trig_inputs();

        for &x in &inputs {
            let (sin_actual, cos_actual) = fastmaths::sincos(x);
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
            let sin_pos = fastmaths::sin(x);
            let sin_neg = fastmaths::sin(-x);
            let cos_pos = fastmaths::cos(x);
            let cos_neg = fastmaths::cos(-x);

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
            let actual = fastmaths::exp(x);
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
            let actual = fastmaths::ln(x);
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
            let sin_actual = fastmaths::sin(x);
            let cos_actual = fastmaths::cos(x);
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
    fn compare_glibc_fastmaths() {
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

            std::println!("| Input | Func | glibc (bits) | fastmaths (bits) | ULP Delta |");
            std::println!("| :--- | :--- | :--- | :--- | :--- |");
            for &x in &test_inputs {
                type CFn = unsafe extern "C" fn(f64) -> f64;
                type RustFn = fn(f64) -> f64;
                type FnSpec = (&'static str, CFn, RustFn);

                let fns: [FnSpec; 4] = [
                    ("exp", *g_exp, fastmaths::exp),
                    ("log", *g_log, fastmaths::ln),
                    ("sin", *g_sin, fastmaths::sin),
                    ("cos", *g_cos, fastmaths::cos),
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
        let lgamma: libloading::Symbol<CFn> = unsafe { lib.get(b"lgamma").unwrap() };
        let tgamma: libloading::Symbol<CFn> = unsafe { lib.get(b"tgamma").unwrap() };
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
            let f = fastmaths::exp(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist exp({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -100.0, 100.0);
            let g = unsafe { exp2(x) };
            let f = fastmaths::exp2(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist exp2({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1.0, 1.0);
            let g = unsafe { expm1(x) };
            let f = fastmaths::expm1(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist expm1({x})"));
        }

        for _ in 0..samples {
            let x = rand_f64_pos(&mut state);
            let g = unsafe { log(x) };
            let f = fastmaths::ln(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist ln({x})"));
        }

        for _ in 0..samples {
            let x = rand_f64_pos(&mut state);
            let g = unsafe { log2(x) };
            let f = fastmaths::log2(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist log2({x})"));
        }

        for _ in 0..samples {
            let x = rand_f64_pos(&mut state);
            let g = unsafe { log10(x) };
            let f = fastmaths::log10(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist log10({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -0.9, 1e6);
            let g = unsafe { log1p(x) };
            let f = fastmaths::log1p(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist log1p({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            if x <= 0.0 && x == x.trunc() {
                continue;
            }
            let g = unsafe { lgamma(x) };
            let f = fastmaths::lgamma(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist lgamma({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            if x <= 0.0 && x == x.trunc() {
                continue;
            }
            let g = unsafe { tgamma(x) };
            let f = fastmaths::tgamma(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist tgamma({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let g = unsafe { sin(x) };
            let f = fastmaths::sin(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist sin({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let g = unsafe { cos(x) };
            let f = fastmaths::cos(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist cos({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let g = unsafe { tan(x) };
            let f = fastmaths::tan(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist tan({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1.0, 1.0);
            let g = unsafe { asin(x) };
            let f = fastmaths::asin(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist asin({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1.0, 1.0);
            let g = unsafe { acos(x) };
            let f = fastmaths::acos(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist acos({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let g = unsafe { atan(x) };
            let f = fastmaths::atan(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist atan({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            let g = unsafe { sinh(x) };
            let f = fastmaths::sinh(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist sinh({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            let g = unsafe { cosh(x) };
            let f = fastmaths::cosh(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist cosh({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            let g = unsafe { tanh(x) };
            let f = fastmaths::tanh(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist tanh({x})"));
        }

        for _ in 0..samples {
            let y = rand_range(&mut state, -1e6, 1e6);
            let x = rand_range(&mut state, -1e6, 1e6);
            if x == 0.0 && y == 0.0 {
                continue;
            }
            let g = unsafe { atan2(y, x) };
            let f = fastmaths::atan2(y, x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist atan2({y},{x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e200, 1e200);
            let y = rand_range(&mut state, -1e200, 1e200);
            let g = unsafe { hypot(x, y) };
            let f = fastmaths::hypot(x, y);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist hypot({x},{y})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let mut y = rand_range(&mut state, 1e-6, 1e6);
            if rand_u64(&mut state) & 1 == 0 {
                y = -y;
            }
            let g = unsafe { fmod(x, y) };
            let f = fastmaths::fmod(x, y);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist fmod({x},{y})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let mut y = rand_range(&mut state, 1e-6, 1e6);
            if rand_u64(&mut state) & 1 == 0 {
                y = -y;
            }
            let g = unsafe { remainder(x, y) };
            let f = fastmaths::remainder(x, y);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist remainder({x},{y})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, 0.1, 10.0);
            let y = rand_range(&mut state, -10.0, 10.0);
            let g = unsafe { pow(x, y) };
            let f = fastmaths::pow(x, y);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist pow({x},{y})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, 0.0, 1e300);
            let g = unsafe { sqrt(x) };
            let f = fastmaths::sqrt(x);
            assert_ulp_eq_glibc(f, g, 1.0, &format!("glibc dist sqrt({x})"));
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e300, 1e300);
            let g = unsafe { cbrt(x) };
            let f = fastmaths::cbrt(x);
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
        let lgamma: libloading::Symbol<CFn> = unsafe { lib.get(b"lgamma").unwrap() };
        let tgamma: libloading::Symbol<CFn> = unsafe { lib.get(b"tgamma").unwrap() };
        let atan: libloading::Symbol<CFn> = unsafe { lib.get(b"atan").unwrap() };
        let atan2: libloading::Symbol<CFn2> = unsafe { lib.get(b"atan2").unwrap() };
        let asin: libloading::Symbol<CFn> = unsafe { lib.get(b"asin").unwrap() };
        let acos: libloading::Symbol<CFn> = unsafe { lib.get(b"acos").unwrap() };
        let sinh: libloading::Symbol<CFn> = unsafe { lib.get(b"sinh").unwrap() };
        let cosh: libloading::Symbol<CFn> = unsafe { lib.get(b"cosh").unwrap() };
        let tanh: libloading::Symbol<CFn> = unsafe { lib.get(b"tanh").unwrap() };
        let asinh: libloading::Symbol<CFn> = unsafe { lib.get(b"asinh").unwrap() };
        let acosh: libloading::Symbol<CFn> = unsafe { lib.get(b"acosh").unwrap() };
        let atanh: libloading::Symbol<CFn> = unsafe { lib.get(b"atanh").unwrap() };
        let erf: libloading::Symbol<CFn> = unsafe { lib.get(b"erf").unwrap() };
        let erfc: libloading::Symbol<CFn> = unsafe { lib.get(b"erfc").unwrap() };
        let exp10: libloading::Symbol<CFn> = unsafe { lib.get(b"exp10").unwrap() };
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
            let f = fastmaths::exp(x);
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
            let f = fastmaths::ln(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr ln ulp {ulp_f} > 1 at {x}");
            assert!(ulp_f <= ulp_g, "fast ln ulp {ulp_f} > glibc {ulp_g} at {x}");
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e6, 1e6);
            let mp = mpfr_sin_f64(x);
            let g = unsafe { sin(x) };
            let f = fastmaths::sin(x);
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
            let f = fastmaths::cos(x);
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
            let f = fastmaths::tan(x);
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
            let f = fastmaths::exp2(x);
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
            let f = fastmaths::expm1(x);
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
            let f = fastmaths::log2(x);
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
            let f = fastmaths::log10(x);
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
            let f = fastmaths::log1p(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr log1p ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast log1p ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            if x <= 0.0 && x == x.trunc() {
                continue;
            }
            let mp = mpfr_lgamma_f64(x);
            let g = unsafe { lgamma(x) };
            let f = fastmaths::lgamma(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr lgamma ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast lgamma ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            if x <= 0.0 && x == x.trunc() {
                continue;
            }
            let mp = mpfr_tgamma_f64(x);
            let g = unsafe { tgamma(x) };
            let f = fastmaths::tgamma(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr tgamma ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast tgamma ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, 0.1, 10.0);
            let y = rand_range(&mut state, -10.0, 10.0);
            let mp = mpfr_pow_f64(x, y);
            let g = unsafe { pow(x, y) };
            let f = fastmaths::pow(x, y);
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
            let f = fastmaths::atan(x);
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
            let f = fastmaths::asin(x);
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
            let f = fastmaths::acos(x);
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
            let f = fastmaths::atan2(y, x);
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
            let f = fastmaths::sinh(x);
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
            let f = fastmaths::cosh(x);
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
            let f = fastmaths::tanh(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr tanh ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast tanh ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -20.0, 20.0);
            let mp = mpfr_asinh_f64(x);
            let g = unsafe { asinh(x) };
            let f = fastmaths::asinh(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr asinh ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast asinh ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, 1.0, 1e6);
            let mp = mpfr_acosh_f64(x);
            let g = unsafe { acosh(x) };
            let f = fastmaths::acosh(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr acosh ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast acosh ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -0.99, 0.99);
            let mp = mpfr_atanh_f64(x);
            let g = unsafe { atanh(x) };
            let f = fastmaths::atanh(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr atanh ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast atanh ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -3.0, 3.0);
            let mp = mpfr_erf_f64(x);
            let g = unsafe { erf(x) };
            let f = fastmaths::erf(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr erf ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast erf ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -3.0, 3.0);
            let mp = mpfr_erfc_f64(x);
            let g = unsafe { erfc(x) };
            let f = fastmaths::erfc(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr erfc ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast erfc ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -50.0, 50.0);
            let mp = mpfr_exp10_f64(x);
            let g = unsafe { exp10(x) };
            let f = fastmaths::exp10(x);
            let ulp_f = ulp_error(f, mp);
            let ulp_g = ulp_error(g, mp);
            assert!(ulp_f <= 1.0, "mpfr exp10 ulp {ulp_f} > 1 at {x}");
            assert!(
                ulp_f <= ulp_g,
                "fast exp10 ulp {ulp_f} > glibc {ulp_g} at {x}"
            );
        }

        for _ in 0..samples {
            let x = rand_range(&mut state, -1e200, 1e200);
            let y = rand_range(&mut state, -1e200, 1e200);
            let mp = mpfr_hypot_f64(x, y);
            let g = unsafe { hypot(x, y) };
            let f = fastmaths::hypot(x, y);
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
            let f = fastmaths::fmod(x, y);
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
            let f = fastmaths::remainder(x, y);
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
            let f = fastmaths::sqrt(x);
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
            let f = fastmaths::cbrt(x);
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
    use proptest::strategy::BoxedStrategy;

    const F64_EXP_BIAS: i32 = 1023;
    const F64_EXP_MIN: i32 = -1022;
    const F64_EXP_MAX: i32 = 1023;
    const F64_MANTISSA_MASK: u64 = (1u64 << 52) - 1;
    const NEAR_ONE_MAX_POW: u32 = 52;
    const TINY_MAX_POW: u32 = 1074;
    const EXP_OVERFLOW: f64 = 709.782_712_893_384;
    const EXP_UNDERFLOW_TO_ZERO: f64 = -745.133_219_101_941_1;
    const SINH_OVERFLOW: f64 = 710.475_860_073_943_9;
    const LN_TABLE_BITS: u32 = 7;
    const LN_TABLE_N: u64 = 1u64 << LN_TABLE_BITS;
    const LN_TABLE_OFF: u64 = 0x3fe6_0000_0000_0000u64;
    const LN_NEAR_ONE_LO: u64 = 0x3fee_0000_0000_0000u64;
    const LN_NEAR_ONE_HI: u64 = 0x3ff1_0900_0000_0000u64;
    const LOG10_NEAR1_BOUND: f64 = 0.4;
    const LN2_DIV_N: f64 = core::f64::consts::LN_2 / 128.0;
    const HYPOT_LARGE_VAL: f64 = f64::from_bits(0x5fe0_0000_0000_0000); // 2^511
    const HYPOT_TINY_VAL: f64 = f64::from_bits(0x2340_0000_0000_0000); // 2^-459
    const HYPOT_EPS: f64 = f64::from_bits(0x3c90_0000_0000_0000); // 2^-54

    fn normal_f64_with_exp(min_exp: i32, max_exp: i32) -> BoxedStrategy<f64> {
        let min = min_exp.max(F64_EXP_MIN);
        let max = max_exp.min(F64_EXP_MAX);
        (any::<bool>(), min..=max, any::<u64>())
            .prop_map(|(neg, exp, mant)| {
                let sign = if neg { 1u64 << 63 } else { 0 };
                let exp_bits = ((exp + F64_EXP_BIAS) as u64) << 52;
                let mant_bits = mant & F64_MANTISSA_MASK;
                f64::from_bits(sign | exp_bits | mant_bits)
            })
            .boxed()
    }

    fn normal_pos_f64_with_exp(min_exp: i32, max_exp: i32) -> BoxedStrategy<f64> {
        let min = min_exp.max(F64_EXP_MIN);
        let max = max_exp.min(F64_EXP_MAX);
        (min..=max, any::<u64>())
            .prop_map(|(exp, mant)| {
                let exp_bits = ((exp + F64_EXP_BIAS) as u64) << 52;
                let mant_bits = mant & F64_MANTISSA_MASK;
                f64::from_bits(exp_bits | mant_bits)
            })
            .boxed()
    }

    fn subnormal_f64() -> BoxedStrategy<f64> {
        (any::<bool>(), 1u64..(1u64 << 52))
            .prop_map(|(neg, mant)| {
                let sign = if neg { 1u64 << 63 } else { 0 };
                f64::from_bits(sign | mant)
            })
            .boxed()
    }

    fn subnormal_pos_f64() -> BoxedStrategy<f64> {
        (1u64..(1u64 << 52)).prop_map(f64::from_bits).boxed()
    }

    fn tiny_positive() -> BoxedStrategy<f64> {
        (1u32..=TINY_MAX_POW)
            .prop_map(|k| 2.0f64.powi(-(k as i32)))
            .boxed()
    }

    fn tiny_signed() -> BoxedStrategy<f64> {
        (1u32..=TINY_MAX_POW, any::<bool>())
            .prop_map(|(k, neg)| {
                let x = 2.0f64.powi(-(k as i32));
                if neg { -x } else { x }
            })
            .boxed()
    }

    fn near_one_signed_open() -> BoxedStrategy<f64> {
        (1u32..=NEAR_ONE_MAX_POW, any::<bool>())
            .prop_map(|(k, neg)| {
                let delta = 2.0f64.powi(-(k as i32));
                let x = 1.0 - delta;
                if neg { -x } else { x }
            })
            .boxed()
    }

    fn near_one_above() -> BoxedStrategy<f64> {
        (1u32..=NEAR_ONE_MAX_POW)
            .prop_map(|k| 1.0 + 2.0f64.powi(-(k as i32)))
            .boxed()
    }

    fn near_one_both() -> BoxedStrategy<f64> {
        (1u32..=NEAR_ONE_MAX_POW, any::<bool>())
            .prop_map(|(k, up)| {
                let delta = 2.0f64.powi(-(k as i32));
                if up { 1.0 + delta } else { 1.0 - delta }
            })
            .boxed()
    }

    fn near_minus_one_open() -> BoxedStrategy<f64> {
        (1u32..=NEAR_ONE_MAX_POW)
            .prop_map(|k| -1.0 + 2.0f64.powi(-(k as i32)))
            .boxed()
    }

    fn ulp_steps(value: f64, max_steps: u32, toward_positive: bool) -> BoxedStrategy<f64> {
        (0u32..=max_steps)
            .prop_map(move |steps| {
                let mut v = value;
                for _ in 0..steps {
                    v = if toward_positive {
                        v.next_up()
                    } else {
                        v.next_down()
                    };
                }
                v
            })
            .boxed()
    }

    fn ulp_steps_exclusive(
        value: f64,
        max_steps: u32,
        toward_positive: bool,
    ) -> BoxedStrategy<f64> {
        (1u32..=max_steps)
            .prop_map(move |steps| {
                let mut v = value;
                for _ in 0..steps {
                    v = if toward_positive {
                        v.next_up()
                    } else {
                        v.next_down()
                    };
                }
                v
            })
            .boxed()
    }

    fn around(value: f64, max_steps: u32) -> BoxedStrategy<f64> {
        prop_oneof![
            1 => Just(value),
            1 => ulp_steps(value, max_steps, true),
            1 => ulp_steps(value, max_steps, false),
        ]
        .boxed()
    }

    fn around_signed(value: f64, max_steps: u32) -> BoxedStrategy<f64> {
        prop_oneof![
            1 => around(value, max_steps),
            1 => around(-value, max_steps),
        ]
        .boxed()
    }

    fn around_signed_below(value: f64, max_steps: u32) -> BoxedStrategy<f64> {
        prop_oneof![
            1 => ulp_steps_exclusive(value, max_steps, false),
            1 => ulp_steps_exclusive(-value, max_steps, true),
        ]
        .boxed()
    }

    fn clamp_below_sinh_overflow(x: f64) -> f64 {
        if !x.is_finite() || x.abs() >= SINH_OVERFLOW {
            let edge = SINH_OVERFLOW.next_down();
            return if x.is_sign_negative() { -edge } else { edge };
        }
        x
    }

    fn tiny_signed_below_tiny() -> BoxedStrategy<f64> {
        (29u32..=TINY_MAX_POW, any::<bool>())
            .prop_map(|(k, neg)| {
                let x = 2.0f64.powi(-(k as i32));
                if neg { -x } else { x }
            })
            .boxed()
    }

    fn tagged_f64<S>(label: &'static str, strat: S) -> BoxedStrategy<(f64, &'static str)>
    where
        S: Strategy<Value = f64> + 'static,
    {
        strat.prop_map(move |x| (x, label)).boxed()
    }

    fn tagged_pair<S>(label: &'static str, strat: S) -> BoxedStrategy<((f64, f64), &'static str)>
    where
        S: Strategy<Value = (f64, f64)> + 'static,
    {
        strat.prop_map(move |v| (v, label)).boxed()
    }

    fn exp_table_boundary_inputs() -> BoxedStrategy<f64> {
        let m_range = -150_000i32..=150_000i32;
        (m_range, any::<bool>())
            .prop_flat_map(|(m, up)| {
                let base = (m as f64) * LN2_DIV_N;
                ulp_steps(base, 32, up)
            })
            .boxed()
    }

    fn ln_table_boundary_inputs() -> BoxedStrategy<f64> {
        (0u32..(LN_TABLE_N as u32), -32i32..=32i32)
            .prop_map(|(i, offset)| {
                let base = LN_TABLE_OFF + ((i as u64) << (52 - LN_TABLE_BITS));
                let ix = if offset >= 0 {
                    base + (offset as u64)
                } else {
                    base - ((-offset) as u64)
                };
                f64::from_bits(ix)
            })
            .boxed()
    }

    fn ln_near_one_boundary_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            2 => around(f64::from_bits(LN_NEAR_ONE_LO), 256),
            2 => around(f64::from_bits(LN_NEAR_ONE_HI), 256),
            2 => around(1.0, 256),
        ]
        .boxed()
    }

    fn log10_near_one_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            2 => around(1.0 - LOG10_NEAR1_BOUND, 256),
            2 => around(1.0, 256),
            2 => around(1.0 + LOG10_NEAR1_BOUND, 256),
        ]
        .boxed()
    }

    fn log_subnormal_boundary_inputs() -> BoxedStrategy<f64> {
        let min_sub = f64::from_bits(1);
        let min_pos = f64::MIN_POSITIVE;
        prop_oneof![
            1 => Just(min_sub),
            1 => Just(min_sub.next_up()),
            1 => Just(min_pos),
            1 => Just(min_pos.next_down()),
        ]
        .boxed()
    }

    fn pow_near_integer_exponent() -> BoxedStrategy<f64> {
        (-1000i32..=1000i32, 1u32..=NEAR_ONE_MAX_POW, any::<bool>())
            .prop_map(|(n, k, up)| {
                let delta = 2.0f64.powi(-(k as i32));
                let base = n as f64;
                if up { base + delta } else { base - delta }
            })
            .boxed()
    }

    fn pow_exp_boundary_inputs() -> BoxedStrategy<(f64, f64)> {
        let delta = (20u32..=40u32).prop_map(|k| 2.0f64.powi(-(k as i32)));
        let target = prop_oneof![
            1 => (EXP_OVERFLOW * 0.99)..(EXP_OVERFLOW * 1.01),
            1 => (EXP_UNDERFLOW_TO_ZERO * 1.01)..(EXP_UNDERFLOW_TO_ZERO * 0.99),
        ];
        (delta, target, any::<bool>())
            .prop_map(|(d, t, above_one)| {
                let x = if above_one { 1.0 + d } else { 1.0 - d };
                let lx = x.ln();
                let y = if lx == 0.0 { 0.0 } else { t / lx };
                (x, y)
            })
            .boxed()
    }

    fn atan2_axis_inputs() -> BoxedStrategy<(f64, f64)> {
        prop_oneof![
            3 => (
                prop_oneof![
                    2 => tiny_signed(),
                    2 => subnormal_f64(),
                    1 => Just(0.0),
                    1 => Just(-0.0),
                ],
                normal_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
            ),
            3 => (
                normal_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
                prop_oneof![
                    2 => tiny_signed(),
                    2 => subnormal_f64(),
                    1 => Just(0.0),
                    1 => Just(-0.0),
                ],
            ),
        ]
        .boxed()
    }

    fn atan2_ratio_stress_inputs() -> BoxedStrategy<(f64, f64)> {
        prop_oneof![
            2 => (normal_f64_with_exp(200, 600), normal_f64_with_exp(-20, 20)),
            2 => (normal_f64_with_exp(-20, 20), normal_f64_with_exp(200, 600)),
        ]
        .boxed()
    }

    fn hypot_threshold_inputs() -> BoxedStrategy<(f64, f64)> {
        let large_boundary =
            (any::<bool>(), 0u32..=32u32, any::<bool>()).prop_map(|(neg, steps, above)| {
                let mut v = HYPOT_LARGE_VAL;
                for _ in 0..steps {
                    v = if above { v.next_up() } else { v.next_down() };
                }
                let v = if neg { -v } else { v };
                let ay = v.abs()
                    * HYPOT_EPS
                    * if above {
                        1.0 + 2.0f64.powi(-20)
                    } else {
                        1.0 - 2.0f64.powi(-20)
                    };
                (v, ay.copysign(v))
            });
        let tiny_boundary =
            (any::<bool>(), 0u32..=32u32, any::<bool>()).prop_map(|(neg, steps, above)| {
                let mut v = HYPOT_TINY_VAL;
                for _ in 0..steps {
                    v = if above { v.next_up() } else { v.next_down() };
                }
                let v = if neg { -v } else { v };
                let ax = (v.abs() / HYPOT_EPS)
                    * if above {
                        1.0 + 2.0f64.powi(-20)
                    } else {
                        1.0 - 2.0f64.powi(-20)
                    };
                (ax.copysign(v), v)
            });
        let subnormals = (subnormal_f64(), subnormal_f64());
        prop_oneof![
            3 => large_boundary,
            3 => tiny_boundary,
            2 => subnormals,
        ]
        .boxed()
    }

    fn pow2_exact() -> BoxedStrategy<f64> {
        (F64_EXP_MIN..=F64_EXP_MAX)
            .prop_map(|k| 2.0f64.powi(k))
            .boxed()
    }

    fn pow2_neighbors() -> BoxedStrategy<f64> {
        (F64_EXP_MIN..=F64_EXP_MAX, 0u32..=8u32, any::<bool>())
            .prop_map(|(k, steps, up)| {
                let mut v = 2.0f64.powi(k);
                for _ in 0..steps {
                    v = if up { v.next_up() } else { v.next_down() };
                }
                v
            })
            .boxed()
    }

    fn pow10_exact() -> BoxedStrategy<f64> {
        (-308i32..=308i32).prop_map(|k| 10.0f64.powi(k)).boxed()
    }

    fn wide_signed_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            2 => subnormal_f64(),
            2 => tiny_signed(),
            6 => normal_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
            1 => Just(0.0),
            1 => Just(-0.0),
        ]
        .boxed()
    }

    fn range_with_edges(min: f64, max: f64) -> BoxedStrategy<f64> {
        prop_oneof![
            6 => min..max,
            2 => ulp_steps(min, 256, true),
            2 => ulp_steps(max, 256, false),
            2 => tiny_signed(),
            1 => Just(0.0),
            1 => Just(-0.0),
        ]
        .boxed()
    }

    fn ptest_exp_inputs() -> BoxedStrategy<f64> {
        let wide = (EXP_UNDERFLOW_TO_ZERO - 5.0)..(EXP_OVERFLOW + 5.0);
        let mid = -50.0..50.0_f64;
        let boundaries = prop_oneof![
            2 => around(EXP_OVERFLOW, 256),
            2 => around(EXP_UNDERFLOW_TO_ZERO, 256),
        ];
        prop_oneof![
            3 => wide,
            2 => mid,
            2 => tiny_signed(),
            2 => exp_table_boundary_inputs(),
            1 => boundaries,
        ]
        .boxed()
    }

    fn ptest_exp2_inputs() -> BoxedStrategy<f64> {
        let wide = -1074.0..1024.0_f64;
        let mid = -20.0..20.0_f64;
        prop_oneof![
            4 => wide,
            2 => mid,
            2 => tiny_signed(),
            1 => ulp_steps(-1074.0, 256, true),
            1 => ulp_steps(1024.0, 256, false),
        ]
        .boxed()
    }

    fn ptest_expm1_inputs() -> BoxedStrategy<f64> {
        let mid = -50.0..50.0_f64;
        prop_oneof![
            5 => mid,
            3 => tiny_signed(),
            1 => ulp_steps(-50.0, 256, true),
            1 => ulp_steps(50.0, 256, false),
        ]
        .boxed()
    }

    fn ptest_exp10_inputs() -> BoxedStrategy<f64> {
        let wide = -308.0..308.0_f64;
        let mid = -10.0..10.0_f64;
        prop_oneof![
            4 => wide,
            2 => mid,
            2 => tiny_signed(),
            1 => ulp_steps(-308.0, 256, true),
            1 => ulp_steps(308.0, 256, false),
        ]
        .boxed()
    }

    fn ptest_ln_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            2 => subnormal_pos_f64(),
            2 => log_subnormal_boundary_inputs(),
            3 => tiny_positive(),
            3 => near_one_both(),
            2 => ln_near_one_boundary_inputs(),
            2 => pow2_exact(),
            2 => pow10_exact(),
            2 => ln_table_boundary_inputs(),
            6 => normal_pos_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
        ]
        .boxed()
    }

    fn ptest_log2_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            2 => subnormal_pos_f64(),
            2 => log_subnormal_boundary_inputs(),
            3 => tiny_positive(),
            3 => near_one_both(),
            4 => pow2_neighbors(),
            2 => ln_table_boundary_inputs(),
            4 => normal_pos_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
        ]
        .boxed()
    }

    fn ptest_log10_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            2 => subnormal_pos_f64(),
            2 => log_subnormal_boundary_inputs(),
            3 => tiny_positive(),
            3 => near_one_both(),
            3 => log10_near_one_inputs(),
            3 => pow10_exact(),
            2 => ln_table_boundary_inputs(),
            4 => normal_pos_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
        ]
        .boxed()
    }

    fn ptest_log1p_inputs() -> BoxedStrategy<f64> {
        let mid = -0.9..0.9_f64;
        prop_oneof![
            4 => mid,
            2 => near_minus_one_open(),
            2 => tiny_signed(),
            2 => normal_pos_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
            1 => ulp_steps(-1.0, 256, true),
        ]
        .boxed()
    }

    fn ptest_trig_inputs() -> BoxedStrategy<f64> {
        let mid = -1.0e6..1.0e6_f64;
        let large = -1.0e20..1.0e20_f64;
        let near_half_pi = (
            -1_000_000i32..=1_000_000i32,
            1u32..=NEAR_ONE_MAX_POW,
            any::<bool>(),
        )
            .prop_map(|(k, p, sign)| {
                let delta = 2.0f64.powi(-(p as i32));
                let base = (k as f64) * (PI / 2.0);
                if sign { base + delta } else { base - delta }
            });
        prop_oneof![
            4 => mid,
            2 => large,
            2 => near_half_pi,
            1 => tiny_signed(),
            1 => Just(0.0),
        ]
        .boxed()
    }

    fn ptest_tan_inputs() -> BoxedStrategy<f64> {
        let mid = -1.0e6..1.0e6_f64;
        let near_singular = (
            -300_000i32..=300_000i32,
            1u32..=NEAR_ONE_MAX_POW,
            any::<bool>(),
        )
            .prop_map(|(k, p, sign)| {
                let delta = 2.0f64.powi(-(p as i32));
                let base = (k as f64) * PI + FRAC_PI_2;
                if sign { base + delta } else { base - delta }
            });
        prop_oneof![
            4 => mid,
            2 => near_singular,
            2 => tiny_signed(),
            1 => ulp_steps(-1.0e6, 256, true),
            1 => ulp_steps(1.0e6, 256, false),
        ]
        .boxed()
    }

    fn ptest_atan_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            5 => -1.0e6..1.0e6_f64,
            2 => tiny_signed(),
            2 => normal_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
            1 => ulp_steps(-1.0e6, 256, true),
            1 => ulp_steps(1.0e6, 256, false),
        ]
        .boxed()
    }

    fn unit_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            4 => -1.0..1.0_f64,
            2 => near_one_signed_open(),
            2 => tiny_signed(),
            1 => Just(1.0),
            1 => Just(-1.0),
        ]
        .boxed()
    }

    fn ptest_atan2_inputs() -> BoxedStrategy<(f64, f64)> {
        let mid = (-1.0e6..1.0e6_f64, -1.0e6..1.0e6_f64);
        let axes = prop_oneof![
            1 => (tiny_signed(), Just(0.0)),
            1 => (Just(0.0), tiny_signed()),
            1 => (tiny_signed(), tiny_signed()),
        ];
        let wide = (
            normal_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
            normal_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
        );
        let specials = proptest::sample::select(vec![
            (0.0, 0.0),
            (-0.0, 0.0),
            (0.0, -0.0),
            (-0.0, -0.0),
            (f64::INFINITY, 1.0),
            (f64::NEG_INFINITY, 1.0),
            (1.0, f64::INFINITY),
            (1.0, f64::NEG_INFINITY),
            (f64::INFINITY, f64::INFINITY),
            (f64::INFINITY, f64::NEG_INFINITY),
            (f64::NEG_INFINITY, f64::INFINITY),
            (f64::NEG_INFINITY, f64::NEG_INFINITY),
        ]);
        prop_oneof![
            4 => mid,
            2 => axes,
            2 => wide,
            2 => atan2_axis_inputs(),
            1 => atan2_ratio_stress_inputs(),
            1 => specials,
        ]
        .boxed()
    }

    fn hypot_arg_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            4 => -1.0e200..1.0e200_f64,
            2 => tiny_signed(),
            2 => normal_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
            1 => Just(0.0),
            1 => Just(-0.0),
        ]
        .boxed()
    }

    fn ptest_hypot_inputs() -> BoxedStrategy<(f64, f64)> {
        prop_oneof![
            4 => (hypot_arg_inputs(), hypot_arg_inputs()),
            2 => hypot_threshold_inputs(),
        ]
        .boxed()
    }

    fn ptest_roundtrip_pos_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            4 => normal_pos_f64_with_exp(-1, 1),
            2 => near_one_both(),
        ]
        .boxed()
    }

    fn ptest_roundtrip_exp_inputs() -> BoxedStrategy<f64> {
        let neg = -1.0..-ROUNDTRIP_EXP_MIN_ABS;
        let pos = ROUNDTRIP_EXP_MIN_ABS..1.0_f64;
        prop_oneof![
            4 => neg,
            4 => pos,
            1 => around(-1.0, 256),
            1 => Just(-ROUNDTRIP_EXP_MIN_ABS),
            1 => ulp_steps_exclusive(-ROUNDTRIP_EXP_MIN_ABS, 256, false),
            1 => Just(ROUNDTRIP_EXP_MIN_ABS),
            1 => ulp_steps_exclusive(ROUNDTRIP_EXP_MIN_ABS, 256, true),
            1 => around(1.0, 256),
        ]
        .boxed()
    }

    fn ptest_hypot_similar_inputs() -> BoxedStrategy<(f64, f64)> {
        (
            normal_pos_f64_with_exp(-100, 100),
            0.5..2.0_f64,
            any::<bool>(),
            any::<bool>(),
        )
            .prop_map(|(base, ratio, neg_x, neg_y)| {
                let x = if neg_x { -base } else { base };
                let y = if neg_y { -(base * ratio) } else { base * ratio };
                (x, y)
            })
            .boxed()
    }

    fn ptest_sinh_inputs() -> BoxedStrategy<f64> {
        let mid = range_with_edges(-100.0, 100.0);
        let thresholds = prop_oneof![
            2 => around_signed(3.725_290_298_461_914e-09, 256),
            2 => around_signed(0.5, 256),
            2 => around_signed(1.0, 256),
            2 => around_signed(22.0, 256),
            2 => around_signed(EXP_OVERFLOW, 256),
            2 => around_signed(SINH_OVERFLOW, 256),
        ];
        let tiny = prop_oneof![
            2 => tiny_signed(),
            2 => subnormal_f64(),
            1 => Just(0.0),
            1 => Just(-0.0),
        ];
        let wide = prop_oneof![
            2 => normal_f64_with_exp(-20, 20),
            1 => normal_f64_with_exp(21, 200),
        ];
        prop_oneof![
            4 => mid,
            4 => thresholds,
            2 => tiny,
            1 => wide,
        ]
        .boxed()
    }

    fn ptest_sinh_nonneg_inputs() -> BoxedStrategy<f64> {
        let mid = range_with_edges(0.0, 100.0);
        let thresholds = prop_oneof![
            2 => around(3.725_290_298_461_914e-09, 256),
            2 => around(0.5, 256),
            2 => around(1.0, 256),
            2 => around(22.0, 256),
            2 => around(EXP_OVERFLOW, 256),
            2 => around_signed_below(SINH_OVERFLOW, 256),
        ];
        let tiny = prop_oneof![
            2 => tiny_positive(),
            2 => subnormal_pos_f64(),
            1 => Just(0.0),
        ];
        let wide = prop_oneof![
            2 => normal_pos_f64_with_exp(-20, 9),
            1 => normal_pos_f64_with_exp(10, 10),
        ];
        prop_oneof![
            4 => mid,
            4 => thresholds,
            2 => tiny,
            1 => wide,
        ]
        .boxed()
    }

    fn ptest_sinh_below_overflow_inputs() -> BoxedStrategy<f64> {
        let mid = range_with_edges(-100.0, 100.0);
        let thresholds = prop_oneof![
            2 => around_signed(3.725_290_298_461_914e-09, 256),
            2 => around_signed(0.5, 256),
            2 => around_signed(1.0, 256),
            2 => around_signed(22.0, 256),
            2 => around_signed(EXP_OVERFLOW, 256),
            2 => around_signed_below(SINH_OVERFLOW, 512),
        ];
        let tiny = prop_oneof![
            2 => tiny_signed(),
            2 => subnormal_f64(),
            1 => Just(0.0),
            1 => Just(-0.0),
        ];
        let wide = (any::<bool>(), normal_pos_f64_with_exp(-20, 9))
            .prop_map(|(neg, val)| if neg { -val } else { val })
            .boxed();
        prop_oneof![
            4 => mid,
            4 => thresholds,
            2 => tiny,
            1 => wide,
        ]
        .prop_map(clamp_below_sinh_overflow)
        .boxed()
    }
    fn ptest_cosh_inputs() -> BoxedStrategy<f64> {
        range_with_edges(-100.0, 100.0)
    }

    fn ptest_tanh_inputs() -> BoxedStrategy<f64> {
        range_with_edges(-20.0, 20.0)
    }

    fn ptest_pow_inputs() -> BoxedStrategy<(f64, f64)> {
        let base = prop_oneof![
            4 => -10.0..10.0_f64,
            2 => tiny_signed(),
            2 => near_one_signed_open(),
            1 => subnormal_f64(),
            1 => Just(-1.0),
            1 => Just(1.0),
            1 => Just(0.0),
        ];
        let exp = prop_oneof![
            3 => -10.0..10.0_f64,
            2 => (-1000i32..=1000i32).prop_map(|k| k as f64),
            2 => pow_near_integer_exponent(),
            1 => (-10i32..=10i32).prop_map(|k| k as f64 + 0.5),
            1 => tiny_signed(),
        ];
        (base, exp).boxed()
    }

    fn ptest_sinh_threshold_inputs() -> BoxedStrategy<(f64, &'static str)> {
        prop_oneof![
            2 => tagged_f64("tiny", around_signed(3.725_290_298_461_914e-09, 512)),
            2 => tagged_f64("half", around_signed(0.5, 512)),
            2 => tagged_f64("one", around_signed(1.0, 512)),
            2 => tagged_f64("small", around_signed(22.0, 512)),
            2 => tagged_f64("exp_hi", around_signed(EXP_OVERFLOW, 512)),
            2 => tagged_f64("sinh_overflow", around_signed(SINH_OVERFLOW, 512)),
        ]
        .boxed()
    }

    fn ptest_exp_threshold_inputs() -> BoxedStrategy<(f64, &'static str)> {
        prop_oneof![
            3 => tagged_f64("table_boundary", exp_table_boundary_inputs()),
            2 => tagged_f64("overflow", around(EXP_OVERFLOW, 512)),
            2 => tagged_f64("underflow", around(EXP_UNDERFLOW_TO_ZERO, 512)),
            1 => tagged_f64("tiny", tiny_signed()),
        ]
        .boxed()
    }

    fn ptest_ln_threshold_inputs() -> BoxedStrategy<(f64, &'static str)> {
        prop_oneof![
            3 => tagged_f64("table_boundary", ln_table_boundary_inputs()),
            2 => tagged_f64("near_one", ln_near_one_boundary_inputs()),
            2 => tagged_f64("pow2_neighbors", pow2_neighbors()),
            1 => tagged_f64("subnormal_boundary", log_subnormal_boundary_inputs()),
        ]
        .boxed()
    }

    fn ptest_log10_threshold_inputs() -> BoxedStrategy<(f64, &'static str)> {
        prop_oneof![
            3 => tagged_f64("near_one", log10_near_one_inputs()),
            2 => tagged_f64("k_branch", normal_pos_f64_with_exp(-1, 0)),
            2 => tagged_f64("subnormal_boundary", log_subnormal_boundary_inputs()),
            1 => tagged_f64("table_boundary", ln_table_boundary_inputs()),
        ]
        .boxed()
    }

    fn ptest_pow_threshold_inputs() -> BoxedStrategy<((f64, f64), &'static str)> {
        let neg_base_int = (-10.0..-0.0_f64, (-100i32..=100i32).prop_map(|k| k as f64));
        let near_one_base = prop_oneof![
            2 => near_one_both(),
            2 => 0.5..2.0_f64,
            1 => tiny_positive(),
        ];
        prop_oneof![
            3 => tagged_pair("near_integer_exp", (near_one_base, pow_near_integer_exponent())),
            2 => tagged_pair("neg_base_int", neg_base_int),
            2 => tagged_pair("exp_boundary", pow_exp_boundary_inputs()),
        ]
        .boxed()
    }

    fn ptest_atan2_threshold_inputs() -> BoxedStrategy<((f64, f64), &'static str)> {
        let specials = proptest::sample::select(vec![
            (0.0, 0.0),
            (-0.0, 0.0),
            (0.0, -0.0),
            (-0.0, -0.0),
            (f64::INFINITY, 1.0),
            (f64::NEG_INFINITY, 1.0),
            (1.0, f64::INFINITY),
            (1.0, f64::NEG_INFINITY),
            (f64::INFINITY, f64::INFINITY),
            (f64::INFINITY, f64::NEG_INFINITY),
            (f64::NEG_INFINITY, f64::INFINITY),
            (f64::NEG_INFINITY, f64::NEG_INFINITY),
        ]);
        prop_oneof![
            3 => tagged_pair("axis", atan2_axis_inputs()),
            2 => tagged_pair("ratio_stress", atan2_ratio_stress_inputs()),
            1 => tagged_pair("specials", specials),
        ]
        .boxed()
    }

    fn ptest_hypot_threshold_inputs() -> BoxedStrategy<((f64, f64), &'static str)> {
        let mixed = (subnormal_f64(), normal_f64_with_exp(-20, 20));
        prop_oneof![
            3 => tagged_pair("thresholds", hypot_threshold_inputs()),
            1 => tagged_pair("mixed_subnormal", mixed),
        ]
        .boxed()
    }

    fn ptest_sqrt_inputs() -> BoxedStrategy<f64> {
        range_with_edges(-1.0e300, 1.0e300)
    }

    fn ptest_cbrt_inputs() -> BoxedStrategy<f64> {
        range_with_edges(-1.0e300, 1.0e300)
    }

    fn nonzero_divisor_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            4 => 1.0e-6..1.0e6_f64,
            4 => -1.0e6..-1.0e-6_f64,
            2 => tiny_signed(),
            2 => normal_f64_with_exp(-20, 20),
        ]
        .boxed()
    }

    fn ptest_fmod_inputs() -> BoxedStrategy<(f64, f64)> {
        (range_with_edges(-1.0e6, 1.0e6), nonzero_divisor_inputs()).boxed()
    }

    fn ptest_remainder_inputs() -> BoxedStrategy<(f64, f64)> {
        (range_with_edges(-1.0e6, 1.0e6), nonzero_divisor_inputs()).boxed()
    }

    fn ptest_asinh_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            4 => -1.0e20..1.0e20_f64,
            2 => tiny_signed(),
            2 => normal_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
            1 => ulp_steps(-1.0e20, 256, true),
            1 => ulp_steps(1.0e20, 256, false),
        ]
        .boxed()
    }

    fn ptest_acosh_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            4 => 1.0..1.0e20_f64,
            2 => near_one_above(),
            2 => normal_pos_f64_with_exp(0, F64_EXP_MAX),
            1 => Just(1.0),
            1 => ulp_steps(1.0, 256, true),
        ]
        .boxed()
    }

    fn ptest_atanh_inputs() -> BoxedStrategy<f64> {
        let mid = -0.999_999..0.999_999_f64;
        prop_oneof![
            5 => mid,
            3 => near_one_signed_open(),
            2 => tiny_signed(),
            1 => Just(0.0),
        ]
        .boxed()
    }

    fn ptest_erf_inputs() -> BoxedStrategy<f64> {
        range_with_edges(-6.0, 6.0)
    }

    fn gamma_inputs() -> BoxedStrategy<f64> {
        let pos = f64::MIN_POSITIVE..20.0_f64;
        let neg_non_int =
            (-20i32..=-1i32, 1u32..=NEAR_ONE_MAX_POW, any::<bool>()).prop_map(|(n, k, sign)| {
                let delta = 2.0f64.powi(-(k as i32));
                let base = n as f64;
                let mut x = if sign { base + delta } else { base - delta };
                if x <= 0.0 && x == x.trunc() {
                    x = nextafter_reference(x, f64::INFINITY);
                }
                x
            });
        let near_int =
            (-20i32..=20i32, 1u32..=NEAR_ONE_MAX_POW, any::<bool>()).prop_map(|(n, k, sign)| {
                let delta = 2.0f64.powi(-(k as i32));
                let base = n as f64;
                let mut x = if sign { base + delta } else { base - delta };
                if x <= 0.0 && x == x.trunc() {
                    x = nextafter_reference(x, f64::INFINITY);
                }
                x
            });
        prop_oneof![
            4 => pos,
            3 => neg_non_int,
            2 => near_int,
            1 => near_one_both(),
            1 => tiny_signed(),
        ]
        .boxed()
    }

    fn ptest_logb_inputs() -> BoxedStrategy<f64> {
        prop_oneof![
            2 => tiny_signed(),
            2 => subnormal_f64(),
            4 => normal_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX),
            1 => Just(f64::MIN_POSITIVE),
            1 => Just(-f64::MIN_POSITIVE),
        ]
        .boxed()
    }

    fn ptest_rounding_inputs() -> BoxedStrategy<f64> {
        let mid = -1.0e6..1.0e6_f64;
        let near_int = (
            -1_000_000i64..=1_000_000i64,
            1u32..=NEAR_ONE_MAX_POW,
            any::<bool>(),
        )
            .prop_map(|(n, k, sign)| {
                let delta = 2.0f64.powi(-(k as i32));
                let base = n as f64;
                if sign { base + delta } else { base - delta }
            });
        let half_int = (
            -1_000_000i64..=1_000_000i64,
            1u32..=NEAR_ONE_MAX_POW,
            any::<bool>(),
        )
            .prop_map(|(n, k, sign)| {
                let delta = 2.0f64.powi(-(k as i32));
                let base = n as f64 + 0.5;
                if sign { base + delta } else { base - delta }
            });
        prop_oneof![
            4 => mid,
            2 => near_int,
            2 => half_int,
            1 => tiny_signed(),
            1 => normal_f64_with_exp(20, F64_EXP_MAX),
        ]
        .boxed()
    }

    fn ptest_fma_inputs() -> BoxedStrategy<(f64, f64, f64)> {
        let mid = (-1.0e3..1.0e3_f64, -1.0e3..1.0e3_f64, -1.0e3..1.0e3_f64);
        let tiny = (tiny_signed(), tiny_signed(), tiny_signed());
        let cancellation = (
            -1.0e3..1.0e3_f64,
            -1.0e3..1.0e3_f64,
            1u32..=NEAR_ONE_MAX_POW,
            any::<bool>(),
        )
            .prop_map(|(x, y, k, sign)| {
                let delta = 2.0f64.powi(-(k as i32));
                let z = -(x * y) + if sign { delta } else { -delta };
                (x, y, z)
            });
        prop_oneof![
            4 => mid,
            2 => tiny,
            2 => cancellation,
        ]
        .boxed()
    }

    fn frexp_inputs() -> BoxedStrategy<f64> {
        wide_signed_inputs()
    }

    fn scalbn_x_inputs() -> BoxedStrategy<f64> {
        wide_signed_inputs()
    }

    fn scalbn_n_inputs() -> BoxedStrategy<i32> {
        prop_oneof![
            4 => -1000i32..1000i32,
            2 => -10i32..10i32,
            1 => Just(0i32),
            1 => Just(1000i32),
            1 => Just(-1000i32),
        ]
        .boxed()
    }

    fn scalbln_n_inputs() -> BoxedStrategy<i64> {
        prop_oneof![
            4 => -1000i64..1000i64,
            2 => -10i64..10i64,
            1 => Just(0i64),
            1 => Just(1000i64),
            1 => Just(-1000i64),
        ]
        .boxed()
    }

    fn ptest_remquo_inputs() -> BoxedStrategy<(f64, f64)> {
        (range_with_edges(-1.0e6, 1.0e6), nonzero_divisor_inputs()).boxed()
    }

    fn ptest_fdim_inputs() -> BoxedStrategy<(f64, f64)> {
        (
            range_with_edges(-1.0e6, 1.0e6),
            range_with_edges(-1.0e6, 1.0e6),
        )
            .boxed()
    }

    fn ptest_nextafter_inputs() -> BoxedStrategy<(f64, f64)> {
        prop_oneof![
            4 => (wide_signed_inputs(), wide_signed_inputs()),
            2 => (tiny_signed(), Just(0.0)),
            2 => (Just(0.0), tiny_signed()),
            1 => (subnormal_f64(), subnormal_f64()),
            1 => (normal_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX), normal_f64_with_exp(F64_EXP_MIN, F64_EXP_MAX)),
        ]
        .boxed()
    }
    proptest! {
        #[test]
        fn ptest_exp_special(x in proptest::sample::select(exp_special_inputs())) {
            let actual = fastmaths::exp(x);
            assert_ulp_eq_exp(actual, x, PROPTEST_ULP_TOL, &format!("exp special({x})"));
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_exp(x in ptest_exp_inputs()) {
            let actual = fastmaths::exp(x);
            assert_ulp_eq_exp(actual, x, PROPTEST_ULP_TOL, &format!("exp({x})"));
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_exp_thresholds((x, bucket) in ptest_exp_threshold_inputs()) {
            let actual = fastmaths::exp(x);
            assert_ulp_eq_exp(actual, x, PROPTEST_ULP_TOL, &format_case("exp", x, bucket));
        }

        #[test]
        fn ptest_exp_ln_roundtrip(x in ptest_roundtrip_pos_inputs()) {
            prop_assume!(x.is_finite() && x > 0.0);
            let y = fastmaths::ln(x);
            prop_assume!(y.is_finite());
            prop_assume!(y.abs() <= 20.0);
            let r = fastmaths::exp(y);
            prop_assume!(r.is_normal());
            assert_ulp_eq(r, x, COMPOSED_ULP_TOL, &format_case("exp_ln_roundtrip", x, "compose"));
        }

        #[test]
        fn ptest_ln_exp_roundtrip(x in ptest_roundtrip_exp_inputs()) {
            prop_assume!(x.is_finite());
            prop_assume!(x.abs() >= ROUNDTRIP_EXP_MIN_ABS);
            prop_assume!(x.abs() <= 20.0);
            let y = fastmaths::exp(x);
            prop_assume!(y.is_normal());
            let r = fastmaths::ln(y);
            assert_ulp_eq(r, x, COMPOSED_ULP_TOL, &format_case("ln_exp_roundtrip", x, "compose"));
        }

        #[test]
        fn ptest_ln(x in ptest_ln_inputs()) {
            if x.is_finite() && x > 0.0 {
                let actual = fastmaths::ln(x);
                let expected = ln_reference(x);
                assert_ulp_eq(
                    actual,
                    expected,
                    PROPTEST_ULP_TOL,
                    &format!("ln({x})"),
                );
            }
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_ln_thresholds((x, bucket) in ptest_ln_threshold_inputs()) {
            if x.is_finite() && x > 0.0 {
                let actual = fastmaths::ln(x);
                let expected = ln_reference(x);
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format_case("ln", x, bucket));
            }
        }

        #[test]
        fn ptest_sin(x in ptest_trig_inputs()) {
            let actual = fastmaths::sin(x);
            let expected = sin_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format!("sin({x})"),
            );
        }

        #[test]
        fn ptest_cos(x in ptest_trig_inputs()) {
            let actual = fastmaths::cos(x);
            let expected = cos_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format!("cos({x})"),
            );
        }

        #[test]
        fn ptest_sincos(x in ptest_trig_inputs()) {
            let (s_actual, c_actual) = fastmaths::sincos(x);
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
        fn ptest_exp2(x in ptest_exp2_inputs()) {
            let actual = fastmaths::exp2(x);
            let expected = exp2_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("exp2({x})"));
        }

        #[test]
        fn ptest_expm1(x in ptest_expm1_inputs()) {
            let actual = fastmaths::expm1(x);
            let expected = expm1_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format!("expm1({x})"),
            );
        }

        #[test]
        fn ptest_log2(x in ptest_log2_inputs()) {
            if x.is_finite() && x > 0.0 {
                let actual = fastmaths::log2(x);
                let expected = log2_reference(x);
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("log2({x})"));
            }
        }

        #[test]
        fn ptest_log2_roundtrip(x in ptest_roundtrip_pos_inputs()) {
            prop_assume!(x.is_finite() && x > 0.0);
            let y = fastmaths::log2(x);
            prop_assume!(y.is_finite());
            let r = fastmaths::exp2(y);
            prop_assume!(r.is_normal());
            assert_ulp_eq(r, x, COMPOSED_ULP_TOL, &format_case("exp2_log2_roundtrip", x, "compose"));
        }

        #[test]
        fn ptest_log10(x in ptest_log10_inputs()) {
            if x.is_finite() && x > 0.0 {
                let actual = fastmaths::log10(x);
                let expected = log10_reference(x);
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("log10({x})"));
            }
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_log10_thresholds((x, bucket) in ptest_log10_threshold_inputs()) {
            if x.is_finite() && x > 0.0 {
                let actual = fastmaths::log10(x);
                let expected = log10_reference(x);
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format_case("log10", x, bucket));
            }
        }

        #[test]
        fn ptest_log1p(x in ptest_log1p_inputs()) {
            let actual = fastmaths::log1p(x);
            let expected = log1p_reference(x);
            if expected.is_nan() {
                assert!(actual.is_nan(), "log1p({x}) expected NaN");
            } else {
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("log1p({x})"));
            }
        }

        #[test]
        fn ptest_tan(x in ptest_tan_inputs()) {
            let actual = fastmaths::tan(x);
            let expected = tan_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format!("tan({x})"),
            );
        }

        #[test]
        fn ptest_atan(x in ptest_atan_inputs()) {
            let actual = fastmaths::atan(x);
            let expected = atan_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("atan({x})"));
        }

        #[test]
        fn ptest_asin(x in unit_inputs()) {
            let actual = fastmaths::asin(x);
            let expected = asin_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("asin({x})"));
        }

        #[test]
        fn ptest_acos(x in unit_inputs()) {
            let actual = fastmaths::acos(x);
            let expected = acos_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("acos({x})"));
        }

        #[test]
        fn ptest_atan2((y, x) in ptest_atan2_inputs()) {
            let actual = fastmaths::atan2(y, x);
            let expected = atan2_reference(y, x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format_case2("atan2", y, x, "proptest"),
            );
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_atan2_thresholds(((y, x), bucket) in ptest_atan2_threshold_inputs()) {
            let actual = fastmaths::atan2(y, x);
            let expected = atan2_reference(y, x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format_case2("atan2", y, x, bucket),
            );
        }

        #[test]
        fn ptest_hypot((x, y) in ptest_hypot_inputs()) {
            let actual = fastmaths::hypot(x, y);
            let expected = hypot_reference(x, y);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format_case2("hypot", x, y, "proptest"),
            );
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_hypot_thresholds(((x, y), bucket) in ptest_hypot_threshold_inputs()) {
            let actual = fastmaths::hypot(x, y);
            let expected = hypot_reference(x, y);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format_case2("hypot", x, y, bucket),
            );
        }

        #[test]
        fn ptest_hypot_sqrt_identity((x, y) in ptest_hypot_similar_inputs()) {
            prop_assume!(x.is_finite() && y.is_finite());
            prop_assume!(x.abs() <= 1.0e150 && y.abs() <= 1.0e150);
            let sum = x * x + y * y;
            prop_assume!(sum.is_finite());
            let expected = fastmaths::sqrt(sum);
            let actual = fastmaths::hypot(x, y);
            assert_ulp_eq(
                actual,
                expected,
                COMPOSED_ULP_TOL,
                &format_case2("hypot_sqrt", x, y, "compose"),
            );
        }

        #[test]
        fn ptest_hypot_scaling_identity((x, y) in ptest_hypot_similar_inputs(), k in -10i32..=10i32) {
            prop_assume!(x.is_finite() && y.is_finite());
            let sx = fastmaths::scalbn(x, k);
            let sy = fastmaths::scalbn(y, k);
            prop_assume!(sx.is_finite() && sy.is_finite());
            let actual = fastmaths::hypot(sx, sy);
            let base = fastmaths::hypot(x, y);
            let expected = fastmaths::scalbn(base, k);
            prop_assume!(actual.is_finite() && expected.is_finite());
            assert_ulp_eq(
                actual,
                expected,
                COMPOSED_ULP_TOL,
                &format_case2("hypot_scaling", sx, sy, "compose"),
            );
        }

        #[test]
        fn ptest_sinh(x in ptest_sinh_inputs()) {
            let actual = fastmaths::sinh(x);
            let expected = sinh_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format_case("sinh", x, "proptest"),
            );
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_sinh_thresholds((x, bucket) in ptest_sinh_threshold_inputs()) {
            let actual = fastmaths::sinh(x);
            let expected = sinh_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format_case("sinh", x, bucket));
        }

        #[test]
        fn ptest_sinh_odd_symmetry(x in ptest_sinh_inputs()) {
            prop_assume!(x.is_finite());
            let a = fastmaths::sinh(x);
            let b = fastmaths::sinh(-x);
            prop_assume!(!a.is_nan() && !b.is_nan());
            assert_eq!(
                b.to_bits(),
                (-a).to_bits(),
                "{}",
                format_case("sinh_odd", x, "property")
            );
        }

        #[test]
        fn ptest_sinh_monotonic_nonneg(x in ptest_sinh_nonneg_inputs()) {
            let x2 = x.next_up();
            prop_assume!(x2.is_finite());
            let s1 = fastmaths::sinh(x);
            let s2 = fastmaths::sinh(x2);
            prop_assume!(!s1.is_nan() && !s2.is_nan());
            assert!(
                s2 >= s1,
                "{}",
                format_case("sinh_monotonic", x, "property")
            );
        }

        #[test]
        fn ptest_sinh_no_premature_overflow(x in ptest_sinh_below_overflow_inputs()) {
            let s = fastmaths::sinh(x);
            assert!(
                s.is_finite(),
                "{}",
                format_case("sinh_no_overflow", x, "property")
            );
        }

        #[test]
        fn ptest_sinh_tiny_exact(x in tiny_signed_below_tiny()) {
            let s = fastmaths::sinh(x);
            assert_eq!(
                s.to_bits(),
                x.to_bits(),
                "{}",
                format_case("sinh_tiny_exact", x, "property")
            );
        }

        #[test]
        fn ptest_cosh(x in ptest_cosh_inputs()) {
            let actual = fastmaths::cosh(x);
            let expected = cosh_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("cosh({x})"));
        }

        #[test]
        fn ptest_tanh(x in ptest_tanh_inputs()) {
            let actual = fastmaths::tanh(x);
            let expected = tanh_reference(x);
            assert_ulp_eq(actual, expected, TANH_ULP_TOL, &format!("tanh({x})"));
        }

        #[test]
        fn ptest_pow((x, y) in ptest_pow_inputs()) {
            let actual = fastmaths::pow(x, y);
            let expected = pow_reference(x, y);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format_case2("pow", x, y, "proptest"),
            );
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_pow_thresholds(((x, y), bucket) in ptest_pow_threshold_inputs()) {
            let actual = fastmaths::pow(x, y);
            let expected = pow_reference(x, y);
            if expected.is_nan() {
                assert!(actual.is_nan(), "{}", format_case2("pow", x, y, bucket));
            } else {
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format_case2("pow", x, y, bucket));
            }
        }

        #[test]
        fn ptest_pow_exp_ln_roundtrip((x, y) in (ptest_roundtrip_pos_inputs(), -2.0..2.0_f64)) {
            prop_assume!(x.is_finite() && x > 0.0 && y.is_finite());
            let lx = fastmaths::ln(x);
            prop_assume!(lx.is_finite());
            let t = lx * y;
            prop_assume!(t.is_finite());
            prop_assume!(t.abs() <= 20.0);
            let expected = fastmaths::exp(t);
            prop_assume!(expected.is_normal());
            let actual = fastmaths::pow(x, y);
            prop_assume!(actual.is_finite());
            assert_ulp_eq(
                actual,
                expected,
                COMPOSED_ULP_TOL,
                &format_case2("pow_exp_ln_roundtrip", x, y, "compose"),
            );
        }

        #[test]
        fn ptest_sqrt(x in ptest_sqrt_inputs()) {
            let actual = fastmaths::sqrt(x);
            let expected = sqrt_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("sqrt({x})"));
        }

        #[test]
        fn ptest_cbrt(x in ptest_cbrt_inputs()) {
            let actual = fastmaths::cbrt(x);
            let expected = cbrt_reference(x);
            assert_ulp_eq(
                actual,
                expected,
                PROPTEST_ULP_TOL,
                &format!("cbrt({x})"),
            );
        }

        #[test]
        fn ptest_fmod((x, y) in ptest_fmod_inputs()) {
            let actual = fastmaths::fmod(x, y);
            let expected = fmod_reference(x, y);
            if expected.is_nan() {
                assert!(actual.is_nan(), "fmod({x},{y}) expected NaN");
            } else {
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("fmod({x},{y})"));
            }
        }

        #[test]
        fn ptest_remainder((x, y) in ptest_remainder_inputs()) {
            let actual = fastmaths::remainder(x, y);
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

        #[test]
        fn ptest_asinh(x in ptest_asinh_inputs()) {
            let actual = fastmaths::asinh(x);
            let expected = asinh_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("asinh({x})"));
        }

        #[test]
        fn ptest_acosh(x in ptest_acosh_inputs()) {
            let actual = fastmaths::acosh(x);
            let expected = acosh_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("acosh({x})"));
        }

        #[test]
        fn ptest_atanh(x in ptest_atanh_inputs()) {
            let actual = fastmaths::atanh(x);
            let expected = atanh_reference(x);
            assert_ulp_eq(actual, expected, ATANH_ULP_TOL, &format!("atanh({x})"));
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_erf(x in ptest_erf_inputs()) {
            let actual = fastmaths::erf(x);
            let expected = erf_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("erf({x})"));
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_erfc(x in ptest_erf_inputs()) {
            let actual = fastmaths::erfc(x);
            let expected = erfc_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("erfc({x})"));
        }

        #[test]
        fn ptest_exp10(x in ptest_exp10_inputs()) {
            let actual = fastmaths::exp10(x);
            let expected = exp10_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("exp10({x})"));
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_lgamma(x in gamma_inputs()) {
            prop_assume!(!(x <= 0.0 && x == x.trunc()));
            let actual = fastmaths::lgamma(x);
            let expected = lgamma_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("lgamma({x})"));
        }

        #[cfg(feature = "mpfr")]
        #[test]
        fn ptest_tgamma(x in gamma_inputs()) {
            prop_assume!(!(x <= 0.0 && x == x.trunc()));
            let actual = fastmaths::tgamma(x);
            let expected = tgamma_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("tgamma({x})"));
        }

        #[test]
        fn ptest_logb(x in ptest_logb_inputs()) {
            if x != 0.0 {
                let actual = fastmaths::logb(x);
                let expected = logb_reference(x);
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("logb({x})"));
            }
        }

        #[test]
        fn ptest_ilogb(x in ptest_logb_inputs()) {
            if x != 0.0 {
                let actual = fastmaths::ilogb(x);
                let expected = ilogb_reference(x);
                assert_eq!(actual, expected, "ilogb({x}) expected {expected}, got {actual}");
            }
        }

        #[test]
        fn ptest_modf(x in ptest_rounding_inputs()) {
            let (frac, int) = fastmaths::modf(x);
            let (frac_e, int_e) = modf_reference(x);
            assert_ulp_eq(frac, frac_e, PROPTEST_ULP_TOL, &format!("modf frac({x})"));
            assert_ulp_eq(int, int_e, PROPTEST_ULP_TOL, &format!("modf int({x})"));
        }

        #[test]
        fn ptest_floor(x in ptest_rounding_inputs()) {
            let actual = fastmaths::floor(x);
            let expected = floor_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("floor({x})"));
        }

        #[test]
        fn ptest_ceil(x in ptest_rounding_inputs()) {
            let actual = fastmaths::ceil(x);
            let expected = ceil_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("ceil({x})"));
        }

        #[test]
        fn ptest_trunc(x in ptest_rounding_inputs()) {
            let actual = fastmaths::trunc(x);
            let expected = trunc_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("trunc({x})"));
        }

        #[test]
        fn ptest_round(x in ptest_rounding_inputs()) {
            let actual = fastmaths::round(x);
            let expected = round_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("round({x})"));
        }

        #[test]
        fn ptest_rint(x in ptest_rounding_inputs()) {
            let actual = fastmaths::rint(x);
            let expected = rint_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("rint({x})"));
        }

        #[test]
        fn ptest_nearbyint(x in ptest_rounding_inputs()) {
            let actual = fastmaths::nearbyint(x);
            let expected = nearbyint_reference(x);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("nearbyint({x})"));
        }

        #[test]
        fn ptest_lrint(x in ptest_rounding_inputs()) {
            let actual = fastmaths::lrint(x);
            let expected = lrint_reference(x);
            assert_eq!(actual, expected, "lrint({x})");
        }

        #[test]
        fn ptest_llrint(x in ptest_rounding_inputs()) {
            let actual = fastmaths::llrint(x);
            let expected = llrint_reference(x);
            assert_eq!(actual, expected, "llrint({x})");
        }

        #[test]
        fn ptest_lround(x in ptest_rounding_inputs()) {
            let actual = fastmaths::lround(x);
            let expected = lround_reference(x);
            assert_eq!(actual, expected, "lround({x})");
        }

        #[test]
        fn ptest_llround(x in ptest_rounding_inputs()) {
            let actual = fastmaths::llround(x);
            let expected = llround_reference(x);
            assert_eq!(actual, expected, "llround({x})");
        }

        #[test]
        fn ptest_fma((x, y, z) in ptest_fma_inputs()) {
            let actual = fastmaths::fma(x, y, z);
            let expected = fma_reference(x, y, z);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("fma({x},{y},{z})"));
        }

        #[test]
        fn ptest_frexp(x in frexp_inputs()) {
            let (m_a, e_a) = fastmaths::frexp(x);
            let (m_e, e_e) = frexp_reference(x);
            assert_ulp_eq(m_a, m_e, PROPTEST_ULP_TOL, &format!("frexp({x}) mantissa"));
            assert_eq!(e_a, e_e, "frexp({x}) exponent");
        }

        #[test]
        fn ptest_scalbn(x in scalbn_x_inputs(), n in scalbn_n_inputs()) {
            let actual = fastmaths::scalbn(x, n);
            let expected = scalbn_reference(x, n);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("scalbn({x},{n})"));
        }

        #[test]
        fn ptest_scalbln(x in scalbn_x_inputs(), n in scalbln_n_inputs()) {
            let actual = fastmaths::scalbln(x, n);
            let expected = scalbln_reference(x, n);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("scalbln({x},{n})"));
        }

        #[test]
        fn ptest_remquo((x, y) in ptest_remquo_inputs()) {
            prop_assume!(y != 0.0);
            let (actual_r, actual_q) = fastmaths::remquo(x, y);
            let (expected_r, expected_q) = remquo_reference(x, y);
            assert_ulp_eq(actual_r, expected_r, PROPTEST_ULP_TOL, &format!("remquo({x},{y})"));
            assert_eq!(actual_q, expected_q, "remquo({x},{y}) quotient");
        }

        #[test]
        fn ptest_fdim((x, y) in ptest_fdim_inputs()) {
            let actual = fastmaths::fdim(x, y);
            let expected = fdim_reference(x, y);
            assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("fdim({x},{y})"));
        }

        #[test]
        fn ptest_fmax((x, y) in ptest_fdim_inputs()) {
            let actual = fastmaths::fmax(x, y);
            let expected = fmax_reference(x, y);
            if actual == 0.0 && expected == 0.0 {
                assert_eq!(actual.to_bits(), expected.to_bits(), "fmax sign mismatch");
            } else {
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("fmax({x},{y})"));
            }
        }

        #[test]
        fn ptest_fmin((x, y) in ptest_fdim_inputs()) {
            let actual = fastmaths::fmin(x, y);
            let expected = fmin_reference(x, y);
            if actual == 0.0 && expected == 0.0 {
                assert_eq!(actual.to_bits(), expected.to_bits(), "fmin sign mismatch");
            } else {
                assert_ulp_eq(actual, expected, PROPTEST_ULP_TOL, &format!("fmin({x},{y})"));
            }
        }

        #[test]
        fn ptest_nextafter((x, y) in ptest_nextafter_inputs()) {
            let actual = fastmaths::nextafter(x, y);
            let expected = nextafter_reference(x, y);
            if actual.is_nan() {
                assert!(expected.is_nan(), "nextafter expected NaN");
            } else {
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "nextafter({x},{y}) expected {expected}, got {actual}"
                );
            }
        }
    }
}

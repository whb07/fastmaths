//! remquo(x,y) implementation.
//!
//! Computes IEEE remainder and low bits of the quotient. Uses exponent alignment
//! and subtractive reduction with correct tie handling, mirroring glibc behavior.

use super::{f64_from_bits, fma_internal, fmod, rint};

const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;
const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;

#[inline(always)]
fn sig_exp(bits: u64) -> (u64, i32) {
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

#[inline(always)]
fn pow2_mod(mut exp: u32, modulus: u128) -> u128 {
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

#[inline(always)]
fn quotient_mod8(ax_bits: u64, ay_bits: u64) -> i32 {
    let (sigx, ex) = sig_exp(ax_bits);
    let (sigy, ey) = sig_exp(ay_bits);
    if sigx == 0 || sigy == 0 {
        return 0;
    }
    let shift = ex - ey;
    if shift >= 0 {
        let d = sigy as u128;
        let d16 = d << 4;
        let pow2_mod_d = pow2_mod(shift as u32, d);
        let pow2_mod_d16 = pow2_mod(shift as u32, d16);
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

#[inline(always)]
pub fn remquo(x: f64, y: f64) -> (f64, i32) {
    let mut hx = x.to_bits();
    let hy = y.to_bits();
    let sx = hx & SIGN_MASK;
    let qs = sx ^ (hy & SIGN_MASK);
    let mut ax = hx & 0x7fff_ffff_ffff_ffffu64;
    let ay = hy & 0x7fff_ffff_ffff_ffffu64;
    let ax0 = ax;
    let ay0 = ay;

    if ay == 0 || ax >= EXP_MASK || ay > EXP_MASK {
        return (f64::NAN, 0);
    }

    if ax != 0 && ay != 0 {
        let sign_y = if (hy & SIGN_MASK) != 0 { -1i64 } else { 1i64 };
        let ay_f = f64_from_bits(ay);
        let t = x / y;
        if t.is_finite() && t.abs() < 4_503_599_627_370_496.0 {
            let n = rint(t);
            let mut n_i = n as i64;
            let mut r = fma_internal(-n, y, x);
            let ayh = 0.5 * ay_f;
            if r > ayh {
                r -= (sign_y as f64) * y;
                n_i += sign_y;
            } else if r < -ayh {
                r += (sign_y as f64) * y;
                n_i -= sign_y;
            }
            if r == ayh || r == -ayh {
                let want_odd = (quotient_mod8(ax0, ay0) & 1) != 0;
                let is_odd = (n_i.abs() & 1) != 0;
                if want_odd != is_odd {
                    let delta = if r > 0.0 { sign_y } else { -sign_y };
                    r -= (delta as f64) * y;
                    n_i += delta;
                }
            }
            let mut q = (n_i.abs() & 0x7) as i32;
            if n_i < 0 {
                q = -q;
            }
            if r == 0.0 {
                r = if sx != 0 { -0.0 } else { 0.0 };
            }
            return (r, q);
        }
    }

    if ay <= 0x7fbf_ffff_ffff_ffffu64 {
        let r = fmod(x, 8.0 * y);
        hx = r.to_bits();
        ax = hx & 0x7fff_ffff_ffff_ffffu64;
    }

    if ax == ay {
        let quo = if qs != 0 { -1 } else { 1 };
        return (0.0 * x, quo);
    }

    let mut xr = f64_from_bits(ax);
    let yr = f64_from_bits(ay);
    if ay <= 0x7fcf_ffff_ffff_ffffu64 && xr >= 4.0 * yr {
        xr -= 4.0 * yr;
    }
    if ay <= 0x7fdf_ffff_ffff_ffffu64 && xr >= 2.0 * yr {
        xr -= 2.0 * yr;
    }

    if ay < 0x0020_0000_0000_0000u64 {
        if xr + xr > yr {
            xr -= yr;
            if xr + xr >= yr {
                xr -= yr;
            }
        }
    } else {
        let y_half = 0.5 * yr;
        if xr > y_half {
            xr -= yr;
            if xr >= y_half {
                xr -= yr;
            }
        }
    }

    let mut quo = quotient_mod8(ax0, ay0);
    if qs != 0 {
        quo = -quo;
    }
    if xr == 0.0 {
        xr = 0.0;
    }
    if sx != 0 {
        xr = -xr;
    }
    (xr, quo)
}

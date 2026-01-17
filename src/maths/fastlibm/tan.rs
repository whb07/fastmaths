//! tan(x) implementation.
//!
//! Performs argument reduction to [-pi/4, pi/4] and evaluates an odd rational
//! polynomial. Quadrant selection and reciprocal identities handle |x|>pi/4.
//! Constants and tables follow fdlibm/glibc conventions.

use super::lo_word;
use super::trig::branred;
use super::utan_tables::*;

const CN: f64 = 134217729.0; // 1 + 2^27
const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn fma_f64(a: f64, b: f64, c: f64) -> f64 {
    use core::arch::x86_64::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
    _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)))
}

#[inline(always)]
fn fma_opt<const USE_FMA: bool>(a: f64, b: f64, c: f64) -> f64 {
    if USE_FMA {
        #[cfg(target_arch = "x86_64")]
        {
            // SAFETY: caller guards via CPUID and target_feature.
            return unsafe { fma_f64(a, b, c) };
        }
    }
    a * b + c
}

#[inline(always)]
fn eadd(x: f64, y: f64) -> (f64, f64) {
    let z = x + y;
    let zz = if x.abs() > y.abs() {
        (x - z) + y
    } else {
        (y - z) + x
    };
    (z, zz)
}

#[inline(always)]
fn mul12<const USE_FMA: bool>(x: f64, y: f64) -> (f64, f64) {
    let z = x * y;
    if USE_FMA {
        let zz = fma_opt::<true>(x, y, -z);
        return (z, zz);
    }
    let p = CN * x;
    let hx = (x - p) + p;
    let tx = x - hx;
    let p = CN * y;
    let hy = (y - p) + p;
    let ty = y - hy;
    let zz = (((hx * hy - z) + hx * ty) + tx * hy) + tx * ty;
    (z, zz)
}

#[inline(always)]
fn div2<const USE_FMA: bool>(x: f64, xx: f64, y: f64, yy: f64) -> (f64, f64) {
    let c = x / y;
    let (u, uu) = mul12::<USE_FMA>(c, y);
    let cc = (((x - u) - uu) + xx - c * yy) / y;
    let z = c + cc;
    let zz = (c - z) + cc;
    (z, zz)
}

#[inline(always)]
fn tan_poly<const USE_FMA: bool>(a: f64, da: f64, n: i32) -> f64 {
    let a2 = a * a;
    let mut t2 = fma_opt::<USE_FMA>(a2, UTAN_D11, UTAN_D9);
    t2 = fma_opt::<USE_FMA>(a2, t2, UTAN_D7);
    t2 = fma_opt::<USE_FMA>(a2, t2, UTAN_D5);
    t2 = fma_opt::<USE_FMA>(a2, t2, UTAN_D3);
    let a3 = a * a2;
    t2 = fma_opt::<USE_FMA>(a3, t2, da);

    if n != 0 {
        let (b, db) = eadd(a, t2);
        let (c, dc) = div2::<USE_FMA>(1.0, 0.0, b, db);
        let y = c + dc;
        -y
    } else {
        a + t2
    }
}

#[inline(always)]
fn tan_table<const USE_FMA: bool>(ya: f64, yya: f64, sy: f64, n: i32) -> f64 {
    let i = (UTAN_MFFTNHF + 256.0 * ya) as i32;
    let idx = i as usize;
    let z0 = ya - UTAN_XFG[idx][0];
    let z = z0 + yya;
    let z2 = z * z;
    let poly = fma_opt::<USE_FMA>(z2, UTAN_E1, UTAN_E0);
    let pz = fma_opt::<USE_FMA>(z * z2, poly, z);
    let fi = UTAN_XFG[idx][1];
    let gi = UTAN_XFG[idx][2];

    if n != 0 {
        let t2 = pz * (fi + gi) / (fi + pz);
        let y = gi - t2;
        -sy * y
    } else {
        let t2 = pz * (gi + fi) / (gi - pz);
        let y = fi + t2;
        sy * y
    }
}

#[inline(always)]
fn tan_impl<const USE_FMA: bool>(x: f64) -> f64 {
    let xb = x.to_bits();
    if (xb & EXP_MASK) == EXP_MASK {
        return f64::NAN;
    }

    let sign = (xb & SIGN_MASK) != 0;
    let w = f64::from_bits(xb & !SIGN_MASK);
    if w <= UTAN_G1 {
        return x;
    }

    if w <= UTAN_G2 {
        let x2 = x * x;
        let mut t2 = fma_opt::<USE_FMA>(x2, UTAN_D11, UTAN_D9);
        t2 = fma_opt::<USE_FMA>(x2, t2, UTAN_D7);
        t2 = fma_opt::<USE_FMA>(x2, t2, UTAN_D5);
        t2 = fma_opt::<USE_FMA>(x2, t2, UTAN_D3);
        let x3 = x * x2;
        return fma_opt::<USE_FMA>(x3, t2, x);
    }

    if w <= UTAN_G3 {
        let i = (UTAN_MFFTNHF + 256.0 * w) as i32;
        let idx = i as usize;
        let z = w - UTAN_XFG[idx][0];
        let z2 = z * z;
        let poly = fma_opt::<USE_FMA>(z2, UTAN_E1, UTAN_E0);
        let pz = fma_opt::<USE_FMA>(z * z2, poly, z);
        let fi = UTAN_XFG[idx][1];
        let gi = UTAN_XFG[idx][2];
        let t2 = pz * (gi + fi) / (gi - pz);
        let y = fi + t2;
        return if sign { -y } else { y };
    }

    if w <= UTAN_G4 {
        let t = x * UTAN_HPINV + UTAN_TOINT;
        let xn = t - UTAN_TOINT;
        let n = (lo_word(t) & 0x0000_0001) as i32;
        let t1 = (x - xn * UTAN_MP1) - xn * UTAN_MP2;
        let mut da = xn * UTAN_MP3;
        let a = t1 - da;
        da = (t1 - a) - da;

        let (ya, yya, sy) = if a < 0.0 {
            (-a, -da, -1.0)
        } else {
            (a, da, 1.0)
        };

        if ya <= UTAN_GY2 {
            return tan_poly::<USE_FMA>(a, da, n);
        }

        return tan_table::<USE_FMA>(ya, yya, sy, n);
    }

    if w <= UTAN_G5 {
        let t = x * UTAN_HPINV + UTAN_TOINT;
        let xn = t - UTAN_TOINT;
        let n = (lo_word(t) & 0x0000_0001) as i32;
        let t1 = (x - xn * UTAN_MP1) - xn * UTAN_MP2;
        let mut da = xn * UTAN_PP3;
        let t = t1 - da;
        da = (t1 - t) - da;
        let t1 = xn * UTAN_PP4;
        let a = t - t1;
        let da = ((t - a) - t1) + da;
        let (a, da) = eadd(a, da);

        let (ya, yya, sy) = if a < 0.0 {
            (-a, -da, -1.0)
        } else {
            (a, da, 1.0)
        };

        if ya <= UTAN_GY2 {
            return tan_poly::<USE_FMA>(a, da, n);
        }

        return tan_table::<USE_FMA>(ya, yya, sy, n);
    }

    let (n_full, a_raw, da_raw) = branred(x);
    let n = n_full & 0x0000_0001;
    let (a, da) = eadd(a_raw, da_raw);

    let (ya, yya, sy) = if a < 0.0 {
        (-a, -da, -1.0)
    } else {
        (a, da, 1.0)
    };

    if ya <= UTAN_GY2 {
        return tan_poly::<USE_FMA>(a, da, n);
    }

    tan_table::<USE_FMA>(ya, yya, sy, n)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn tan_fma(x: f64) -> f64 {
    tan_impl::<true>(x)
}

#[inline(always)]
fn tan_generic(x: f64) -> f64 {
    tan_impl::<false>(x)
}

#[inline(always)]
pub fn tan(x: f64) -> f64 {
    #[cfg(target_arch = "x86_64")]
    if super::cpu_has_fma() {
        // SAFETY: guarded by CPUID.
        return unsafe { tan_fma(x) };
    }
    tan_generic(x)
}

use super::trig::branred;
use super::utan_tables::*;
use super::{hi_word, lo_word};

const CN: f64 = 134217729.0; // 1 + 2^27

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn fma_f64(a: f64, b: f64, c: f64) -> f64 {
    use core::arch::x86_64::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
    _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)))
}

#[inline(always)]
fn fma_fast(a: f64, b: f64, c: f64) -> f64 {
    #[cfg(target_arch = "x86_64")]
    if super::cpu_has_fma() {
        // SAFETY: guarded by CPUID.
        return unsafe { fma_f64(a, b, c) };
    }
    a * b + c
}

#[inline]
fn eadd(x: f64, y: f64) -> (f64, f64) {
    let z = x + y;
    let zz = if x.abs() > y.abs() {
        (x - z) + y
    } else {
        (y - z) + x
    };
    (z, zz)
}

#[inline]
fn mul12(x: f64, y: f64) -> (f64, f64) {
    let z = x * y;
    if super::cpu_has_fma() {
        let zz = fma_fast(x, y, -z);
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

#[inline]
fn div2(x: f64, xx: f64, y: f64, yy: f64) -> (f64, f64) {
    let c = x / y;
    let (u, uu) = mul12(c, y);
    let cc = (((x - u) - uu) + xx - c * yy) / y;
    let z = c + cc;
    let zz = (c - z) + cc;
    (z, zz)
}

#[inline]
fn tan_poly(a: f64, da: f64, n: i32) -> f64 {
    let a2 = a * a;
    let mut t2 = UTAN_D9 + a2 * UTAN_D11;
    t2 = UTAN_D7 + a2 * t2;
    t2 = UTAN_D5 + a2 * t2;
    t2 = UTAN_D3 + a2 * t2;
    t2 = da + a * a2 * t2;

    if n != 0 {
        let (b, db) = eadd(a, t2);
        let (c, dc) = div2(1.0, 0.0, b, db);
        let y = c + dc;
        -y
    } else {
        a + t2
    }
}

#[inline]
fn tan_table(ya: f64, yya: f64, sy: f64, n: i32) -> f64 {
    let i = (UTAN_MFFTNHF + 256.0 * ya) as i32;
    let idx = i as usize;
    let z0 = ya - UTAN_XFG[idx][0];
    let z = z0 + yya;
    let z2 = z * z;
    let pz = z + z * z2 * (UTAN_E0 + z2 * UTAN_E1);
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

#[inline]
pub fn tan(x: f64) -> f64 {
    let ux = hi_word(x);
    if (ux & 0x7ff0_0000) == 0x7ff0_0000 {
        return f64::NAN;
    }

    let w = x.abs();
    if w <= UTAN_G1 {
        return x;
    }

    if w <= UTAN_G2 {
        let x2 = x * x;
        let mut t2 = UTAN_D9 + x2 * UTAN_D11;
        t2 = UTAN_D7 + x2 * t2;
        t2 = UTAN_D5 + x2 * t2;
        t2 = UTAN_D3 + x2 * t2;
        t2 *= x * x2;
        return x + t2;
    }

    if w <= UTAN_G3 {
        let i = (UTAN_MFFTNHF + 256.0 * w) as i32;
        let idx = i as usize;
        let z = w - UTAN_XFG[idx][0];
        let z2 = z * z;
        let pz = z + z * z2 * (UTAN_E0 + z2 * UTAN_E1);
        let fi = UTAN_XFG[idx][1];
        let gi = UTAN_XFG[idx][2];
        let t2 = pz * (gi + fi) / (gi - pz);
        let y = fi + t2;
        return if x < 0.0 { -y } else { y };
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
            return tan_poly(a, da, n);
        }

        return tan_table(ya, yya, sy, n);
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
            return tan_poly(a, da, n);
        }

        return tan_table(ya, yya, sy, n);
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
        return tan_poly(a, da, n);
    }

    tan_table(ya, yya, sy, n)
}

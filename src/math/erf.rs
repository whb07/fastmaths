//! erf/erfc implementation (core-math/glibc port).
//!
//! Piecewise approximations with table-driven argument reduction and
//! double-double arithmetic for correctly rounded results (<= 1 ULP).
//! Constants and tables are sourced from glibc/core-math (see glibc/).

use super::{copysign, floor, fma_internal, ldexp, rint};
use super::{erf_data, erfc_data};

const SIGN_MASK: u64 = 0x8000_0000_0000_0000;

#[inline(always)]
fn fma(a: f64, b: f64, c: f64) -> f64 {
    fma_internal(a, b, c)
}

#[inline(always)]
fn asuint64(x: f64) -> u64 {
    x.to_bits()
}

#[inline(always)]
fn asdouble(x: u64) -> f64 {
    f64::from_bits(x)
}

// Add a + b, such that *hi + *lo approximates a + b (|a| >= |b|).
#[inline(always)]
fn fast_two_sum(hi: &mut f64, lo: &mut f64, a: f64, b: f64) {
    let s = a + b;
    let e = s - a;
    *hi = s;
    *lo = b - e;
}

#[inline(always)]
fn two_sum(hi: &mut f64, lo: &mut f64, a: f64, b: f64) {
    let s = a + b;
    let aa = s - b;
    let bb = s - aa;
    let da = a - aa;
    let db = b - bb;
    *hi = s;
    *lo = da + db;
}

// Multiply exactly a and b, such that *hi + *lo = a * b.
#[inline(always)]
fn a_mul(hi: &mut f64, lo: &mut f64, a: f64, b: f64) {
    *hi = a * b;
    *lo = fma(a, b, -*hi);
}

// Multiply a double with a double-double: a * (bh + bl).
#[inline(always)]
fn s_mul(hi: &mut f64, lo: &mut f64, a: f64, bh: f64, bl: f64) {
    a_mul(hi, lo, a, bh);
    *lo = fma(a, bl, *lo);
}

// Returns (ah + al) * (bh + bl) - (al * bl).
#[inline(always)]
fn d_mul(hi: &mut f64, lo: &mut f64, ah: f64, al: f64, bh: f64, bl: f64) {
    a_mul(hi, lo, ah, bh);
    *lo = fma(ah, bl, *lo);
    *lo = fma(al, bh, *lo);
}

// Add a + (bh + bl), assuming |a| >= |bh|.
#[inline(always)]
fn fast_sum(hi: &mut f64, lo: &mut f64, a: f64, bh: f64, bl: f64) {
    fast_two_sum(hi, lo, a, bh);
    *lo += bl;
}

#[inline(always)]
fn roundeven_finite(x: f64) -> f64 {
    rint(x)
}

// Assuming 0 <= z <= f64::from_bits(0x4017afb48dc96626), put in h+l an approximation of erf(z).
// Return err: |(h + l)/erf(z) - 1| < err*|h+l|.
fn cr_erf_fast(h: &mut f64, l: &mut f64, mut z: f64) -> f64 {
    if z < 0.0625 {
        let mut z2h = 0.0;
        let mut z2l = 0.0;
        a_mul(&mut z2h, &mut z2l, z, z);
        let z4 = z2h * z2h;
        let c9 = fma(erf_data::C0[7], z2h, erf_data::C0[6]);
        let c5 = fma(erf_data::C0[5], z2h, erf_data::C0[4]);
        let c5 = fma(c9, z4, c5);
        let mut th = 0.0;
        let mut tl = 0.0;
        a_mul(&mut th, &mut tl, z2h, c5);
        fast_two_sum(h, l, erf_data::C0[2], th);
        *l += tl + erf_data::C0[3];
        let h_copy = *h;
        a_mul(&mut th, &mut tl, z2h, *h);
        tl += fma(z2h, *l, erf_data::C0[1]);
        fast_two_sum(h, l, erf_data::C0[0], th);
        *l += fma(z2l, h_copy, tl);
        a_mul(h, &mut tl, *h, z);
        *l = fma(*l, z, tl);
        return f64::from_bits(0x3ba7800000000000);
    }

    let v = floor(16.0 * z);
    let i = (16.0 * z) as usize;
    z = (z - 0.03125) - 0.0625 * v;
    let c = &erf_data::C[i - 1];
    let z2 = z * z;
    let z4 = z2 * z2;
    let c9 = fma(c[12], z, c[11]);
    let c7 = fma(c[10], z, c[9]);
    let c5 = fma(c[8], z, c[7]);
    let mut c3h = 0.0;
    let mut c3l = 0.0;
    fast_two_sum(&mut c3h, &mut c3l, c[5], z * c[6]);
    let c7 = fma(c9, z2, c7);
    let mut tl = 0.0;
    let c3h0 = c3h;
    fast_two_sum(&mut c3h, &mut tl, c3h0, c5 * z2);
    c3l += tl;
    let c3h1 = c3h;
    fast_two_sum(&mut c3h, &mut tl, c3h1, c7 * z4);
    c3l += tl;
    let mut th = 0.0;
    let mut tl2 = 0.0;
    a_mul(&mut th, &mut tl2, z, c3h);
    let mut c2h = 0.0;
    let mut c2l = 0.0;
    fast_two_sum(&mut c2h, &mut c2l, c[4], th);
    c2l += fma(z, c3l, tl2);
    a_mul(&mut th, &mut tl2, z, c2h);
    fast_two_sum(h, l, c[2], th);
    *l += tl2 + fma(z, c2l, c[3]);
    a_mul(&mut th, &mut tl2, z, *h);
    tl2 = fma(z, *l, tl2);
    fast_two_sum(h, l, c[0], th);
    *l += tl2 + c[1];
    f64::from_bits(0x3ba1100000000000)
}

// For |z| < 1/8, assuming z >= 2^-61, thus no underflow can occur.
fn cr_erf_accurate_tiny(h: &mut f64, l: &mut f64, z: f64, exceptions: bool) {
    if exceptions {
        let mut i = 0usize;
        let mut j = erf_data::EXCEPTIONS_TINY.len();
        while i + 1 < j {
            let k = (i + j) / 2;
            if erf_data::EXCEPTIONS_TINY[k][0] <= z {
                i = k;
            } else {
                j = k;
            }
        }
        if z == erf_data::EXCEPTIONS_TINY[i][0] {
            *h = erf_data::EXCEPTIONS_TINY[i][1];
            *l = erf_data::EXCEPTIONS_TINY[i][2];
            return;
        }
    }

    let z2 = z * z;
    let mut th = 0.0;
    let mut tl = 0.0;
    *h = erf_data::P[21 / 2 + 4];
    for a in (13..=19).rev().step_by(2) {
        *h = fma(*h, z2, erf_data::P[a / 2 + 4]);
    }
    *l = 0.0;
    for a in (9..=11).rev().step_by(2) {
        a_mul(&mut th, &mut tl, *h, z);
        tl = fma(*l, z, tl);
        a_mul(h, l, th, z);
        *l = fma(tl, z, *l);
        fast_two_sum(h, &mut tl, erf_data::P[a / 2 + 4], *h);
        *l += tl;
    }
    for a in (1..=7).rev().step_by(2) {
        a_mul(&mut th, &mut tl, *h, z);
        tl = fma(*l, z, tl);
        a_mul(h, l, th, z);
        *l = fma(tl, z, *l);
        fast_two_sum(h, &mut tl, erf_data::P[a - 1], *h);
        *l += erf_data::P[a] + tl;
    }
    a_mul(h, &mut tl, *h, z);
    *l = fma(*l, z, tl);
}

// Accurate erf approximation for 0 <= z <= f64::from_bits(0x4017afb48dc96626).
#[cold]
#[inline(never)]
fn cr_erf_accurate(h: &mut f64, l: &mut f64, z: f64) {
    for ex in erf_data::EXCEPTIONS.iter() {
        if z == ex[0] {
            *h = ex[1];
            *l = ex[2];
            return;
        }
    }
    if z < 0.125 {
        cr_erf_accurate_tiny(h, l, z, true);
        return;
    }
    let v = floor(8.0 * z);
    let i = (8.0 * z) as usize;
    let zz = (z - 0.0625) - 0.125 * v;
    let p = &erf_data::C2[i - 1];
    *h = p[26];
    for j in (11..=17).rev() {
        *h = fma(*h, zz, p[8 + j]);
    }
    *l = 0.0;
    let mut th = 0.0;
    let mut tl = 0.0;
    for j in (8..=10).rev() {
        a_mul(&mut th, &mut tl, *h, zz);
        tl = fma(*l, zz, tl);
        two_sum(h, l, p[8 + j], th);
        *l += tl;
    }
    for j in (0..=7).rev() {
        a_mul(&mut th, &mut tl, *h, zz);
        tl = fma(*l, zz, tl);
        two_sum(h, l, p[2 * j], th);
        *l += p[2 * j + 1] + tl;
    }
}

#[inline(always)]
pub fn erf(x: f64) -> f64 {
    let z = x.abs();
    let ux = asuint64(z);
    if ux > 0x4017_afb4_8dc9_6626u64 {
        let os = copysign(1.0, x);
        if ux > 0x7ff0_0000_0000_0000u64 {
            return x + x;
        }
        if ux == 0x7ff0_0000_0000_0000u64 {
            return os;
        }
        return os - f64::from_bits(0x3c90000000000000) * os;
    }
    if z < f64::from_bits(0x3c20000000000000) {
        if x == 0.0 {
            return x;
        }
        let y = erf_data::C0[0] * x;
        let sx = x * f64::from_bits(0x4690000000000000);
        let mut h = 0.0;
        let mut l = 0.0;
        a_mul(&mut h, &mut l, erf_data::C0[0], sx);
        l = fma(erf_data::C0[1], sx, l);
        l += h - y * f64::from_bits(0x4690000000000000);
        return fma(l, f64::from_bits(0x3950000000000000), y);
    }

    let mut h = 0.0;
    let mut l = 0.0;
    let err = cr_erf_fast(&mut h, &mut l, z);
    let t = asuint64(x);
    let mut u = asuint64(h);
    let mut v = asuint64(l);
    u ^= t & SIGN_MASK;
    v ^= t & SIGN_MASK;
    let uf = asdouble(u);
    let vf = asdouble(v);
    let left = uf + fma(err, -uf, vf);
    let right = uf + fma(err, uf, vf);
    if left == right {
        return left;
    }
    cr_erf_accurate(&mut h, &mut l, z);
    if x >= 0.0 { h + l } else { (-h) + (-l) }
}

// Approximation for exp(x) with x = xh + xl, small magnitude.
#[inline(always)]
fn q_1(hi: &mut f64, lo: &mut f64, zh: f64, zl: f64) {
    let z = zh + zl;
    let mut q = fma(erf_data::Q_1[4], zh, erf_data::Q_1[3]);
    q = fma(q, z, erf_data::Q_1[2]);
    fast_two_sum(hi, lo, erf_data::Q_1[1], q * z);
    d_mul(hi, lo, zh, zl, *hi, *lo);
    fast_sum(hi, lo, erf_data::Q_1[0], *hi, *lo);
}

// Approximation of exp(x) where x = xh + xl.
#[inline(always)]
fn exp_1(hi: &mut f64, lo: &mut f64, xh: f64, xl: f64) {
    const INVLOG2: f64 = f64::from_bits(0x40b71547652b82fe);
    let k = roundeven_finite(xh * INVLOG2);
    const LOG2H: f64 = f64::from_bits(0x3f262e42fefa39ef);
    const LOG2L: f64 = f64::from_bits(0x3bbabc9e3b39803f);
    let mut kh = 0.0;
    let mut kl = 0.0;
    s_mul(&mut kh, &mut kl, k, LOG2H, LOG2L);
    let mut yh = 0.0;
    let mut yl = 0.0;
    fast_two_sum(&mut yh, &mut yl, xh - kh, xl);
    yl -= kl;
    let k_i = k as i64;
    let m = (k_i >> 12) + 0x3ff;
    let i2 = ((k_i >> 6) & 0x3f) as usize;
    let i1 = (k_i & 0x3f) as usize;
    let t1h = erf_data::T1[i2][0];
    let t1l = erf_data::T1[i2][1];
    let t2h = erf_data::T2[i1][0];
    let t2l = erf_data::T2[i1][1];
    d_mul(hi, lo, t2h, t2l, t1h, t1l);
    let mut qh = 0.0;
    let mut ql = 0.0;
    q_1(&mut qh, &mut ql, yh, yl);
    d_mul(hi, lo, *hi, *lo, qh, ql);
    let df = asdouble((m as u64) << 52);
    *hi *= df;
    *lo *= df;
}

// Put in 2^e*(h+l) an approximation of exp(xh+xl) for -742 <= xh+xl <= -2.92.
#[cold]
#[inline(never)]
fn exp_accurate(h: &mut f64, l: &mut f64, e: &mut i32, xh: f64, xl: f64) {
    const INVLOG2: f64 = f64::from_bits(0x3ff71547652b82fe);
    let k = roundeven_finite(xh * INVLOG2) as i32;
    const LOG2H: f64 = f64::from_bits(0x3fe62e42fefa39ef);
    const LOG2L: f64 = f64::from_bits(0x3c7abc9e3b398000);
    const LOG2T: f64 = f64::from_bits(0x398f97b57a079a19);
    let mut yh = fma(-(k as f64), LOG2H, xh);
    let mut th = 0.0;
    let mut tl = 0.0;
    two_sum(&mut th, &mut tl, -(k as f64) * LOG2L, xl);
    let mut yl = 0.0;
    let yh0 = yh;
    fast_two_sum(&mut yh, &mut yl, yh0, th);
    yl = fma(-(k as f64), LOG2T, yl + tl);
    *h = erfc_data::E2[19 + 8];
    for i in (16..=18).rev() {
        *h = fma(*h, yh, erfc_data::E2[i + 8]);
    }
    a_mul(&mut th, &mut tl, *h, yh);
    tl = fma(*h, yl, tl);
    fast_two_sum(h, l, erfc_data::E2[15 + 8], th);
    *l += tl;
    for i in (8..=14).rev() {
        a_mul(&mut th, &mut tl, *h, yh);
        tl = fma(*h, yl, tl);
        tl = fma(*l, yh, tl);
        fast_two_sum(h, l, erfc_data::E2[i + 8], th);
        *l += tl;
    }
    for i in (0..=7).rev() {
        a_mul(&mut th, &mut tl, *h, yh);
        tl = fma(*h, yl, tl);
        tl = fma(*l, yh, tl);
        fast_two_sum(h, l, erfc_data::E2[2 * i], th);
        *l += tl + erfc_data::E2[2 * i + 1];
    }
    *e = k;
}

// Fast asymptotic erfc for large x.
#[inline(never)]
fn erfc_asympt_fast(h: &mut f64, l: &mut f64, x: f64) -> f64 {
    if x >= f64::from_bits(0x4039db1bb14e15ca) {
        *h = 0.0;
        *l = 0.0;
        return 1.0;
    }
    let mut eh = 0.0;
    let mut el = 0.0;
    let mut uh = 0.0;
    let mut ul = 0.0;
    a_mul(&mut uh, &mut ul, x, x);
    exp_1(&mut eh, &mut el, -uh, -ul);
    let yh = 1.0 / x;
    let yl = yh * fma(-x, yh, 1.0);
    const THRESHOLD: [f64; 6] = [
        f64::from_bits(0x3fbd500000000000),
        f64::from_bits(0x3fc59da6ca291ba6),
        f64::from_bits(0x3fcbc00000000000),
        f64::from_bits(0x3fd0c00000000000),
        f64::from_bits(0x3fd3800000000000),
        f64::from_bits(0x3fd6300000000000),
    ];
    let mut i = 0usize;
    while i < THRESHOLD.len() && yh > THRESHOLD[i] {
        i += 1;
    }
    let p = &erfc_data::T[i];
    a_mul(&mut uh, &mut ul, yh, yh);
    ul = fma(2.0 * yh, yl, ul);
    let mut zh = p[12];
    zh = fma(zh, uh, p[11]);
    zh = fma(zh, uh, p[10]);
    s_mul(h, l, zh, uh, ul);
    let mut zl = 0.0;
    fast_two_sum(&mut zh, &mut zl, p[9], *h);
    zl += *l;
    for j in (3usize..=15).rev().step_by(2) {
        d_mul(h, l, zh, zl, uh, ul);
        fast_two_sum(&mut zh, &mut zl, p[j.div_ceil(2)], *h);
        zl += *l;
    }
    d_mul(h, l, zh, zl, uh, ul);
    fast_two_sum(&mut zh, &mut zl, p[0], *h);
    zl += *l + p[1];
    d_mul(&mut uh, &mut ul, zh, zl, yh, yl);
    d_mul(h, l, uh, ul, eh, el);
    if *h >= f64::from_bits(0x044151b9a3fdd5c9) {
        return f64::from_bits(0x3bbd900000000000) * *h;
    }
    f64::from_bits(0x0010000000000000)
}

// Fast path for erfc, returning absolute error bound.
fn cr_erfc_fast(h: &mut f64, l: &mut f64, x: f64) -> f64 {
    if x < 0.0 {
        let mut err = cr_erf_fast(h, l, -x);
        err *= *h;
        let mut t = 0.0;
        fast_two_sum(h, &mut t, 1.0, *h);
        *l += t;
        return err + f64::from_bits(0x3994000000000000);
    }
    const THRESHOLD1: f64 = f64::from_bits(0x400713786d9c7c09);
    if x <= THRESHOLD1 {
        let mut err = cr_erf_fast(h, l, x);
        err *= *h;
        let mut t = 0.0;
        fast_two_sum(h, &mut t, 1.0, -*h);
        *l = t - *l;
        if x >= f64::from_bits(0x3fde861fbb24c00a) {
            return err;
        }
        return err + f64::from_bits(0x3974000000000000);
    }
    erfc_asympt_fast(h, l, x)
}

// Accurate asymptotic erfc for larger x.
#[cold]
#[inline(never)]
fn erfc_asympt_accurate(x: f64) -> f64 {
    for ex in erfc_data::EXCEPTIONS.iter() {
        if x == ex[0] {
            return ex[1] + ex[2];
        }
    }
    if x == f64::from_bits(0x403a8f7bfbd15495) {
        return fma(
            f64::from_bits(0x0000000000000001),
            -0.25,
            f64::from_bits(0x000667bd620fd95b),
        );
    }
    let mut h = 0.0;
    let mut l = 0.0;
    let mut eh = 0.0;
    let mut el = 0.0;
    let mut uh = 0.0;
    let mut ul = 0.0;
    a_mul(&mut uh, &mut ul, x, x);
    let mut e = 0i32;
    exp_accurate(&mut eh, &mut el, &mut e, -uh, -ul);
    let yh = 1.0 / x;
    let yl = yh * fma(-x, yh, 1.0);
    const THRESHOLD: [f64; 10] = [
        f64::from_bits(0x3fb4500000000000),
        f64::from_bits(0x3fbe000000000000),
        f64::from_bits(0x3fc3f00000000000),
        f64::from_bits(0x3fc9500000000000),
        f64::from_bits(0x3fcf500000000000),
        f64::from_bits(0x3fd3100000000000),
        f64::from_bits(0x3fd7100000000000),
        f64::from_bits(0x3fdbc00000000000),
        f64::from_bits(0x3fe0b00000000000),
        f64::from_bits(0x3fe3000000000000),
    ];
    let mut i = 0usize;
    while i < THRESHOLD.len() && yh > THRESHOLD[i] {
        i += 1;
    }
    let p = &erfc_data::TACC[i];
    a_mul(&mut uh, &mut ul, yh, yh);
    ul = fma(2.0 * yh, yl, ul);
    let mut zh = p[14 + 6 + i];
    let mut zl = 0.0;
    let mut th = 0.0;
    let mut tl = 0.0;
    let top = 27 + 2 * i;
    let mut j = top as i32;
    while j >= 13 {
        a_mul(&mut th, &mut tl, zh, uh);
        tl = fma(zh, ul, tl);
        tl = fma(zl, uh, tl);
        two_sum(&mut zh, &mut zl, p[((j - 1) / 2) as usize + 6], th);
        zl += tl;
        j -= 2;
    }
    let mut j2 = 11;
    while j2 >= 1 {
        a_mul(&mut th, &mut tl, zh, uh);
        tl = fma(zh, ul, tl);
        tl = fma(zl, uh, tl);
        two_sum(&mut zh, &mut zl, p[(j2 - 1) as usize], th);
        zl += tl + p[j2 as usize];
        j2 -= 2;
    }
    a_mul(&mut uh, &mut ul, zh, yh);
    ul = fma(zh, yl, ul);
    ul = fma(zl, yh, ul);
    let uh0 = uh;
    let ul0 = ul;
    fast_two_sum(&mut uh, &mut ul, uh0, ul0);
    a_mul(&mut h, &mut l, uh, eh);
    l = fma(uh, el, l);
    l = fma(ul, eh, l);
    let mut res = ldexp(h + l, e);
    if res < f64::from_bits(0x0010000000000000) {
        let mut corr = h - ldexp(res, -e);
        corr += l;
        res += ldexp(corr, e);
    }
    res
}

#[cold]
#[inline(never)]
fn cr_erfc_accurate(x: f64) -> f64 {
    let mut h = 0.0;
    let mut l = 0.0;
    let mut t = 0.0;
    if x < 0.0 {
        for ex in erfc_data::EXCEPTIONS_ACCURATE.iter() {
            if x == ex[0] {
                return ex[1] + ex[2];
            }
        }
        cr_erf_accurate(&mut h, &mut l, -x);
        let h0 = h;
        fast_two_sum(&mut h, &mut t, 1.0, h0);
        l += t;
        return h + l;
    }
    if x <= f64::from_bits(0x3ffb59ffb450828c) {
        for ex in erfc_data::EXCEPTIONS_ACCURATE_2.iter() {
            if x == ex[0] {
                return ex[1] + ex[2];
            }
        }
        cr_erf_accurate(&mut h, &mut l, x);
        let h0 = h;
        fast_two_sum(&mut h, &mut t, 1.0, -h0);
        l = t - l;
        return h + l;
    }
    erfc_asympt_accurate(x)
}

#[inline(always)]
pub fn erfc(x: f64) -> f64 {
    let t = asuint64(x);
    let at = t & 0x7fff_ffff_ffff_ffffu64;
    if t >= 0x8000_0000_0000_0000u64 {
        if t >= 0xc017_744f_8f74_e94bu64 {
            if t >= 0xfff0_0000_0000_0000u64 {
                if t == 0xfff0_0000_0000_0000u64 {
                    return 2.0;
                }
                return x + x;
            }
            return 2.0 - f64::from_bits(0x3c90000000000000);
        }
        if f64::from_bits(0xbc9c5bf891b4ef6a) <= x {
            return fma(-x, f64::from_bits(0x3c90000000000000), 1.0);
        }
    } else {
        if at >= 0x403b_39dc_41e4_8bfdu64 {
            if at >= 0x7ff0_0000_0000_0000u64 {
                if at == 0x7ff0_0000_0000_0000u64 {
                    return 0.0;
                }
                return x + x;
            }
            return f64::from_bits(0x0000000000000001) * 0.25;
        }
        if x <= f64::from_bits(0x3c8c5bf891b4ef6a) {
            return fma(-x, f64::from_bits(0x3c90000000000000), 1.0);
        }
    }
    let mut h = 0.0;
    let mut l = 0.0;
    let err = cr_erfc_fast(&mut h, &mut l, x);
    let left = h + (l - err);
    let right = h + (l + err);
    if left == right {
        return left;
    }
    cr_erfc_accurate(x)
}

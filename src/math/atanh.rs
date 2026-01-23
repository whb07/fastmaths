//! atanh(x) implementation.
//!
//! Small |x| uses the core-math polynomial with double-double correction.
//! For the remaining range, use log1p-difference or log-ratio to keep
//! <= 1 ULP accuracy with good throughput.

use super::atanh_data::{CH, CL};
use super::{asdouble, copysign, fasttwosum, fma_wrap};
use super::{log::ln, log1p};

const TINY_BITS: u64 = 0x3e4d_12ed_0af1_a27f; // 0x1.d12ed0af1a27fp-27
const SMALL_BITS: u64 = 0x3fd0_0000_0000_0000; // 0x1p-2
const ONE_BITS: u64 = 0x3ff0_0000_0000_0000;
const INF_BITS: u64 = 0x7ff0_0000_0000_0000;
const LOG1P_BOUND: f64 = 0.5;

const SMALL_C: [f64; 9] = [
    f64::from_bits(0x3fc9_9999_9999_999a),
    f64::from_bits(0x3fc2_4924_9249_2244),
    f64::from_bits(0x3fbc_71c7_1c79_715f),
    f64::from_bits(0x3fb7_45d1_6f77_7723),
    f64::from_bits(0x3fb3_b13c_a417_4634),
    f64::from_bits(0x3fb1_10c9_7249_89bd),
    f64::from_bits(0x3fae_2d17_608a_5b2e),
    f64::from_bits(0x3faa_0b56_308c_ba0b),
    f64::from_bits(0x3faf_b634_1208_ad2e),
];

const SMALL_T: f64 = f64::from_bits(0x3c75_5555_5555_5555); // 0x1.5555555555555p-56
const SMALL_BASE: f64 = f64::from_bits(0x3fd5_5555_5555_5555); // 0x1.5555555555555p-2
const EPS_X4: f64 = f64::from_bits(0x3cad_0000_0000_0000); // 0x1.dp-53
const EPS_TINY: f64 = f64::from_bits(0x3980_0000_0000_0000); // 0x1p-103

#[inline(always)]
fn asuint64(x: f64) -> u64 {
    x.to_bits()
}

#[inline(always)]
fn muldd_acc(xh: f64, xl: f64, ch: f64, cl: f64, l: &mut f64) -> f64 {
    let ahlh = ch * xl;
    let alhh = cl * xh;
    let ahhh = ch * xh;
    let mut ahhl = fma_wrap(ch, xh, -ahhh);
    ahhl += alhh + ahlh;
    let chh = ahhh + ahhl;
    *l = (ahhh - chh) + ahhl;
    chh
}

#[inline(always)]
fn muldd_acc2(xh: f64, xl: f64, ch: f64, cl: f64, l: &mut f64) -> f64 {
    let ahlh = ch * xl;
    let alhh = cl * xh;
    let ahhh = ch * xh;
    let mut ahhl = fma_wrap(ch, xh, -ahhh);
    ahhl += alhh + ahlh;
    fasttwosum(ahhh, ahhl, l)
}

#[inline(always)]
fn mulddd3(xh: f64, xl: f64, ch: f64, l: &mut f64) -> f64 {
    let hh = xh * ch;
    *l = fma_wrap(ch, xh, -hh) + xl * ch;
    hh
}

#[inline(always)]
fn polydd3(xh: f64, xl: f64, n: usize, c: &[[f64; 2]], l: &mut f64) -> f64 {
    let mut i = n - 1;
    let mut cl = 0.0;
    let mut ch = fasttwosum(c[i][0], *l, &mut cl);
    cl += c[i][1];
    while i > 0 {
        i -= 1;
        ch = muldd_acc2(xh, xl, ch, cl, &mut cl);
        let mut tl = 0.0;
        let th = fasttwosum(c[i][0], ch, &mut tl);
        ch = th;
        cl += tl + c[i][1];
    }
    *l = cl;
    ch
}

#[inline(always)]
fn as_atanh_zero(x: f64) -> f64 {
    let x2 = x * x;
    let x2l = fma_wrap(x, x, -x2);
    let y2 = x2 * (CL[0] + x2 * (CL[1] + x2 * (CL[2] + x2 * (CL[3] + x2 * (CL[4])))));
    let mut y2v = y2;
    let mut y1 = polydd3(x2, x2l, CH.len(), CH, &mut y2v);
    y1 = mulddd3(y1, y2v, x, &mut y2v);
    y1 = muldd_acc2(x2, x2l, y1, y2v, &mut y2v);
    let mut y1e = 0.0;
    let y0 = fasttwosum(x, y1, &mut y1e);
    let mut y2e = 0.0;
    let mut y1h = fasttwosum(y1e, y2v, &mut y2e);
    let mut t = asuint64(y1h);
    if (t & (!0u64 >> 12)) == 0 {
        let w = asuint64(y2e);
        if ((w ^ t) >> 63) != 0 {
            t = t.wrapping_sub(1);
        } else {
            t = t.wrapping_add(1);
        }
        y1h = asdouble(t);
    }
    y0 + y1h
}

#[inline(always)]
pub fn atanh(x: f64) -> f64 {
    let ax = x.abs();
    let aix = asuint64(ax);
    if aix >= ONE_BITS {
        if aix == ONE_BITS {
            return if x.is_sign_negative() {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }
        if aix > INF_BITS {
            return x + x;
        }
        return f64::NAN;
    }

    if aix < SMALL_BITS {
        if aix < TINY_BITS {
            return fma_wrap(x, f64::from_bits(0x3c80_0000_0000_0000), x);
        }
        let x2 = x * x;
        let dx2 = fma_wrap(x, x, -x2);
        let x4 = x2 * x2;
        let x3 = x2 * x;
        let x8 = x4 * x4;
        let dx3 = fma_wrap(x2, x, -x3) + dx2 * x;
        let p = (SMALL_C[0] + x2 * SMALL_C[1])
            + x4 * (SMALL_C[2] + x2 * SMALL_C[3])
            + x8 * ((SMALL_C[4] + x2 * SMALL_C[5])
                + x4 * (SMALL_C[6] + x2 * SMALL_C[7])
                + x8 * SMALL_C[8]);
        let t = fma_wrap(x2, p, SMALL_T);
        let mut pl = 0.0;
        let mut ph = fasttwosum(SMALL_BASE, t, &mut pl);
        ph = muldd_acc(ph, pl, x3, dx3, &mut pl);
        let mut tl = 0.0;
        ph = fasttwosum(x, ph, &mut tl);
        pl += tl;
        let eps = x * (x4 * EPS_X4 + EPS_TINY);
        let lb = ph + (pl - eps);
        let ub = ph + (pl + eps);
        if lb == ub {
            return lb;
        }
        return as_atanh_zero(x);
    }

    if ax <= LOG1P_BOUND {
        let y = 0.5 * (log1p(ax) - log1p(-ax));
        return copysign(y, x);
    }

    let ratio = (1.0 + ax) / (1.0 - ax);
    let y = 0.5 * ln(ratio);
    copysign(y, x)
}

#![allow(
    clippy::collapsible_if,
    clippy::needless_late_init,
    clippy::needless_range_loop,
    dead_code
)]

use super::sincos_tab::SINCOS_TAB;
use super::{floor_f64, hi_word, lo_word, scalbn_internal, with_hi_lo};

// ========= fdlibm/glibc-grade sin/cos =========
//
// Implements:
// - __ieee754_rem_pio2 (argument reduction)
// - __kernel_rem_pio2 (Payne–Hanek core, 24-bit limbs)
// - __kernel_sin / __kernel_cos (accurate kernels with tail)

const TWO_OVER_PI: [u32; 66] = [
    0xa2f983u32,
    0x6e4e44u32,
    0x1529fcu32,
    0x2757d1u32,
    0xf534ddu32,
    0xc0db62u32,
    0x95993cu32,
    0x439041u32,
    0xfe5163u32,
    0xabdebbu32,
    0xc561b7u32,
    0x246e3au32,
    0x424dd2u32,
    0xe00649u32,
    0x2eea09u32,
    0xd1921cu32,
    0xfe1debu32,
    0x1cb129u32,
    0xa73ee8u32,
    0x8235f5u32,
    0x2ebb44u32,
    0x84e99cu32,
    0x7026b4u32,
    0x5f7e41u32,
    0x3991d6u32,
    0x398353u32,
    0x39f49cu32,
    0x845f8bu32,
    0xbdf928u32,
    0x3b1ff8u32,
    0x97ffdeu32,
    0x05980fu32,
    0xef2f11u32,
    0x8b5a0au32,
    0x6d1f6du32,
    0x367ecfu32,
    0x27cb09u32,
    0xb74f46u32,
    0x3f669eu32,
    0x5fea2du32,
    0x7527bau32,
    0xc7ebe5u32,
    0xf17b3du32,
    0x0739f7u32,
    0x8a5292u32,
    0xea6bfbu32,
    0x5fb11fu32,
    0x8d5d08u32,
    0x560330u32,
    0x46fc7bu32,
    0x6babf0u32,
    0xcfbc20u32,
    0x9af436u32,
    0x1da9e3u32,
    0x91615eu32,
    0xe61b08u32,
    0x659985u32,
    0x5f14a0u32,
    0x68408du32,
    0xffd880u32,
    0x4d7327u32,
    0x310606u32,
    0x1556cau32,
    0x73a8c9u32,
    0x60e27bu32,
    0xc08c6bu32,
];

const NPIO2_HW: [u32; 32] = [
    0x3ff921fbu32,
    0x400921fbu32,
    0x4012d97cu32,
    0x401921fbu32,
    0x401f6a7au32,
    0x4022d97cu32,
    0x4025fdbbu32,
    0x402921fbu32,
    0x402c463au32,
    0x402f6a7au32,
    0x4031475cu32,
    0x4032d97cu32,
    0x40346b9cu32,
    0x4035fdbbu32,
    0x40378fdbu32,
    0x403921fbu32,
    0x403ab41bu32,
    0x403c463au32,
    0x403dd85au32,
    0x403f6a7au32,
    0x40407e4cu32,
    0x4041475cu32,
    0x4042106cu32,
    0x4042d97cu32,
    0x4043a28cu32,
    0x40446b9cu32,
    0x404534acu32,
    0x4045fdbbu32,
    0x4046c6cbu32,
    0x40478fdbu32,
    0x404858ebu32,
    0x404921fbu32,
];

// invpio2: 53 bits of 2/pi, and split pi/2 pieces (fdlibm)
const HALF: f64 = 5.00000000000000000000e-01;
const TWO24: f64 = 1.67772160000000000000e+07; // 2^24
const INVPIO2: f64 = f64::from_bits(0x3fe4_5f30_6dc9_c883);

const PIO2_1: f64 = 1.57079632673412561417e+00; // 0x3FF921FB54400000
const PIO2_1T: f64 = 6.07710050650619224932e-11; // 0x3DD0B4611A626331
const PIO2_2: f64 = 6.07710050630396597660e-11; // 0x3DD0B4611A600000
const PIO2_2T: f64 = 2.02226624879595063154e-21; // 0x3BA3198A2E037073
const PIO2_3: f64 = 2.02226624871116645580e-21; // 0x3BA3198A2E000000
const PIO2_3T: f64 = 8.47842766036889956997e-32; // 0x397B839A252049C1

// __kernel_sin coefficients (fdlibm)
const KS_HALF: f64 = 5.00000000000000000000e-01;
const S1: f64 = -1.66666666666666324348e-01;
const S2: f64 = 8.33333333332248946124e-03;
const S3: f64 = -1.98412698298579493134e-04;
const S4: f64 = 2.75573137070700676789e-06;
const S5: f64 = -2.50507602534068634195e-08;
const S6: f64 = 1.58969099521155010221e-10;

// __kernel_cos coefficients (fdlibm)
const KC_ONE: f64 = 1.00000000000000000000e+00;
const C1: f64 = 4.16666666666666019037e-02;
const C2: f64 = -1.38888888888741095749e-03;
const C3: f64 = 2.48015872894767294178e-05;
const C4: f64 = -2.75573143513906633035e-07;
const C5: f64 = 2.08757232129817482790e-09;
const C6: f64 = -1.13596475577881948265e-11;

// ========= IBM Accurate Mathematical Library sin/cos =========
const IBM_SN3: f64 = -1.66666666666664880952546298448555E-01;
const IBM_SN5: f64 = 8.33333214285722277379541354343671E-03;
const IBM_CS2: f64 = 4.99999999999999999999950396842453E-01;
const IBM_CS4: f64 = -4.16666666666664434524222570944589E-02;
const IBM_CS6: f64 = 1.38888874007937613028114285595617E-03;

const IBM_S1: f64 = f64::from_bits(0xbfc5_5555_5555_5555);
const IBM_S2: f64 = f64::from_bits(0x3f81_1111_1111_0ece);
const IBM_S3: f64 = f64::from_bits(0xbf2a_01a0_19db_08b8);
const IBM_S4: f64 = f64::from_bits(0x3ec7_1de2_7b9a_7ed9);
const IBM_S5: f64 = f64::from_bits(0xbe5a_ddff_c2fc_df59);

const IBM_BIG: f64 = f64::from_bits(0x42c8_0000_0000_0000);
const IBM_HP0: f64 = f64::from_bits(0x3ff9_21fb_5444_2d18);
const IBM_HP1: f64 = f64::from_bits(0x3c91_a626_3314_5c07);
const IBM_MP1: f64 = f64::from_bits(0x3ff9_21fb_5800_0000);
const IBM_MP2: f64 = f64::from_bits(0xbe4d_de97_3c00_0000);
const IBM_PP3: f64 = f64::from_bits(0xbc8c_b3b3_9800_0000);
const IBM_PP4: f64 = f64::from_bits(0xbacd_747f_23e3_2ed7);
const IBM_HPINV: f64 = f64::from_bits(0x3fe4_5f30_6dc9_c883);
const IBM_TOINT: f64 = f64::from_bits(0x4338_0000_0000_0000);
// 2^27 + 1 (exact). Used by branred.c to split x into two numbers.
const IBM_SPLIT: f64 = f64::from_bits(0x41a0_0000_0200_0000);
const SIGN_BIT: u64 = 0x8000_0000_0000_0000u64;

const BR_T576: f64 = f64::from_bits(0x63f0_0000_0000_0000);
const BR_TM600: f64 = f64::from_bits(0x1a70_0000_0000_0000);
const BR_TM24: f64 = f64::from_bits(0x3e70_0000_0000_0000);
const BR_BIG: f64 = f64::from_bits(0x4338_0000_0000_0000);
const BR_BIG1: f64 = f64::from_bits(0x4358_0000_0000_0000);
const BR_MP2: f64 = f64::from_bits(0xbe4d_de97_4000_0000);

const TOVERP: [f64; 75] = [
    10680707.0, 7228996.0, 1387004.0, 2578385.0, 16069853.0, 12639074.0, 9804092.0, 4427841.0,
    16666979.0, 11263675.0, 12935607.0, 2387514.0, 4345298.0, 14681673.0, 3074569.0, 13734428.0,
    16653803.0, 1880361.0, 10960616.0, 8533493.0, 3062596.0, 8710556.0, 7349940.0, 6258241.0,
    3772886.0, 3769171.0, 3798172.0, 8675211.0, 12450088.0, 3874808.0, 9961438.0, 366607.0,
    15675153.0, 9132554.0, 7151469.0, 3571407.0, 2607881.0, 12013382.0, 4155038.0, 6285869.0,
    7677882.0, 13102053.0, 15825725.0, 473591.0, 9065106.0, 15363067.0, 6271263.0, 9264392.0,
    5636912.0, 4652155.0, 7056368.0, 13614112.0, 10155062.0, 1944035.0, 9527646.0, 15080200.0,
    6658437.0, 6231200.0, 6832269.0, 16767104.0, 5075751.0, 3212806.0, 1398474.0, 7579849.0,
    6349435.0, 12618859.0, 4703257.0, 12806093.0, 14477321.0, 2786137.0, 12875403.0, 9837734.0,
    14528324.0, 13719321.0, 343717.0,
];

#[inline(always)]
pub(crate) fn kernel_sin(x: f64, y: f64, iy: i32) -> f64 {
    let ix = hi_word(x) & 0x7fff_ffff;
    if ix < 0x3e40_0000 {
        if (x as i32) == 0 {
            return x;
        }
    }
    let z = x * x;
    let v = z * x;
    let r = S2 + z * (S3 + z * (S4 + z * (S5 + z * S6)));
    if iy == 0 {
        x + v * (S1 + z * r)
    } else {
        x - ((z * (KS_HALF * y - v * r) - y) - v * S1)
    }
}

#[inline(always)]
pub(crate) fn kernel_cos(x: f64, y: f64) -> f64 {
    let ix = hi_word(x) & 0x7fff_ffff;
    if ix < 0x3e40_0000 {
        if (x as i32) == 0 {
            return KC_ONE;
        }
    }
    let z = x * x;
    let r = z * (C1 + z * (C2 + z * (C3 + z * (C4 + z * (C5 + z * C6)))));
    if ix < 0x3fd3_3333 {
        KC_ONE - (0.5 * z - (z * r - x * y))
    } else {
        let qx: f64;
        if ix > 0x3fe9_0000 {
            qx = 0.28125;
        } else {
            // qx = x/4 with low 32 bits cleared
            let hi = ix - 0x0020_0000;
            qx = with_hi_lo(hi, 0);
        }
        let hz = 0.5 * z - qx;
        let a = KC_ONE - qx;
        a - (hz - (z * r - x * y))
    }
}

// ---- __kernel_rem_pio2 (Payne–Hanek core) ----
const INIT_JK: [i32; 4] = [2, 3, 4, 6];
const PIO2_CHUNKS: [f64; 8] = [
    1.57079625129699707031e+00,
    7.54978941586159635335e-08,
    5.39030252995776476554e-15,
    3.28200341580791294123e-22,
    1.27065575308067607349e-29,
    1.22933308981111328932e-36,
    2.73370053816464559624e-44,
    2.16741683877804819444e-51,
];
const KR_ZERO: f64 = 0.0;
const KR_ONE: f64 = 1.0;
const KR_TWO24: f64 = 1.67772160000000000000e+07;
const KR_TWON24: f64 = 5.96046447753906250000e-08;

#[inline(always)]
fn kernel_rem_pio2(x: &[f64; 3], y: &mut [f64; 2], e0: i32, nx: i32, prec: i32) -> i32 {
    let mut iq = [0i32; 20];
    let mut f = [0f64; 20];
    let mut fq = [0f64; 20];
    let mut q = [0f64; 20];

    let jk = INIT_JK[prec as usize];
    let jp = jk;

    let jx = nx - 1;
    let mut jv = (e0 - 3) / 24;
    if jv < 0 {
        jv = 0;
    }
    let mut q0 = e0 - 24 * (jv + 1);

    let mut j = jv - jx;
    let m = jx + jk;
    let jx_us = jx as usize;

    for i in 0..=(m as usize) {
        f[i] = if j < 0 {
            KR_ZERO
        } else {
            TWO_OVER_PI[j as usize] as f64
        };
        j += 1;
    }

    for i in 0..=(jk as usize) {
        let mut fw = 0.0;
        for jj in 0..=jx_us {
            fw += x[jj] * f[(jx + (i as i32) - (jj as i32)) as usize];
        }
        q[i] = fw;
    }

    let mut jz = jk;

    'recompute: loop {
        // distill q[] into iq[] reversingly
        let mut z = q[jz as usize];
        let mut i = 0;
        let mut jj = jz;
        while jj > 0 {
            let fw = ((KR_TWON24 * z) as i32) as f64;
            iq[i] = (z - KR_TWO24 * fw) as i32;
            z = q[(jj - 1) as usize] + fw;
            i += 1;
            jj -= 1;
        }

        // compute n
        z = scalbn_internal(z, q0);
        z -= 8.0 * floor_f64(z * 0.125);
        let mut n = z as i32;
        z -= n as f64;

        let mut ih = 0;
        if q0 > 0 {
            let i2 = iq[(jz - 1) as usize] >> (24 - q0);
            n += i2;
            iq[(jz - 1) as usize] -= i2 << (24 - q0);
            ih = iq[(jz - 1) as usize] >> (23 - q0);
        } else if q0 == 0 {
            ih = iq[(jz - 1) as usize] >> 23;
        } else if z >= 0.5 {
            ih = 2;
        }

        if ih > 0 {
            n += 1;
            let mut carry = 0;
            for i in 0..(jz as usize) {
                let jv = iq[i];
                if carry == 0 {
                    if jv != 0 {
                        carry = 1;
                        iq[i] = 0x1_000000 - jv;
                    }
                } else {
                    iq[i] = 0x0ffffff - jv;
                }
            }
            if q0 > 0 {
                match q0 {
                    1 => iq[(jz - 1) as usize] &= 0x7fffff,
                    2 => iq[(jz - 1) as usize] &= 0x3fffff,
                    _ => {}
                }
            }
            if ih == 2 {
                z = KR_ONE - z;
                if carry != 0 {
                    z -= scalbn_internal(KR_ONE, q0);
                }
            }
        }

        // check if recomputation needed
        if z == KR_ZERO {
            let mut jacc = 0;
            for i in ((jk as usize)..=(jz as usize - 1)).rev() {
                jacc |= iq[i];
            }
            if jacc == 0 {
                // need recomputation
                let mut k = 1;
                while iq[(jk - k) as usize] == 0 {
                    k += 1;
                }
                for ii in (jz + 1)..=(jz + k) {
                    let idx = (jv + ii) as usize;
                    f[(jx + ii) as usize] = TWO_OVER_PI[idx] as f64;
                    let mut fw = 0.0;
                    for jj in 0..=jx_us {
                        fw += x[jj] * f[(jx + ii - (jj as i32)) as usize];
                    }
                    q[ii as usize] = fw;
                }
                jz += k;
                continue 'recompute;
            }
        }

        // chop off zero terms
        if z == 0.0 {
            jz -= 1;
            q0 -= 24;
            while iq[jz as usize] == 0 {
                jz -= 1;
                q0 -= 24;
            }
        } else {
            z = scalbn_internal(z, -q0);
            if z >= KR_TWO24 {
                let fw = ((KR_TWON24 * z) as i32) as f64;
                iq[jz as usize] = (z - KR_TWO24 * fw) as i32;
                jz += 1;
                q0 += 24;
                iq[jz as usize] = fw as i32;
            } else {
                iq[jz as usize] = z as i32;
            }
        }

        // convert iq[] chunks to q[]
        let mut fw = scalbn_internal(KR_ONE, q0);
        for i in (0..=(jz as usize)).rev() {
            q[i] = fw * (iq[i] as f64);
            fw *= KR_TWON24;
        }

        // compute PIo2 * q
        for i in (0..=(jz as usize)).rev() {
            let mut fw2 = 0.0;
            let mut k = 0usize;
            while k <= (jp as usize) && k <= (jz as usize - i) {
                fw2 += PIO2_CHUNKS[k] * q[i + k];
                k += 1;
            }
            fq[jz as usize - i] = fw2;
        }

        // compress fq into y (prec 2 -> 2 terms)
        let mut fw3 = 0.0;
        for i in (0..=(jz as usize)).rev() {
            fw3 += fq[i];
        }
        y[0] = if ih == 0 { fw3 } else { -fw3 };
        let mut fw4 = fq[0] - fw3;
        for i in 1..=(jz as usize) {
            fw4 += fq[i];
        }
        y[1] = if ih == 0 { fw4 } else { -fw4 };

        return n & 7;
    }
}

// ---- __ieee754_rem_pio2 wrapper ----
#[inline(always)]
fn rem_pio2(x: f64) -> (i32, f64, f64) {
    let hx = hi_word(x) as i32;
    let ix = (hx & 0x7fff_ffff) as u32;

    // |x| <= pi/4
    if ix <= 0x3fe9_21fbu32 {
        return (0, x, 0.0);
    }

    // |x| < 3pi/4: special n=±1
    if ix < 0x4002_d97cu32 {
        if hx > 0 {
            let z = x - PIO2_1;
            if ix != 0x3ff9_21fbu32 {
                let y0 = z - PIO2_1T;
                let y1 = (z - y0) - PIO2_1T;
                return (1, y0, y1);
            } else {
                let z2 = z - PIO2_2;
                let y0 = z2 - PIO2_2T;
                let y1 = (z2 - y0) - PIO2_2T;
                return (1, y0, y1);
            }
        } else {
            let z = x + PIO2_1;
            if ix != 0x3ff9_21fbu32 {
                let y0 = z + PIO2_1T;
                let y1 = (z - y0) + PIO2_1T;
                return (-1, y0, y1);
            } else {
                let z2 = z + PIO2_2;
                let y0 = z2 + PIO2_2T;
                let y1 = (z2 - y0) + PIO2_2T;
                return (-1, y0, y1);
            }
        }
    }

    // medium: |x| <= 2^19*(pi/2)
    if ix <= 0x4139_21fbu32 {
        let t = x.abs();
        #[cfg(target_feature = "fma")]
        let n = (super::fma_internal(t, INVPIO2, HALF)) as i32;
        #[cfg(not(target_feature = "fma"))]
        let n = (t * INVPIO2 + HALF) as i32;
        let fn_ = n as f64;

        let mut r = t - fn_ * PIO2_1;
        let mut w = fn_ * PIO2_1T;

        let mut y0 = r - w;

        if n < 32 && ix != NPIO2_HW[(n - 1) as usize] {
            // quick path
        } else {
            let j = (ix >> 20) as i32;
            let mut i = j - (((hi_word(y0) >> 20) & 0x7ff) as i32);
            if i > 16 {
                // 2nd iteration
                let t2 = r;
                w = fn_ * PIO2_2;
                r = t2 - w;
                w = fn_ * PIO2_2T - ((t2 - r) - w);
                y0 = r - w;
                i = j - (((hi_word(y0) >> 20) & 0x7ff) as i32);
                if i > 49 {
                    // 3rd iteration
                    let t3 = r;
                    w = fn_ * PIO2_3;
                    r = t3 - w;
                    w = fn_ * PIO2_3T - ((t3 - r) - w);
                    y0 = r - w;
                }
            }
        }

        let y1 = (r - y0) - w;
        if hx < 0 {
            return (-n, -y0, -y1);
        } else {
            return (n, y0, y1);
        }
    }

    // large: Payne–Hanek core
    if ix >= 0x7ff0_0000u32 {
        return (0, f64::NAN, f64::NAN);
    }

    let e0 = ((ix >> 20) as i32) - 1046; // ilogb(z)-23

    // set exponent so z is scaled
    let hi = ix - ((e0 as u32) << 20);
    let mut z = with_hi_lo(hi, lo_word(x));

    let mut tx = [0.0f64; 3];
    for i in 0..2 {
        tx[i] = (z as i32) as f64;
        z = (z - tx[i]) * TWO24;
    }
    tx[2] = z;

    let mut nx = 3;
    while nx > 0 && tx[nx - 1] == 0.0 {
        nx -= 1;
    }
    if nx == 0 {
        return (0, 0.0, 0.0);
    }

    let mut yy = [0.0f64; 2];
    let n = kernel_rem_pio2(&tx, &mut yy, e0, nx as i32, 2);
    if hx < 0 {
        (-(n as i32), -yy[0], -yy[1])
    } else {
        (n as i32, yy[0], yy[1])
    }
}

#[inline(always)]
fn taylor_sin(xx: f64, x: f64, dx: f64) -> f64 {
    let poly = (((IBM_S5 * xx + IBM_S4) * xx + IBM_S3) * xx + IBM_S2) * xx + IBM_S1;
    let t = (poly * x - 0.5 * dx) * xx + dx;
    x + t
}

#[inline(always)]
fn sincos_table_lookup(u: f64) -> (f64, f64, f64, f64) {
    let idx = (lo_word(u) as usize) << 2;
    let sn = SINCOS_TAB[idx];
    let ssn = SINCOS_TAB[idx + 1];
    let cs = SINCOS_TAB[idx + 2];
    let ccs = SINCOS_TAB[idx + 3];
    (sn, ssn, cs, ccs)
}

#[inline(always)]
fn do_cos(mut x: f64, mut dx: f64) -> f64 {
    if x < 0.0 {
        dx = -dx;
    }
    let u = IBM_BIG + x.abs();
    x = x.abs() - (u - IBM_BIG) + dx;

    let xx = x * x;
    let s = x + x * xx * (IBM_SN3 + xx * IBM_SN5);
    let c = xx * (IBM_CS2 + xx * (IBM_CS4 + xx * IBM_CS6));
    let (sn, ssn, cs, ccs) = sincos_table_lookup(u);
    let cor = (ccs - s * ssn - cs * c) - sn * s;
    cs + cor
}

#[inline(always)]
fn do_sin(mut x: f64, mut dx: f64) -> f64 {
    let sign = x.to_bits() & SIGN_BIT;
    let absx = x.abs();
    if absx < 0.126 {
        if sign != 0 {
            dx = -dx;
        }
        let res = taylor_sin(absx * absx, absx, dx);
        return f64::from_bits(res.to_bits() ^ sign);
    }
    if sign != 0 {
        dx = -dx;
    }
    let u = IBM_BIG + absx;
    x = absx - (u - IBM_BIG);

    let xx = x * x;
    let s = x + (dx + x * xx * (IBM_SN3 + xx * IBM_SN5));
    let c = x * dx + xx * (IBM_CS2 + xx * (IBM_CS4 + xx * IBM_CS6));
    let (sn, ssn, cs, ccs) = sincos_table_lookup(u);
    let cor = (ssn + s * ccs - sn * c) + cs * s;
    let res = sn + cor;
    f64::from_bits(res.to_bits() ^ sign)
}

// ---- FMA-specialized do_sin/do_cos (for glibc-like ifunc parity) ----

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn fma_f64(a: f64, b: f64, c: f64) -> f64 {
    use core::arch::x86_64::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
    _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)))
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn do_cos_fma(mut x: f64, mut dx: f64) -> f64 {
    if x < 0.0 {
        dx = -dx;
    }
    let absx = x.abs();
    let u = IBM_BIG + absx;
    x = absx - (u - IBM_BIG) + dx;

    let xx = x * x;
    let t = unsafe { fma_f64(xx, IBM_SN5, IBM_SN3) };
    let s = unsafe { fma_f64(x * xx, t, x) };
    let t = unsafe { fma_f64(xx, IBM_CS6, IBM_CS4) };
    let t = unsafe { fma_f64(xx, t, IBM_CS2) };
    let c = xx * t;

    let (sn, ssn, cs, ccs) = sincos_table_lookup(u);
    let cor = unsafe { fma_f64(-sn, s, fma_f64(-cs, c, fma_f64(-s, ssn, ccs))) };
    cs + cor
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn do_sin_fma(mut x: f64, mut dx: f64) -> f64 {
    let sign = x.to_bits() & SIGN_BIT;
    let absx = x.abs();
    if absx < 0.126 {
        if sign != 0 {
            dx = -dx;
        }
        let xx = absx * absx;
        let poly = unsafe {
            fma_f64(
                fma_f64(fma_f64(fma_f64(IBM_S5, xx, IBM_S4), xx, IBM_S3), xx, IBM_S2),
                xx,
                IBM_S1,
            )
        };
        let t = unsafe { fma_f64(fma_f64(poly, absx, -0.5 * dx), xx, dx) };
        let res = absx + t;
        return f64::from_bits(res.to_bits() ^ sign);
    }
    if sign != 0 {
        dx = -dx;
    }
    let u = IBM_BIG + absx;
    x = absx - (u - IBM_BIG);

    let xx = x * x;
    let t = unsafe { fma_f64(xx, IBM_SN5, IBM_SN3) };
    let s = x + unsafe { fma_f64(x * xx, t, dx) };

    let t = unsafe { fma_f64(xx, IBM_CS6, IBM_CS4) };
    let t = unsafe { fma_f64(xx, t, IBM_CS2) };
    let v = xx * t;
    let c = unsafe { fma_f64(x, dx, v) };

    let (sn, ssn, cs, ccs) = sincos_table_lookup(u);
    let cor = unsafe { fma_f64(cs, s, fma_f64(-sn, c, fma_f64(s, ccs, ssn))) };
    let res = sn + cor;
    f64::from_bits(res.to_bits() ^ sign)
}

#[cfg(not(target_arch = "x86_64"))]
#[inline(always)]
unsafe fn do_cos_fma(x: f64, dx: f64) -> f64 {
    let _ = (x, dx);
    unreachable!()
}

#[cfg(not(target_arch = "x86_64"))]
#[inline(always)]
unsafe fn do_sin_fma(x: f64, dx: f64) -> f64 {
    let _ = (x, dx);
    unreachable!()
}

#[inline(always)]
fn do_sincos(a: f64, da: f64, n: i32) -> f64 {
    let mut retval = if (n & 1) != 0 {
        do_cos(a, da)
    } else {
        do_sin(a, da)
    };
    if (n & 2) != 0 {
        retval = -retval;
    }
    retval
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn do_sincos_fma(a: f64, da: f64, n: i32) -> f64 {
    let mut retval = if (n & 1) != 0 {
        unsafe { do_cos_fma(a, da) }
    } else {
        unsafe { do_sin_fma(a, da) }
    };
    if (n & 2) != 0 {
        retval = -retval;
    }
    retval
}

#[inline(always)]
fn reduce_sincos(x: f64) -> (i32, f64, f64) {
    let t = x * IBM_HPINV + IBM_TOINT;
    let xn = t - IBM_TOINT;
    let n = (lo_word(t) & 3) as i32;

    let y = (x - xn * IBM_MP1) - xn * IBM_MP2;

    let t1 = xn * IBM_PP3;
    let t2 = y - t1;
    let mut db = (y - t2) - t1;

    let t1 = xn * IBM_PP4;
    let b = t2 - t1;
    db += (t2 - b) - t1;

    (n, b, db)
}

#[inline(always)]
fn adjust_gor(gor: f64, k: i32) -> f64 {
    let mut hi = hi_word(gor);
    hi = hi.wrapping_sub(((k * 24) as u32) << 20);
    let lo = lo_word(gor);
    f64::from_bits(((hi as u64) << 32) | (lo as u64))
}

#[inline(always)]
pub(super) fn branred(x: f64) -> (i32, f64, f64) {
    let mut r = [0.0f64; 6];
    let mut sum = 0.0;

    let x = x * BR_TM600;
    let t = x * IBM_SPLIT;
    let x1 = t - (t - x);
    let x2 = x - x1;

    let mut k = (((hi_word(x1) >> 20) & 2047) as i32 - 450) / 24;
    if k < 0 {
        k = 0;
    }
    let mut gor = adjust_gor(BR_T576, k);
    for i in 0..6 {
        r[i] = x1 * TOVERP[k as usize + i] * gor;
        gor *= BR_TM24;
    }
    for i in 0..3 {
        let s = (r[i] + BR_BIG) - BR_BIG;
        sum += s;
        r[i] -= s;
    }
    let mut t = 0.0;
    for i in 0..6 {
        t += r[5 - i];
    }
    let mut bb = (((((r[0] - t) + r[1]) + r[2]) + r[3]) + r[4]) + r[5];
    let s = (t + BR_BIG) - BR_BIG;
    sum += s;
    t -= s;
    let mut b = t + bb;
    bb += t - b;
    let s = (sum + BR_BIG1) - BR_BIG1;
    sum -= s;

    let b1 = b;
    let bb1 = bb;
    let sum1 = sum;
    sum = 0.0;

    let mut k = (((hi_word(x2) >> 20) & 2047) as i32 - 450) / 24;
    if k < 0 {
        k = 0;
    }
    let mut gor = adjust_gor(BR_T576, k);
    for i in 0..6 {
        r[i] = x2 * TOVERP[k as usize + i] * gor;
        gor *= BR_TM24;
    }
    for i in 0..3 {
        let s = (r[i] + BR_BIG) - BR_BIG;
        sum += s;
        r[i] -= s;
    }
    t = 0.0;
    for i in 0..6 {
        t += r[5 - i];
    }
    bb = (((((r[0] - t) + r[1]) + r[2]) + r[3]) + r[4]) + r[5];
    let s = (t + BR_BIG) - BR_BIG;
    sum += s;
    t -= s;
    b = t + bb;
    bb += t - b;
    let s = (sum + BR_BIG1) - BR_BIG1;
    sum -= s;

    let b2 = b;
    let bb2 = bb;
    let sum2 = sum;

    let mut sum = sum1 + sum2;
    let mut b = b1 + b2;
    let bb = if b1.abs() > b2.abs() {
        (b1 - b) + b2
    } else {
        (b2 - b) + b1
    };

    if b > 0.5 {
        b -= 1.0;
        sum += 1.0;
    } else if b < -0.5 {
        b += 1.0;
        sum -= 1.0;
    }

    let s = b + (bb + bb1 + bb2);
    let t = ((b - s) + bb) + (bb1 + bb2);
    let b = s * IBM_SPLIT;
    let t1 = b - (b - s);
    let t2 = s - t1;
    let b = s * IBM_HP0;
    let bb = (((t1 * IBM_MP1 - b) + t1 * BR_MP2) + t2 * IBM_MP1)
        + (t2 * BR_MP2 + s * IBM_HP1 + t * IBM_HP0);
    let s = b + bb;
    let t = (b - s) + bb;
    let n = (sum as i32) & 3;
    (n, s, t)
}

#[inline(always)]
fn sin_generic(x: f64) -> f64 {
    let k = hi_word(x) & 0x7fff_ffff;
    if k < 0x3e50_0000 {
        return x;
    }
    if k < 0x3feb_6000 {
        return do_sin(x, 0.0);
    }
    if k < 0x4003_68fd {
        let t = IBM_HP0 - x.abs();
        let v = do_cos(t, IBM_HP1);
        return if x.is_sign_negative() { -v } else { v };
    }
    if k < 0x4199_21fb {
        let (n, a, da) = reduce_sincos(x);
        return do_sincos(a, da, n);
    }
    if k < 0x7ff0_0000 {
        let (n, a, da) = branred(x);
        return do_sincos(a, da, n);
    }
    x * 0.0
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn sin_fma(x: f64) -> f64 {
    let k = hi_word(x) & 0x7fff_ffff;
    if k < 0x3e50_0000 {
        return x;
    }
    if k < 0x3feb_6000 {
        return unsafe { do_sin_fma(x, 0.0) };
    }
    if k < 0x4003_68fd {
        let t = IBM_HP0 - x.abs();
        let v = unsafe { do_cos_fma(t, IBM_HP1) };
        return if x.is_sign_negative() { -v } else { v };
    }
    if k < 0x4199_21fb {
        let (n, a, da) = reduce_sincos(x);
        return unsafe { do_sincos_fma(a, da, n) };
    }
    if k < 0x7ff0_0000 {
        let (n, a, da) = branred(x);
        return unsafe { do_sincos_fma(a, da, n) };
    }
    x * 0.0
}

#[inline(always)]
pub(super) fn sin(x: f64) -> f64 {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        if super::cpu_has_fma() {
            return sin_fma(x);
        }
    }
    sin_generic(x)
}

#[inline(always)]
fn cos_generic(x: f64) -> f64 {
    let k = hi_word(x) & 0x7fff_ffff;
    if k < 0x3e40_0000 {
        return 1.0;
    }
    if k < 0x3feb_6000 {
        return do_cos(x, 0.0);
    }
    if k < 0x4003_68fd {
        let y = IBM_HP0 - x.abs();
        let a = y + IBM_HP1;
        let da = (y - a) + IBM_HP1;
        return do_sin(a, da);
    }
    if k < 0x4199_21fb {
        let (n, a, da) = reduce_sincos(x);
        return do_sincos(a, da, n + 1);
    }
    if k < 0x7ff0_0000 {
        let (n, a, da) = branred(x);
        return do_sincos(a, da, n + 1);
    }
    x * 0.0
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn cos_fma(x: f64) -> f64 {
    let k = hi_word(x) & 0x7fff_ffff;
    if k < 0x3e40_0000 {
        return 1.0;
    }
    if k < 0x3feb_6000 {
        return unsafe { do_cos_fma(x, 0.0) };
    }
    if k < 0x4003_68fd {
        let y = IBM_HP0 - x.abs();
        let a = y + IBM_HP1;
        let da = (y - a) + IBM_HP1;
        return unsafe { do_sin_fma(a, da) };
    }
    if k < 0x4199_21fb {
        let (n, a, da) = reduce_sincos(x);
        return unsafe { do_sincos_fma(a, da, n + 1) };
    }
    if k < 0x7ff0_0000 {
        let (n, a, da) = branred(x);
        return unsafe { do_sincos_fma(a, da, n + 1) };
    }
    x * 0.0
}

#[inline(always)]
pub(super) fn cos(x: f64) -> f64 {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        if super::cpu_has_fma() {
            return cos_fma(x);
        }
    }
    cos_generic(x)
}

#[inline(always)]
fn sincos_generic(x: f64) -> (f64, f64) {
    let k = hi_word(x) & 0x7fff_ffff;
    if k < 0x3e40_0000 {
        return (x, 1.0);
    }
    if k < 0x3e50_0000 {
        return (x, 1.0);
    }
    if k < 0x3feb_6000 {
        return (do_sin(x, 0.0), do_cos(x, 0.0));
    }
    if k < 0x4003_68fd {
        let t = IBM_HP0 - x.abs();
        let sin_v = do_cos(t, IBM_HP1);
        let sin_v = if x.is_sign_negative() { -sin_v } else { sin_v };
        let a = t + IBM_HP1;
        let da = (t - a) + IBM_HP1;
        let cos_v = do_sin(a, da);
        return (sin_v, cos_v);
    }
    if k < 0x4199_21fb {
        let (n, a, da) = reduce_sincos(x);
        return (do_sincos(a, da, n), do_sincos(a, da, n + 1));
    }
    if k < 0x7ff0_0000 {
        let (n, a, da) = branred(x);
        return (do_sincos(a, da, n), do_sincos(a, da, n + 1));
    }
    let nan = x * 0.0;
    (nan, nan)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn sincos_fma(x: f64) -> (f64, f64) {
    let k = hi_word(x) & 0x7fff_ffff;
    if k < 0x3e40_0000 {
        return (x, 1.0);
    }
    if k < 0x3e50_0000 {
        return (x, 1.0);
    }
    if k < 0x3feb_6000 {
        return (unsafe { do_sin_fma(x, 0.0) }, unsafe { do_cos_fma(x, 0.0) });
    }
    if k < 0x4003_68fd {
        let t = IBM_HP0 - x.abs();
        let sin_v = unsafe { do_cos_fma(t, IBM_HP1) };
        let sin_v = if x.is_sign_negative() { -sin_v } else { sin_v };
        let a = t + IBM_HP1;
        let da = (t - a) + IBM_HP1;
        let cos_v = unsafe { do_sin_fma(a, da) };
        return (sin_v, cos_v);
    }
    if k < 0x4199_21fb {
        let (n, a, da) = reduce_sincos(x);
        return (unsafe { do_sincos_fma(a, da, n) }, unsafe {
            do_sincos_fma(a, da, n + 1)
        });
    }
    if k < 0x7ff0_0000 {
        let (n, a, da) = branred(x);
        return (unsafe { do_sincos_fma(a, da, n) }, unsafe {
            do_sincos_fma(a, da, n + 1)
        });
    }
    let nan = x * 0.0;
    (nan, nan)
}

#[inline(always)]
pub fn sincos(x: f64) -> (f64, f64) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        if super::cpu_has_fma() {
            return sincos_fma(x);
        }
    }
    sincos_generic(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sin_cos_identity() {
        for i in 0..1000 {
            let x = (i as f64) * 0.1;
            let s = sin(x);
            let c = cos(x);
            let identity = s * s + c * c;
            assert!(
                (identity - 1.0).abs() < 1e-15,
                "Identity failed for x={x}: got {identity}"
            );
        }
    }
}

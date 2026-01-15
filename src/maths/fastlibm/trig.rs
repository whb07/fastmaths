#![allow(
    clippy::collapsible_if,
    clippy::needless_late_init,
    clippy::needless_range_loop
)]

use super::{
    f64_to_bits, floor_f64, hi_word, is_inf_bits, is_nan_bits, lo_word, scalbn, with_hi_lo,
};

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

#[inline(always)]
fn kernel_sin(x: f64, y: f64, iy: i32) -> f64 {
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
fn kernel_cos(x: f64, y: f64) -> f64 {
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
        z = scalbn(z, q0);
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
                    z -= scalbn(KR_ONE, q0);
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
            z = scalbn(z, -q0);
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
        let mut fw = scalbn(KR_ONE, q0);
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
        let n = (t.mul_add(INVPIO2, HALF)) as i32;
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

#[inline(always)]
pub(super) fn sin(x: f64) -> f64 {
    let ux = f64_to_bits(x);
    if is_nan_bits(ux) || is_inf_bits(ux) {
        return f64::NAN;
    }

    let (n, y0, y1) = rem_pio2(x);
    match n & 3 {
        0 => kernel_sin(y0, y1, 1),
        1 => kernel_cos(y0, y1),
        2 => -kernel_sin(y0, y1, 1),
        _ => -kernel_cos(y0, y1),
    }
}

#[inline(always)]
pub(super) fn cos(x: f64) -> f64 {
    let ux = f64_to_bits(x);
    if is_nan_bits(ux) || is_inf_bits(ux) {
        return f64::NAN;
    }

    let (n, y0, y1) = rem_pio2(x);
    match n & 3 {
        0 => kernel_cos(y0, y1),
        1 => -kernel_sin(y0, y1, 1),
        2 => -kernel_cos(y0, y1),
        _ => kernel_sin(y0, y1, 1),
    }
}

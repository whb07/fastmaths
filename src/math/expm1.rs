//! expm1(x) implementation.
//!
//! Uses a dedicated small-x polynomial to avoid catastrophic cancellation and
//! switches to exp(x)-1 for larger magnitudes. Polynomial degree and constants
//! match fdlibm/glibc style minimax fits.

use super::{LN2_HI, LN2_LO, hi_word, lo_word, with_hi_lo};

const ONE: f64 = 1.0;
const HUGE: f64 = 1.0e300;
const TINY: f64 = 1.0e-300;
const O_THRESHOLD: f64 = f64::from_bits(0x4086_2e42_fefa_39ef);
const INVLN2: f64 = f64::from_bits(0x3ff7_1547_652b_82fe);
const Q: [f64; 6] = [
    1.0,
    f64::from_bits(0xbfa1_1111_1111_10f4),
    f64::from_bits(0x3f5a_01a0_19fe_5585),
    f64::from_bits(0xbf14_ce19_9eaa_dbb7),
    f64::from_bits(0x3ed0_cfca_86e6_5239),
    f64::from_bits(0xbe8a_fdb7_6e09_c32d),
];

#[inline(always)]
fn set_high_word(x: f64, hi: u32) -> f64 {
    with_hi_lo(hi, lo_word(x))
}

#[inline(always)]
fn expm1_generic(mut x: f64) -> f64 {
    if x == 0.0 {
        return x;
    }
    let mut hx = hi_word(x);
    let xsb = hx & 0x8000_0000;
    hx &= 0x7fff_ffff;

    if hx >= 0x4043_687a {
        if hx >= 0x4086_2e42 {
            if hx >= 0x7ff0_0000 {
                let low = lo_word(x);
                if ((hx & 0xfffff) | low) != 0 {
                    return x + x;
                }
                return if xsb == 0 { x } else { -1.0 };
            }
            if x > O_THRESHOLD {
                return HUGE * HUGE;
            }
        }
        if xsb != 0 {
            return TINY - ONE;
        }
    }
    let mut k = 0i32;
    let mut c = 0.0;

    if hx > 0x3fd6_2e42 {
        let (hi, lo) = if hx < 0x3ff0_a2b2 {
            if xsb == 0 {
                k = 1;
                (x - LN2_HI, LN2_LO)
            } else {
                k = -1;
                (x + LN2_HI, -LN2_LO)
            }
        } else {
            k = (INVLN2 * x + if xsb == 0 { 0.5 } else { -0.5 }) as i32;
            let t = k as f64;
            (x - t * LN2_HI, t * LN2_LO)
        };
        x = hi - lo;
        c = (hi - x) - lo;
    } else if hx < 0x3c90_0000 {
        let t = HUGE + x;
        return x - (t - (HUGE + x));
    }

    let hfx = 0.5 * x;
    let hxs = x * hfx;
    // Keep evaluation order aligned with fdlibm/glibc to avoid
    // cross-libm rounding deltas from fused operations.
    let r1 = ONE + hxs * Q[1];
    let h2 = hxs * hxs;
    let r2 = Q[2] + hxs * Q[3];
    let h4 = h2 * h2;
    let r3 = Q[4] + hxs * Q[5];
    let r1 = (r1 + h2 * r2) + h4 * r3;
    let t = 3.0 - r1 * hfx;
    let mut e = hxs * ((r1 - t) / (6.0 - x * t));

    if k == 0 {
        return x - (x * e - hxs);
    }

    e = (x * (e - c) - c) - hxs;
    if k == -1 {
        return 0.5 * (x - e) - 0.5;
    }
    if k == 1 {
        if x < -0.25 {
            return -2.0 * (e - (x + 0.5));
        }
        return ONE + 2.0 * (x - e);
    }

    if k <= -2 || k > 56 {
        let mut y = ONE - (e - x);
        let high = hi_word(y).wrapping_add((k as u32) << 20);
        y = set_high_word(y, high);
        return y - ONE;
    }

    if k < 20 {
        let t = set_high_word(ONE, 0x3ff0_0000 - (0x0020_0000 >> k));
        let mut y = t - (e - x);
        let high = hi_word(y).wrapping_add((k as u32) << 20);
        y = set_high_word(y, high);
        return y;
    }

    let t = set_high_word(ONE, (0x3ffu32.wrapping_sub(k as u32)) << 20);
    let mut y = x - (e + t);
    y += ONE;
    let high = hi_word(y).wrapping_add((k as u32) << 20);
    set_high_word(y, high)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn fma_f64(a: f64, b: f64, c: f64) -> f64 {
    use core::arch::x86_64::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
    _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)))
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn expm1_fma(mut x: f64) -> f64 {
    if x == 0.0 {
        return x;
    }
    let mut hx = hi_word(x);
    let xsb = hx & 0x8000_0000;
    hx &= 0x7fff_ffff;

    if hx >= 0x4043_687a {
        if hx >= 0x4086_2e42 {
            if hx >= 0x7ff0_0000 {
                let low = lo_word(x);
                if ((hx & 0xfffff) | low) != 0 {
                    return x + x;
                }
                return if xsb == 0 { x } else { -1.0 };
            }
            if x > O_THRESHOLD {
                return HUGE * HUGE;
            }
        }
        if xsb != 0 {
            return TINY - ONE;
        }
    }
    let mut k = 0i32;
    let mut c = 0.0;

    if hx > 0x3fd6_2e42 {
        let (hi, lo) = if hx < 0x3ff0_a2b2 {
            if xsb == 0 {
                k = 1;
                (x - LN2_HI, LN2_LO)
            } else {
                k = -1;
                (x + LN2_HI, -LN2_LO)
            }
        } else {
            k = (INVLN2 * x + if xsb == 0 { 0.5 } else { -0.5 }) as i32;
            let t = k as f64;
            (x - t * LN2_HI, t * LN2_LO)
        };
        x = hi - lo;
        c = (hi - x) - lo;
    } else if hx < 0x3c90_0000 {
        let t = HUGE + x;
        return x - (t - (HUGE + x));
    }

    let hfx = 0.5 * x;
    let hxs = x * hfx;
    let r1 = unsafe { fma_f64(hxs, Q[1], ONE) };
    let h2 = hxs * hxs;
    let r2 = unsafe { fma_f64(hxs, Q[3], Q[2]) };
    let h4 = h2 * h2;
    let r3 = unsafe { fma_f64(hxs, Q[5], Q[4]) };
    let r1 = unsafe { fma_f64(h2, r2, r1) };
    let r1 = unsafe { fma_f64(h4, r3, r1) };
    let t = unsafe { fma_f64(r1, -hfx, 3.0) };
    let mut e = hxs * ((r1 - t) / (6.0 - x * t));

    if k == 0 {
        return x - (x * e - hxs);
    }

    e = (x * (e - c) - c) - hxs;
    if k == -1 {
        return 0.5 * (x - e) - 0.5;
    }
    if k == 1 {
        if x < -0.25 {
            return -2.0 * (e - (x + 0.5));
        }
        return ONE + 2.0 * (x - e);
    }

    if k <= -2 || k > 56 {
        let mut y = ONE - (e - x);
        let high = hi_word(y).wrapping_add((k as u32) << 20);
        y = set_high_word(y, high);
        return y - ONE;
    }

    if k < 20 {
        let t = set_high_word(ONE, 0x3ff0_0000 - (0x0020_0000 >> k));
        let mut y = t - (e - x);
        let high = hi_word(y).wrapping_add((k as u32) << 20);
        y = set_high_word(y, high);
        return y;
    }

    let t = set_high_word(ONE, (0x3ffu32.wrapping_sub(k as u32)) << 20);
    let mut y = x - (e + t);
    y += ONE;
    let high = hi_word(y).wrapping_add((k as u32) << 20);
    set_high_word(y, high)
}

#[inline(always)]
pub fn expm1(x: f64) -> f64 {
    #[cfg(target_arch = "x86_64")]
    if super::cpu_has_fma() {
        // SAFETY: guarded by CPUID.
        return unsafe { expm1_fma(x) };
    }
    expm1_generic(x)
}

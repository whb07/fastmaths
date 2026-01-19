//! erf/erfc implementation.
//!
//! Piecewise rational approximations for |x| in several ranges (small, medium,
//! large). Uses compensated exp tails to preserve accuracy in erfc tails. Coefficients
//! are fdlibm-derived minimax fits.
use super::{exp, expm1, f64_to_bits, hi_word, with_hi_lo};

const ERX: f64 = 8.450_629_115_104_675_292_97e-01;
const EFX8: f64 = 1.027_033_336_764_100_690_53e+00;
const PP0: f64 = 1.283_791_670_955_125_585_61e-01;
const PP1: f64 = -3.250_421_072_470_014_993_70e-01;
const PP2: f64 = -2.848_174_957_559_851_047_66e-02;
const PP3: f64 = -5.770_270_296_489_441_591_57e-03;
const PP4: f64 = -2.376_301_665_665_016_260_84e-05;
const QQ1: f64 = 3.979_172_239_591_553_528_19e-01;
const QQ2: f64 = 6.502_224_998_876_729_444_85e-02;
const QQ3: f64 = 5.081_306_281_875_765_627_76e-03;
const QQ4: f64 = 1.324_947_380_043_216_445_26e-04;
const QQ5: f64 = -3.960_228_278_775_368_123_20e-06;
const PA0: f64 = -2.362_118_560_752_659_440_77e-03;
const PA1: f64 = 4.148_561_186_837_483_316_66e-01;
const PA2: f64 = -3.722_078_760_357_013_238_47e-01;
const PA3: f64 = 3.183_466_199_011_617_536_74e-01;
const PA4: f64 = -1.108_946_942_823_966_774_76e-01;
const PA5: f64 = 3.547_830_432_561_823_593_71e-02;
const PA6: f64 = -2.166_375_594_868_790_843_00e-03;
const QA1: f64 = 1.064_208_804_008_442_282_86e-01;
const QA2: f64 = 5.403_979_177_021_710_489_37e-01;
const QA3: f64 = 7.182_865_441_419_626_628_68e-02;
const QA4: f64 = 1.261_712_198_087_616_421_12e-01;
const QA5: f64 = 1.363_708_391_202_905_073_62e-02;
const QA6: f64 = 1.198_449_984_679_910_741_70e-02;
const RA0: f64 = -9.864_944_034_847_148_227_05e-03;
const RA1: f64 = -6.938_585_727_071_817_643_72e-01;
const RA2: f64 = -1.055_862_622_532_329_098_14e+01;
const RA3: f64 = -6.237_533_245_032_600_603_96e+01;
const RA4: f64 = -1.623_966_694_625_734_703_55e+02;
const RA5: f64 = -1.846_050_929_067_110_359_94e+02;
const RA6: f64 = -8.128_743_550_630_659_342_46e+01;
const RA7: f64 = -9.814_329_344_169_145_485_92e+00;
const SA1: f64 = 1.965_127_166_743_925_712_92e+01;
const SA2: f64 = 1.376_577_541_435_190_426_00e+02;
const SA3: f64 = 4.345_658_774_752_292_288_21e+02;
const SA4: f64 = 6.453_872_717_332_678_803_36e+02;
const SA5: f64 = 4.290_081_400_275_678_333_86e+02;
const SA6: f64 = 1.086_350_055_417_794_351_34e+02;
const SA7: f64 = 6.570_249_770_319_281_701_35e+00;
const SA8: f64 = -6.042_441_521_485_809_874_38e-02;
const RB0: f64 = -9.864_942_924_700_099_285_97e-03;
const RB1: f64 = -7.992_832_376_805_230_065_74e-01;
const RB2: f64 = -1.775_795_491_775_475_198_89e+01;
const RB3: f64 = -1.606_363_848_558_219_160_62e+02;
const RB4: f64 = -6.375_664_433_683_896_277_22e+02;
const RB5: f64 = -1.025_095_131_611_077_249_54e+03;
const RB6: f64 = -4.835_191_916_086_513_970_19e+02;
const SB1: f64 = 3.033_806_074_348_245_829_24e+01;
const SB2: f64 = 3.257_925_129_965_739_188_26e+02;
const SB3: f64 = 1.536_729_586_084_436_959_94e+03;
const SB4: f64 = 3.199_858_219_508_595_539_08e+03;
const SB5: f64 = 2.553_050_406_433_164_425_83e+03;
const SB6: f64 = 4.745_285_412_069_553_672_15e+02;
const SB7: f64 = -2.244_095_244_658_581_833_62e+01;
const SPLIT: f64 = 134_217_729.0; // 2^27 + 1

#[inline(always)]
fn fabs(x: f64) -> f64 {
    f64::from_bits(f64_to_bits(x) & 0x7fff_ffff_ffff_ffffu64)
}

#[inline(always)]
fn fma(x: f64, y: f64, z: f64) -> f64 {
    super::fma_internal(x, y, z)
}

#[inline(always)]
fn recip_refine(d: f64) -> f64 {
    let r0 = 1.0 / d;
    let err0 = fma(d, r0, -1.0);
    let r1 = r0 - r0 * err0;
    let err1 = fma(d, r1, -1.0);
    r1 - r1 * err1
}

#[inline(always)]
fn div_refine(n: f64, d: f64) -> f64 {
    let r0 = n / d;
    let p = d * r0;
    let e = fma(d, r0, -p);
    let rem = (n - p) - e;
    r0 + rem / d
}

#[inline(always)]
fn div_refine2(n: f64, d: f64) -> f64 {
    let r0 = div_refine(n, d);
    let p = d * r0;
    let e = fma(d, r0, -p);
    let rem = (n - p) - e;
    r0 + rem / d
}

#[inline(always)]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    (s, err)
}

#[inline(always)]
fn split(a: f64) -> (f64, f64) {
    let t = SPLIT * a;
    let hi = t - (t - a);
    let lo = a - hi;
    (hi, lo)
}

#[inline(always)]
fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    let (ah, al) = split(a);
    let (bh, bl) = split(b);
    let err = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
    (p, err)
}

#[inline(always)]
fn dd_mul_scalar(ah: f64, al: f64, b: f64) -> (f64, f64) {
    let (p, pe) = two_prod(ah, b);
    let s = pe + al * b;
    two_sum(p, s)
}

#[inline(always)]
fn dd_add_scalar(ah: f64, al: f64, b: f64) -> (f64, f64) {
    let (s, e) = two_sum(ah, b);
    two_sum(s, e + al)
}

#[inline(always)]
fn div_dd(nh: f64, nl: f64, dh: f64, dl: f64) -> f64 {
    let r0 = nh / dh;
    let p = dh * r0;
    let e1 = fma(dh, r0, -p);
    let rem = ((nh - p) - e1) + (nl - r0 * dl);
    r0 + rem / dh
}

#[inline(always)]
fn muldd2(xh: f64, xl: f64, ch: f64, cl: f64, l: &mut f64) -> f64 {
    let ahhh = ch * xh;
    *l = (ch * xl + cl * xh) + fma(ch, xh, -ahhh);
    ahhh
}

#[inline(always)]
fn exp_dd_parts(h: f64, l: f64) -> (f64, f64) {
    let eh = exp(h);
    if !eh.is_finite() || l == 0.0 {
        return (eh, 0.0);
    }
    let em1 = expm1(l);
    let ch = 1.0 + em1;
    let cl = em1 - (ch - 1.0);
    let mut pl = 0.0;
    let ph = muldd2(eh, 0.0, ch, cl, &mut pl);
    (ph, pl)
}

#[inline(always)]
fn with_set_low_word(x: f64, low: u32) -> f64 {
    with_hi_lo(hi_word(x), low)
}

fn erfc1(x: f64) -> f64 {
    let s = fabs(x) - 1.0;
    let (mut ph, mut pl) = (PA6, 0.0);
    (ph, pl) = dd_mul_scalar(ph, pl, s);
    (ph, pl) = dd_add_scalar(ph, pl, PA5);
    (ph, pl) = dd_mul_scalar(ph, pl, s);
    (ph, pl) = dd_add_scalar(ph, pl, PA4);
    (ph, pl) = dd_mul_scalar(ph, pl, s);
    (ph, pl) = dd_add_scalar(ph, pl, PA3);
    (ph, pl) = dd_mul_scalar(ph, pl, s);
    (ph, pl) = dd_add_scalar(ph, pl, PA2);
    (ph, pl) = dd_mul_scalar(ph, pl, s);
    (ph, pl) = dd_add_scalar(ph, pl, PA1);
    (ph, pl) = dd_mul_scalar(ph, pl, s);
    (ph, pl) = dd_add_scalar(ph, pl, PA0);
    let (mut qh, mut ql) = (QA6, 0.0);
    (qh, ql) = dd_mul_scalar(qh, ql, s);
    (qh, ql) = dd_add_scalar(qh, ql, QA5);
    (qh, ql) = dd_mul_scalar(qh, ql, s);
    (qh, ql) = dd_add_scalar(qh, ql, QA4);
    (qh, ql) = dd_mul_scalar(qh, ql, s);
    (qh, ql) = dd_add_scalar(qh, ql, QA3);
    (qh, ql) = dd_mul_scalar(qh, ql, s);
    (qh, ql) = dd_add_scalar(qh, ql, QA2);
    (qh, ql) = dd_mul_scalar(qh, ql, s);
    (qh, ql) = dd_add_scalar(qh, ql, QA1);
    (qh, ql) = dd_mul_scalar(qh, ql, s);
    (qh, ql) = dd_add_scalar(qh, ql, 1.0);
    let y = div_dd(ph, pl, qh, ql);
    let t = 1.0 - ERX;
    let r = t - y;
    let err = (t - r) - y;
    r + err
}

fn erfc2(ix: u32, mut x: f64) -> f64 {
    if ix < 0x3ff3_d70a {
        return erfc1(x);
    }
    x = fabs(x);
    let xx = x * x;
    let s = recip_refine(xx);
    let (rh, rl, sh, sl) = if ix < 0x4006_db6d {
        let (mut rh, mut rl) = (RA7, 0.0);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RA6);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RA5);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RA4);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RA3);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RA2);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RA1);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RA0);

        let (mut sh, mut sl) = (SA8, 0.0);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SA7);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SA6);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SA5);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SA4);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SA3);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SA2);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SA1);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, 1.0);
        (rh, rl, sh, sl)
    } else {
        let (mut rh, mut rl) = (RB6, 0.0);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RB5);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RB4);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RB3);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RB2);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RB1);
        (rh, rl) = dd_mul_scalar(rh, rl, s);
        (rh, rl) = dd_add_scalar(rh, rl, RB0);

        let (mut sh, mut sl) = (SB7, 0.0);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SB6);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SB5);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SB4);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SB3);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SB2);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, SB1);
        (sh, sl) = dd_mul_scalar(sh, sl, s);
        (sh, sl) = dd_add_scalar(sh, sl, 1.0);
        (rh, rl, sh, sl)
    };
    let z = with_set_low_word(x, 0);
    let base = fma(-z, z, -0.5625);
    let r_div = div_dd(rh, rl, sh, sl);
    let tail = fma(z - x, z + x, r_div);
    let sum = base + tail;
    let sum_tail = (base - sum) + tail;
    let (eh, el) = exp_dd_parts(sum, sum_tail);
    div_dd(eh, el, x, 0.0)
}

#[inline(always)]
pub fn erf(x: f64) -> f64 {
    let mut ix = hi_word(x);
    let sign = (ix >> 31) as usize;
    ix &= 0x7fff_ffff;

    if ix >= 0x7ff0_0000 {
        return 1.0 - 2.0 * (sign as f64) + 1.0 / x;
    }

    if ix < 0x3feb_0000 {
        if ix < 0x3e30_0000 {
            return 0.125 * (8.0 * x + EFX8 * x);
        }
        let z = x * x;
        let mut r = PP4;
        r = fma(z, r, PP3);
        r = fma(z, r, PP2);
        r = fma(z, r, PP1);
        r = fma(z, r, PP0);
        let mut s = QQ5;
        s = fma(z, s, QQ4);
        s = fma(z, s, QQ3);
        s = fma(z, s, QQ2);
        s = fma(z, s, QQ1);
        s = fma(z, s, 1.0);
        let y = r / s;
        return x + x * y;
    }

    let y = if ix < 0x4018_0000 {
        1.0 - erfc2(ix, x)
    } else {
        let x1p_1022 = f64::from_bits(0x0010_0000_0000_0000u64);
        1.0 - x1p_1022
    };

    if sign != 0 { -y } else { y }
}

#[inline(always)]
pub fn erfc(x: f64) -> f64 {
    let mut ix = hi_word(x);
    let sign = (ix >> 31) as usize;
    ix &= 0x7fff_ffff;

    if ix >= 0x7ff0_0000 {
        return 2.0 * (sign as f64) + 1.0 / x;
    }

    if ix < 0x3feb_0000 {
        if ix < 0x3c70_0000 {
            return 1.0 - x;
        }
        let z = x * x;
        let mut r = PP4;
        r = fma(z, r, PP3);
        r = fma(z, r, PP2);
        r = fma(z, r, PP1);
        r = fma(z, r, PP0);
        let mut s = QQ5;
        s = fma(z, s, QQ4);
        s = fma(z, s, QQ3);
        s = fma(z, s, QQ2);
        s = fma(z, s, QQ1);
        s = fma(z, s, 1.0);
        let y = r / s;
        let t = fma(x, y, x);
        let t_err = fma(x, y, x - t);
        if sign != 0 || ix < 0x3fd0_0000 {
            return (1.0 - t) - t_err;
        }
        return (0.5 - (t - 0.5)) - t_err;
    }

    if ix < 0x403c_0000 {
        if sign != 0 {
            return 2.0 - erfc2(ix, x);
        }
        return erfc2(ix, x);
    }

    let x1p_1022 = f64::from_bits(0x0010_0000_0000_0000u64);
    if sign != 0 {
        2.0 - x1p_1022
    } else {
        x1p_1022 * x1p_1022
    }
}

//! pow(x,y) implementation.
//!
//! Handles special cases (negative bases, integer exponents, NaN/Inf) and reduces
//! to exp(y*log(x)) for the general case. Uses ln/exp cores with split constants
//! for accuracy.

use super::{LN2_HI, LN2_LO, f64_from_bits, f64_to_bits, fma_internal, ln};

const POW_LOG_TABLE_BITS: u32 = 7;
const POW_LOG_N: u64 = 1u64 << POW_LOG_TABLE_BITS;
const POW_LOG_OFF: u64 = 0x3fe6_9555_0000_0000u64;
const SIGN_BIT: u64 = 0x8000_0000_0000_0000u64;

const POW_OVERFLOW_LOG: f64 = 709.782_712_893_384;
const POWI_EXP_CUTOFF: i64 = 64;
const POW_LOG_A: [f64; 7] = [
    f64::from_bits(0xbfe0_0000_0000_0000),
    f64::from_bits(0xbfe5_5555_5555_5560),
    f64::from_bits(0x3fe0_0000_0000_0006),
    f64::from_bits(0x3fe9_9999_9959_554e),
    f64::from_bits(0xbfe5_5555_5529_a47a),
    f64::from_bits(0xbff2_495b_9b48_45e9),
    f64::from_bits(0x3ff0_002b_8b26_3fc3),
];

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn fma_f64(a: f64, b: f64, c: f64) -> f64 {
    use core::arch::x86_64::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
    _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)))
}

const POW_LOG_INVC: [u64; 128] = [
    0x3ff6a00000000000u64,
    0x3ff6800000000000u64,
    0x3ff6600000000000u64,
    0x3ff6400000000000u64,
    0x3ff6200000000000u64,
    0x3ff6000000000000u64,
    0x3ff5e00000000000u64,
    0x3ff5c00000000000u64,
    0x3ff5a00000000000u64,
    0x3ff5800000000000u64,
    0x3ff5600000000000u64,
    0x3ff5600000000000u64,
    0x3ff5400000000000u64,
    0x3ff5200000000000u64,
    0x3ff5000000000000u64,
    0x3ff4e00000000000u64,
    0x3ff4c00000000000u64,
    0x3ff4a00000000000u64,
    0x3ff4a00000000000u64,
    0x3ff4800000000000u64,
    0x3ff4600000000000u64,
    0x3ff4400000000000u64,
    0x3ff4200000000000u64,
    0x3ff4000000000000u64,
    0x3ff4000000000000u64,
    0x3ff3e00000000000u64,
    0x3ff3c00000000000u64,
    0x3ff3a00000000000u64,
    0x3ff3a00000000000u64,
    0x3ff3800000000000u64,
    0x3ff3600000000000u64,
    0x3ff3400000000000u64,
    0x3ff3400000000000u64,
    0x3ff3200000000000u64,
    0x3ff3000000000000u64,
    0x3ff3000000000000u64,
    0x3ff2e00000000000u64,
    0x3ff2c00000000000u64,
    0x3ff2c00000000000u64,
    0x3ff2a00000000000u64,
    0x3ff2800000000000u64,
    0x3ff2600000000000u64,
    0x3ff2600000000000u64,
    0x3ff2400000000000u64,
    0x3ff2400000000000u64,
    0x3ff2200000000000u64,
    0x3ff2000000000000u64,
    0x3ff2000000000000u64,
    0x3ff1e00000000000u64,
    0x3ff1c00000000000u64,
    0x3ff1c00000000000u64,
    0x3ff1a00000000000u64,
    0x3ff1a00000000000u64,
    0x3ff1800000000000u64,
    0x3ff1600000000000u64,
    0x3ff1600000000000u64,
    0x3ff1400000000000u64,
    0x3ff1400000000000u64,
    0x3ff1200000000000u64,
    0x3ff1000000000000u64,
    0x3ff1000000000000u64,
    0x3ff0e00000000000u64,
    0x3ff0e00000000000u64,
    0x3ff0c00000000000u64,
    0x3ff0c00000000000u64,
    0x3ff0a00000000000u64,
    0x3ff0a00000000000u64,
    0x3ff0800000000000u64,
    0x3ff0800000000000u64,
    0x3ff0600000000000u64,
    0x3ff0400000000000u64,
    0x3ff0400000000000u64,
    0x3ff0200000000000u64,
    0x3ff0200000000000u64,
    0x3ff0000000000000u64,
    0x3ff0000000000000u64,
    0x3fefc00000000000u64,
    0x3fef800000000000u64,
    0x3fef400000000000u64,
    0x3fef000000000000u64,
    0x3feec00000000000u64,
    0x3fee800000000000u64,
    0x3fee400000000000u64,
    0x3fee200000000000u64,
    0x3fede00000000000u64,
    0x3feda00000000000u64,
    0x3fed600000000000u64,
    0x3fed400000000000u64,
    0x3fed000000000000u64,
    0x3fecc00000000000u64,
    0x3feca00000000000u64,
    0x3fec600000000000u64,
    0x3fec400000000000u64,
    0x3fec000000000000u64,
    0x3febe00000000000u64,
    0x3feba00000000000u64,
    0x3feb800000000000u64,
    0x3feb400000000000u64,
    0x3feb200000000000u64,
    0x3feae00000000000u64,
    0x3feac00000000000u64,
    0x3feaa00000000000u64,
    0x3fea600000000000u64,
    0x3fea400000000000u64,
    0x3fea000000000000u64,
    0x3fe9e00000000000u64,
    0x3fe9c00000000000u64,
    0x3fe9a00000000000u64,
    0x3fe9600000000000u64,
    0x3fe9400000000000u64,
    0x3fe9200000000000u64,
    0x3fe9000000000000u64,
    0x3fe8c00000000000u64,
    0x3fe8a00000000000u64,
    0x3fe8800000000000u64,
    0x3fe8600000000000u64,
    0x3fe8400000000000u64,
    0x3fe8200000000000u64,
    0x3fe7e00000000000u64,
    0x3fe7c00000000000u64,
    0x3fe7a00000000000u64,
    0x3fe7800000000000u64,
    0x3fe7600000000000u64,
    0x3fe7400000000000u64,
    0x3fe7200000000000u64,
    0x3fe7000000000000u64,
    0x3fe6e00000000000u64,
    0x3fe6c00000000000u64,
];

const POW_LOG_LOGC: [u64; 128] = [
    0xbfd62c82f2b9c800u64,
    0xbfd5d1bdbf580800u64,
    0xbfd5767717455800u64,
    0xbfd51aad872df800u64,
    0xbfd4be5f95777800u64,
    0xbfd4618bc21c6000u64,
    0xbfd404308686a800u64,
    0xbfd3a64c55694800u64,
    0xbfd347dd9a988000u64,
    0xbfd2e8e2bae12000u64,
    0xbfd2895a13de8800u64,
    0xbfd2895a13de8800u64,
    0xbfd22941fbcf7800u64,
    0xbfd1c898c1699800u64,
    0xbfd1675cababa800u64,
    0xbfd1058bf9ae4800u64,
    0xbfd0a324e2739000u64,
    0xbfd0402594b4d000u64,
    0xbfd0402594b4d000u64,
    0xbfcfb9186d5e4000u64,
    0xbfcef0adcbdc6000u64,
    0xbfce27076e2af000u64,
    0xbfcd5c216b4fc000u64,
    0xbfcc8ff7c79aa000u64,
    0xbfcc8ff7c79aa000u64,
    0xbfcbc286742d9000u64,
    0xbfcaf3c94e80c000u64,
    0xbfca23bc1fe2b000u64,
    0xbfca23bc1fe2b000u64,
    0xbfc9525a9cf45000u64,
    0xbfc87fa06520d000u64,
    0xbfc7ab890210e000u64,
    0xbfc7ab890210e000u64,
    0xbfc6d60fe719d000u64,
    0xbfc5ff3070a79000u64,
    0xbfc5ff3070a79000u64,
    0xbfc526e5e3a1b000u64,
    0xbfc44d2b6ccb8000u64,
    0xbfc44d2b6ccb8000u64,
    0xbfc371fc201e9000u64,
    0xbfc29552f81ff000u64,
    0xbfc1b72ad52f6000u64,
    0xbfc1b72ad52f6000u64,
    0xbfc0d77e7cd09000u64,
    0xbfc0d77e7cd09000u64,
    0xbfbfec9131dbe000u64,
    0xbfbe27076e2b0000u64,
    0xbfbe27076e2b0000u64,
    0xbfbc5e548f5bc000u64,
    0xbfba926d3a4ae000u64,
    0xbfba926d3a4ae000u64,
    0xbfb8c345d631a000u64,
    0xbfb8c345d631a000u64,
    0xbfb6f0d28ae56000u64,
    0xbfb51b073f062000u64,
    0xbfb51b073f062000u64,
    0xbfb341d7961be000u64,
    0xbfb341d7961be000u64,
    0xbfb16536eea38000u64,
    0xbfaf0a30c0118000u64,
    0xbfaf0a30c0118000u64,
    0xbfab42dd71198000u64,
    0xbfab42dd71198000u64,
    0xbfa77458f632c000u64,
    0xbfa77458f632c000u64,
    0xbfa39e87b9fec000u64,
    0xbfa39e87b9fec000u64,
    0xbf9f829b0e780000u64,
    0xbf9f829b0e780000u64,
    0xbf97b91b07d58000u64,
    0xbf8fc0a8b0fc0000u64,
    0xbf8fc0a8b0fc0000u64,
    0xbf7fe02a6b100000u64,
    0xbf7fe02a6b100000u64,
    0x0000000000000000u64,
    0x0000000000000000u64,
    0x3f80101575890000u64,
    0x3f90205658938000u64,
    0x3f98492528c90000u64,
    0x3fa0415d89e74000u64,
    0x3fa466aed42e0000u64,
    0x3fa894aa149fc000u64,
    0x3faccb73cdddc000u64,
    0x3faeea31c006c000u64,
    0x3fb1973bd1466000u64,
    0x3fb3bdf5a7d1e000u64,
    0x3fb5e95a4d97a000u64,
    0x3fb700d30aeac000u64,
    0x3fb9335e5d594000u64,
    0x3fbb6ac88dad6000u64,
    0x3fbc885801bc4000u64,
    0x3fbec739830a2000u64,
    0x3fbfe89139dbe000u64,
    0x3fc1178e8227e000u64,
    0x3fc1aa2b7e23f000u64,
    0x3fc2d1610c868000u64,
    0x3fc365fcb0159000u64,
    0x3fc4913d8333b000u64,
    0x3fc527e5e4a1b000u64,
    0x3fc6574ebe8c1000u64,
    0x3fc6f0128b757000u64,
    0x3fc7898d85445000u64,
    0x3fc8beafeb390000u64,
    0x3fc95a5adcf70000u64,
    0x3fca93ed3c8ae000u64,
    0x3fcb31d8575bd000u64,
    0x3fcbd087383be000u64,
    0x3fcc6ffbc6f01000u64,
    0x3fcdb13db0d49000u64,
    0x3fce530effe71000u64,
    0x3fcef5ade4dd0000u64,
    0x3fcf991c6cb3b000u64,
    0x3fd07138604d5800u64,
    0x3fd0c42d67616000u64,
    0x3fd1178e8227e800u64,
    0x3fd16b5ccbacf800u64,
    0x3fd1bf99635a6800u64,
    0x3fd214456d0eb800u64,
    0x3fd2bef07cdc9000u64,
    0x3fd314f1e1d36000u64,
    0x3fd36b6776be1000u64,
    0x3fd3c25277333000u64,
    0x3fd419b423d5e800u64,
    0x3fd4718dc271c800u64,
    0x3fd4c9e09e173000u64,
    0x3fd522ae0738a000u64,
    0x3fd57bf753c8d000u64,
    0x3fd5d5bddf596000u64,
];

const POW_LOG_LOGCTAIL: [u64; 128] = [
    0x3cfab42428375680u64,
    0xbd1ca508d8e0f720u64,
    0xbd2362a4d5b6506du64,
    0xbce684e49eb067d5u64,
    0xbd041b6993293ee0u64,
    0x3d13d82f484c84ccu64,
    0x3cdc42f3ed820b3au64,
    0x3d20b1c686519460u64,
    0x3d25594dd4c58092u64,
    0x3d267b1e99b72bd8u64,
    0x3d15ca14b6cfb03fu64,
    0x3d15ca14b6cfb03fu64,
    0xbd165a242853da76u64,
    0xbd1fafbc68e75404u64,
    0x3d1f1fc63382a8f0u64,
    0xbd26a8c4fd055a66u64,
    0xbd0c6bee7ef4030eu64,
    0xbcf036b89ef42d7fu64,
    0xbcf036b89ef42d7fu64,
    0x3d0d572aab993c87u64,
    0x3d2b26b79c86af24u64,
    0xbd172f4f543fff10u64,
    0x3d21ba91bbca681bu64,
    0x3d27794f689f8434u64,
    0x3d27794f689f8434u64,
    0x3d194eb0318bb78fu64,
    0x3cba4e633fcd9066u64,
    0xbd258c64dc46c1eau64,
    0xbd258c64dc46c1eau64,
    0xbd2ad1d904c1d4e3u64,
    0x3d2bbdbf7fdbfa09u64,
    0x3d2bdb9072534a58u64,
    0x3d2bdb9072534a58u64,
    0xbd10e46aa3b2e266u64,
    0xbd1e9e439f105039u64,
    0xbd1e9e439f105039u64,
    0xbd20de8b90075b8fu64,
    0x3d170cc16135783cu64,
    0x3d170cc16135783cu64,
    0x3cf178864d27543au64,
    0xbd248d301771c408u64,
    0xbd2e80a41811a396u64,
    0xbd2e80a41811a396u64,
    0x3d0a699688e85bf4u64,
    0x3d0a699688e85bf4u64,
    0xbd2575545ca333f2u64,
    0x3d2a342c2af0003cu64,
    0x3d2a342c2af0003cu64,
    0xbd1d0c57585fbe06u64,
    0x3d253935e85baac8u64,
    0x3d253935e85baac8u64,
    0x3d137c294d2f5668u64,
    0x3d137c294d2f5668u64,
    0xbd269737c93373dau64,
    0x3d1f025b61c65e57u64,
    0x3d1f025b61c65e57u64,
    0x3d2c5edaccf913dfu64,
    0x3d2c5edaccf913dfu64,
    0x3d147c5e768fa309u64,
    0x3d2d599e83368e91u64,
    0x3d2d599e83368e91u64,
    0x3d1c827ae5d6704cu64,
    0x3d1c827ae5d6704cu64,
    0xbd2cfc4634f2a1eeu64,
    0xbd2cfc4634f2a1eeu64,
    0x3cf502b7f526feaau64,
    0x3cf502b7f526feaau64,
    0xbd2980267c7e09e4u64,
    0xbd2980267c7e09e4u64,
    0xbd288d5493faa639u64,
    0xbcdf1e7cf6d3a69cu64,
    0xbcdf1e7cf6d3a69cu64,
    0xbd19e23f0dda40e4u64,
    0xbd19e23f0dda40e4u64,
    0x0000000000000000u64,
    0x0000000000000000u64,
    0xbd10c76b999d2be8u64,
    0xbd23dc5b06e2f7d2u64,
    0xbd2aa0ba325a0c34u64,
    0x3d0111c05cf1d753u64,
    0xbd2c167375bdfd28u64,
    0xbd197995d05a267du64,
    0xbd1a68f247d82807u64,
    0xbd0e113e4fc93b7bu64,
    0xbd25325d560d9e9bu64,
    0x3d2cc85ea5db4ed7u64,
    0xbd2c69063c5d1d1eu64,
    0x3cec1e8da99ded32u64,
    0x3d23115c3abd47dau64,
    0xbd1390802bf768e5u64,
    0x3d2646d1c65aacd3u64,
    0xbd2dc068afe645e0u64,
    0xbd2534d64fa10afdu64,
    0x3d21ef78ce2d07f2u64,
    0x3d2ca78e44389934u64,
    0x3d039d6ccb81b4a1u64,
    0x3cc62fa8234b7289u64,
    0x3d25837954fdb678u64,
    0x3d2633e8e5697dc7u64,
    0x3d19cf8b2c3c2e78u64,
    0xbd25118de59c21e1u64,
    0xbd1c661070914305u64,
    0xbd073d54aae92cd1u64,
    0x3d07f22858a0ff6fu64,
    0xbd28724350562169u64,
    0xbd0c358d4eace1aau64,
    0xbd2d4bc4595412b6u64,
    0xbcf1ec72c5962bd2u64,
    0xbd2aff2af715b035u64,
    0x3cc212276041f430u64,
    0xbcca211565bb8e11u64,
    0x3d1bcbecca0cdf30u64,
    0x3cf89cdb16ed4e91u64,
    0x3d27188b163ceae9u64,
    0xbd2c210e63a5f01cu64,
    0x3d2b9acdf7a51681u64,
    0x3d2ca6ed5147bdb7u64,
    0x3d0a87deba46baeau64,
    0x3d2a9cfa4a5004f4u64,
    0xbd28e27ad3213cb8u64,
    0x3d116ecdb0f177c8u64,
    0x3d183b54b606bd5cu64,
    0x3d08e436ec90e09du64,
    0xbd2f27ce0967d675u64,
    0xbd2e20891b0ad8a4u64,
    0x3d2ebe708164c759u64,
    0x3d1fadedee5d40efu64,
    0xbd0a0b2a08a465dcu64,
];

#[inline]
fn classify_integer(x: f64) -> (bool, bool) {
    if !x.is_finite() {
        return (false, false);
    }
    let ux = f64_to_bits(x) & 0x7fff_ffff_ffff_ffffu64;
    let exp = ((ux >> 52) & 0x7ff) as i32;
    if exp < 1023 {
        return (x == 0.0, false);
    }
    if exp > 1075 {
        return (true, false);
    }
    let frac_bits = 52 - (exp - 1023);
    if frac_bits == 0 {
        let odd = (ux & 1) != 0;
        return (true, odd);
    }
    let mask = (1u64 << frac_bits) - 1;
    if (ux & mask) != 0 {
        return (false, false);
    }
    let odd = ((ux >> frac_bits) & 1) != 0;
    (true, odd)
}

#[inline]
fn dd_mul(a_hi: f64, a_lo: f64, b_hi: f64, b_lo: f64) -> (f64, f64) {
    let p = a_hi * b_hi;
    if !p.is_finite() || !a_hi.is_finite() || !b_hi.is_finite() {
        return (p, 0.0);
    }
    let e = fma_internal(a_hi, b_hi, -p) + (a_hi * b_lo + a_lo * b_hi) + (a_lo * b_lo);
    let hi = p + e;
    let lo = (p - hi) + e;
    (hi, lo)
}

#[inline]
fn dd_recip(a_hi: f64, a_lo: f64) -> f64 {
    if a_hi == 0.0 {
        if a_lo == 0.0 {
            return f64::INFINITY;
        }
        return 1.0 / a_lo;
    }
    // Newton refinement of reciprocal using double-double product.
    let mut r = 1.0 / a_hi;
    // For very small a_hi, the naive reciprocal can overflow to inf and the
    // Newton step would turn into inf + (-inf) = NaN. In that regime the true
    // reciprocal overflows too (a_lo is only a tiny correction), so return +inf.
    if !r.is_finite() {
        return f64::INFINITY;
    }
    for _ in 0..2 {
        let (p_hi, p_lo) = dd_mul(a_hi, a_lo, r, 0.0);
        let mut err = fma_internal(-p_hi, 1.0, 1.0);
        err -= p_lo;
        r = r + r * err;
    }
    r
}

fn powi(base: f64, mut exp: i64) -> f64 {
    if exp == 0 {
        return 1.0;
    }
    let neg = exp < 0;
    if neg {
        exp = -exp;
    }
    let mut acc_hi = 1.0;
    let mut acc_lo = 0.0;
    let mut b_hi = base;
    let mut b_lo = 0.0;
    let mut e = exp as u64;
    while e != 0 {
        if (e & 1) != 0 {
            let (hi, lo) = dd_mul(acc_hi, acc_lo, b_hi, b_lo);
            acc_hi = hi;
            acc_lo = lo;
        }
        let (hi, lo) = dd_mul(b_hi, b_lo, b_hi, b_lo);
        b_hi = hi;
        b_lo = lo;
        e >>= 1;
    }
    if !neg {
        return acc_hi + acc_lo;
    }

    dd_recip(acc_hi, acc_lo)
}

#[inline]
fn log_inline_generic(ix: u64, tail: &mut f64) -> f64 {
    let tmp = ix.wrapping_sub(POW_LOG_OFF);
    let i = ((tmp >> (52 - POW_LOG_TABLE_BITS)) & (POW_LOG_N - 1)) as usize;
    let k = ((tmp as i64) >> 52) as i32;
    let iz = ix.wrapping_sub(tmp & (0x0fffu64 << 52));
    let z = f64_from_bits(iz);
    let kd = k as f64;

    let invc = f64_from_bits(POW_LOG_INVC[i]);
    let logc = f64_from_bits(POW_LOG_LOGC[i]);
    let logctail = f64_from_bits(POW_LOG_LOGCTAIL[i]);

    let zhi = f64_from_bits((iz.wrapping_add(1u64 << 31)) & (!0u64 << 32));
    let zlo = z - zhi;
    let rhi = zhi * invc - 1.0;
    let rlo = zlo * invc;
    let r = rhi + rlo;

    let t1 = kd * LN2_HI + logc;
    let t2 = t1 + r;
    let lo1 = kd * LN2_LO + logctail;
    let lo2 = t1 - t2 + r;

    let ar = POW_LOG_A[0] * r;
    let ar2 = r * ar;
    let ar3 = r * ar2;

    let arhi = POW_LOG_A[0] * rhi;
    let arhi2 = rhi * arhi;
    let hi = t2 + arhi2;
    let lo3 = rlo * (ar + arhi);
    let lo4 = t2 - hi + arhi2;

    let p = ar3
        * (POW_LOG_A[1]
            + r * POW_LOG_A[2]
            + ar2 * (POW_LOG_A[3] + r * POW_LOG_A[4] + ar2 * (POW_LOG_A[5] + r * POW_LOG_A[6])));
    let lo = lo1 + lo2 + lo3 + lo4 + p;
    let y = hi + lo;
    *tail = hi - y + lo;
    y
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn log_inline_fma(ix: u64, tail: &mut f64) -> f64 {
    let tmp = ix.wrapping_sub(POW_LOG_OFF);
    let i = ((tmp >> (52 - POW_LOG_TABLE_BITS)) & (POW_LOG_N - 1)) as usize;
    let k = ((tmp as i64) >> 52) as i32;
    let iz = ix.wrapping_sub(tmp & (0x0fffu64 << 52));
    let z = f64_from_bits(iz);
    let kd = k as f64;

    let invc = f64_from_bits(POW_LOG_INVC[i]);
    let logc = f64_from_bits(POW_LOG_LOGC[i]);
    let logctail = f64_from_bits(POW_LOG_LOGCTAIL[i]);

    let r = unsafe { fma_f64(z, invc, -1.0) };

    let t1 = kd * LN2_HI + logc;
    let t2 = t1 + r;
    let lo1 = kd * LN2_LO + logctail;
    let lo2 = t1 - t2 + r;

    let ar = POW_LOG_A[0] * r;
    let ar2 = r * ar;
    let ar3 = r * ar2;

    let hi = t2 + ar2;
    let lo3 = unsafe { fma_f64(ar, r, -ar2) };
    let lo4 = t2 - hi + ar2;

    let p = ar3
        * (POW_LOG_A[1]
            + r * POW_LOG_A[2]
            + ar2 * (POW_LOG_A[3] + r * POW_LOG_A[4] + ar2 * (POW_LOG_A[5] + r * POW_LOG_A[6])));
    let lo = lo1 + lo2 + lo3 + lo4 + p;
    let y = hi + lo;
    *tail = hi - y + lo;
    y
}

#[inline]
fn mul_log_generic(y: f64, hi: f64, lo: f64) -> (f64, f64) {
    let yhi = f64_from_bits(f64_to_bits(y) & (!0u64 << 27));
    let ylo = y - yhi;
    let lhi = f64_from_bits(f64_to_bits(hi) & (!0u64 << 27));
    let llo = hi - lhi + lo;
    let ehi = yhi * lhi;
    let elo = ylo * lhi + y * llo;
    (ehi, elo)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn mul_log_fma(y: f64, hi: f64, lo: f64) -> (f64, f64) {
    let ehi = y * hi;
    let elo = y * lo + unsafe { fma_f64(y, hi, -ehi) };
    (ehi, elo)
}

#[inline]
fn pow_exp_generic(x: f64, y: f64) -> f64 {
    let mut ix = f64_to_bits(x);
    if (ix & 0x7ff0_0000_0000_0000u64) == 0 {
        let xn = x * f64_from_bits(0x4330_0000_0000_0000u64); // 2^52
        ix = f64_to_bits(xn);
        ix &= 0x7fff_ffff_ffff_ffffu64;
        ix = ix.wrapping_sub(52u64 << 52);
    } else {
        ix &= 0x7fff_ffff_ffff_ffffu64;
    }

    let mut lo = 0.0;
    let hi = log_inline_generic(ix, &mut lo);
    let (ehi, elo) = mul_log_generic(y, hi, lo);
    super::exp::exp_with_tail_generic(ehi, elo)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
unsafe fn pow_exp_fma(x: f64, y: f64) -> f64 {
    let mut ix = f64_to_bits(x);
    if (ix & 0x7ff0_0000_0000_0000u64) == 0 {
        let xn = x * f64_from_bits(0x4330_0000_0000_0000u64); // 2^52
        ix = f64_to_bits(xn);
        ix &= 0x7fff_ffff_ffff_ffffu64;
        ix = ix.wrapping_sub(52u64 << 52);
    } else {
        ix &= 0x7fff_ffff_ffff_ffffu64;
    }

    let mut lo = 0.0;
    let hi = unsafe { log_inline_fma(ix, &mut lo) };
    let (ehi, elo) = unsafe { mul_log_fma(y, hi, lo) };
    unsafe { super::exp::exp_with_tail_fma(ehi, elo) }
}

#[inline]
fn pow_exp(x: f64, y: f64) -> f64 {
    #[cfg(target_arch = "x86_64")]
    if super::cpu_has_fma() {
        // SAFETY: guarded by CPUID.
        return unsafe { pow_exp_fma(x, y) };
    }
    pow_exp_generic(x, y)
}

#[inline]
fn apply_sign(x: f64, neg: bool) -> f64 {
    if neg {
        f64_from_bits(x.to_bits() | SIGN_BIT)
    } else {
        x
    }
}

#[inline]
pub fn pow(x: f64, y: f64) -> f64 {
    if y == 0.0 {
        return 1.0;
    }
    if x == 1.0 {
        return 1.0;
    }
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }
    if y.is_infinite() {
        let ax = x.abs();
        if ax == 1.0 {
            return f64::NAN;
        }
        if ax > 1.0 {
            return if y.is_sign_positive() {
                f64::INFINITY
            } else {
                0.0
            };
        }
        return if y.is_sign_positive() {
            0.0
        } else {
            f64::INFINITY
        };
    }
    if x == 0.0 {
        let (y_is_int, y_is_odd) = classify_integer(y);
        if y.is_sign_positive() {
            return if y_is_int && x.is_sign_negative() && y_is_odd {
                -0.0
            } else {
                0.0
            };
        }
        return if y_is_int && x.is_sign_negative() && y_is_odd {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }
    if x.is_infinite() {
        let (_, y_is_odd) = classify_integer(y);
        return if y.is_sign_positive() {
            if x.is_sign_negative() && y_is_odd {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            }
        } else if x.is_sign_negative() && y_is_odd {
            -0.0
        } else {
            0.0
        };
    }

    if x > 0.0 {
        return pow_exp(x, y);
    }

    let (y_is_int, y_is_odd) = classify_integer(y);
    if !y_is_int {
        return f64::NAN;
    }
    let sign_neg = y_is_odd;
    let ax = -x;
    let y_abs = y.abs();
    if y_abs < (1u64 << 53) as f64 {
        let yi = y as i64;
        if yi.abs() > POWI_EXP_CUTOFF {
            let res = pow_exp(ax, y);
            return apply_sign(res, sign_neg);
        }
        if yi < 0 && ax > 1.0 {
            let log_ax = ln(ax);
            if log_ax * (-y) > POW_OVERFLOW_LOG {
                let res = pow_exp(ax, y);
                return apply_sign(res, sign_neg);
            }
        }
        let res = powi(ax, yi);
        return apply_sign(res, sign_neg);
    }

    let res = pow_exp(ax, y);
    apply_sign(res, sign_neg)
}

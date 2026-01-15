#![allow(non_camel_case_types)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::unusual_byte_groupings)]

mod cos;
mod exp;
mod log;
mod sin;
mod trig;

pub use cos::cos;
pub use exp::exp;
pub use log::ln;
pub use sin::sin;

// ========= bit helpers =========

#[inline(always)]
fn f64_from_bits(u: u64) -> f64 {
    f64::from_bits(u)
}
#[inline(always)]
fn f64_to_bits(x: f64) -> u64 {
    x.to_bits()
}

#[inline(always)]
fn hi_word(x: f64) -> u32 {
    (f64_to_bits(x) >> 32) as u32
}
#[inline(always)]
fn lo_word(x: f64) -> u32 {
    (f64_to_bits(x) & 0xffff_ffffu64) as u32
}
#[inline(always)]
fn with_hi_lo(hi: u32, lo: u32) -> f64 {
    f64_from_bits(((hi as u64) << 32) | (lo as u64))
}

#[inline(always)]
fn get_exp_bits(u: u64) -> i32 {
    ((u >> 52) & 0x7ff) as i32
}

#[inline(always)]
fn is_nan_bits(u: u64) -> bool {
    (u & 0x7ff0_0000_0000_0000u64) == 0x7ff0_0000_0000_0000u64
        && (u & 0x000f_ffff_ffff_ffffu64) != 0
}
#[inline(always)]
fn is_inf_bits(u: u64) -> bool {
    (u & 0x7fff_ffff_ffff_ffffu64) == 0x7ff0_0000_0000_0000u64
}

/// scalbn(x, n): multiply by 2^n without calling any libm.
#[inline(always)]
fn scalbn(mut x: f64, n: i32) -> f64 {
    let ux = f64_to_bits(x);
    let e = get_exp_bits(ux);
    if e == 0 {
        if x == 0.0 {
            return x;
        }
        // normalize
        x *= f64_from_bits(0x4350_0000_0000_0000u64); // 2^54
        let uy = f64_to_bits(x);
        let ey = get_exp_bits(uy) - 54;
        let ne = ey + n;
        if ne <= 0 {
            // underflow/subnormal
            return x
                * f64_from_bits(((ne + 1023 + 54) as u64) << 52)
                * f64_from_bits(0x3c90_0000_0000_0000u64);
        }
        return f64_from_bits((uy & 0x800f_ffff_ffff_ffffu64) | ((ne as u64) << 52));
    }
    if e == 0x7ff {
        return x;
    }
    let ne = e + n;
    if ne <= 0 {
        if ne <= -52 {
            return 0.0 * x;
        }
        let mant = (ux & 0x000f_ffff_ffff_ffffu64) | 0x0010_0000_0000_0000u64;
        let shift = (1 - ne) as u32;
        let sub = mant >> shift;
        return f64_from_bits((ux & 0x8000_0000_0000_0000u64) | sub);
    }
    if ne >= 0x7ff {
        return if x.is_sign_negative() {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }
    f64_from_bits((ux & 0x800f_ffff_ffff_ffffu64) | ((ne as u64) << 52))
}

/// floor(x) implemented via bit manipulation (no libm).
#[inline(always)]
fn floor_f64(x: f64) -> f64 {
    let u = f64_to_bits(x);
    let sx = u >> 63;
    let e = ((u >> 52) & 0x7ff) as i32;
    if e == 0x7ff {
        return x;
    } // NaN/Inf
    if e == 0 {
        // |x| < 2^-1022
        return if sx == 1 && (u << 1) != 0 { -1.0 } else { 0.0 };
    }
    let j0 = e - 1023;
    if j0 < 0 {
        // |x| < 1
        return if sx == 1 && x != 0.0 { -1.0 } else { 0.0 };
    }
    if j0 >= 52 {
        return x;
    }
    let mask = (1u64 << (52 - j0)) - 1;
    if (u & mask) == 0 {
        return x;
    }
    let mut ui = u & !mask;
    if sx == 1 {
        // negative: floor moves away from zero
        ui = ui.wrapping_add(1u64 << (52 - j0));
    }
    f64_from_bits(ui)
}

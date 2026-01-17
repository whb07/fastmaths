use super::{remainder, rint};

const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;

#[inline(always)]
pub fn remquo(x: f64, y: f64) -> (f64, i32) {
    let hx = x.to_bits();
    let hy = y.to_bits();
    let ax = hx & 0x7fff_ffff_ffff_ffffu64;
    let ay = hy & 0x7fff_ffff_ffff_ffffu64;

    if ay == 0 || ax >= EXP_MASK || ay > EXP_MASK {
        return (f64::NAN, 0);
    }

    let r = remainder(x, y);
    if r.is_nan() {
        return (r, 0);
    }

    let q = rint((x - r) / y);
    let q_int = q as i64;
    let q_abs = if q_int < 0 {
        q_int.wrapping_neg() as u64
    } else {
        q_int as u64
    };
    let mut quo = (q_abs & 0x7) as i32;
    if q_int < 0 {
        quo = -quo;
    }
    (r, quo)
}

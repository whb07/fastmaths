use super::fmod;

const SIGN_MASK: u64 = 0x8000_0000_0000_0000u64;
const EXP_MASK: u64 = 0x7ff0_0000_0000_0000u64;
const MAX_HALF: u64 = 0x7fe0_0000_0000_0000u64;

#[inline(always)]
pub fn remainder(x: f64, y: f64) -> f64 {
    let hx = x.to_bits();
    let hy = y.to_bits();
    let sx = hx >> 63;

    let ax_bits = hx & !SIGN_MASK;
    let ay_bits = hy & !SIGN_MASK;

    if ax_bits >= EXP_MASK || ay_bits > EXP_MASK || ay_bits == 0 {
        return f64::NAN;
    }

    let ax = x.abs();
    let ay = y.abs();
    if ax <= ay * 0.5 {
        return x;
    }
    let mut r = if ay_bits < MAX_HALF {
        let t = fmod(x, y + y).abs();
        if t + t > ay {
            let mut v = t - ay;
            if v + v >= ay {
                v -= ay;
            } else if v == 0.0 {
                v = 0.0;
            }
            v
        } else {
            t
        }
    } else {
        let mut v = x.abs();
        let y_half = ay * 0.5;
        if v > y_half {
            v -= ay;
            if v >= y_half {
                v -= ay;
            } else if v == 0.0 {
                v = 0.0;
            }
        }
        v
    };

    if sx != 0 {
        r = -r;
    }
    r
}

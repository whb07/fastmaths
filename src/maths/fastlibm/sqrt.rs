use super::scalbn_internal;

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "sse2")]
unsafe fn sqrt_sse(x: f64) -> f64 {
    use core::arch::x86_64::{_mm_cvtsd_f64, _mm_set_sd, _mm_sqrt_sd};
    _mm_cvtsd_f64(_mm_sqrt_sd(_mm_set_sd(0.0), _mm_set_sd(x)))
}

#[inline]
fn next_up(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }
    if x == 0.0 {
        return f64::from_bits(1);
    }
    let mut ux = x.to_bits();
    if x.is_sign_negative() {
        ux -= 1;
    } else {
        ux += 1;
    }
    f64::from_bits(ux)
}

#[inline]
fn next_down(x: f64) -> f64 {
    if x.is_nan() || x == f64::NEG_INFINITY {
        return x;
    }
    if x == 0.0 {
        return -f64::from_bits(1);
    }
    let mut ux = x.to_bits();
    if x.is_sign_negative() {
        ux += 1;
    } else {
        ux -= 1;
    }
    f64::from_bits(ux)
}

#[inline]
fn sqrt_fallback(x: f64) -> f64 {
    let mut ax = x;
    let mut scale = 0;
    let mut ux = ax.to_bits();
    if ((ux >> 52) & 0x7ff) == 0 {
        // Normalize subnormals.
        ax = scalbn_internal(ax, 54);
        scale = -27;
        ux = ax.to_bits();
    }

    let mut y = f64::from_bits((ux >> 1) + 0x1ff8_0000_0000_0000);
    y = 0.5 * (y + ax / y);
    y = 0.5 * (y + ax / y);
    y = 0.5 * (y + ax / y);
    y = 0.5 * (y + ax / y);
    y = 0.5 * (y + ax / y);
    y = 0.5 * (y + ax / y);

    if scale != 0 {
        y = scalbn_internal(y, scale);
    }

    let y2 = y * y;
    if y2 < x {
        let y_next = next_up(y);
        if y_next * y_next <= x {
            y = y_next;
        }
    } else if y2 > x {
        let y_prev = next_down(y);
        if y_prev * y_prev >= x {
            y = y_prev;
        }
    }

    y
}

#[inline]
pub fn sqrt(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == 0.0 {
        return x;
    }
    if x.is_infinite() {
        return f64::INFINITY;
    }
    if x < 0.0 {
        return f64::NAN;
    }

    #[cfg(target_arch = "x86_64")]
    unsafe {
        sqrt_sse(x)
    }

    #[cfg(not(target_arch = "x86_64"))]
    {
        sqrt_fallback(x)
    }
}

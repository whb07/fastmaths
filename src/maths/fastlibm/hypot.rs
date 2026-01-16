use super::{scalbn, sqrt};

#[inline]
pub fn hypot(x: f64, y: f64) -> f64 {
    let mut ax = x.abs();
    let mut ay = y.abs();

    if ax.is_nan() || ay.is_nan() {
        return f64::NAN;
    }
    if ax.is_infinite() || ay.is_infinite() {
        return f64::INFINITY;
    }

    if ax < ay {
        core::mem::swap(&mut ax, &mut ay);
    }
    if ax == 0.0 {
        return 0.0;
    }

    let exp = ((ax.to_bits() >> 52) & 0x7ff) as i32 - 1023;
    if exp > 600 {
        let sx = scalbn(ax, -600);
        let sy = scalbn(ay, -600);
        return scalbn(sqrt(sx * sx + sy * sy), 600);
    }
    if exp < -500 {
        let sx = scalbn(ax, 600);
        let sy = scalbn(ay, 600);
        return scalbn(sqrt(sx * sx + sy * sy), -600);
    }

    sqrt(ax * ax + ay * ay)
}

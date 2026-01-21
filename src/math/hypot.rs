//! hypot(x,y) implementation.
//!
//! Scales inputs to avoid overflow/underflow, orders by magnitude, then computes
//! sqrt(x^2+y^2) with guarded squaring. This follows fdlibm-style stable formulas.

use super::{cpu_has_fma, fma_internal, sqrt};

const SCALE: f64 = f64::from_bits(0x1a70_0000_0000_0000); // 2^-600
const LARGE_VAL: f64 = f64::from_bits(0x5fe0_0000_0000_0000); // 2^511
const TINY_VAL: f64 = f64::from_bits(0x2340_0000_0000_0000); // 2^-459
const EPS: f64 = f64::from_bits(0x3c90_0000_0000_0000); // 2^-54

#[inline(always)]
fn kernel(ax: f64, ay: f64) -> f64 {
    if cpu_has_fma() {
        let t1 = ay + ay;
        let t2 = ax - ay;
        if t1 >= ax {
            return sqrt(fma_internal(t1, ax, t2 * t2));
        }
        return sqrt(fma_internal(ax, ax, ay * ay));
    }

    let mut h = sqrt(ax * ax + ay * ay);
    let (t1, t2) = if h <= 2.0 * ay {
        let delta = h - ay;
        let t1 = ax * (2.0 * delta - ax);
        let t2 = (delta - 2.0 * (ax - ay)) * delta;
        (t1, t2)
    } else {
        let delta = h - ax;
        let t1 = 2.0 * delta * (ax - 2.0 * ay);
        let t2 = (4.0 * delta - ay) * ay + delta * delta;
        (t1, t2)
    };
    h -= (t1 + t2) / (2.0 * h);
    h
}

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
    if ax > LARGE_VAL {
        if ay <= ax * EPS {
            return ax + ay;
        }
        return kernel(ax * SCALE, ay * SCALE) / SCALE;
    }
    if ay < TINY_VAL {
        if ax >= ay / EPS {
            return ax + ay;
        }
        let res = kernel(ax / SCALE, ay / SCALE) * SCALE;
        return res;
    }
    if ay <= ax * EPS {
        return ax + ay;
    }

    kernel(ax, ay)
}

//! sin(x) implementation.
//!
//! Thin wrapper around the shared trig range reducer in trig.rs. Uses kernel
//! polynomials on |x|≤pi/4 with Payne–Hanek reduction for huge arguments.

#[inline(always)]
pub fn sin(x: f64) -> f64 {
    super::trig::sin(x)
}

//! cos(x) implementation.
//!
//! Thin wrapper around the shared trig range reducer in trig.rs. Uses kernel
//! polynomials on |x|≤pi/4 with Payne–Hanek reduction for huge arguments.

#[inline(always)]
pub fn cos(x: f64) -> f64 {
    super::trig::cos(x)
}

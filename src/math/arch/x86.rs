#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "fma")]
pub(crate) unsafe fn fma_hw(a: f64, b: f64, c: f64) -> f64 {
    use core::arch::x86_64::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
    _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)))
}

#[cfg(target_arch = "x86")]
#[target_feature(enable = "fma")]
pub(crate) unsafe fn fma_hw(a: f64, b: f64, c: f64) -> f64 {
    use core::arch::x86::{_mm_cvtsd_f64, _mm_fmadd_sd, _mm_set_sd};
    _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)))
}

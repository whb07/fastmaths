// Architecture-specific helpers (e.g., FMA).

#[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
mod x86;
#[cfg(target_arch = "aarch64")]
mod aarch64;

#[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
pub(crate) use x86::fma_hw;

#[cfg(target_arch = "aarch64")]
pub(crate) use aarch64::fma_hw;

#[inline(always)]
pub(crate) unsafe fn fma_hw(a: f64, b: f64, c: f64) -> f64 {
    let out: f64;
    unsafe {
        core::arch::asm!(
            "fmadd {out:d}, {a:d}, {b:d}, {c:d}",
            out = out(vreg) out,
            a = in(vreg) a,
            b = in(vreg) b,
            c = in(vreg) c,
            options(pure, nomem, nostack)
        );
    }
    out
}

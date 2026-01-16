use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_ilogb};

fn fast_ilogb(x: f64) -> f64 {
    fastlibm::ilogb(x) as f64
}

fn glibc_ilogb_f64(x: f64) -> f64 {
    glibc_ilogb(x) as f64
}

fn bench_ilogb(c: &mut Criterion) {
    let inputs = [1e-300, 1e-10, 1e-6, 1.0, 2.0, 1024.0, 1e20];
    let common = gen_range(1024, -100.0, 100.0, 0x1d01);

    let mut group = c.benchmark_group("ilogb/smoke");
    bench_inputs(&mut group, &inputs, fast_ilogb, glibc_ilogb_f64);
    group.finish();

    let mut group = c.benchmark_group("ilogb/common");
    bench_inputs(&mut group, &common, fast_ilogb, glibc_ilogb_f64);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_ilogb(&mut c);
    c.final_summary();
}

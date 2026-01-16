use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_acos};

fn bench_acos(c: &mut Criterion) {
    let inputs = [-1.0, -0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9, 1.0];
    let common = gen_range(1024, -1.0, 1.0, 0x246b);

    let mut group = c.benchmark_group("acos/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::acos, glibc_acos);
    group.finish();

    let mut group = c.benchmark_group("acos/common");
    bench_inputs(&mut group, &common, fastlibm::acos, glibc_acos);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_acos(&mut c);
    c.final_summary();
}

use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_asin};

fn bench_asin(c: &mut Criterion) {
    let inputs = [-1.0, -0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9, 1.0];
    let common = gen_range(1024, -1.0, 1.0, 0x246a);

    let mut group = c.benchmark_group("asin/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::asin, glibc_asin);
    group.finish();

    let mut group = c.benchmark_group("asin/common");
    bench_inputs(&mut group, &common, fastlibm::asin, glibc_asin);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_asin(&mut c);
    c.final_summary();
}

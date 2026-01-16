use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_exp10};

fn bench_exp10(c: &mut Criterion) {
    let inputs = [-10.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 2.0, 10.0, 100.0];
    let common = gen_range(1024, -5.0, 5.0, 0x1601);
    let wide = gen_range(1024, -100.0, 100.0, 0x1602);

    let mut group = c.benchmark_group("exp10/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::exp10, glibc_exp10);
    group.finish();

    let mut group = c.benchmark_group("exp10/common");
    bench_inputs(&mut group, &common, fastlibm::exp10, glibc_exp10);
    group.finish();

    let mut group = c.benchmark_group("exp10/wide");
    bench_inputs(&mut group, &wide, fastlibm::exp10, glibc_exp10);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_exp10(&mut c);
    c.final_summary();
}

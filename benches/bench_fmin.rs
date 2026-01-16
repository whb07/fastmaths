use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs2, configure_criterion, gen_pairs, glibc_fmin};

fn bench_fmin(c: &mut Criterion) {
    let inputs = [(0.0, -0.0), (1.0, 2.0), (-1.0, -2.0), (5.3, 2.1)];
    let pairs = gen_pairs(1024, -100.0, 100.0, 0x1a01);

    let mut group = c.benchmark_group("fmin/smoke");
    bench_inputs2(&mut group, &inputs, fastlibm::fmin, glibc_fmin);
    group.finish();

    let mut group = c.benchmark_group("fmin/common");
    bench_inputs2(&mut group, &pairs, fastlibm::fmin, glibc_fmin);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_fmin(&mut c);
    c.final_summary();
}

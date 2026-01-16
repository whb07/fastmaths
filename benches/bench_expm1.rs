use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, glibc_expm1};

fn bench_expm1(c: &mut Criterion) {
    let inputs = [-50.0, -10.0, -1.0, -1e-6, 0.0, 1e-6, 0.5, 1.0, 10.0, 50.0];
    let mut group = c.benchmark_group("expm1/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::expm1, glibc_expm1);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_expm1(&mut c);
    c.final_summary();
}

use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, glibc_atan};

fn bench_atan(c: &mut Criterion) {
    let inputs = [-1e6, -10.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 10.0, 1e6];
    let mut group = c.benchmark_group("atan/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::atan, glibc_atan);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_atan(&mut c);
    c.final_summary();
}

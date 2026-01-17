use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_ceil};

fn bench_ceil(c: &mut Criterion) {
    let inputs = [
        -1e6, -100.5, -10.9, -1.1, -1.0, -0.9, -0.5, -0.0, 0.0, 0.5, 0.9, 1.0, 1.1, 10.1, 100.5,
        1e6,
    ];
    let common = gen_range(2048, -1e6, 1e6, 0x1213);

    let mut group = c.benchmark_group("ceil/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::ceil, glibc_ceil);
    group.finish();

    let mut group = c.benchmark_group("ceil/common");
    bench_inputs(&mut group, &common, fastlibm::ceil, glibc_ceil);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_ceil(&mut c);
    c.final_summary();
}

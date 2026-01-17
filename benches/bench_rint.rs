use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_rint};

fn bench_rint(c: &mut Criterion) {
    let inputs = [
        -1e6, -100.5, -10.5, -1.5, -1.1, -1.0, -0.5, -0.0, 0.0, 0.5, 1.0, 1.1, 1.5, 10.5, 100.5,
        1e6,
    ];
    let common = gen_range(2048, -1e6, 1e6, 0x1819);

    let mut group = c.benchmark_group("rint/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::rint, glibc_rint);
    group.finish();

    let mut group = c.benchmark_group("rint/common");
    bench_inputs(&mut group, &common, fastlibm::rint, glibc_rint);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_rint(&mut c);
    c.final_summary();
}

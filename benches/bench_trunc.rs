use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_trunc};

fn bench_trunc(c: &mut Criterion) {
    let inputs = [
        -1e6, -100.9, -10.5, -1.1, -1.0, -0.9, -0.5, -0.0, 0.0, 0.5, 0.9, 1.0, 1.1, 10.5, 100.9,
        1e6,
    ];
    let common = gen_range(2048, -1e6, 1e6, 0x1415);

    let mut group = c.benchmark_group("trunc/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::trunc, glibc_trunc);
    group.finish();

    let mut group = c.benchmark_group("trunc/common");
    bench_inputs(&mut group, &common, fastlibm::trunc, glibc_trunc);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_trunc(&mut c);
    c.final_summary();
}

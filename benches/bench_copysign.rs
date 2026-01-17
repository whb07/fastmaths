use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs2, configure_criterion, gen_pairs, glibc_copysign};

fn bench_copysign(c: &mut Criterion) {
    let inputs = [
        (1.0, 1.0),
        (1.0, -1.0),
        (-1.0, 1.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (-0.0, 1.0),
        (1e-300, -1.0),
        (-1e-300, 1.0),
        (1e6, -1e6),
    ];
    let common = gen_pairs(2048, -1e6, 1e6, 0x2425);

    let mut group = c.benchmark_group("copysign/smoke");
    bench_inputs2(&mut group, &inputs, fastlibm::copysign, glibc_copysign);
    group.finish();

    let mut group = c.benchmark_group("copysign/common");
    bench_inputs2(&mut group, &common, fastlibm::copysign, glibc_copysign);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_copysign(&mut c);
    c.final_summary();
}

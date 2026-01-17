use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs2, configure_criterion, gen_pairs, glibc_pow};

fn bench_pow(c: &mut Criterion) {
    let inputs = [
        (2.0, 3.0),
        (2.0, -3.0),
        (10.0, 0.5),
        (0.5, 2.0),
        (-2.0, 3.0),
        (-2.0, 4.0),
        (1.1, 10.0),
        (1e-2, 3.0),
        (1e2, -2.0),
    ];
    let common = gen_pairs(1024, -10.0, 10.0, 0x1234);

    let mut group = c.benchmark_group("pow/smoke");
    bench_inputs2(&mut group, &inputs, fastlibm::pow, glibc_pow);
    group.finish();

    let mut group = c.benchmark_group("pow/common");
    bench_inputs2(&mut group, &common, fastlibm::pow, glibc_pow);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_pow(&mut c);
    c.final_summary();
}

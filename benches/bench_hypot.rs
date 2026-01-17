use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs2, configure_criterion, gen_pairs, glibc_hypot};

fn bench_hypot(c: &mut Criterion) {
    let inputs = [
        (3.0, 4.0),
        (1e-300, 1e-300),
        (1e-100, 1e-100),
        (1e-6, 1e-6),
        (1.0, 1.0),
        (1e6, 1e6),
        (1e100, 1e100),
    ];
    let wide = gen_pairs(1024, -1e6, 1e6, 0x2718);
    let huge = gen_pairs(1024, -1e300, 1e300, 0x3141_5926);

    let mut group = c.benchmark_group("hypot/smoke");
    bench_inputs2(&mut group, &inputs, fastlibm::hypot, glibc_hypot);
    group.finish();

    let mut group = c.benchmark_group("hypot/wide");
    bench_inputs2(&mut group, &wide, fastlibm::hypot, glibc_hypot);
    group.finish();

    let mut group = c.benchmark_group("hypot/huge");
    bench_inputs2(&mut group, &huge, fastlibm::hypot, glibc_hypot);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_hypot(&mut c);
    c.final_summary();
}

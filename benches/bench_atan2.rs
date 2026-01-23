use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs2, configure_criterion, gen_pairs, glibc_atan2};

fn bench_atan2(c: &mut Criterion) {
    let inputs = [
        (0.0, 1.0),
        (1.0, 0.0),
        (-1.0, 0.0),
        (1.0, 1.0),
        (-1.0, 1.0),
        (1.0, -1.0),
        (-1.0, -1.0),
        (1e-6, 1.0),
        (1.0, 1e-6),
        (1e6, 1.0),
        (1.0, 1e6),
    ];
    let common = gen_pairs(1024, -1e3, 1e3, 0x3141);
    let wide = gen_pairs(1024, -1e6, 1e6, 0x5926);

    let mut group = c.benchmark_group("atan2/smoke");
    bench_inputs2(&mut group, &inputs, fastmaths::atan2, glibc_atan2);
    group.finish();

    let mut group = c.benchmark_group("atan2/common");
    bench_inputs2(&mut group, &common, fastmaths::atan2, glibc_atan2);
    group.finish();

    let mut group = c.benchmark_group("atan2/wide");
    bench_inputs2(&mut group, &wide, fastmaths::atan2, glibc_atan2);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_atan2(&mut c);
    c.final_summary();
}

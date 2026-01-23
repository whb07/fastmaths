use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_round};

fn bench_round(c: &mut Criterion) {
    let inputs = [
        -1e6, -100.5, -10.5, -1.5, -1.1, -1.0, -0.5, -0.0, 0.0, 0.5, 1.0, 1.1, 1.5, 10.5, 100.5,
        1e6,
    ];
    let common = gen_range(2048, -1e6, 1e6, 0x1617);

    let mut group = c.benchmark_group("round/smoke");
    bench_inputs(&mut group, &inputs, fastmaths::round, glibc_round);
    group.finish();

    let mut group = c.benchmark_group("round/common");
    bench_inputs(&mut group, &common, fastmaths::round, glibc_round);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_round(&mut c);
    c.final_summary();
}

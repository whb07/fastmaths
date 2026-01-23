use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs_i64, configure_criterion, gen_range, glibc_lround};

fn bench_lround(c: &mut Criterion) {
    let inputs = [
        -1e6, -100.5, -10.5, -1.5, -1.1, -1.0, -0.5, -0.0, 0.0, 0.5, 1.0, 1.1, 1.5, 10.5, 100.5,
        1e6,
    ];
    let common = gen_range(2048, -1e6, 1e6, 0x2021);

    let mut group = c.benchmark_group("lround/smoke");
    bench_inputs_i64(&mut group, &inputs, fastmaths::lround, glibc_lround);
    group.finish();

    let mut group = c.benchmark_group("lround/common");
    bench_inputs_i64(&mut group, &common, fastmaths::lround, glibc_lround);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_lround(&mut c);
    c.final_summary();
}

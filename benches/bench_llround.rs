use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs_i64, configure_criterion, gen_range, glibc_llround};

fn bench_llround(c: &mut Criterion) {
    let inputs = [
        -1e6, -100.5, -10.5, -1.5, -1.1, -1.0, -0.5, -0.0, 0.0, 0.5, 1.0, 1.1, 1.5, 10.5, 100.5,
        1e6,
    ];
    let common = gen_range(2048, -1e6, 1e6, 0x2223);

    let mut group = c.benchmark_group("llround/smoke");
    bench_inputs_i64(&mut group, &inputs, fastmaths::llround, glibc_llround);
    group.finish();

    let mut group = c.benchmark_group("llround/common");
    bench_inputs_i64(&mut group, &common, fastmaths::llround, glibc_llround);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_llround(&mut c);
    c.final_summary();
}

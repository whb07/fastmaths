use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs_i64, configure_criterion, gen_range, glibc_llrint};

fn bench_llrint(c: &mut Criterion) {
    let inputs = [
        -1e6, -100.5, -10.5, -1.5, -1.1, -1.0, -0.5, -0.0, 0.0, 0.5, 1.0, 1.1, 1.5, 10.5, 100.5,
        1e6,
    ];
    let common = gen_range(2048, -1e6, 1e6, 0x1e1f);

    let mut group = c.benchmark_group("llrint/smoke");
    bench_inputs_i64(&mut group, &inputs, fastmaths::llrint, glibc_llrint);
    group.finish();

    let mut group = c.benchmark_group("llrint/common");
    bench_inputs_i64(&mut group, &common, fastmaths::llrint, glibc_llrint);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_llrint(&mut c);
    c.final_summary();
}

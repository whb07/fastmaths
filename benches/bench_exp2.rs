use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, glibc_exp2};

fn bench_exp2(c: &mut Criterion) {
    let inputs = [
        -1074.0, -100.0, -10.0, -1.0, -1e-6, 0.0, 1e-6, 0.5, 1.0, 2.0, 10.0, 100.0, 1023.0,
    ];
    let mut group = c.benchmark_group("exp2/smoke");
    bench_inputs(&mut group, &inputs, fastmaths::exp2, glibc_exp2);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_exp2(&mut c);
    c.final_summary();
}

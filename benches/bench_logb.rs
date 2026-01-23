use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_logb};

fn bench_logb(c: &mut Criterion) {
    let inputs = [1e-300, 1e-10, 1e-6, 1.0, 2.0, 1024.0, 1e20];
    let common = gen_range(1024, -100.0, 100.0, 0x1c01);

    let mut group = c.benchmark_group("logb/smoke");
    bench_inputs(&mut group, &inputs, fastmaths::logb, glibc_logb);
    group.finish();

    let mut group = c.benchmark_group("logb/common");
    bench_inputs(&mut group, &common, fastmaths::logb, glibc_logb);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_logb(&mut c);
    c.final_summary();
}

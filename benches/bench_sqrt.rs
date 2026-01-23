use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_sqrt};

fn bench_sqrt(c: &mut Criterion) {
    let inputs = [0.0, 1e-300, 1e-6, 0.5, 1.0, 2.0, 10.0, 1e6, 1e300];
    let common = gen_range(1024, 0.0, 1e6, 0x4242);
    let huge = gen_range(1024, 0.0, 1e300, 0x7777);

    let mut group = c.benchmark_group("sqrt/smoke");
    bench_inputs(&mut group, &inputs, fastmaths::sqrt, glibc_sqrt);
    group.finish();

    let mut group = c.benchmark_group("sqrt/common");
    bench_inputs(&mut group, &common, fastmaths::sqrt, glibc_sqrt);
    group.finish();

    let mut group = c.benchmark_group("sqrt/huge");
    bench_inputs(&mut group, &huge, fastmaths::sqrt, glibc_sqrt);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_sqrt(&mut c);
    c.final_summary();
}

use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_cbrt};

fn bench_cbrt(c: &mut Criterion) {
    let inputs = [0.0, -1e-300, -1.0, -8.0, 1.0, 8.0, 1e-6, 1e6, 1e300];
    let common = gen_range(1024, -1e6, 1e6, 0x2468_ace0);
    let huge = gen_range(1024, -1e300, 1e300, 0x1357_9bdf);

    let mut group = c.benchmark_group("cbrt/smoke");
    bench_inputs(&mut group, &inputs, fastmaths::cbrt, glibc_cbrt);
    group.finish();

    let mut group = c.benchmark_group("cbrt/common");
    bench_inputs(&mut group, &common, fastmaths::cbrt, glibc_cbrt);
    group.finish();

    let mut group = c.benchmark_group("cbrt/huge");
    bench_inputs(&mut group, &huge, fastmaths::cbrt, glibc_cbrt);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_cbrt(&mut c);
    c.final_summary();
}

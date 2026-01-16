use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_asinh};

fn bench_asinh(c: &mut Criterion) {
    let inputs = [-100.0, -10.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 10.0, 100.0];
    let common = gen_range(1024, -5.0, 5.0, 0x1101);
    let wide = gen_range(1024, -100.0, 100.0, 0x1102);

    let mut group = c.benchmark_group("asinh/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::asinh, glibc_asinh);
    group.finish();

    let mut group = c.benchmark_group("asinh/common");
    bench_inputs(&mut group, &common, fastlibm::asinh, glibc_asinh);
    group.finish();

    let mut group = c.benchmark_group("asinh/wide");
    bench_inputs(&mut group, &wide, fastlibm::asinh, glibc_asinh);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_asinh(&mut c);
    c.final_summary();
}

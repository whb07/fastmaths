use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_cosh};

fn bench_cosh(c: &mut Criterion) {
    let inputs = [-20.0, -10.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 10.0, 20.0];
    let common = gen_range(1024, -5.0, 5.0, 0x357c);
    let wide = gen_range(1024, -20.0, 20.0, 0x357d);

    let mut group = c.benchmark_group("cosh/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::cosh, glibc_cosh);
    group.finish();

    let mut group = c.benchmark_group("cosh/common");
    bench_inputs(&mut group, &common, fastlibm::cosh, glibc_cosh);
    group.finish();

    let mut group = c.benchmark_group("cosh/wide");
    bench_inputs(&mut group, &wide, fastlibm::cosh, glibc_cosh);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_cosh(&mut c);
    c.final_summary();
}

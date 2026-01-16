use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_erfc};

fn bench_erfc(c: &mut Criterion) {
    let inputs = [-6.0, -3.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 3.0, 6.0];
    let common = gen_range(1024, -3.0, 3.0, 0x1501);
    let wide = gen_range(1024, -6.0, 6.0, 0x1502);

    let mut group = c.benchmark_group("erfc/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::erfc, glibc_erfc);
    group.finish();

    let mut group = c.benchmark_group("erfc/common");
    bench_inputs(&mut group, &common, fastlibm::erfc, glibc_erfc);
    group.finish();

    let mut group = c.benchmark_group("erfc/wide");
    bench_inputs(&mut group, &wide, fastlibm::erfc, glibc_erfc);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_erfc(&mut c);
    c.final_summary();
}

use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_atanh};

fn bench_atanh(c: &mut Criterion) {
    let inputs = [-0.99, -0.5, -1e-6, 0.0, 1e-6, 0.5, 0.99];
    let common = gen_range(1024, -0.9, 0.9, 0x1301);
    let tight = gen_range(1024, -0.999, 0.999, 0x1302);

    let mut group = c.benchmark_group("atanh/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::atanh, glibc_atanh);
    group.finish();

    let mut group = c.benchmark_group("atanh/common");
    bench_inputs(&mut group, &common, fastlibm::atanh, glibc_atanh);
    group.finish();

    let mut group = c.benchmark_group("atanh/tight");
    bench_inputs(&mut group, &tight, fastlibm::atanh, glibc_atanh);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_atanh(&mut c);
    c.final_summary();
}

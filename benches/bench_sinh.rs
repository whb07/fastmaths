use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_sinh};

fn bench_sinh(c: &mut Criterion) {
    let inputs = [-20.0, -10.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 10.0, 20.0];
    let common = gen_range(1024, -5.0, 5.0, 0x357a);
    let wide = gen_range(1024, -20.0, 20.0, 0x357b);

    let mut group = c.benchmark_group("sinh/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::sinh, glibc_sinh);
    group.finish();

    let mut group = c.benchmark_group("sinh/common");
    bench_inputs(&mut group, &common, fastlibm::sinh, glibc_sinh);
    group.finish();

    let mut group = c.benchmark_group("sinh/wide");
    bench_inputs(&mut group, &wide, fastlibm::sinh, glibc_sinh);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_sinh(&mut c);
    c.final_summary();
}

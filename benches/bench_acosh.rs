use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_acosh};

fn bench_acosh(c: &mut Criterion) {
    let inputs = [1.0, 1.0 + 1e-12, 1.1, 1.5, 2.0, 10.0, 100.0, 1e6];
    let common = gen_range(1024, 1.0, 10.0, 0x1201);
    let wide = gen_range(1024, 1.0, 1e6, 0x1202);

    let mut group = c.benchmark_group("acosh/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::acosh, glibc_acosh);
    group.finish();

    let mut group = c.benchmark_group("acosh/common");
    bench_inputs(&mut group, &common, fastlibm::acosh, glibc_acosh);
    group.finish();

    let mut group = c.benchmark_group("acosh/wide");
    bench_inputs(&mut group, &wide, fastlibm::acosh, glibc_acosh);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_acosh(&mut c);
    c.final_summary();
}

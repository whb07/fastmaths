use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_log1p};

fn bench_log1p(c: &mut Criterion) {
    let inputs = [-0.9, -0.5, -0.1, -1e-6, 0.0, 1e-6, 0.1, 1.0, 10.0, 1e6];
    let common = gen_range(1024, -0.9, 10.0, 0x135a);

    let mut group = c.benchmark_group("log1p/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::log1p, glibc_log1p);
    group.finish();

    let mut group = c.benchmark_group("log1p/common");
    bench_inputs(&mut group, &common, fastlibm::log1p, glibc_log1p);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_log1p(&mut c);
    c.final_summary();
}

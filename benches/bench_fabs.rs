use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_fabs};

fn bench_fabs(c: &mut Criterion) {
    let inputs = [-1e6, -100.5, -1.0, -0.0, 0.0, 0.5, 1.0, 100.5, 1e6];
    let common = gen_range(2048, -1e6, 1e6, 0x2627);

    let mut group = c.benchmark_group("fabs/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::fabs, glibc_fabs);
    group.finish();

    let mut group = c.benchmark_group("fabs/common");
    bench_inputs(&mut group, &common, fastlibm::fabs, glibc_fabs);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_fabs(&mut c);
    c.final_summary();
}

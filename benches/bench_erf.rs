use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_erf};

fn bench_erf(c: &mut Criterion) {
    let inputs = [-6.0, -3.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 3.0, 6.0];
    let common = gen_range(1024, -3.0, 3.0, 0x1401);
    let wide = gen_range(1024, -6.0, 6.0, 0x1402);

    let mut group = c.benchmark_group("erf/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::erf, glibc_erf);
    group.finish();

    let mut group = c.benchmark_group("erf/common");
    bench_inputs(&mut group, &common, fastlibm::erf, glibc_erf);
    group.finish();

    let mut group = c.benchmark_group("erf/wide");
    bench_inputs(&mut group, &wide, fastlibm::erf, glibc_erf);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_erf(&mut c);
    c.final_summary();
}
